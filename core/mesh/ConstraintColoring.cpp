// core/mesh/ConstraintColoring.cpp
// 约束着色的实现：并行的极大独立集迭代（Luby 风格），权重是边索引的确定性哈希。
#include "core/mesh/ConstraintColoring.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>

#include "core/math/Parallel.h"
#include "core/mesh/Mesh.h"

namespace pd {

namespace {

constexpr uint32_t kUncolored = 0xFFFFFFFFu;

/// 确定性优先级：边索引 → 32 位哈希（splitmix64 的最后一步）。
///
/// **绝对不能用 `rand()`**：着色结果必须只由拓扑决定，否则每次运行的颜色→边映射都不同，
/// 求和顺序随之改变，`docs/plan.md` M2 要的"与线程数无关、逐位可复现"就无从谈起。
/// 用"索引的哈希"还有一个好处：它不依赖遍历顺序，也不依赖线程数。
inline uint32_t priorityOf(uint32_t edgeIndex) {
  uint64_t x = static_cast<uint64_t>(edgeIndex) * 0x9E3779B97F4A7C15ULL + 0x2545F4914F6CDD1DULL;
  x ^= x >> 30;
  x *= 0xBF58476D1CE4E5B9ULL;
  x ^= x >> 27;
  x *= 0x94D049BB133111EBULL;
  x ^= x >> 31;
  return static_cast<uint32_t>(x >> 32);
}

/// 顶点 → 关联边（CSR，拓扑级）。
struct VertexEdgeAdjacency {
  std::vector<uint32_t> start;  ///< 大小 vertexCount+1
  std::vector<uint32_t> edges;
};

VertexEdgeAdjacency buildAdjacency(const Mesh& mesh) {
  const int nv = mesh.vertexCount();
  const int ne = mesh.edgeCount();
  VertexEdgeAdjacency adj;
  adj.start.assign(static_cast<std::size_t>(nv) + 1, 0);

  for (const Edge& e : mesh.edges) {
    adj.start[static_cast<std::size_t>(e.a) + 1] += 1;
    adj.start[static_cast<std::size_t>(e.b) + 1] += 1;
  }
  for (int v = 0; v < nv; ++v) {
    adj.start[static_cast<std::size_t>(v) + 1] += adj.start[static_cast<std::size_t>(v)];
  }
  adj.edges.assign(static_cast<std::size_t>(2) * static_cast<std::size_t>(ne), 0);
  std::vector<uint32_t> cursor(adj.start.begin(), adj.start.end() - 1);
  for (int i = 0; i < ne; ++i) {
    const Edge& e = mesh.edges[static_cast<std::size_t>(i)];
    adj.edges[cursor[static_cast<std::size_t>(e.a)]++] = static_cast<uint32_t>(i);
    adj.edges[cursor[static_cast<std::size_t>(e.b)]++] = static_cast<uint32_t>(i);
  }
  return adj;
}

/// 边 e 的两个端点上的关联边（用于冲突判定）。
struct NeighborScanner {
  const Mesh& mesh;
  const VertexEdgeAdjacency& adj;

  /// 遍历 e 的所有邻居（跳过自身）。
  template <typename Fn>
  void forEach(uint32_t e, Fn&& fn) const {
    const Edge& edge = mesh.edges[e];
    const uint32_t endpoints[2] = {static_cast<uint32_t>(edge.a), static_cast<uint32_t>(edge.b)};
    for (const uint32_t v : endpoints) {
      for (uint32_t k = adj.start[v]; k < adj.start[v + 1]; ++k) {
        const uint32_t f = adj.edges[k];
        if (f != e) fn(f);
      }
    }
  }
};

/// (优先级, 边索引) 是否严格优于 (优先级, 边索引) —— 用全序消除并列，避免死循环。
inline bool better(uint32_t ea, uint32_t eb, const std::vector<uint32_t>& prio) {
  return prio[ea] > prio[eb] || (prio[ea] == prio[eb] && ea < eb);
}

}  // namespace

void buildConstraintColoring(const Mesh& mesh, ConstraintColoring& out) {
  const int ne = mesh.edgeCount();
  out.colorOf.assign(static_cast<std::size_t>(ne), kUncolored);
  out.start.clear();
  out.order.clear();
  out.colorCount = 0;
  out.built = true;
  if (ne <= 0) {
    out.start.push_back(0);
    return;
  }

  const VertexEdgeAdjacency adj = buildAdjacency(mesh);
  const NeighborScanner scan{mesh, adj};

  std::vector<uint32_t> prio(static_cast<std::size_t>(ne));
  for (int i = 0; i < ne; ++i) prio[static_cast<std::size_t>(i)] = priorityOf(static_cast<uint32_t>(i));

  // 未着色的边，始终保持**索引升序** —— "结果与线程数无关"的一部分：
  // 每轮的处理顺序、以及同色内的边序，都由索引唯一决定。
  std::vector<uint32_t> remaining(static_cast<std::size_t>(ne));
  for (int i = 0; i < ne; ++i) remaining[static_cast<std::size_t>(i)] = static_cast<uint32_t>(i);
  std::vector<uint32_t> nextRemaining;
  std::vector<uint32_t> candidate(static_cast<std::size_t>(ne), kUncolored);

  int maxColor = -1;
  bool firstRound = true;
  while (!remaining.empty()) {
    // ---- 并行的一趟：为每条剩余边判断"本轮能否着色"，并算出候选颜色 ----
    // 只读 colorOf（本轮之前的状态），逐元素写 candidate ⇒ 无竞争、与线程数无关。
    parallelFor(remaining.size(), [&](std::size_t b, std::size_t e) {
      for (std::size_t k = b; k < e; ++k) {
        const uint32_t ed = remaining[k];
        candidate[ed] = kUncolored;

        bool allNeighborsColored = true;
        bool isLocalMax = true;
        scan.forEach(ed, [&](uint32_t f) {
          if (out.colorOf[f] == kUncolored) {
            allNeighborsColored = false;
            if (better(f, ed, prio)) isLocalMax = false;
          }
        });

        // 可着色的两个条件（Jones–Plassmann）：
        //   ① 邻居都已着色 —— 现在就可以取"最小可用颜色"；
        //   ② 我在剩余集合里是局部极大 —— 说明我的邻居都比我差，先给我颜色不会挡别人。
        // 同一轮里被着色的边两两不相邻（相邻两条不可能同时满足 ① 或 ②），
        // 所以"最小可用颜色"只需避开**上一轮及更早**的邻居颜色即可。
        if (!allNeighborsColored && !isLocalMax) continue;

        // 最小可用颜色：扫一遍已着色邻居用掉的颜色。
        // 不用位掩码是为了不受 64 色上限约束（不规则网格的 Δ 可能很大）；
        // 内层最多被扫描 Δ 次 × 颜色数，实测量级可忽略。
        uint32_t pick = 0;
        bool clash = true;
        while (clash) {
          clash = false;
          scan.forEach(ed, [&](uint32_t f) {
            if (!clash && out.colorOf[f] == pick) {
              clash = true;
              ++pick;
            }
          });
        }
        candidate[ed] = pick;
      }
    });

    // ---- 应用（串行 O(E)：一次性成本，不值得为它引入并行归约） ----
    bool anyColored = false;
    nextRemaining.clear();
    nextRemaining.reserve(remaining.size());
    for (const uint32_t ed : remaining) {
      const uint32_t c = candidate[ed];
      if (c != kUncolored) {
        out.colorOf[ed] = c;
        maxColor = std::max(maxColor, static_cast<int>(c));
        anyColored = true;
      } else {
        nextRemaining.push_back(ed);
      }
    }
    if (!anyColored) {
      // 理论上不可达：剩余集合里优先级最大的那条边必然是局部极大。
      // 真走到这里说明着色逻辑被改坏了，宁可停下也不要产出"不恰当的着色"。
      std::fprintf(stderr,
                   "[buildConstraintColoring] 某一轮没有着色任何边（剩余 %zu 条）。\n"
                   "  这不可能发生 —— 剩余集合中优先级最大者必是局部极大。请检查判定逻辑。\n",
                   remaining.size());
      std::fflush(stderr);
      std::abort();
    }
    remaining.swap(nextRemaining);
    firstRound = false;
  }
  (void)firstRound;

  out.colorCount = maxColor + 1;
  out.start.assign(static_cast<std::size_t>(out.colorCount) + 1, 0);

  // 计数排序把边按颜色分组（同色内仍是索引升序）。
  // 注意中间计数用 start 的 1..C 位置，避免再开一个数组。
  for (const uint32_t c : out.colorOf) {
    out.start[static_cast<std::size_t>(c) + 1] += 1;
  }
  for (int c = 0; c < out.colorCount; ++c) {
    out.start[static_cast<std::size_t>(c) + 1] += out.start[static_cast<std::size_t>(c)];
  }
  out.order.assign(static_cast<std::size_t>(ne), 0);
  std::vector<uint32_t> cursor(out.start.begin(), out.start.end() - 1);
  for (int i = 0; i < ne; ++i) {
    const uint32_t c = out.colorOf[static_cast<std::size_t>(i)];
    out.order[cursor[static_cast<std::size_t>(c)]++] = static_cast<uint32_t>(i);
  }

  // 每个颜色类内部**按最小端点顶点排序**（并列时按边索引，保证确定性）。
  //
  // 为什么必须这么做（实测踩到的坑）：散射按颜色分轮执行时，同色边会**并发写 b**。
  // 如果同色内还是边序号序，线程拿到的分块散布在整条边表上 ⇒ 它们写的顶点也是散布的，
  // 相邻顶点常落在同一条 64 字节 cache line 里 ⇒ 该 line 在核心之间来回搬（假共享）。
  // 40×40（dim=4800，仅 600 条 cache line）实测散射比归约方案慢 1.87×，主因就是它。
  // 按顶点序排列后，每个线程的分块对应一段**连续顶点区**，
  // 只有区段边界的少量 line 还会争抢（200×200 上本来就赢 5×，加上这条更稳）。
  //
  // 这不改变结果：同一顶点在每一色里最多被写一次，同色内的先后顺序与浮点结果无关。
  for (int c = 0; c < out.colorCount; ++c) {
    const auto first = out.order.begin() + static_cast<std::ptrdiff_t>(out.start[static_cast<std::size_t>(c)]);
    const auto last = out.order.begin() + static_cast<std::ptrdiff_t>(out.start[static_cast<std::size_t>(c) + 1]);
    std::sort(first, last, [&mesh](uint32_t x, uint32_t y) {
      const uint32_t mx = static_cast<uint32_t>(std::min(mesh.edges[x].a, mesh.edges[x].b));
      const uint32_t my = static_cast<uint32_t>(std::min(mesh.edges[y].a, mesh.edges[y].b));
      return mx != my ? mx < my : x < y;
    });
  }
}

}  // namespace pd
