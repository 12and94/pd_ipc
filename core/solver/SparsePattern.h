// core/solver/SparsePattern.h
// 稀疏结构的最小描述：一个"3x3 非零块"的位置（顶点 i 的行、顶点 j 的列）。
//
// 单独放在这里，避免 Mesh 与求解器接口互相包含（求解器只需要结构，不需要网格内容）。
#pragma once

namespace pd {

struct BlockEntry {
  int row = -1;  ///< 块行 = 顶点索引
  int col = -1;  ///< 块列 = 顶点索引
};

}  // namespace pd
