// _perf/check_report.js —— 在 Node 里真正跑一遍报告里的渲染代码，检查生成的 SVG 几何是否正常。
// 用法：node _perf/check_report.js
// 为什么需要它：报告是"自包含 HTML + 客户端 JS 渲染"，用浏览器肉眼只能看到结果；
// 这里用最小 DOM 桩把 renderStack() 的真代码跑一遍，直接量条形的右端是否越界、
// 最大段是不是全局回代 —— 之前的单位混用 bug（累计 ms 与 ms/子步同刻度）就是这么被发现的。
const fs = require('fs');
const path = require('path');

const file = path.join(path.resolve(__dirname, '..'), 'build', '_perf', 'perf-report.html');
const html = fs.readFileSync(file, 'utf8');
const src = html.match(/<script>([\s\S]*)<\/script>/)[1];

const store = {};
global.document = {
  getElementById: (id) => ({
    set innerHTML(v) { store[id] = v; },
    get innerHTML() { return store[id] || ''; },
  }),
};
new Function(src)();

const svg = store.stack || '';
const rects = [...svg.matchAll(/<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)"/g)]
  .map((m) => ({ x: +m[1], w: +m[3] }));
const labelW = 132, axisW = 1180 - 132 - 210;
const right = Math.max(...rects.map((r) => r.x + r.w));
const reds = (svg.match(/#d62828/g) || []).length;

console.log('段数          :', rects.length);
console.log('最右端点      :', right.toFixed(1), '/ 轴右端', labelW + axisW,
            right > labelW + axisW + 0.5 ? '  ← ⚠ 溢出' : '  ✓ 未越界');
console.log('红色段（回代）:', reds, '(期望 20 = 4 工况 × 5 线程)');
console.log('组标尺文本    :', (svg.match(/标尺：最长条/g) || []).length, '  “最大：”标注:', (svg.match(/最大：/g) || []).length);
console.log('一次性标注    :', svg.includes('组装+分解（一次性）'));

// ---------------------------------------------------------------------------
// 新增板块（剪切 / 弯曲 / 采样变体 / 惯性右端）的自检。
// 这些板块同样是"客户端 JS 渲染"，肉眼看浏览器只能看到结果；这里用同一个 DOM 桩把它们
// 也跑一遍，并抽几个关键数字做**反查**（例如隔行采样的应变必须与 off 几乎相同 ——
// 那是"退化"的量化特征，解析错列就会立刻体现出来）。
// ---------------------------------------------------------------------------
const NEW_IDS = ['shbars', 'shtable', 'bdbars', 'bdtable', 'varbars', 'vartable', 'inerttable'];
const empty = NEW_IDS.filter((id) => !store[id] || store[id].length < 40);
console.log('\n新板块渲染    :', NEW_IDS.length + ' 个 id，非空 '
            + (NEW_IDS.length - empty.length) + '/' + NEW_IDS.length
            + (empty.length ? ('  ← ⚠ 空：' + empty.join(', ')) : '  ✓'));

// 条形右端的上限 = 左标签宽 + 条宽上限。
//   · pairs()      （剪切/弯曲）: W=1180, labelW=200, padR=150 ⇒ 条宽上限 830 ⇒ 右端 1030
//   · simpleBars() （采样变体）  : W=1180, labelW=200, padR=190 ⇒ 条宽上限 790 ⇒ 右端  990
const axisNew = 200 + (1180 - 200 - 190);
const axisPair = 200 + (1180 - 200 - 150);
for (const [id, axis] of [['shbars', axisPair], ['bdbars', axisPair], ['varbars', axisNew]]) {
  const s = store[id] || '';
  const rs = [...s.matchAll(/<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)"/g)]
    .map((m) => ({ x: +m[1], w: +m[3] }));
  const r = rs.length ? Math.max(...rs.map((x) => x.x + x.w)) : 0;
  console.log('  ' + id.padEnd(10) + ': ' + String(rs.length).padStart(2) + ' 条，最右 '
              + r.toFixed(1) + ' / 轴右端 ' + axis + (r <= axis + 0.5 ? '  ✓' : '  ← ⚠ 溢出'));
}
const cell = '<td[^>]*>([^<]*)<\\/td>';
const rowRe = (label) => new RegExp('<tr><td>' + label + '<\\/td>' + cell + cell + cell + cell);

// 剪切：投影+散射 应当明显变大（那条通路整条都参与）
const sh = (store.shtable || '').match(/<tr><td>投影\+散射<\/td><td>([\d.]+)<\/td><td>([\d.]+)<\/td><td[^>]*>([+-][\d.]+)%/);
console.log('  剪切 投影+散射:', sh ? (sh[1] + ' → ' + sh[2] + ' ms (' + sh[3] + '%)  '
  + (Number(sh[3]) > 20 ? '✓ 明显变大' : '← ⚠ 变化过小')) : '← ⚠ 没抓到');
// 弯曲：单子步 与 全局回代 都该涨；投影+散射 几乎不变（结构性收益）
const bd = (store.bdtable || '').match(/<tr><td>单子步\(ms\)（每子步）<\/td><td>([\d.]+)<\/td><td>([\d.]+)<\/td><td[^>]*>([+-][\d.]+)%/);
const bdPS = (store.bdtable || '').match(/<tr><td>投影\+散射<\/td><td>([\d.]+)<\/td><td>([\d.]+)<\/td><td[^>]*>([+-][\d.]+)%/);
console.log('  弯曲 单子步  :', bd ? (bd[1] + ' → ' + bd[2] + ' ms (' + bd[3] + '%)') : '← ⚠ 没抓到');
console.log('  弯曲 投影+散射:', bdPS ? (bdPS[1] + ' → ' + bdPS[2] + ' ms (' + bdPS[3] + '%)  '
  + (Math.abs(Number(bdPS[3])) < 10 ? '✓ 几乎不变（弯曲无局部步/散射）' : '← ⚠ 超出预期')) : '← ⚠ 没抓到');

// 采样变体：14 个配置、6 个"退化"判定、隔行应变 ≈ off
const vt = store.vartable || '';
const cfgCount = (vt.match(/<tr>/g) || []).length - 1;
const degrade = (vt.match(/class="bad">退化/g) || []).length;   // 只数判定单元格（脚注里也出现过"退化"二字）
// 变体表每行 10 格：名字 刚度 stencil 单子步 vsOff L nnz 因子 nnz 应变max 弹性能 结论
// ⇒ 跳过 7 格（cell 自带一个捕获组），第 7 个捕获组就是"应变 max"。
// ⚠️ 不能用 "跳过 7 格后再 ([\\d.]+)"：跳过 7 格后光标停在<下一个 <td> 上，那个模式永远匹配不上。
const strainOf = (name) => {
  const esc = name.replace(/[.*+?^${}()|[\]\\/]/g, '\\$&');
  const m = vt.match(new RegExp('<tr><td>' + esc + '<\\/td>' + cell + cell + cell + cell + cell + cell + cell));
  return m ? Number(m[7]) : null;
};
const offStrain = strainOf('off');
const s2Strain = strainOf('std/s2');
console.log('  采样变体 行数:', cfgCount, '(期望 14)   标"退化"的行:', degrade, '(期望 6)');
console.log('  隔行退化反查  :', offStrain && s2Strain
  ? ('off ' + offStrain + ' vs std/s2 ' + s2Strain + '  相对差 '
     + (100 * Math.abs(s2Strain - offStrain) / offStrain).toFixed(1) + '%  '
     + (Math.abs(s2Strain - offStrain) / offStrain < 0.02 ? '✓ 几乎相同（退化成立）' : '← ⚠ 差异过大'))
  : '← ⚠ 没抓到');

// 惯性右端：3 档线程、1 线程那档必须有"理论差值 = 0"的标注
const it = store.inerttable || '';
console.log('  惯性右端 行数:', (it.match(/<tr>/g) || []).length - 1, '(期望 3)   '
            + (it.includes('机制失效档') ? '✓ 含 1 线程对照说明' : '← ⚠ 缺少 1 线程对照说明'));

// ---------------------------------------------------------------------------
// 四套约束集的阶段构成（§十）：3 个工况 × 4 套配置。做两类反查：
//  ① 结构故事：只加弯曲时"投影+散射"占比几乎不动、"全局回代"占比明显上跳；
//  ② 反直觉结论：同一工况下"剪切+弯曲"的 L nnz 更多、但**因子 nnz 更少**（确定性指标，
//     不是时间噪声）—— 这正是"两者不比只加弯曲慢"的原因。
// ---------------------------------------------------------------------------
const csRows = [];
for (const chunk of (store.cstable || '').split('<table>').slice(1)) {
  for (const c of chunk.matchAll(/<tr>([\s\S]*?)<\/tr>/g)) {
    const tds = [...c[1].matchAll(/<td[^>]*>([\s\S]*?)<\/td>/g)].map((m) => m[1]);
    if (tds.length < 10) continue;
    // 注意：没开弯曲的配置，stencil 列是 '—'（不是数字）⇒ num 必须容错
    const num = (s) => { const m = String(s).match(/[\d.]+/); return m ? Number(m[0]) : 0; };
    const share = (s) => Number(String(s).replace(/[\s\S]*<span[^>]*>/, '').match(/[\d.]+/)[0]);
    // 列序号要按 STAGES 的顺序数：0..6 是元信息，7=predict 8=projScatter 9=assemble
    // 10=solve 11=unpackScan 12=velocity 13=other（数错一列就会得到离谱的"占比变化"）
    csRows.push({
      name: tds[0], edges: num(tds[1]), stencils: num(tds[2]), lNnz: num(tds[3]),
      factorNnz: num(tds[4]), total: num(tds[5]),
      ms: { projScatter: num(tds[8]), solve: num(tds[10]) },
      share: { projScatter: share(tds[8]), solve: share(tds[10]) },
    });
  }
}
console.log('\n约束集对照 行数:', csRows.length, '(期望 12 = 3 工况 × 4 套)');
let okCS = csRows.length === 12;
for (let g = 0; g < 3; ++g) {
  const grp = csRows.slice(g * 4, g * 4 + 4);
  const none = grp.find((r) => r.name.indexOf('默认') === 0);
  const bend = grp.find((r) => r.name === '只加弯曲');
  const both = grp.find((r) => r.name.indexOf('剪切 + 弯曲') === 0);
  if (!none || !bend || !both) { okCS = false; continue; }
  const dPS = bend.share.projScatter - none.share.projScatter;
  const dSolve = bend.share.solve - none.share.solve;
  // 结构故事要看**绝对耗时比**，不能看占比：总时间涨了，占比自然被稀释
  // （40×40 上投影+散射 0.301→0.308 只有 +2 %，但占比掉了 4.6 点就是因为分母变大）。
  const rPS = bend.ms.projScatter / none.ms.projScatter;
  const rSolve = bend.ms.solve / none.ms.solve;
  const fDrop = 100 * (both.factorNnz - bend.factorNnz) / bend.factorNnz;
  const lGrow = 100 * (both.lNnz - bend.lNnz) / bend.lNnz;
  // 判据用"占比变化"而不是"绝对倍数"：绝对倍数随规模变（投影+散射 1.02× → 1.61×，
  // 因为 100×100 上整条迭代都更受访存压力影响），而"回代占比显著上跳、投影+散射占比不涨"
  // 这三档都成立、也是这条结论真正要表达的。
  const pass = dSolve > 5 && dPS <= 0.5 && fDrop < 0 && lGrow > 0;
  if (!pass) okCS = false;
  console.log('  工况' + (g + 1) + '：只加弯曲 ⇒ 投影+散射 ' + rPS.toFixed(2) + '×、全局回代 '
              + rSolve.toFixed(2) + '×（占比 ' + (dPS >= 0 ? '+' : '') + dPS.toFixed(1) + '/+'
              + dSolve.toFixed(1) + ' 点）；两者一起时 L nnz +' + lGrow.toFixed(1) + '% 而因子 nnz '
              + fDrop.toFixed(1) + '%  ' + (pass ? '✓' : '← ⚠'));
}
const csNote = (store.cstable || '').includes('因子 nnz 反而更少');
console.log('  反直觉结论标注:', csNote ? '✓ 已写入表格脚注' : '← ⚠ 缺失');
const csAxis = 150 + (1180 - 150 - 250);   // renderConstraintSets 的轴右端
const csRects = [...(store.csbars || '').matchAll(/<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)"/g)]
  .map((m) => ({ x: +m[1], w: +m[3] }));
const csRight = csRects.length ? Math.max(...csRects.map((r) => r.x + r.w)) : 0;
console.log('  csbars      :', csRects.length, '条，最右', csRight.toFixed(1), '/ 轴右端', csAxis,
            csRight <= csAxis + 0.5 ? '  ✓' : '  ← ⚠ 溢出');
okCS = okCS && csNote && csRects.length >= 12 && csRight <= csAxis + 0.5;

// ---------------------------------------------------------------------------
// 采集时间（2026-09-24 新增）：每个小节标题下方应有一行"采集时间 …"。
// 存在的理由：不同小节的数据来自**不同会话**，只有同一小节内部才可比。
// 条数必须与小节数对上（§一 关键结论是混合来源、§十四 是数据来源表，这两节不插）。
// ---------------------------------------------------------------------------
const stamps = (html.match(/class="sub stamp"/g) || []).length;
const hasRule = html.includes('组内可比、组间不可比');
console.log('\n采集时间行    :', stamps, '(期望 12 = §二…§十三)',
            hasRule ? '  ✓ 含"组内可比、组间不可比"' : '  ← ⚠ 缺页首口径句');

const ok = rects.length > 0 && right <= labelW + axisW + 0.5 && reds === 20
  && empty.length === 0 && cfgCount === 14 && degrade === 6 && it.includes('机制失效档')
  && offStrain && s2Strain && Math.abs(s2Strain - offStrain) / offStrain < 0.02
  && sh && Number(sh[3]) > 20 && bdPS && Math.abs(Number(bdPS[3])) < 10 && okCS
  && stamps === 12 && hasRule;
console.log(ok ? '结论：几何自检通过' : '结论：**自检未通过**');
process.exit(ok ? 0 : 1);
