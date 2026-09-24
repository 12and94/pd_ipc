// _perf/make_report.js
// 从 build/_perf/logs 与 build/_baseline 的基准日志里解析分阶段数据，生成自包含的性能看板。
// 用法： node _perf/make_report.js      （产物写到 build/_perf/perf-report.html）
const fs = require('fs');
const path = require('path');

const ROOT = path.resolve(__dirname, '..');
const LOGDIR = path.join(ROOT, 'build', '_perf', 'logs');
const BASE = path.join(ROOT, 'build', '_baseline');

const STAGES = [
  { key: 'predict', label: '预测', color: '#8ecae6' },
  { key: 'projScatter', label: '投影+散射（融合）', color: '#fb8500' },
  { key: 'assemble', label: '组装+分解（一次性）', color: '#6b705c' },
  { key: 'solve', label: '全局回代（三分量并行）', color: '#d62828' },
  { key: 'unpackScan', label: '解包+判据（含残差）', color: '#9d4edd' },
  { key: 'velocity', label: '速度更新', color: '#588157' },
  { key: 'other', label: '其它', color: '#adb5bd' }
];
const LABEL2KEY = {
  '预测': 'predict', '局部步': 'localStep', '散射': 'scatter', '投影+散射': 'projScatter',
  '组装+分解': 'assemble', '全局回代': 'solve', '解包+判据': 'unpackScan', '速度更新': 'velocity',
  '其它': 'other'
};

const fmt = (v, d = 3) => (Math.round(v * 10 ** d) / 10 ** d).toFixed(d);   // Node 侧的控制台摘要用

function readText(p) {  if (!fs.existsSync(p)) return null;
  const buf = fs.readFileSync(p);
  // PowerShell 的 `*>` 重定向在这里写出的是 UTF-16LE（带 BOM），要按字节判断，
  // 否则中文标签会全部乱码、正则一条都匹配不上。
  if (buf.length >= 2 && buf[0] === 0xff && buf[1] === 0xfe) {
    return buf.toString('utf16le').replace(/^\uFEFF/, '');
  }
  return buf.toString('utf8').replace(/^\uFEFF/, '');
}

// 解析一次 pd_bench 运行
function parseBench(file) {
  const txt = readText(file);
  if (!txt) return null;
  const out = { file: path.basename(file), stages: {}, total: null, avgIter: null, threads: null, vertices: null, edges: null };
  const head = txt.match(/顶点\s+(\d+)\s+约束\s+(\d+).*线程\s+(\d+)/);
  if (head) { out.vertices = +head[1]; out.edges = +head[2]; out.threads = +head[3]; }
  const tot = txt.match(/子步数\s+(\d+).*单子步\s+([\d.]+)\s+ms\s+平均迭代\s+([\d.]+)/);
  if (tot) { out.steps = +tot[1]; out.total = +tot[2]; out.avgIter = +tot[3]; }
  for (const line of txt.split('\n')) {
    const m = line.match(/^\s*([^\s]+)\s+([\d.]+)\s+([\d.]+)%/);
    if (m && LABEL2KEY[m[1]]) out.stages[LABEL2KEY[m[1]]] = +m[2];
  }
  return out;
}

// 取多次重复里 total 最小的一次
function bestOf(files) {
  let best = null;
  for (const f of files) {
    const r = parseBench(f);
    if (r && r.total !== null && (!best || r.total < best.total)) best = r;
  }
  return best;
}

function reps(dir, prefix, n) {
  return Array.from({ length: n }, (_, i) => path.join(dir, `${prefix}_${i + 1}.txt`));
}

// ---------------------------------------------------------------- 线程扫描
const WORKLOADS = [
  { id: '40x40_it40', label: '40×40 · iters 40 · 残差 0.3', note: 'PD 迭代撞满上限的工况' },
  { id: '60x60_it2', label: '60×60 · iters 2（交互档）', note: '接近查看器稳态' },
  { id: '100x100_it10', label: '100×100 · iters 10', note: '1 万顶点' },
  { id: '200x200_it10', label: '200×200 · iters 10', note: '4 万顶点 / 7.96 万约束' }
];
const THREADS = [1, 2, 4, 8, 18];
const SWEEP_REPS = [1, 2, 3];   // 与 refresh_report.ps1 的 -SweepReps 默认值一致（逐轮交替，跨轮取最小）
const sweep = WORKLOADS.map(w => ({
  ...w,
  series: THREADS.map(t => bestOf(SWEEP_REPS.map(r => path.join(LOGDIR, `${w.id}_t${t}_r${r}.txt`))))
                  .filter(Boolean)
}));

// ---------------------------------------------------------------- 成对对照（同一会话）
function pair(label, aFiles, bFiles, aName, bName, threads) {
  const a = bestOf(aFiles), b = bestOf(bFiles);
  if (!a || !b) return null;
  return { label, threads, a: { name: aName, total: a.total, file: a.file, stages: a.stages },
           b: { name: bName, total: b.total, file: b.file, stages: b.stages } };
}
const history = [
  pair('40×40/iters40 · 18 线程', [path.join(BASE, 'p1_40x40_it40_t18_O.txt')], [path.join(BASE, 'p1_40x40_it40_t18_N.txt')],
       '改造前（Phase 0）', 'Phase 1：一个子步一个区域', 18),
  pair('40×40/iters40 · 4 线程', reps(BASE, 'c_40x40_it40', 3), reps(BASE, 'c_40x40_it40_p2', 3),
       'Phase 1', 'Phase 2：约束图着色', 4),
  pair('40×40/iters40 · 同二进制', reps(BASE, 'd_40x40_it40', 3), reps(BASE, 'd_40x40_it40_18', 3),
       '当前默认（4 线程）', '旧默认（18 线程）', '4 vs 18'),
  pair('60×60/iters2 交互档 · 同二进制', reps(BASE, 'd_60x60_it2', 3), reps(BASE, 'd_60x60_it2_18', 3),
       '当前默认（4 线程）', '旧默认（18 线程）', '4 vs 18'),
  pair('Phase 4b：全局步回代自己实现 · 40×40/iters40 · 4 线程',
       reps(BASE, 'sol_40x40_it40_eig', 3), reps(BASE, 'sol_40x40_it40_new', 3),
       'Eigen solve', '自己实现回代（scatter）', 4),
  pair('Phase 4b：全局步回代自己实现 · 200×200/iters10 · 4 线程',
       reps(BASE, 'sol_200x200_it10_eig', 3), reps(BASE, 'sol_200x200_it10_new', 3),
       'Eigen solve', '自己实现回代（scatter）', 4)
].filter(Boolean);

// ---------------------------------------------------------------- 查看器时间线
function parseViewer(file, name) {
  const txt = readText(file);
  if (!txt) return null;
  const pts = [];
  let sec = 0;
  for (const line of txt.split('\n')) {
    const m = line.match(/^fps\s+([\d.]+)\s+帧\s+([\d.]+)\s+ms\s+物理\s+([\d.]+)\s+ms\s+迭代\s+(\d+).*残差\s+([\d.eE+-]+)/);
    if (m) { pts.push({ t: sec++, fps: +m[1], frame: +m[2], physics: +m[3], iters: +m[4], residual: +m[5] }); }
  }
  return pts.length ? { name, pts } : null;
}
const viewer = [
  parseViewer(path.join(BASE, 'viewer_long_03.txt'), '改造前（Phase 0 代码、18 线程）'),
  parseViewer(path.join(BASE, 'viewer_4t.txt'), 'Phase 2 + 默认 4 线程'),
  parseViewer(path.join(BASE, 'viewer_p3.txt'), 'Phase 3（当前）')
].filter(Boolean);

// ---------------------------------------------------------------- Phase 3（残差并行化）
// 数据来自 _perf/ab_residual.ps1 的日志（固定列宽、由同一脚本生成）：
//   用例 基线 新版 变化% 基线迭代 新版迭代 基线散射 新版散射 基线解包 新版解包 基线回代 新版回代
// 列顺序变了就要同步改这里（脚本与报告是一对）。
function parsePhase3(file) {
  const txt = readText(file);
  if (!txt) return [];
  const rows = [];
  const re = /^(\S+)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)%\s+(\S+)\s+(\S+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$/;
  for (const line of txt.split('\n')) {
    const m = line.match(re);
    if (!m) continue;
    rows.push({ name: m[1], baseTotal: +m[2], newTotal: +m[3], delta: +m[4],
                baseIter: +m[5], newIter: +m[6],
                baseScatter: +m[7], newScatter: +m[8],
                baseUnpack: +m[9], newUnpack: +m[10],
                baseSolve: +m[11], newSolve: +m[12] });
  }
  return rows;
}
const phase3 = parsePhase3(path.join(BASE, 'phase3_ab.log'));
// Phase 4c（全局步按 x/y/z 三分量并行）的 A/B 日志**故意沿用同一列格式**，
// 所以上面这个解析器可以原样复用 —— 脚本 _perf/refresh_report.ps1 与这一段是一对，
// 列顺序变了要同步改两处。
const phase4c = parsePhase3(path.join(BASE, 'phase4c_ab.log'));
// Phase 4d（投影与散射融合）的 A/B：同样是上面那个列格式（"散射"那一列放的是
// "投影+散射"：基线 = 局部步+散射 之和，新版 = 融合后的桶）。
const phase4d = parsePhase3(path.join(BASE, 'phase4d_ab.log'));

// ---------------------------------------------------------------- 附加约束（剪切 / 弯曲）
// 数据来自 _perf/ab_shear.ps1 与 ab_bend.ps1 的"两值对照"日志格式：
//   取最小值对照：
//   项            关闭        开启        变化
//   单子步(ms)    2.967       3.481       +17.3%
//   ...
//   off : 约束 = 3120  L nnz = 23496  因子 nnz = 66027  末次残差 = ...  最终应变 max = ...  弹性能 = ...
// 这类日志是"开/关成对、交替 N 次取最小"，语义与 Phase 3/4c/4d 不同（这里**没有"谁更好"**：
// 弯曲/剪切是可选物理，代价本来就是正的），所以单独解析、单独渲染。
function parseTwoWay(file) {
  const txt = readText(file);
  if (!txt) return null;
  const out = { rows: [], off: {}, on: {}, title: '' };
  const lines = txt.split('\n');
  out.title = (lines[0] || '').replace(/^=+\s*/, '').replace(/\s*=+$/, '');
  let inTable = false;
  for (const line of lines) {
    if (line.includes('取最小值对照')) { inTable = true; continue; }
    if (inTable) {
      const m = line.match(/^\s*(.+?)\s{2,}([\d.]+)\s+([\d.]+)\s+([+-][\d.]+)%\s*$/);
      if (m) { out.rows.push({ name: m[1].trim(), off: +m[2], on: +m[3], pct: +m[4] }); }
      if (/^\s*$/.test(line) && out.rows.length) inTable = false;
    }
    const s = line.match(/^(off|on)\s*:\s*(.+)$/);
    if (s) {
      const t = s[2];
      const tag = s[1] === 'off' ? out.off : out.on;
      const pick = (re) => { const mm = t.match(re); return mm ? +mm[1] : null; };
      tag.constraints = pick(/约束\s*=\s*(\d+)/);
      tag.edgesAdded = pick(/新增对角边\s*=\s*(\d+)/);
      tag.stencils = pick(/stencil\s*=\s*(\d+)/);
      tag.lNnz = pick(/L nnz\s*=\s*(\d+)/);
      tag.factorNnz = pick(/因子 nnz\s*=\s*(\d+)/);
      tag.residual = t.match(/末次残差\s*=\s*([\d.eE+-]+)/)?.[1] ?? null;
      tag.strainMax = t.match(/最终应变 max\s*=\s*([\d.eE+-]+)/)?.[1] ?? null;
      tag.energy = t.match(/弹性能\s*=\s*([\d.eE+-]+)/)?.[1] ?? null;
    }
  }
  return out.rows.length ? out : null;
}
const shearAB = parseTwoWay(path.join(BASE, 'shear_ab.log'));
const bendAB = parseTwoWay(path.join(BASE, 'bend_ab.log'));

// ---------------------------------------------------------------- 弯曲 stencil 采样变体
// 数据来自 _perf/ab_bend_variants.ps1（多配置 × 同轮交替取最小）。两张表：
//   取最小值对照：配置 刚度 stencil 单子步 vs off 回代 投影散射 组装 解包 应变max 应变avg 弹性能
//   规模与峰值：  配置 L nnz 因子 nnz 全程峰值应变 末次残差
// 列之间用 2 个以上空格分隔（固定列宽），配置名里可能含**一个**空格（如 "std/s2 x2k"）。
function parseVariants(file) {
  const txt = readText(file);
  if (!txt) return null;
  const out = { configs: [], byName: {} };
  const lines = txt.split('\n');
  let mode = null;
  for (const line of lines) {
    if (line.startsWith('取最小值对照')) { mode = 'min'; continue; }
    if (line.startsWith('规模与峰值')) { mode = 'scale'; continue; }
    if (/^\s*$/.test(line) || line.startsWith('配置') || line.startsWith('各轮')) continue;
    const parts = line.trim().split(/\s{2,}/);
    if (!parts.length) continue;
    if (mode === 'min' && parts.length === 12) {
      const c = { name: parts[0], k: parts[1], stencils: +parts[2], perStep: +parts[3],
                  vsOff: parts[4], solve: +parts[5], projScatter: +parts[6], assemble: +parts[7],
                  unpack: +parts[8], strainMax: +parts[9], strainAvg: +parts[10], elastic: +parts[11] };
      out.configs.push(c);
      out.byName[c.name] = c;
    } else if (mode === 'scale' && parts.length === 5) {
      const c = out.byName[parts[0]];
      if (c) { c.lNnz = +parts[1]; c.factorNnz = +parts[2]; c.peakStrain = +parts[3]; c.residual = parts[4]; }
    }
  }
  return out.configs.length ? out : null;
}
const variants = parseVariants(path.join(BASE, 'bend_variants_ab.log'));

// ---------------------------------------------------------------- 惯性右端并行化（三方，结论：否）
// 数据来自 _perf/ab_inertial.ps1 的第三方对照（旧-single / 控制-single / 新-for）；
// 关键是"1 线程"那一档：单线程时 omp for 与串行等价，理论差值必须为 0 —— 它测出 1.7 %，
// 说明那点差别来自代码生成/布局，不是并行。详见 docs/perf.md §14。
function parseInertial(file) {
  const txt = readText(file);
  if (!txt) return null;
  const out = [];
  let cur = null;
  for (const line of txt.split('\n')) {
    const h = line.match(/^---\s*(\d+)\s*线程\s*---/);
    if (h) { cur = { threads: +h[1], rows: [] }; out.push(cur); continue; }
    if (!cur) continue;
    const r = line.match(/^\s*(\S+)\s+各轮\s+.*→\s*最小\s+([\d.]+)\s*ms/);
    if (r) { cur.rows.push({ name: r[1], ms: +r[2] }); continue; }
    const f = line.match(/新-for vs 旧-single：\s*([+-][\d.]+)%\s+新-for vs 控制-single：\s*([+-][\d.]+)%/);
    if (f) { cur.newVsOld = +f[1]; cur.newVsCtrl = +f[2]; continue; }
    const b = line.match(/旧-single vs 控制-single.*?：\s*([+-][\d.]+)%/);
    if (b) { cur.floor = +b[1]; }
  }
  return out.length ? out : null;
}
const inertial = parseInertial(path.join(BASE, 'inertial_ab3.log'));

// ---------------------------------------------------------------- 四套约束集的阶段构成
// 数据来自 _perf/constraint_sets.ps1：同一工况下跑 4 套约束集
// （默认 / 只剪切 / 只弯曲 / 两者都开），每套 2 次取单子步最小值。
// 日志是**原样的 pd_bench 输出**，所以这里直接复用 parseBench（连解析器都不用新写）。
// ⚠️ 四套的**物理不同**（约束集不同）：本节只比"时间与构成"，绝不比应变/能量。
const CS_DIR = path.join(ROOT, 'build', '_perf', 'logs_cs');
const CS_CFGS = [['none', '默认（距离 + pin）'], ['shear', '只加剪切'], ['bend', '只加弯曲'],
                 ['both', '剪切 + 弯曲']];
const CS_WORKLOADS = [
  { id: '40x40_it40', label: '40×40 · iters 40 · 300 子步 · 4 线程' },
  { id: '60x60_it2', label: '60×60 · iters 2（交互档）· 4 线程' },
  { id: '100x100_it10', label: '100×100 · iters 10 · 4 线程' }
];
function logMeta(file) {
  const txt = readText(file);
  if (!txt) return { stencils: 0, lNnz: 0, factorNnz: 0 };
  const pick = (re) => { const m = txt.match(re); return m ? +m[1] : 0; };
  return {
    stencils: pick(/弯曲约束：新增 stencil\s+(\d+)/),   // 只在开启弯曲时打印
    lNnz: pick(/nnz\s+(\d+)（结构）/),
    factorNnz: pick(/因子 nnz\s+(\d+)/)
  };
}
const constraintSets = CS_WORKLOADS.map((w) => ({
  id: w.id, label: w.label,
  rows: CS_CFGS.map(([cid, label]) => {
    const files = [1, 2, 3].map((r) => path.join(CS_DIR, `${w.id}__${cid}_r${r}.txt`));
    const best = bestOf(files);
    if (!best) return null;
    const metas = files.map(logMeta);
    const pickMax = (key) => metas.reduce((a, b) => Math.max(a, b[key]), 0);
    return { cfg: cid, label, stages: best.stages, steps: best.steps, total: best.total,
             avgIter: best.avgIter, edges: best.edges, stencils: pickMax('stencils'),
             lNnz: pickMax('lNnz'), factorNnz: pickMax('factorNnz') };
  }).filter(Boolean)
})).filter((w) => w.rows.length);

// ---------------------------------------------------------------- 采集时间
// 看板里的数字全部来自日志文件，而**不同小节的文件来自不同会话**（有的还是几天前采的）。
// 把每个小节的数据文件最后写入时间印在该节标题下方，配合页首那句"组内可比、组间不可比"，
// 免得有人把两个小节横着比（跨会话的绝对 ms 本来就不可比 —— docs/perf.md §0 的纪律）。
const fmtStamp = (d) => {
  const p = (n) => String(n).padStart(2, '0');
  return `${d.getFullYear()}-${p(d.getMonth() + 1)}-${p(d.getDate())} ${p(d.getHours())}:${p(d.getMinutes())}`;
};
function stamp(files) {
  const ts = files.map((f) => { try { return fs.statSync(f).mtime; } catch { return null; } })
                   .filter(Boolean).sort((a, b) => a - b);
  if (!ts.length) return '数据缺失';
  const a = fmtStamp(ts[0]), b = fmtStamp(ts[ts.length - 1]);
  return a === b ? a : `${a} … ${b}`;
}
const sweepFiles = WORKLOADS.flatMap((w) => THREADS.flatMap((t) => SWEEP_REPS.map((r) => path.join(LOGDIR, `${w.id}_t${t}_r${r}.txt`))));
const csFiles = fs.existsSync(CS_DIR)
  ? fs.readdirSync(CS_DIR).filter((f) => f.endsWith('.txt')).map((f) => path.join(CS_DIR, f)) : [];
const SECTION_SRC = {
  '二': sweepFiles, '三': sweepFiles,
  '四': history.flatMap((h) => [path.join(BASE, h.a.file), path.join(BASE, h.b.file)]),
  '五': [path.join(BASE, 'phase3_ab.log')],
  '六': [path.join(BASE, 'phase4c_ab.log')],
  '七': [path.join(BASE, 'phase4d_ab.log')],
  '八': [path.join(BASE, 'shear_ab.log')],
  '九': [path.join(BASE, 'bend_ab.log')],
  '十': csFiles,
  '十一': [path.join(BASE, 'bend_variants_ab.log')],
  '十二': [path.join(BASE, 'inertial_ab3.log')],
  '十三': [path.join(BASE, 'viewer_long_03.txt'), path.join(BASE, 'viewer_4t.txt'), path.join(BASE, 'viewer_p3.txt')]
};
// 单点实现：在每个 <h2>（形如 "三、线程数扫描"）后面插一行采集时间。
// §一（关键结论，混合来源）与 §十四（数据来源表）不插。
function withStamps(h) {
  return h.replace(/(<h2>([一二三四五六七八九十]+)、[^<]*<\/h2>)/g, (m, whole, num) => {
    const files = SECTION_SRC[num];
    if (!files) return whole;
    return whole + `<div class="sub stamp">采集时间 ${stamp(files)}（该节数据文件的最后写入时间）</div>`;
  });
}

// ---------------------------------------------------------------- HTML
// KPI 卡片里的 Phase 4c 数字**从数据里算**（不要手写）：同一轮会话的另一次测量与
// docs/perf.md §9 的 −42.3 % 会差零点几个百分点（噪声），手写就会和下面的表自相矛盾。
const p4kpi = (phase4c.find(c => c.name === '40x40_it40_t4') || {}).delta;
const p4kpiText = (p4kpi === undefined) ? '—' : ('−' + Math.abs(p4kpi).toFixed(1) + '%');

const DATA = { sweep, history, viewer, phase3, phase4c, phase4d, shearAB, bendAB, variants, inertial,
               constraintSets, stages: STAGES,
               generated: new Date().toISOString().slice(0, 16).replace('T', ' ') };

const rawHtml = `<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>pd_ipc 性能看板 · 分函数层级</title>
<style>
  :root { color-scheme: dark; }
  body { margin:0; padding:24px 28px 60px; background:#12141a; color:#e6e8ee;
         font:14px/1.6 "Segoe UI","Microsoft YaHei",system-ui,sans-serif; }
  h1 { font-size:22px; margin:0 0 4px; }
  h2 { font-size:17px; margin:34px 0 6px; padding-bottom:6px; border-bottom:1px solid #2a2f3a; }
  .sub { color:#98a2b3; font-size:13px; margin-bottom:18px; }
  .stamp { margin:2px 0 10px; font-size:12.5px; color:#7d879a; }
  .card { background:#171a21; border:1px solid #232833; border-radius:10px; padding:16px 18px; margin:14px 0 22px; }
  .row { display:flex; gap:22px; flex-wrap:wrap; align-items:flex-start; }
  .legend { display:flex; gap:14px; flex-wrap:wrap; font-size:12.5px; color:#b7bfcc; margin:4px 0 12px; }
  .legend i { display:inline-block; width:11px; height:11px; border-radius:2px; margin-right:5px; vertical-align:-1px; }
  table { border-collapse:collapse; font-size:13px; }
  th,td { border-bottom:1px solid #232833; padding:5px 12px 5px 0; text-align:right; white-space:nowrap; }
  th:first-child, td:first-child { text-align:left; }
  th { color:#98a2b3; font-weight:600; }
  .tot { font-variant-numeric:tabular-nums; }
  .best { color:#7ee787; font-weight:700; }
  .bad { color:#ff9f9f; }
  .warn { background:#2a2118; border-left:3px solid #fb8500; padding:10px 14px; border-radius:6px; font-size:13px; color:#e8d3b8; }
  .kpi { display:flex; gap:28px; flex-wrap:wrap; }
  .kpi div { min-width:150px; }
  .kpi b { display:block; font-size:24px; color:#7ee787; font-variant-numeric:tabular-nums; }
  .kpi span { color:#98a2b3; font-size:12.5px; }
  .toggle { background:#222735; border:1px solid #333a4a; color:#dfe3ea; border-radius:6px; padding:6px 12px; cursor:pointer; font-size:13px; }
  .toggle.on { background:#2f6feb; border-color:#2f6feb; color:#fff; }
  svg text { fill:#c7ced9; font-size:11px; }
  code { background:#20242e; padding:1px 5px; border-radius:4px; font-size:12.5px; }
</style></head><body>
<h1>pd_ipc 性能看板 · 分函数层级</h1>
<div class="sub">数据来自 <code>build/_perf/logs</code> 与 <code>build/_baseline</code> 的基准日志；生成于 ${DATA.generated}。
口径：固定工况、同一轮会话内交替重复 <b>2–5 次取最小值</b>、MSVC Release、<code>/openmp</code>（vcomp）。
hover 图形可看精确数值。<br>
⚠️ <b>组内可比、组间不可比</b>：各小节的数据来自<b>不同会话</b>（每节标题下方印了该节数据文件的采集时间），
所以只有同一小节内部的成对/成组数字才能横着比 —— 跨节比绝对 ms 没有意义（docs/perf.md §0 的纪律）。</div>

<h2>一、关键结论</h2>
<div class="card kpi">
  <div><b>${p4kpiText}</b><span>Phase 4c：全局步按 x/y/z<br>三分量并行 · 40×40/iters40 · 4 线程</span></div>
  <div><b>−26%</b><span>40×40/iters40：默认线程数<br>由 18 改成 4 的收益（Phase 4a）</span></div>
  <div><b>−13.6%</b><span>Phase 3 残差并行化<br>40×40/iters40 · 18 线程</span></div>
  <div><b>−18.2%</b><span>Phase 4b 全局步回代自己实现<br>200×200/iters10 · 4 线程</span></div>
  <div><b>0.0%</b><span>Phase 3/4b/4c 的物理输出变化<br>（求和顺序未变 ⇒ 逐位相同）</span></div>
  <div><b>0.55 ms</b><span>60×60 交互档稳态物理耗时<br>（Phase 4c 后重测；改造前 1.9 ms）</span></div>
  <div><b>16/16</b><span>散射跨线程 + 重复运行<br>位级相同的次数</span></div>
</div>

<h2>二、分函数层级：每个阶段各花多少</h2>
<div class="card">
  <div class="row">
    <button class="toggle on" id="modeAbs" onclick="setMode('abs')">绝对值（ms）</button>
    <button class="toggle" id="modePct" onclick="setMode('pct')">占比（%）</button>
  </div>
  <div class="sub" style="margin:0 0 10px">
    一行 = 一个线程数；条形 = <b>单子步的墙钟时间</b>（阶段计时桶的累计值 ÷ 子步数），按阶段堆叠；
    右端是"单子步总计"和"占比最大的那个阶段"。<b>每个工况各自归一化</b>（标尺写在每组右上角），
    所以长度只在同一组内可比，不要跨工况比较。悬停任一段可看精确数值。<br>
    注意 <span style="color:#6b705c">■ 组装+分解</span> 是<b>一次性</b>成本（只在 stamp 变化时发生）：
    这里按"累计 ÷ 子步数"折算，因此子步数少的工况（例如 200×200 只有 8 个子步）里它会明显偏大。
  </div>
  <div class="legend" id="legend"></div>
  <div id="stack"></div>
</div>

<h2>三、线程数扫描</h2>
<div class="card"><div class="sub" style="margin:0 0 10px">
  每条曲线一个工况，绿点 = 该工况<b>本轮日志里</b>最小的线程档。<br>
  ✅ 本轮扫描是<b>逐轮交替</b>（轮次在最外层，每轮把 5 个线程档背靠背跑完，跨轮取最小值）
  ⇒ 会话内的机器漂移对各档是公平的，"绿点"第一次有可信的口径。<br>
  ⚠️ 但它仍然只是"这一轮日志里最小"：各档之间差几个百分点就是噪声；而且<b>组内可比、组间不可比</b>
  （别的节来自别的会话）。默认线程数是<b>机器相关参数</b>（<code>min(hw−2, 4)</code>，本机 = 4），
  换机器或换量级必须重标定 —— 见 <code>docs/open-issues.md</code> §0 与 <code>HANDOFF.md</code> §4。
</div><div id="lines"></div></div>

<h2>四、优化历程（同一会话内的成对对照）</h2>
<div class="card"><div id="hist"></div></div>

<h2>五、Phase 3：残差并行化（同一会话交替 A/B，4 次取最小值）</h2>
<div class="card">
  <div class="sub" style="margin:0 0 10px">
    基线 = 着色散射 + <b>串行残差</b>；新版 = 着色散射 + <b>逐顶点 gather 残差</b>（与解包趟融合）。
    上图：单子步总耗时；下图：<code>解包+判据</code> 阶段（残差所在的那一桶，含串行趟与 barrier 等待）。
    1 线程那一例走的是<b>同一条串行路径</b>，所以它应当零变化 —— 这正是"没有拖慢 1 线程参照"的证据。</div>
  <div id="p3total"></div>
  <div id="p3unpack"></div>
  <div id="p3table"></div>
</div>

<h2>六、Phase 4c：全局步按 x/y/z 三分量并行（同一会话交替 A/B，5 次取最小值）</h2>
<div class="card">
  <div class="sub" style="margin:0 0 10px">
    基线 = 整趟串行回代（改动前二进制）；新版 = 按 x/y/z 三分量拆分。<br>
    为什么能拆：<code>L</code> 的每一块都是标量 × I₃（惯性 m/h²·I、距离约束 ±κ·I、pin 行 I）
    ⇒ <code>L = Ã ⊗ I₃</code> ⇒ 它的图是 <b>3 个互不相连的连通分量</b>，消元填充不可能跨分量
    ⇒ 3n 个方程本来就是 <b>3 个互不相干的标量系统</b>（免费拥有 3 个独立右端项）。回代是<b>延迟受限</b>的
    （达成带宽只有流式读写的 9–13 %），3 条独立链分给 3 条线程正好把"未决访存请求数"提高约 3 倍；
    <b>零同步、零竞争</b>。<br>
    上图：单子步总耗时；下图：<code>全局回代</code> 阶段。两侧<b>物理输出逐位相同</b>
    （迭代数 / 末次残差 / 应变 / 能量逐字一致）；1 线程那一例没有并行可拿，只有局部性收益。<br>
    ⚠️ <code>全局回代</code> 桶的口径在本轮变了：原来是 <code>omp single</code> 里"唯一执行者自己的
    执行时间"（不含等待），现在是按分量 <code>omp for</code> 之后的<b>阶段墙钟</b>（含等最慢者）
    —— 与散射 / 解包桶口径一致，但跟历史上"回代占 76 %"不能直接比较。
  </div>
  <div id="p4total"></div>
  <div id="p4solve"></div>
  <div id="p4table"></div>
</div>

<h2>七、Phase 4d：投影与散射融合（同一会话交替 A/B，5 次取最小值）</h2>
<div class="card">
  <div class="sub" style="margin:0 0 10px">
    基线 = 两趟（<code>projectInRegion</code> 物化 <code>targets</code> → <code>scatterIntoInRegion</code> 读它累加）；
    新版 = 合成一趟，按颜色就地算 <code>d_c</code> 直接累加，<code>targets</code> 不再物化。<br>
    为什么是"零数值风险"：两条路**共用同一份算术 helper**（<code>projectEdgeValue</code> /
    <code>scatterEdgeValue</code>）⇒ 每个操作与操作数顺序逐字相同 ⇒ 结果<b>逐位相同</b>
    （已用三档工况逐字比对；<code>test_primitives</code> 也有"融合路径 == 两趟路径逐位"的断言）。<br>
    上图：单子步总耗时；下图：<code>投影+散射</code>（基线是『局部步+散射』之和）。<br>
    ⚠️ <b>桶口径</b>：本轮起 <code>局部步</code> 与 <code>散射</code> 合并为 <code>投影+散射</code>，
    引用历史占比时要相加再比。
  </div>
  <div id="p5total"></div>
  <div id="p5scatter"></div>
  <div id="p5table"></div>
</div>

<h2>八、可选物理：剪切约束（同一会话交替 A/B，5 次取最小值）</h2>
<div class="card">
  <div class="sub" style="margin:0 0 10px">
    剪切 = 给规则网格每个四边形加一条<b>对角边</b>（复用距离约束的整套实现：装配 / 着色 / 投影 /
    散射 / 三分量并行回代<b>一个字都没改</b>）。所以它是"多一类约束"的<b>纯代价</b>：
    约束数 3120 → 4641，<code>投影+散射</code>那一桶明显变大。<br>
    开关：<code>SceneConfig::shearStiffness</code> / <code>PD_SHEAR=&lt;刚度&gt;</code>，<b>默认 0 = 关闭</b>
    （关闭时网格/边表/装配矩阵/残差口径与改动前逐字相同）。<br>
    ⚠️ 这里<b>没有"谁更好"</b>：开启是可选物理，代价本来就是正的（红色 = 关闭、绿色 = 开启，只表示两侧的值）。
  </div>
  <div id="shbars"></div>
  <div id="shtable"></div>
</div>

<h2>九、可选物理：线性（中点）弯曲约束（同一会话交替 A/B，5 次取最小值）</h2>
<div class="card">
  <div class="sub" style="margin:0 0 10px">
    弯曲 = 沿每行/每列取连续三个顶点 <code>(a,b,c)</code>，<code>A_c x = x_a − 2x_b + x_c</code>、
    <code>E = (k/2)‖A_c x‖²</code>（等价于"b 与两邻居中点之间一根静止长度 0 的弹簧"）。<br>
    选这个形式（而不是二面角）的理由是四条结构性收益：<code>A_c</code> 常数 ⇒ <b>只分解一次</b>；
    每块是 <code>k·w_p·w_q·I₃</code> ⇒ <b>⊗ 结构与三分量并行保住</b>；流形线性 ⇒ 投影 <code>p_c ≡ 0</code>
    ⇒ <b>局部步与散射一个字都不用改</b>；右端只多一个与位形无关的常向量。<br>
    代价几乎全在<b>全局回代</b>：因子 nnz 66,027 → 182,613（2.77×），而 <code>L</code> 的 nnz 只 +78 %
    ⇒ 涨的是<b>消元填充</b>，不是矩阵本身。<br>
    <b>默认 0 = 关闭</b>（<code>bendStiffness</code> / <code>PD_BEND=&lt;k&gt;</code>）。
  </div>
  <div id="bdbars"></div>
  <div id="bdtable"></div>
</div>

<h2>十、四套约束集的阶段构成：只剪切 / 只弯曲 / 两者都开（同一工况、4 线程）</h2>
<div class="card">
  <div class="sub" style="margin:0 0 10px">
    第二节画的是<b>默认约束集</b>（距离 + pin）。这里把另外三套也采一遍，回答"多出来的时间挪到哪一块了"：<br>
    · <b>只加剪切</b>：时间主要进 <code>投影+散射</code>（剪切走的就是这条完整通路：每条对角边都要投影+散射）；<br>
    · <b>只加弯曲</b>：代价几乎全进 <code>全局回代</code>（弯曲没有局部步/散射）——
      <b>但要按规模分开读</b>：40×40/60×60 上回代涨 2.0×/2.4×，而 <code>投影+散射</code> 几乎不动（+2 %/+7 %）；
      到了 100×100 回代涨 <b>5.9×</b>、连 <code>投影+散射</code> 也跟着涨 61 %（整条迭代的访存压力都变大了）；<br>
    · <b>两者都开</b>：不是简单相加 —— 100×100 上它比"只加弯曲"<b>还快</b>（19.9 vs 26.0 ms/子步），
      原因见下表脚注（<b>填充不是约束数的单调函数</b>）。<br>
    ⚠️ 四套的<b>物理不同</b>（约束集不同）⇒ 这里只比<b>时间与构成</b>，绝不比应变/能量；
    每套跑 3 次取单子步最小值（同一轮会话）。日志：<code>build/_perf/logs_cs/</code>。
  </div>
  <div id="csbars"></div>
  <div id="cstable"></div>
</div>

<h2>十一、弯曲 stencil 的采样变体：一个有效、一个退化（同一会话交替 A/B，5 次取最小值）</h2>
<div class="card">
  <div class="sub" style="margin:0 0 10px">
    把"取哪些 stencil"变成开关（<code>PD_BEND_SAMPLING=standard|rows|cols|checker</code>、
    <code>PD_BEND_STRIDE=1|2|3</code>），因为弯曲的代价主要在<b>消元填充</b>上，少取 stencil 就直接省填充。<br>
    <b>结论</b>：✅ <b>单向（rows/cols）有效</b> —— 因子 nnz −49 %、单子步只 +15.9 %（对照全采样 +75 %），
    但材料变成<b>各向异性</b>；❌ <b>隔行采样（stride ≥ 2）是退化的</b> —— 未被采样的中心顶点可以自由吸收
    全部曲率（把奇数顶点取成两侧偶数顶点的中点即精确满足全部被采样约束）⇒ <b>付了代价却几乎无效果</b>
    （k 加 100 倍，应变只变 0.4 %；对照完整采样 −75.5 %）；
    ⚠️ 棋盘同样规模下比单向<b>更贵也更弱</b>，不推荐。<br>
    下表"单子步"是每子步 ms，"因子 nnz / L nnz"来自同一轮运行；"应变 max"越小说明越挺。
  </div>
  <div id="varbars"></div>
  <div id="vartable"></div>
</div>

<h2>十二、被否掉的优化：惯性右端并行化（三方 A/B，7 次取最小值）</h2>
<div class="card">
  <div class="sub" style="margin:0 0 10px">
    一个子步里惯性右端 <code>bBase = (M/h²)x̂</code> 原本是 <code>omp single</code>
    （1 条线程干、其余在 barrier 空转）。改成 <code>omp for</code>（按顶点切给所有线程）能赚吗？<b>不能</b>。<br>
    <b>决定性证据是"1 线程"那一档</b>：单线程时 <code>omp for</code> 与串行循环是同一件事，理论差值必须为 0，
    可它已经"快"了 1.7 % ⇒ 那点差别来自<b>代码生成 / 布局</b>（顺带还去掉了每子步一次 <code>b.resize()</code>），
    不是并行。去掉这个地板后 4T/18T 只剩 −1.2 % / −0.6 %，落在"同源码两次构建"的地板（0.2–1.4 %）里。<br>
    量级本来也不够：这一步每子步只跑一次、O(n)=1600 顶点的纯读写（串行 3–5 µs），4T 下子步 2630 µs
    ⇒ 即使并行无限快，上限也只有 <b>0.1 %</b> 量级。代码<b>已回退</b>。
  </div>
  <div id="inerttable"></div>
  <div class="sub" style="margin:8px 0 0">
    <b>方法论收获</b>：凡"把某一步并行化 / 换写法"的对照，必须留一个<b>让该机制失效的档位</b>（这里 = 1 线程）
    作为对照；否则 1–2 % 的代码布局抖动会被当成并行收益写进文档。
  </div>
</div>

<h2>十三、交互窗口实测（60×60，每秒钟一个采样点）</h2>
<div class="card">
  <div class="sub" style="margin:0 0 10px">前 ~10 秒是瞬态（迭代用尽、约 20 FPS），之后跨过残差门槛进入稳态。</div>
  <div id="viewer"></div>
</div>

<h2>十四、精确数字与数据来源</h2>
<div class="card" id="tables"></div>

<div class="card warn">
  <b>测量纪律（同样重要）</b><br>
  ① 只比较"两侧物理完全相同"的运行（迭代数 / 末次残差 / 最终应变一致）；<br>
  ② 旧新二进制必须<b>同一轮会话内交替、重复多次取最小值</b> —— 单次运行的波动可达 19 %，
     本轮就出现过"一次测得慢 35 %、重复取最小值后是快 9 %"的情况；<br>
  ③ 换过 OpenMP 开关后必须核对<b>二进制实际链接的运行时</b>：本项目踩过 CMake 缓存残留
     <code>/openmp:llvm</code> 的坑，导致几轮数据失真（详见 docs/perf.md §6.2）；<br>
  ④ 分阶段耗时只作定性参考，占比波动约 ±5 %。
</div>

<script>
const DATA = ${JSON.stringify(DATA)};
const STAGES = DATA.stages;
let MODE = 'abs';

function el(tag, attrs, inner) {
  const a = Object.entries(attrs || {}).map(([k, v]) => \` \${k}="\${v}"\`).join('');
  return inner === undefined ? \`<\${tag}\${a}/>\` : \`<\${tag}\${a}>\${inner}</\${tag}>\`;
}
function fmt(v, d = 3) { return (Math.round(v * 10 ** d) / 10 ** d).toFixed(d); }

function renderLegend() {
  document.getElementById('legend').innerHTML = STAGES
    .map(s => \`<span><i style="background:\${s.color}"></i>\${s.label}</span>\`).join('');
}

function renderStack() {
  const W = 1180, rowH = 30, labelW = 132, padR = 210, axisW = W - labelW - padR;
  let y = 14, out = '';
  for (const w of DATA.sweep) {
    // 阶段计时桶是**累计值**，而 total 是**每子步**：必须先把阶段折算成 ms/子步，
    // 否则两者混在一个刻度里 —— 条形会全部溢出（这张图早期就是这么错的：
    // 屏幕上只剩前两个阶段，看起来像"预测+局部步 占满"，而回代明明是 70%）。
    const rows = w.series.map(function (s) {
      const st = {};
      let sum = 0;
      for (const sg of STAGES) {
        const v = (s.stages[sg.key] || 0) / (s.steps || 1);
        st[sg.key] = v;
        sum += v;
      }
      return { threads: s.threads, total: s.total, st: st, sum: sum };
    });
    const totals = rows.map(function (r) { return r.total; });
    const bestT = Math.min.apply(null, totals);
    const maxAbs = Math.max.apply(null, totals);

    out += el('text', { x: 0, y: y + 12, style: 'fill:#e6e8ee;font-weight:600' }, w.label);
    out += el('text', { x: labelW + axisW + 8, y: y + 12, style: 'fill:#6b7280' },
              '标尺：最长条 = ' + fmt(maxAbs, maxAbs < 10 ? 2 : 1) + ' ms/子步（本组独立归一化）');
    y += 18;

    const axisBottom = y + rows.length * rowH;
    if (MODE === 'abs') {
      for (let k = 0; k <= 4; ++k) {
        const v = maxAbs * k / 4;
        const x = labelW + v * (axisW / maxAbs);
        out += el('line', { x1: fmt(x, 1), y1: y, x2: fmt(x, 1), y2: axisBottom + 4, stroke: '#232833' });
        out += el('text', { x: fmt(x, 1), y: axisBottom + 16, 'text-anchor': 'middle', style: 'fill:#6b7280' },
                  fmt(v, v < 10 ? 2 : 1));
      }
      out += el('text', { x: fmt(labelW + axisW + 8, 1), y: axisBottom + 16, style: 'fill:#6b7280' }, 'ms/子步');
    }

    for (const r of rows) {
      const isBest = Math.abs(r.total - bestT) < 1e-9;
      out += el('text', { x: 0, y: y + rowH / 2 + 4, style: 'fill:' + (isBest ? '#7ee787' : '#c7ced9') },
                r.threads + ' 线程');
      let x = labelW;
      let bigKey = null;
      let bigV = -1;
      for (const sg of STAGES) {
        const v = r.st[sg.key] || 0;
        if (v > bigV) { bigV = v; bigKey = sg.key; }
        const wpx = MODE === 'abs' ? v * (axisW / maxAbs) : 100 * v / r.sum * (axisW / 100);
        if (wpx < 0.4) continue;
        const tip = sg.label + '：' + fmt(v) + ' ms/子步（' + fmt(100 * v / r.sum, 1) + '%）' +
                    (sg.key === 'assemble' ? '　※一次性成本，按累计÷子步数折算' : '');
        out += el('rect', { x: fmt(x, 2), y: y + 3, width: fmt(Math.max(wpx, 0.4), 2), height: rowH - 8,
                            fill: sg.color, rx: 1.5 }, el('title', {}, tip));
        if (wpx > 46) {
          out += el('text', { x: fmt(x + wpx / 2, 1), y: y + rowH / 2 + 4, 'text-anchor': 'middle',
                              style: 'fill:#0b0d12;font-weight:600' }, fmt(v, v < 10 ? 2 : 1));
        }
        x += wpx;
      }
      let bigName = '';
      for (const sg of STAGES) { if (sg.key === bigKey) bigName = sg.label; }
      const rx = labelW + axisW + 8;
      out += el('text', { x: fmt(rx, 1), y: y + rowH / 2 - 1, style: 'fill:' + (isBest ? '#7ee787' : '#e6e8ee') + ';font-weight:600' },
                fmt(r.total, r.total < 10 ? 3 : 1) + ' ms');
      out += el('text', { x: fmt(rx, 1), y: y + rowH / 2 + 13, style: 'fill:#6b7280' },
                '最大：' + bigName + ' ' + fmt(bigV > 0 ? 100 * bigV / r.sum : 0, 0) + '%');
      y += rowH;
    }
    y += 26;
  }
  document.getElementById('stack').innerHTML =
    el('svg', { viewBox: '0 0 ' + W + ' ' + y, width: '100%', height: y }, out);
}
function setMode(m) {
  MODE = m;
  document.getElementById('modeAbs').className = 'toggle' + (m === 'abs' ? ' on' : '');
  document.getElementById('modePct').className = 'toggle' + (m === 'pct' ? ' on' : '');
  renderStack();
}

function renderLines() {
  const W = 1180, H = 260, pad = 52, gapX = 18;
  const n = DATA.sweep.length;
  const cw = (W - pad * 2 - gapX * (n - 1)) / n;
  let out = '';
  DATA.sweep.forEach((w, i) => {
    const xs = pad + i * (cw + gapX);
    const vals = w.series.map(s => s.total);
    const lo = Math.min(...vals) * 0.95, hi = Math.max(...vals) * 1.05;
    const X = t => xs + (Math.log2(t) - Math.log2(1)) / (Math.log2(18) - Math.log2(1)) * cw;
    const Y = v => 30 + (hi - v) / (hi - lo) * (H - 70);
    const best = Math.min(...vals);
    let svg = el('text', { x: xs, y: 16, style: 'fill:#e6e8ee;font-weight:600' }, w.label);
    svg += el('line', { x1: xs, y1: H - 34, x2: xs + cw, y2: H - 34, stroke: '#2a2f3a' });
    const pts = w.series.map(s => [X(s.threads), Y(s.total)]);
    svg += el('polyline', { points: pts.map(p => p.map(v => fmt(v, 1)).join(',')).join(' '),
                            fill: 'none', stroke: '#2f6feb', 'stroke-width': 2 });
    w.series.forEach((s, k) => {
      const isBest = Math.abs(s.total - best) < 1e-9;
      svg += el('circle', { cx: fmt(pts[k][0], 1), cy: fmt(pts[k][1], 1), r: isBest ? 6 : 3.5,
                           fill: isBest ? '#7ee787' : '#8ab4f8' },
                el('title', {}, \`\${s.threads} 线程: \${fmt(s.total)} ms/子步（平均迭代 \${fmt(s.avgIter, 2)}）\`));
      svg += el('text', { x: fmt(pts[k][0], 1), y: fmt(pts[k][1] - 10, 1),
                          'text-anchor': 'middle', style: 'fill:' + (isBest ? '#7ee787' : '#98a2b3') },
                fmt(s.total, s.total < 10 ? 2 : 0));
    });
    svg += el('text', { x: xs, y: H - 14, style: 'fill:#98a2b3' }, '1→2→4→8→18 线程（对数刻度）');
    out += el('g', {}, svg);
  });
  document.getElementById('lines').innerHTML =
    el('svg', { viewBox: \`0 0 \${W} \${H + 20}\`, width: '100%', height: H + 20 }, out);
}

function renderHist() {
  let out = '<table><tr><th>对照（同一会话）</th><th>线程</th><th>改造/旧</th><th>改造/新</th><th>改善</th></tr>';
  for (const h of DATA.history) {
    const imp = 1 - h.b.total / h.a.total;
    out += \`<tr><td>\${h.label}</td><td>\${h.threads}</td>
      <td>\${h.a.name}: \${fmt(h.a.total)} ms</td>
      <td>\${h.b.name}: \${fmt(h.b.total)} ms</td>
      <td class="\${imp > 0 ? 'best' : 'bad'}">\${(imp * 100).toFixed(1)}%</td></tr>\`;
  }
  out += '</table>';
  document.getElementById('hist').innerHTML = out;
}

function renderViewer() {
  const W = 1180, H = 240, pad = 56;
  const all = DATA.viewer.flatMap(s => s.pts);
  const maxT = Math.max(...all.map(p => p.t)), maxY = Math.max(...all.map(p => p.physics)) * 1.05;
  const colors = ['#ff9f9f', '#7ee787'];
  let out = '';
  out += el('line', { x1: pad, y1: H - 30, x2: W - 20, y2: H - 30, stroke: '#2a2f3a' });
  out += el('line', { x1: pad, y1: 20, x2: pad, y2: H - 30, stroke: '#2a2f3a' });
  out += el('text', { x: 6, y: 26, style: 'fill:#98a2b3' }, '物理 ms');
  out += el('text', { x: W - 150, y: H - 10, style: 'fill:#98a2b3' }, '秒');
  DATA.viewer.forEach((s, i) => {
    const pts = s.pts.map(p => [pad + p.t / maxT * (W - pad - 20), 20 + (1 - p.physics / maxY) * (H - 50)]);
    out += el('polyline', { points: pts.map(p => p.map(v => fmt(v, 1)).join(',')).join(' '),
                            fill: 'none', stroke: colors[i], 'stroke-width': 1.8 });
    out += el('text', { x: W - 420, y: 24 + i * 16, style: 'fill:' + colors[i] }, s.name);
  });
  document.getElementById('viewer').innerHTML = el('svg', { viewBox: \`0 0 \${W} \${H}\`, width: '100%', height: H }, out);
}

function renderTables() {
  let out = '';
  for (const w of DATA.sweep) {
    out += \`<h3 style="margin:18px 0 6px;font-size:14px">\${w.label} <span style="color:#98a2b3;font-weight:400">（\${w.note}；单位 ms，括号为占比）</span></h3>\`;
    out += '<table><tr><th>阶段</th>' + w.series.map(s => \`<th>\${s.threads} 线程</th>\`).join('') + '</tr>';
    for (const st of STAGES) {
      out += \`<tr><td>\${st.label}</td>\` + w.series.map(s => {
        const v = s.stages[st.key] || 0;
        return \`<td>\${fmt(v)} <span style="color:#6b7280">(\${fmt(100 * v / s.total, 1)}%)</span></td>\`;
      }).join('') + '</tr>';
    }
    out += '<tr><td><b>合计</b></td>' + w.series.map(s => {
      const isBest = Math.abs(s.total - Math.min(...w.series.map(z => z.total))) < 1e-9;
      return \`<td class="\${isBest ? 'best' : ''}"><b>\${fmt(s.total)}</b></td>\`;
    }).join('') + '</tr>';
    out += \`<tr><td>平均迭代</td>\` + w.series.map(s => \`<td>\${fmt(s.avgIter, 2)}</td>\`).join('') + '</tr>';
    out += \`<tr><td>日志文件</td>\` + w.series.map(s => \`<td style="font-size:11px;color:#6b7280">\${s.file}</td>\`).join('') + '</tr>';
    out += '</table>';
  }
  document.getElementById('tables').innerHTML = out;
}

// 画一组"两两成对"的横条：每例两种颜色（基线 / 新版），长度按给定字段的比例
function pairs(cases, targetId, title, fieldBase, fieldNew, unit) {
  const W = 1180, rowH = 34, gap = 12, labelW = 200, padR = 150;
  let y = 6, out = '';
  out += el('text', { x: 0, y: y + 12, style: 'fill:#e6e8ee;font-weight:600' }, title);
  y += 24;
  let max = 0;
  for (const c of cases) max = Math.max(max, c[fieldBase], c[fieldNew]);
  for (const c of cases) {
    const bv = c[fieldBase], nv = c[fieldNew];
    const bw = Math.max((W - labelW - padR) * bv / max, 0.5);
    const nw = Math.max((W - labelW - padR) * nv / max, 0.5);
    out += el('text', { x: 0, y: y + rowH / 2 + 5, style: 'fill:#c7ced9' }, c.name);
    out += el('rect', { x: labelW, y: y + 1, width: fmt(bw, 2), height: rowH / 2 - 3, fill: '#ff9f9f', rx: 2 },
              el('title', {}, '基线: ' + fmt(bv, 1) + ' ' + unit));
    out += el('rect', { x: labelW, y: y + rowH / 2, width: fmt(nw, 2), height: rowH / 2 - 3, fill: '#7ee787', rx: 2 },
              el('title', {}, '新版: ' + fmt(nv, 1) + ' ' + unit));
    const pct = 100 * (nv - bv) / bv;
    const tip = pct < 0 ? fmt(pct, 1) + '%' : '+' + fmt(pct, 1) + '%';
    out += el('text', { x: fmt(labelW + Math.max(bw, nw) + 8, 1), y: y + rowH / 2 + 5,
                        style: 'fill:' + (pct < 0 ? '#7ee787' : '#ff9f9f') },
              fmt(bv, 1) + ' → ' + fmt(nv, 1) + '  (' + tip + ')');
    y += rowH + gap;
  }
  document.getElementById(targetId).innerHTML =
    el('svg', { viewBox: '0 0 ' + W + ' ' + y, width: '100%', height: y }, out);
}

// 通用的"同一会话 A/B"渲染：总耗时 + 第二个关注的桶 + 明细表（Phase 3 与 Phase 4c 共用）
function renderAB(cases, idTotal, idSecond, idTable, opt) {
  if (!cases.length) {
    document.getElementById(idTotal).innerHTML =
      '<i style="color:#98a2b3">缺少 ' + opt.missing + '（先跑 ' + opt.script + '）</i>';
    return;
  }
  pairs(cases, idTotal, '单子步总耗时（ms）：红=基线 · 绿=新版', 'baseTotal', 'newTotal', 'ms/子步');
  pairs(cases, idSecond, opt.titleSecond, opt.fBase, opt.fNew, 'ms');
  let t = '<table><tr><th>用例</th><th>基线 ms/子步</th><th>新版 ms/子步</th><th>改善</th><th>迭代（基线/新版）</th>'
        + '<th>' + opt.colSecond + '（基线→新版）</th><th>全局回代（基线→新版）</th></tr>';
  for (const c of cases) {
    const imp = 1 - c.newTotal / c.baseTotal;
    t += '<tr><td>' + c.name + '</td><td>' + fmt(c.baseTotal) + '</td><td>' + fmt(c.newTotal) + '</td>'
       + '<td class="' + (imp > 0 ? 'best' : '') + '">' + (imp * 100).toFixed(1) + '%</td>'
       + '<td>' + c.baseIter + ' / ' + c.newIter + '</td>'
       + '<td>' + fmt(c[opt.fBase], 0) + ' → ' + fmt(c[opt.fNew], 0) + '</td>'
       + '<td>' + fmt(c.baseSolve, 0) + ' → ' + fmt(c.newSolve, 0) + '</td></tr>';
  }
  t += '</table><div class="sub" style="margin-top:8px">' + opt.footnote + '</div>';
  document.getElementById(idTable).innerHTML = t;
}

function renderPhase3() {
  renderAB(DATA.phase3, 'p3total', 'p3unpack', 'p3table', {
    missing: 'build/_baseline/phase3_ab.log', script: '_perf/ab_residual.ps1',
    titleSecond: '"解包+判据" 阶段累计耗时（ms；残差就在这一桶）',
    fBase: 'baseUnpack', fNew: 'newUnpack', colSecond: '解包+判据',
    footnote: '迭代数逐例相同 ⇒ 两侧物理一致（残差两个口径逐位相同，放行判据不会分叉）。'
  });
}

function renderPhase4c() {
  renderAB(DATA.phase4c, 'p4total', 'p4solve', 'p4table', {
    missing: 'build/_baseline/phase4c_ab.log', script: '_perf/refresh_report.ps1',
    titleSecond: '"全局回代" 阶段累计耗时（ms；三分量并行所在的那一桶）',
    fBase: 'baseSolve', fNew: 'newSolve', colSecond: '全局回代',
    footnote: '两侧迭代数相同、物理输出逐字相同（末次残差 / 应变 / 能量）⇒ 只是"同样结果、更少时间"；'
            + '1 线程那一例走同一条串行路径，只该有小幅局部性收益。'
  });
}

function renderPhase4d() {
  renderAB(DATA.phase4d, 'p5total', 'p5scatter', 'p5table', {
    missing: 'build/_baseline/phase4d_ab.log', script: '_perf/refresh_report.ps1（A/B 部分）',
    titleSecond: '"投影+散射" 阶段累计耗时（ms；基线为『局部步+散射』之和）',
    fBase: 'baseScatter', fNew: 'newScatter', colSecond: '投影+散射',
    footnote: '两侧迭代数与物理输出逐字相同 ⇒ 只是"同样结果、更少时间"；'
            + '收益随规模递减（大网格那一档本来就是访存受限）。'
  });
}

// ---------------------------------------------------------------- 新增板块的渲染

// 单序列水平条形图（用于"多配置对照"，比较哪一档贵）
function simpleBars(items, targetId, title, unit, fmtV) {
  const W = 1180, rowH = 30, gap = 10, labelW = 200, padR = 190;
  let y = 6, out = '';
  out += el('text', { x: 0, y: y + 12, style: 'fill:#e6e8ee;font-weight:600' }, title);
  y += 24;
  const max = items.reduce((a, b) => Math.max(a, b.value), 0) || 1;
  for (const it of items) {
    const w = Math.max((W - labelW - padR) * it.value / max, 0.5);
    out += el('text', { x: 0, y: y + rowH / 2 + 5, style: 'fill:#c7ced9' }, it.name);
    out += el('rect', { x: labelW, y: y + 2, width: fmt(w, 2), height: rowH - 8,
                        fill: it.color || '#4a9eda', rx: 2 },
              el('title', {}, it.name + ': ' + fmtV(it.value) + ' ' + unit));
    out += el('text', { x: fmt(labelW + w + 8, 1), y: y + rowH / 2 + 5, style: 'fill:#98a2b3' },
              fmtV(it.value) + (it.note ? '   ' + it.note : ''));
    y += rowH + gap;
  }
  document.getElementById(targetId).innerHTML =
    el('svg', { viewBox: '0 0 ' + W + ' ' + y, width: '100%', height: y }, out);
}

// 两值对照（剪切 / 弯曲）：阶段条形 + 明细表 + 规模与物理读数
function renderTwoWay(data, ids, opt) {
  if (!data) {
    document.getElementById(ids.bars).innerHTML =
      '<i style="color:#98a2b3">缺少 ' + opt.missing + '（先跑 ' + opt.script + '）</i>';
    return;
  }
  const bars = data.rows.map(r => ({ name: r.name, baseTotal: r.off, newTotal: r.on }));
  pairs(bars, ids.bars, opt.barTitle, 'baseTotal', 'newTotal', 'ms');
  let t = '<table><tr><th>项（累计 ms）</th><th>关闭</th><th>开启</th><th>变化</th></tr>';
  for (const r of data.rows) {
    const key = r.name.indexOf('单子步') === 0;
    t += '<tr><td>' + r.name + (key ? '（每子步）' : '') + '</td><td>' + fmt(r.off) + '</td><td>'
       + fmt(r.on) + '</td><td class="' + (Math.abs(r.pct) < 1 ? '' : 'bad') + '">'
       + (r.pct > 0 ? '+' : '') + r.pct.toFixed(1) + '%</td></tr>';
  }
  t += '</table>';
  const o = data.off, n = data.on;
  const scale = (k, f) => (o[k] === null || o[k] === undefined) ? '' :
    ('<tr><td>' + f + '</td><td>' + o[k] + '</td><td>' + n[k] + '</td><td class="'
     + (n[k] > o[k] ? 'bad' : 'best') + '">' + (o[k] ? ((100 * (n[k] - o[k]) / o[k]).toFixed(1) + '%') : '—')
     + '</td></tr>');
  t += '<table style="margin-top:14px"><tr><th>规模（同一轮运行）</th><th>关闭</th><th>开启</th><th>变化</th></tr>'
     + scale('constraints', '约束数') + scale('stencils', '弯曲 stencil 数')
     + scale('lNnz', 'L nnz') + scale('factorNnz', '因子 nnz') + '</table>';
  t += '<table style="margin-top:14px"><tr><th>物理读数</th><th>关闭</th><th>开启</th><th>变化</th></tr>'
     + '<tr><td>末次残差 m/s²</td><td>' + o.residual + '</td><td>' + n.residual + '</td><td>—</td></tr>'
     + '<tr><td>最终应变 max</td><td>' + o.strainMax + '</td><td>' + n.strainMax + '</td><td>—</td></tr>'
     + '<tr><td>弹性能</td><td>' + o.energy + '</td><td>' + n.energy + '</td><td>—</td></tr></table>';
  t += '<div class="sub" style="margin-top:8px">' + opt.footnote + '</div>';
  document.getElementById(ids.table).innerHTML = t;
}

function renderShear() {
  renderTwoWay(DATA.shearAB, { bars: 'shbars', table: 'shtable' }, {
    missing: 'build/_baseline/shear_ab.log', script: '_perf/ab_shear.ps1',
    barTitle: '各阶段累计耗时（ms）：红 = 关闭 · 绿 = 开启（这里是"代价"，不是"改善"）',
    footnote: '读法：<code>投影+散射</code>涨得最多（剪切走的就是这条完整通路：每条对角边都要投影+散射），'
            + '而回代只因矩阵略稠而小幅变慢；装配只发生一次（0.3 万子步里摊薄）。'
  });
}

function renderBend() {
  renderTwoWay(DATA.bendAB, { bars: 'bdbars', table: 'bdtable' }, {
    missing: 'build/_baseline/bend_ab.log', script: '_perf/ab_bend.ps1 -Bend 1000',
    barTitle: '各阶段累计耗时（ms）：红 = 关闭 · 绿 = 开启（这里是"代价"，不是"改善"）',
    footnote: '读法：<code>投影+散射</code>几乎没变（弯曲<b>没有</b>局部步/散射 —— 这是选这个形式换来的结构性收益），'
            + '代价几乎全在<b>全局回代</b>上；<code>解包+判据</code>涨是因为残差里必须多算一项弯曲力'
            + '（少算就是假阴性）。'
  });
}

function renderVariants() {
  const V = DATA.variants;
  if (!V) {
    document.getElementById('varbars').innerHTML =
      '<i style="color:#98a2b3">缺少 build/_baseline/bend_variants_ab.log（先跑 _perf/ab_bend_variants.ps1）</i>';
    return;
  }
  const color = (n) => n.indexOf('off') === 0 ? '#adb5bd'
    : (n.indexOf('std/s1') === 0 ? '#d62828'
      : (n.indexOf('rows') === 0 || n.indexOf('cols') === 0 ? '#7ee787'
        : (n.indexOf('checker') === 0 ? '#fb8500' : '#4a9eda')));
  const off = V.byName['off'] ? V.byName['off'].perStep : Math.min(...V.configs.map(c => c.perStep));
  simpleBars(V.configs.map(c => ({ name: c.name, value: c.perStep, color: color(c.name),
                                   note: c.name === 'off' ? '' : c.vsOff })),
             'varbars', '单子步（ms/子步）：灰=不开弯曲 · 绿=单向（可用）· 红=全采样 · 橙=棋盘',
             'ms', v => v.toFixed(2));
  let t = '<table><tr><th>配置</th><th>刚度</th><th>stencil</th><th>单子步 ms</th><th>vs off</th>'
        + '<th>L nnz</th><th>因子 nnz</th><th>应变 max</th><th>弹性能</th><th>结论</th></tr>';
  for (const c of V.configs) {
    let verdict = '—';
    // 判定顺序有讲究：先看"隔行"，再看单向/全采样/棋盘 —— 否则 "checker/s2" 会被
    // 归到"棋盘"那一档，漏掉它其实也是退化的（stride ≥ 2 才是决定性因素）。
    if (c.name === 'off') verdict = '基准（不启用弯曲）';
    else if (c.name.indexOf('/s2') >= 0 || c.name.indexOf('/s3') >= 0) {
      verdict = '<span class="bad">退化：付了代价、几乎无效果</span>';
    } else if (c.name.indexOf('rows') === 0 || c.name.indexOf('cols') === 0) {
      verdict = '<b>有效 · 因子 nnz −49 %</b>，但各向异性';
    } else if (c.name.indexOf('std/s1') === 0) verdict = '<b>各向同性，代价 +75 %</b>';
    else verdict = '比单向更贵也更弱 · 不推荐';
    t += '<tr><td>' + c.name + '</td><td>' + c.k + '</td><td>' + c.stencils + '</td><td>' + fmt(c.perStep)
       + '</td><td>' + c.vsOff + '</td><td>' + (c.lNnz ?? '—') + '</td><td>' + (c.factorNnz ?? '—')
       + '</td><td>' + (c.strainMax !== undefined ? c.strainMax.toFixed(5) : '—') + '</td><td>'
       + (c.elastic !== undefined ? c.elastic.toFixed(4) : '—') + '</td><td>' + verdict + '</td></tr>';
  }
  t += '</table><div class="sub" style="margin-top:8px">'
     + '"应变 max"是同一工况结束时的最大相对应变（越小越挺）；<b>隔行采样</b>那一行的应变与 off 几乎相同'
     + '（0.0578 vs 0.0586），却多付了 +11.5 % 的时间 —— 这就是"退化"的量化含义。详见 docs/perf.md §13。</div>';
  document.getElementById('vartable').innerHTML = t;
}

function renderInertial() {
  const I = DATA.inertial;
  const el2 = document.getElementById('inerttable');
  if (!I) { el2.innerHTML = '<i style="color:#98a2b3">缺少 build/_baseline/inertial_ab3.log</i>'; return; }
  let t = '<table><tr><th>线程</th><th>旧-single</th><th>控制-single<br>（同源码另一次构建）</th><th>新-for</th>'
        + '<th>新 vs 旧</th><th>新 vs 控制</th><th>构建/噪声地板</th></tr>';
  const names = { '1': '1（机制失效档 → 理论差值 = 0）', '4': '4（默认）', '18': '18' };
  for (const b of I) {
    const get = (n) => (b.rows.find(r => r.name === n) || {}).ms;
    const cell = (v) => v === undefined ? '—' : fmt(v);
    const d = (v) => v === undefined ? '—' : (v > 0 ? '+' : '') + v.toFixed(1) + '%';
    t += '<tr><td>' + (names[b.threads] || b.threads) + '</td><td>' + cell(get('旧-single')) + '</td><td>'
       + cell(get('控制-single')) + '</td><td><b>' + cell(get('新-for')) + '</b></td><td class="bad">'
       + d(b.newVsOld) + '</td><td class="bad">' + d(b.newVsCtrl) + '</td><td>' + d(b.floor) + '</td></tr>';
  }
  t += '</table>';
  el2.innerHTML = t;
}

// 四套约束集的阶段构成（堆叠条 + 明细表）。
// 与 §二 同一套口径：阶段桶是**累计值**，必须除以子步数才是 ms/子步；
// 每一组（工况）独立归一化，所以组间不要比长度，只比组内构成。
function renderConstraintSets() {
  const CS = DATA.constraintSets || [];
  if (!CS.length) {
    document.getElementById('csbars').innerHTML =
      '<i style="color:#98a2b3">缺少 build/_perf/logs_cs/*.txt（先跑 _perf/constraint_sets.ps1）</i>';
    return;
  }
  const W = 1180, rowH = 30, labelW = 150, padR = 250, axisW = W - labelW - padR;
  let y = 14, out = '';
  for (const w of CS) {
    const rows = w.rows.map((r) => {
      const st = {}; let sum = 0;
      for (const sg of STAGES) { const v = (r.stages[sg.key] || 0) / (r.steps || 1); st[sg.key] = v; sum += v; }
      return { label: r.label, total: r.total, st, sum, raw: r };
    });
    const base = rows[0].total;
    const maxAbs = Math.max.apply(null, rows.map((r) => r.total));
    out += el('text', { x: 0, y: y + 12, style: 'fill:#e6e8ee;font-weight:600' }, w.label);
    out += el('text', { x: labelW + axisW + 8, y: y + 12, style: 'fill:#6b7280' },
              '标尺：最长条 = ' + fmt(maxAbs, maxAbs < 10 ? 2 : 1) + ' ms/子步（本组独立归一化）');
    y += 18;
    const axisBottom = y + rows.length * rowH;
    for (let k = 0; k <= 4; ++k) {
      const v = maxAbs * k / 4, x = labelW + v * (axisW / maxAbs);
      out += el('line', { x1: fmt(x, 1), y1: y, x2: fmt(x, 1), y2: axisBottom + 4, stroke: '#232833' });
      out += el('text', { x: fmt(x, 1), y: axisBottom + 16, 'text-anchor': 'middle', style: 'fill:#6b7280' },
                fmt(v, maxAbs < 10 ? 2 : 0));
    }
    out += el('text', { x: fmt(labelW + axisW + 8, 1), y: axisBottom + 16, style: 'fill:#6b7280' }, 'ms/子步');
    for (const r of rows) {
      const d = 100 * (r.total - base) / base;
      out += el('text', { x: 0, y: y + rowH / 2 + 4, style: 'fill:' + (r.label.indexOf('默认') === 0 ? '#7ee787' : '#c7ced9') },
                r.label);
      let x = labelW;
      for (const sg of STAGES) {
        const v = r.st[sg.key] || 0;
        const wpx = v * (axisW / maxAbs);
        if (wpx < 0.4) continue;
        const tip = sg.label + '：' + fmt(v) + ' ms/子步（' + fmt(100 * v / r.sum, 1) + '%）' +
                    (sg.key === 'assemble' ? '　※一次性成本，按累计÷子步数折算' : '');
        out += el('rect', { x: fmt(x, 2), y: y + 3, width: fmt(Math.max(wpx, 0.4), 2), height: rowH - 8,
                            fill: sg.color, rx: 1.5 }, el('title', {}, tip));
        if (wpx > 46) {
          out += el('text', { x: fmt(x + wpx / 2, 1), y: y + rowH / 2 + 4, 'text-anchor': 'middle',
                              style: 'fill:#0b0d12;font-weight:600' }, fmt(v, 2));
        }
        x += wpx;
      }
      out += el('text', { x: fmt(labelW + axisW + 8, 1), y: y + rowH / 2 + 4,
                          style: 'fill:' + (d > 1 ? '#ff9f9f' : '#7ee787') },
                fmt(r.total, 2) + ' ms' + (Math.abs(d) < 0.05 ? '' : '  (' + (d > 0 ? '+' : '') + fmt(d, 1) + '%)'));
      y += rowH + 6;
    }
    y += 18;
  }
  document.getElementById('csbars').innerHTML =
    el('svg', { viewBox: '0 0 ' + W + ' ' + y, width: '100%', height: y }, out);

  let t = '';
  for (const w of CS) {
    const base = w.rows[0].total;
    const bendOnly = w.rows.find((r) => r.cfg === 'bend');
    const both = w.rows.find((r) => r.cfg === 'both');
    t += '<div class="sub" style="margin:16px 0 6px"><b>' + w.label + '</b></div>';
    t += '<table><tr><th>约束集</th><th>边/约束</th><th>弯曲 stencil</th><th>L nnz</th><th>因子 nnz</th>'
       + '<th>单子步 ms</th><th>Δ vs 默认</th>'
       + STAGES.map((s) => '<th>' + s.label + '<br>ms / 占比</th>').join('') + '</tr>';
    for (const r of w.rows) {
      const sum = STAGES.reduce((a, s) => a + (r.stages[s.key] || 0) / (r.steps || 1), 0);
      const d = 100 * (r.total - base) / base;
      t += '<tr><td>' + r.label + '</td><td>' + r.edges + '</td><td>' + (r.stencils || '—') + '</td><td>'
         + r.lNnz + '</td><td>' + r.factorNnz + '</td><td>' + fmt(r.total, 3) + '</td><td class="'
         + (d > 1 ? 'bad' : '') + '">' + (r.label.indexOf('默认') === 0 ? '—' : ('+' + fmt(d, 1) + '%')) + '</td>';
      for (const s of STAGES) {
        const v = (r.stages[s.key] || 0) / (r.steps || 1);
        t += '<td>' + fmt(v, 3) + '<br><span class="sub">' + fmt(100 * v / sum, 1) + '%</span></td>';
      }
      t += '</tr>';
    }
    t += '</table>';
    if (bendOnly && both) {
      const fB = bendOnly.factorNnz, fX = both.factorNnz;
      const pct = 100 * (fX - fB) / fB;
      t += '<div class="sub" style="margin-top:8px"><b>⚠️ 反直觉但可复现</b>：'
         + '"剪切 + 弯曲"的 <code>L</code> 非零比"只加弯曲"<b>更多</b>（' + bendOnly.lNnz + ' → ' + both.lNnz
         + '，+' + fmt(100 * (both.lNnz - bendOnly.lNnz) / bendOnly.lNnz, 1) + '%），'
         + '但<b>因子 nnz 反而更少</b>（' + fB + ' → ' + fX + '，' + fmt(pct, 1) + '%）⇒ 回代更快'
         + '（' + fmt(bendOnly.total, 2) + ' → ' + fmt(both.total, 2) + ' ms/子步）。<br>'
         + '原因：<b>填充量由"消元顺序 + 整体结构"共同决定，不是约束数的单调函数</b> —— 加上剪切对角边之后，'
         + 'AMD 给出的排序变了，而新排序对弯曲那部分结构恰好更省填充。'
         + '这既是"多一类约束不一定更慢"的例子，也提示"给含弯曲的矩阵找一个更好的排序"还有空间'
         + '（本项目只做过 AMD + 自然序，没有试过 ND/其他排序）。</div>';
    }
  }
  document.getElementById('cstable').innerHTML = t;
}

renderLegend(); renderStack(); renderLines(); renderHist(); renderPhase3(); renderPhase4c(); renderPhase4d();
renderShear(); renderBend(); renderConstraintSets(); renderVariants(); renderInertial();
renderViewer(); renderTables();
</script>
</body></html>`;

// 采集时间那一行在模板之外插进去（单点实现，见上面的 withStamps）。
const html = withStamps(rawHtml);

const outFile = path.join(ROOT, 'build', '_perf', 'perf-report.html');
fs.writeFileSync(outFile, Buffer.from(html, 'utf8'));   // 显式 UTF-8 无 BOM

// 控制台摘要
console.log('=== 线程扫描（ms/子步，最小值）===');
for (const w of sweep) {
  const cells = w.series.map(s => `${s.threads}T=${fmt(s.total)}`).join('  ');
  const best = w.series.reduce((a, b) => (a.total < b.total ? a : b));
  console.log(`  ${w.id.padEnd(14)} ${cells}   最优=${best.threads} 线程`);
}
console.log('\n=== 成对对照 ===');
for (const h of history) {
  console.log(`  ${h.label} (${h.threads}): ${h.a.name} ${fmt(h.a.total)} → ${h.b.name} ${fmt(h.b.total)}  (${((1 - h.b.total / h.a.total) * 100).toFixed(1)}%)`);
}
console.log('\n=== 查看器时间线 ===');
for (const v of viewer) console.log(`  ${v.name}: ${v.pts.length} 个采样点`);
console.log('\n=== Phase 3（残差并行化，同一会话 A/B）===');
for (const c of phase3) {
  console.log(`  ${c.name.padEnd(18)} ${fmt(c.baseTotal)} → ${fmt(c.newTotal)}  (${c.delta.toFixed(1)}%)  ` +
              `解包+判据 ${fmt(c.baseUnpack, 1)} → ${fmt(c.newUnpack, 1)}  迭代 ${c.baseIter}/${c.newIter}`);
}
console.log('\n=== Phase 4c（全局步按 x/y/z 三分量并行，同一会话 A/B）===');
for (const c of phase4c) {
  console.log(`  ${c.name.padEnd(18)} ${fmt(c.baseTotal)} → ${fmt(c.newTotal)}  (${c.delta.toFixed(1)}%)  ` +
              `全局回代 ${fmt(c.baseSolve, 1)} → ${fmt(c.newSolve, 1)}  迭代 ${c.baseIter}/${c.newIter}`);
}
console.log(`\nHTML: ${outFile}  (${fs.statSync(outFile).size} 字节)`);
