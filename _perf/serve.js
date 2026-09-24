// _perf/serve.js —— 只为看这份性能看板起的静态服务（仅本机回环，端口 8137）。
// 它不碰 DSH 的 Web GUI（3081），只是把**产物目录** _perf/ 暴露成 http://127.0.0.1:8137/。
// 脚本自己在 _perf/ 下（2026-09-24 从 _perf/ 迁入版本控制），服务的是 _perf/。
const http = require('http');
const fs = require('fs');
const path = require('path');

const ROOT = path.join(path.resolve(__dirname, '..'), 'build', '_perf');
const PORT = 8137;

http.createServer((req, res) => {
  const name = decodeURIComponent((req.url || '/').split('?')[0]);
  const rel = name === '/' ? 'perf-report.html' : name.replace(/^\/+/, '');
  const file = path.join(ROOT, rel);
  if (!file.startsWith(ROOT)) { res.writeHead(403); return res.end('forbidden'); }
  fs.readFile(file, (err, buf) => {
    if (err) { res.writeHead(404); return res.end('not found'); }
    const ext = path.extname(file).toLowerCase();
    const type = ext === '.html' ? 'text/html; charset=utf-8'
               : ext === '.json' ? 'application/json; charset=utf-8'
               : 'text/plain; charset=utf-8';
    res.writeHead(200, { 'Content-Type': type });
    res.end(buf);
  });
}).listen(PORT, '127.0.0.1', () => console.log(`serving ${ROOT} at http://127.0.0.1:${PORT}/`));
