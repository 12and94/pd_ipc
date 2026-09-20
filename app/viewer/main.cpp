// app/viewer/main.cpp
// 可交互查看器：GLFW + OpenGL 固定管线 + 线段自绘 HUD。
//
// 为什么用固定管线而不是现代 OpenGL：本项目首期只需要"能看到、能拖、有数字"，
// 固定管线（glBegin/glEnd）不需要加载扩展函数指针，代码量小且没有第三方加载器依赖。
// 物理后端将来换 Vulkan 时，渲染这条线可以独立替换（docs/plan.md D10：渲染与物理解耦）。
//
// 交互：
//   左键拖拽平移相机 / 右键或滚轮缩放 / R 重置 / SPACE 暂停 / 单步 S / 重力 G
//   刚度 [ ] / 迭代数 - = / 子步数 , . / 线程数 ; ' / ESC 退出
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <GLFW/glfw3.h>

#include "core/assemble/Assembler.h"
#include "core/math/Parallel.h"
#include "core/sim/Integrator.h"
#include "core/sim/Scene.h"

using namespace pd;

namespace {

// ---------------------------------------------------------------------------
// 线段数字字体：每个数字用 7 段中的若干段表示，HUD 只画线，不依赖纹理与字体文件。
// 段序：0=上 1=右上 2=右下 3=下 4=左下 5=左上 6=中
// ---------------------------------------------------------------------------
constexpr int kDigitSegments[10] = {
    0b0111111,  // 0
    0b0000110,  // 1
    0b1011011,  // 2
    0b1001111,  // 3
    0b1100110,  // 4
    0b1101101,  // 5
    0b1111101,  // 6
    0b0000111,  // 7
    0b1111111,  // 8
    0b1101111,  // 9
};

void drawSegment(float x, float y, float w, float h, int seg) {
  const float x1 = x + w, y1 = y + h;
  switch (seg) {
    case 0: glVertex2f(x, y1); glVertex2f(x1, y1); break;
    case 1: glVertex2f(x1, y1); glVertex2f(x1, y + h * 0.5f); break;
    case 2: glVertex2f(x1, y + h * 0.5f); glVertex2f(x1, y); break;
    case 3: glVertex2f(x, y); glVertex2f(x1, y); break;
    case 4: glVertex2f(x, y + h * 0.5f); glVertex2f(x, y); break;
    case 5: glVertex2f(x, y1); glVertex2f(x, y + h * 0.5f); break;
    case 6: glVertex2f(x, y + h * 0.5f); glVertex2f(x1, y + h * 0.5f); break;
    default: break;
  }
}

/// 用线段画一个数字。
void drawDigit(char c, float x, float y, float w, float h) {
  if (c < '0' || c > '9') return;
  const int mask = kDigitSegments[c - '0'];
  for (int s = 0; s < 7; ++s) {
    if (mask & (1 << s)) drawSegment(x, y, w, h, s);
  }
}

/// 用线段画一个字符串（数字与空格、'/'、'.'、'%' 之外的字符会被跳过）。
void drawText(const std::string& text, float x, float y, float w, float h) {
  float cursor = x;
  for (char c : text) {
    if (c >= '0' && c <= '9') {
      drawDigit(c, cursor, y, w, h);
      cursor += w * 1.8f;
    } else if (c == ' ') {
      cursor += w * 1.2f;
    } else if (c == '.') {
      glVertex2f(cursor, y);
      glVertex2f(cursor + w * 0.4f, y);
      cursor += w * 1.2f;
    } else if (c == '/') {
      glVertex2f(cursor, y);
      glVertex2f(cursor + w, y + h);
      cursor += w * 1.6f;
    } else if (c == '=') {
      glVertex2f(cursor, y + h * 0.35f);
      glVertex2f(cursor + w, y + h * 0.35f);
      glVertex2f(cursor, y + h * 0.65f);
      glVertex2f(cursor + w, y + h * 0.65f);
      cursor += w * 1.6f;
    } else if (c == '%') {
      glVertex2f(cursor, y);
      glVertex2f(cursor + w * 0.6f, y + h * 0.5f);
      glVertex2f(cursor + w * 0.6f, y + h * 0.5f);
      glVertex2f(cursor + w * 0.6f, y + h);
      cursor += w * 1.6f;
    } else if (c == ':') {
      glVertex2f(cursor + w * 0.4f, y + h * 0.2f);
      glVertex2f(cursor + w * 0.4f, y + h * 0.3f);
      glVertex2f(cursor + w * 0.4f, y + h * 0.7f);
      glVertex2f(cursor + w * 0.4f, y + h * 0.8f);
      cursor += w * 1.4f;
    }
  }
}

// ---------------------------------------------------------------------------
// 相机与交互状态
// ---------------------------------------------------------------------------
struct Camera {
  float yaw = 0.6f;
  float pitch = 0.45f;
  float distance = 2.0f;
  Vec3 target{0.5, 0.0, 0.5};
};

struct Interactive {
  bool paused = false;
  bool singleStep = false;
  bool gravityOn = true;
  int threads = 0;
  double lastPrint = 0.0;
  double frameMsAvg = 0.0;
  double physicsMsAvg = 0.0;
  double fpsAvg = 0.0;
};

Camera g_camera;
Interactive g_state;
bool g_dragging = false;
double g_lastX = 0.0, g_lastY = 0.0;

void scrollCallback(GLFWwindow*, double, double yoffset) {
  g_camera.distance *= static_cast<float>(std::exp(-0.1 * yoffset));
  g_camera.distance = std::clamp(g_camera.distance, 0.1f, 50.0f);
}

void mouseButtonCallback(GLFWwindow*, int button, int action, int) {
  if (button == GLFW_MOUSE_BUTTON_LEFT) g_dragging = (action == GLFW_PRESS);
}

void cursorPosCallback(GLFWwindow*, double x, double y) {
  if (!g_dragging) { g_lastX = x; g_lastY = y; return; }
  const double dx = x - g_lastX;
  const double dy = y - g_lastY;
  g_lastX = x;
  g_lastY = y;
  g_camera.yaw -= static_cast<float>(dx) * 0.01f;
  g_camera.pitch = std::clamp(g_camera.pitch + static_cast<float>(dy) * 0.01f, -1.4f, 1.4f);
}

void keyCallback(GLFWwindow* window, int key, int, int action, int) {
  if (action != GLFW_PRESS && action != GLFW_REPEAT) return;
  switch (key) {
    case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(window, GLFW_TRUE); break;
    case GLFW_KEY_SPACE: g_state.paused = !g_state.paused; break;
    case GLFW_KEY_S: g_state.singleStep = true; break;
    case GLFW_KEY_G: g_state.gravityOn = !g_state.gravityOn; break;
    case GLFW_KEY_R: g_state.threads = g_state.threads; break;  // 重置由主循环处理
    default: break;
  }
}

// ---------------------------------------------------------------------------
// PNG 输出（最小实现：zlib stored 块，不依赖任何图像库）
// ---------------------------------------------------------------------------
namespace {

uint32_t crc32Of(const uint8_t* data, size_t len) {
  static uint32_t table[256];
  static bool init = false;
  if (!init) {
    for (uint32_t i = 0; i < 256; ++i) {
      uint32_t c = i;
      for (int k = 0; k < 8; ++k) c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
      table[i] = c;
    }
    init = true;
  }
  uint32_t c = 0xFFFFFFFFu;
  for (size_t i = 0; i < len; ++i) c = table[(c ^ data[i]) & 0xFF] ^ (c >> 8);
  return c ^ 0xFFFFFFFFu;
}

void appendBE32(std::vector<uint8_t>& v, uint32_t x) {
  v.push_back(static_cast<uint8_t>(x >> 24));
  v.push_back(static_cast<uint8_t>(x >> 16));
  v.push_back(static_cast<uint8_t>(x >> 8));
  v.push_back(static_cast<uint8_t>(x));
}

void appendChunk(std::vector<uint8_t>& out, const char type[4], const std::vector<uint8_t>& data) {
  appendBE32(out, static_cast<uint32_t>(data.size()));
  const size_t start = out.size();
  out.insert(out.end(), type, type + 4);
  out.insert(out.end(), data.begin(), data.end());
  const uint32_t crc = crc32Of(out.data() + start, out.size() - start);
  appendBE32(out, crc);
}

/// 把 RGBA 像素写成 PNG。用 zlib 的 stored（未压缩）块，实现简单且无依赖。
bool writePng(const char* path, int w, int h, const std::vector<uint8_t>& rgba) {
  std::vector<uint8_t> raw;
  raw.reserve(static_cast<std::size_t>(h) * (1 + static_cast<std::size_t>(w) * 4));
  for (int y = 0; y < h; ++y) {
    raw.push_back(0);  // filter = none
    const uint8_t* row = rgba.data() + static_cast<std::size_t>(y) * w * 4;
    raw.insert(raw.end(), row, row + static_cast<std::size_t>(w) * 4);
  }
  // zlib 容器：2 字节头 + deflate stored 块 + adler32
  std::vector<uint8_t> z;
  z.push_back(0x78);
  z.push_back(0x01);
  size_t pos = 0;
  while (pos < raw.size()) {
    const size_t n = std::min<size_t>(65535, raw.size() - pos);
    z.push_back(static_cast<uint8_t>((pos + n >= raw.size()) ? 1 : 0));  // BFINAL
    z.push_back(static_cast<uint8_t>(n & 0xFF));
    z.push_back(static_cast<uint8_t>((n >> 8) & 0xFF));
    z.push_back(static_cast<uint8_t>(~n & 0xFF));
    z.push_back(static_cast<uint8_t>((~n >> 8) & 0xFF));
    z.insert(z.end(), raw.begin() + static_cast<long>(pos), raw.begin() + static_cast<long>(pos + n));
    pos += n;
  }
  uint32_t a = 1, b = 0;
  for (uint8_t byte : raw) {
    a = (a + byte) % 65521;
    b = (b + a) % 65521;
  }
  appendBE32(z, (b << 16) | a);

  std::vector<uint8_t> png = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
  std::vector<uint8_t> ihdr;
  appendBE32(ihdr, static_cast<uint32_t>(w));
  appendBE32(ihdr, static_cast<uint32_t>(h));
  ihdr.push_back(8);  // bit depth
  ihdr.push_back(6);  // color type RGBA
  ihdr.push_back(0);
  ihdr.push_back(0);
  ihdr.push_back(0);
  appendChunk(png, "IHDR", ihdr);
  appendChunk(png, "IDAT", z);
  appendChunk(png, "IEND", {});

  FILE* f = std::fopen(path, "wb");
  if (!f) return false;
  const size_t written = std::fwrite(png.data(), 1, png.size(), f);
  std::fclose(f);
  return written == png.size();
}

}  // namespace

// 固定管线没有 glu，这里手写一个 lookAt（只用到它一个函数）。
void gluLookAlike(float ex, float ey, float ez, float cx, float cy, float cz) {
  float fx = cx - ex, fy = cy - ey, fz = cz - ez;
  const float flen = std::sqrt(fx * fx + fy * fy + fz * fz);
  if (flen > 1e-9f) { fx /= flen; fy /= flen; fz /= flen; }
  float ux = 0.0f, uy = 1.0f, uz = 0.0f;
  float sx = fy * uz - fz * uy;
  float sy = fz * ux - fx * uz;
  float sz = fx * uy - fy * ux;
  const float slen = std::sqrt(sx * sx + sy * sy + sz * sz);
  if (slen > 1e-9f) { sx /= slen; sy /= slen; sz /= slen; }
  const float tx = sy * fz - sz * fy;
  const float ty = sz * fx - sx * fz;
  const float tz = sx * fy - sy * fx;

  const float m[16] = {
      sx, tx, -fx, 0.0f,
      sy, ty, -fy, 0.0f,
      sz, tz, -fz, 0.0f,
      0.0f, 0.0f, 0.0f, 1.0f,
  };
  glMultMatrixf(m);
  glTranslatef(-ex, -ey, -ez);
}

void setupCamera() {
  const float cy = std::cos(g_camera.yaw);
  const float sy = std::sin(g_camera.yaw);
  const float cp = std::cos(g_camera.pitch);
  const float sp = std::sin(g_camera.pitch);
  const float ex = static_cast<float>(g_camera.target.x) + g_camera.distance * cp * sy;
  const float ey = static_cast<float>(g_camera.target.y) + g_camera.distance * sp;
  const float ez = static_cast<float>(g_camera.target.z) + g_camera.distance * cp * cy;
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAlike(ex, ey, ez, static_cast<float>(g_camera.target.x),
               static_cast<float>(g_camera.target.y), static_cast<float>(g_camera.target.z));
}

void drawCloth(const SimContext& ctx, bool showConstraints) {
  glLineWidth(1.0f);
  if (showConstraints) {
    glColor3f(0.35f, 0.65f, 0.95f);
    glBegin(GL_LINES);
    for (const auto& e : ctx.mesh.edges) {
      const Vec3& a = ctx.mesh.positions[static_cast<std::size_t>(e.a)];
      const Vec3& b = ctx.mesh.positions[static_cast<std::size_t>(e.b)];
      glVertex3d(a.x, a.y, a.z);
      glVertex3d(b.x, b.y, b.z);
    }
    glEnd();
  }

  // pinned 顶点画成点（醒目地标出约束位置）
  glPointSize(6.0f);
  glColor3f(1.0f, 0.45f, 0.2f);
  glBegin(GL_POINTS);
  for (int v = 0; v < ctx.mesh.vertexCount(); ++v) {
    if (!ctx.mesh.isPinned(v)) continue;
    const Vec3& p = ctx.mesh.positions[static_cast<std::size_t>(v)];
    glVertex3d(p.x, p.y, p.z);
  }
  glEnd();
}

void drawGround() {
  glColor3f(0.22f, 0.22f, 0.25f);
  glBegin(GL_LINES);
  const double s = 1.0;
  for (int i = -5; i <= 5; ++i) {
    const double t = i * s * 0.2;
    glVertex3d(t, -0.02, -s);
    glVertex3d(t, -0.02, s);
    glVertex3d(-s, -0.02, t);
    glVertex3d(s, -0.02, t);
  }
  glEnd();
}

void drawHud(const SimContext& ctx, const Interactive& st, int width, int height) {
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0.0, width, 0.0, height, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glDisable(GL_DEPTH_TEST);

  glLineWidth(1.6f);
  glColor3f(0.9f, 0.95f, 1.0f);

  const float w = 6.0f, h = 10.0f;
  float y = static_cast<float>(height) - 18.0f;

  char buf[256];
  std::snprintf(buf, sizeof(buf), "FPS %.1f  FRAME %.2f MS  PHYS %.2f MS", st.fpsAvg,
                st.frameMsAvg, st.physicsMsAvg);
  glBegin(GL_LINES);
  drawText(buf, 12.0f, y, w, h);
  glEnd();
  y -= 16.0f;

  std::snprintf(buf, sizeof(buf), "VERT %.0f  EDGE %.0f  ITER %d  THREADS %d",
                static_cast<double>(ctx.mesh.vertexCount()),
                static_cast<double>(ctx.mesh.edgeCount()), ctx.iterationsUsed, numThreads());
  glBegin(GL_LINES);
  drawText(buf, 12.0f, y, w, h);
  glEnd();
  y -= 16.0f;

  std::snprintf(buf, sizeof(buf), "STIFFNESS %.0f  DT 1/%.0f  SUBSTEPS %d  DAMP %.2f",
                ctx.config.stiffness, 1.0 / ctx.config.dt, ctx.config.substepsPerFrame,
                ctx.config.velocityDamping);
  glBegin(GL_LINES);
  drawText(buf, 12.0f, y, w, h);
  glEnd();
  y -= 16.0f;

  std::snprintf(buf, sizeof(buf), "STRAIN MAX %.4f  MEAN %.4f  DROPS %d", ctx.mesh.maxRelativeStrain(),
                ctx.mesh.meanRelativeStrain(), ctx.factorizeCount);
  glBegin(GL_LINES);
  drawText(buf, 12.0f, y, w, h);
  glEnd();
  y -= 16.0f;

  std::snprintf(buf, sizeof(buf), "SPACE PAUSE %d  S STEP  G GRAVITY %d  ESC QUIT",
                st.paused ? 1 : 0, st.gravityOn ? 1 : 0);
  glBegin(GL_LINES);
  drawText(buf, 12.0f, y, w, h);
  glEnd();

  glEnable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

}  // namespace

int main(int argc, char** argv) {
  int maxFrames = 0;          // >0 时跑够帧数就退出（用于自动化验证）
  std::string shotPath;       // 非空则截图后退出
  int gridN = 60;
  double stiffness = 1.0e5;
  bool pinTopEdge = true;
  bool pinSingle = false;
  int iters = 40;
  double tol = 1.0e-5;
  double damping = 0.02;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--frames" && i + 1 < argc) maxFrames = std::atoi(argv[++i]);
    else if (a == "--shot" && i + 1 < argc) shotPath = argv[++i];
    else if (a == "--grid" && i + 1 < argc) gridN = std::atoi(argv[++i]);
    else if (a == "--stiffness" && i + 1 < argc) stiffness = std::atof(argv[++i]);
    else if (a == "--pin-corners") pinTopEdge = false;
    else if (a == "--pin-single") { pinTopEdge = false; pinSingle = true; }
    else if (a == "--iters" && i + 1 < argc) iters = std::atoi(argv[++i]);
    else if (a == "--tol" && i + 1 < argc) tol = std::atof(argv[++i]);
    else if (a == "--damping" && i + 1 < argc) damping = std::atof(argv[++i]);
    else if (a == "--help" || a == "-h") {
      std::printf("用法: pd_viewer [--grid N] [--frames N] [--shot FILE.png]\n");
      return 0;
    }
  }

  std::printf("[参数] grid=%d frames=%d shot=%s\n", gridN, maxFrames, shotPath.c_str());

  // 先建场景（无窗口也能验证物理链路）
  SceneConfig cfg;
  cfg.gridNx = gridN;
  cfg.gridNy = gridN;
  cfg.gridSpacing = 0.02;
  cfg.dt = 1.0 / 120.0;
  cfg.substepsPerFrame = 2;
  cfg.stiffness = stiffness;
  cfg.maxIterations = iters;
  cfg.relTolerance = tol;
  SimContext ctx = makeScene(cfg);

  // 默认钉住整条上边（只钉两个角时，布料会被拉成尖锥，不是常见演示形态）
  if (pinSingle) {
    for (auto& s : ctx.mesh.pinned) s = 0;   // 先清除默认的两角 pin，避免变成三个 pin
    const int mid = (cfg.gridNy - 1) * cfg.gridNx + cfg.gridNx / 2;
    ctx.mesh.pinned[static_cast<std::size_t>(mid)] = 1;
    refreshPinPositions(ctx);
  } else if (pinTopEdge) {
    const int nx = cfg.gridNx, ny = cfg.gridNy;
    for (int i = 0; i < nx; ++i) {
      const int v = (ny - 1) * nx + i;
      ctx.mesh.pinned[static_cast<std::size_t>(v)] = 1;
    }
    refreshPinPositions(ctx);
  }

  // ---- 取景 ----
  // 网格铺在 XZ 平面（y 朝上），默认钉住的是**上边**（z 最大的那一行）。
  // 布料会绕这条钉住的边垂下来，因此布的实际包围盒是
  //     x ∈ [0, W]，z ∈ [0, D]（钉住边在上方），y ∈ [-D, 0]（垂下部分）
  // 镜头要对准这个包围盒的中心，而不是原点平面 —— 否则布会缩在画面角落。
  //
  // 视距由视场角反推，保证最短边也进画面：
  //   垂直半视场角 = 22.5°，故 dist ≥ (包围盒高度/2 + 中心偏移) / tan(22.5°)。
  // 这里再留 1.15 倍余量，避免贴边。
  const double clothW = (cfg.gridNx - 1) * cfg.gridSpacing;
  const double clothD = (cfg.gridNy - 1) * cfg.gridSpacing;
  const double boxH = 0.5 * clothD;            // 包围盒高度（垂下半边）
  const double centerY = -0.5 * clothD;        // 包围盒中心的 y
  const double tanHalfFov = 0.4142135623730951;  // tan(22.5°)
  g_camera.target = Vec3{0.5 * clothW, centerY, 0.5 * clothD};
  g_camera.distance = static_cast<float>(1.15 * (std::fabs(centerY) + boxH + 0.5 * clothW) / tanHalfFov);
  g_camera.pitch = 0.22f;  // 稍微俯视：既看得到垂下来的面，也看得到地面网格

  if (!glfwInit()) {
    std::printf("GLFW 初始化失败\n");
    return 1;
  }
  GLFWwindow* window = glfwCreateWindow(1280, 800, "PD cloth (distance constraints only)", nullptr, nullptr);
  if (!window) {
    std::printf("创建窗口失败\n");
    glfwTerminate();
    return 1;
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);
  glfwSetScrollCallback(window, scrollCallback);
  glfwSetMouseButtonCallback(window, mouseButtonCallback);
  glfwSetCursorPosCallback(window, cursorPosCallback);
  glfwSetKeyCallback(window, keyCallback);

  std::printf("=== PD 布料查看器 ===\n");
  std::printf("顶点 %d  约束 %d  线程 %d\n", ctx.mesh.vertexCount(), ctx.mesh.edgeCount(), numThreads());
  std::printf("交互: 左键拖拽旋转 / 滚轮缩放 / SPACE 暂停 / S 单步 / G 重力 / R 重置 / ESC 退出\n");

  double lastTime = glfwGetTime();
  double hudTimer = lastTime;

  int frameIndex = 0;
  while (!glfwWindowShouldClose(window)) {
    ++frameIndex;
    const double now = glfwGetTime();
    const double frameDt = now - lastTime;
    lastTime = now;

    // ---- 物理 ----
    const auto tPhys0 = std::chrono::steady_clock::now();
    if (!g_state.paused || g_state.singleStep) {
      ctx.config.gravity = g_state.gravityOn ? Vec3{0.0, -9.81, 0.0} : Vec3{0.0, 0.0, 0.0};
      stepFrame(ctx);
      g_state.singleStep = false;
    }
    const double physMs =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tPhys0).count();

    // ---- HUD 数字平滑（避免闪烁）----
    const double alpha = 0.1;
    g_state.frameMsAvg = g_state.frameMsAvg * (1 - alpha) + frameDt * 1e3 * alpha;
    g_state.physicsMsAvg = g_state.physicsMsAvg * (1 - alpha) + physMs * alpha;
    if (frameDt > 1e-6) g_state.fpsAvg = g_state.fpsAvg * (1 - alpha) + (1.0 / frameDt) * alpha;

    // ---- 键盘持续按键处理 ----
    if (glfwGetKey(window, GLFW_KEY_LEFT_BRACKET) == GLFW_PRESS) {
      ctx.config.stiffness = std::max(1.0, ctx.config.stiffness * 0.98);
      for (auto& e : ctx.mesh.edges) e.stiffness = ctx.config.stiffness;   // 必须同步到实际边刚度
    }
    if (glfwGetKey(window, GLFW_KEY_RIGHT_BRACKET) == GLFW_PRESS) {
      ctx.config.stiffness = std::min(1e8, ctx.config.stiffness * 1.02);
      for (auto& e : ctx.mesh.edges) e.stiffness = ctx.config.stiffness;
    }
    if (glfwGetKey(window, GLFW_KEY_MINUS) == GLFW_PRESS) {
      ctx.config.maxIterations = std::max(1, ctx.config.maxIterations - 1);
    }
    if (glfwGetKey(window, GLFW_KEY_EQUAL) == GLFW_PRESS) {
      ctx.config.maxIterations = std::min(200, ctx.config.maxIterations + 1);
    }
    if (glfwGetKey(window, GLFW_KEY_COMMA) == GLFW_PRESS) {
      ctx.config.substepsPerFrame = std::max(1, ctx.config.substepsPerFrame - 1);
    }
    if (glfwGetKey(window, GLFW_KEY_PERIOD) == GLFW_PRESS) {
      ctx.config.substepsPerFrame = std::min(16, ctx.config.substepsPerFrame + 1);
    }
    if (glfwGetKey(window, GLFW_KEY_SEMICOLON) == GLFW_PRESS) {
      setNumThreads(std::max(1, numThreads() - 1));
    }
    if (glfwGetKey(window, GLFW_KEY_APOSTROPHE) == GLFW_PRESS) {
      setNumThreads(numThreads() + 1);
    }
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
      // 重置：位置回到静止形，速度清零
      ctx.mesh.positions = ctx.mesh.restPositions;
      for (auto& v : ctx.mesh.velocities) v = Vec3{};
    }

    // ---- 渲染 ----
    int width = 0, height = 0;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);
    glClearColor(0.09f, 0.10f, 0.12f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    const double aspect = (height > 0) ? static_cast<double>(width) / height : 1.0;
    const double fovDeg = 45.0, nearZ = 0.01, farZ = 100.0;
    const double f = 1.0 / std::tan(fovDeg * 3.14159265358979323846 / 360.0);
    const double proj[16] = {
        f / aspect, 0, 0, 0,
        0, f, 0, 0,
        0, 0, (farZ + nearZ) / (nearZ - farZ), -1,
        0, 0, (2 * farZ * nearZ) / (nearZ - farZ), 0,
    };
    glLoadMatrixd(proj);
    setupCamera();

    drawGround();
    drawCloth(ctx, true);

    drawHud(ctx, g_state, width, height);

    glfwSwapBuffers(window);
    glfwPollEvents();

    // 截图 / 帧数上限（自动化验证用）
    const bool timeToStop = (maxFrames > 0 && frameIndex >= maxFrames);
    if (timeToStop) {
      int fw = 0, fh = 0;
      glfwGetFramebufferSize(window, &fw, &fh);
      std::vector<uint8_t> pixels(static_cast<std::size_t>(fw) * fh * 4);
      glPixelStorei(GL_PACK_ALIGNMENT, 1);
      glReadPixels(0, 0, fw, fh, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
      // OpenGL 是自下而上，PNG 需要自上而下
      std::vector<uint8_t> flipped(pixels.size());
      for (int y = 0; y < fh; ++y) {
        std::memcpy(&flipped[static_cast<std::size_t>(y) * fw * 4],
                    &pixels[static_cast<std::size_t>(fh - 1 - y) * fw * 4],
                    static_cast<std::size_t>(fw) * 4);
      }
      if (!shotPath.empty()) {
        const bool ok = writePng(shotPath.c_str(), fw, fh, flipped);
        std::printf("截图 %s: %s (%dx%d)\n", shotPath.c_str(), ok ? "成功" : "失败", fw, fh);
      }
      std::printf("[自动退出] 已跑 %d 帧，平均 %.1f FPS，物理 %.2f ms/帧，迭代 %d，应变 max %.4f\n",
                  frameIndex, g_state.fpsAvg, g_state.physicsMsAvg, ctx.iterationsUsed,
                  ctx.mesh.maxRelativeStrain());
      break;
    }

    // 每秒在控制台也打一行，便于无 GUI 环境核对
    if (now - hudTimer > 1.0) {
      hudTimer = now;
      std::printf("fps %.1f  帧 %.2f ms  物理 %.2f ms  迭代 %d  应变 max %.4f  mean %.4f  分解 %d 次\n",
                  g_state.fpsAvg, g_state.frameMsAvg, g_state.physicsMsAvg, ctx.iterationsUsed,
                  ctx.mesh.maxRelativeStrain(), ctx.mesh.meanRelativeStrain(), ctx.factorizeCount);
    }
  }

  glfwDestroyWindow(window);
  glfwTerminate();
  // 退出原因诊断：区分"跑到 --frames 限额"（timeToStop）与"窗口被关闭"。
  // 注：stdout 重定向到文件时是全缓冲，必须显式 fflush，否则日志看起来是空的。
  std::printf("[退出] 共渲染 %d 帧，窗口关闭事件=%s\n", frameIndex,
              (maxFrames > 0 && frameIndex >= maxFrames) ? "否（达到 --frames 限额）" : "是");
  std::fflush(stdout);
  return 0;
}
