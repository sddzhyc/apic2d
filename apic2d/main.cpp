#include <cfloat>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "array2_utils.h"
#include "fluidsim.h"

// #define RENDER // Uncomment this line to use the RenderWidget for rendering

#ifdef RENDER
#include "RenderWidget.h"
#else
    #include "gluvi.h"
    #include "openglutils.h"

    #include "backends/imgui_impl_glut.h"
    #include "backends/imgui_impl_opengl2.h"
    #include "imgui.h"
#endif

using namespace std;

// Try changing the grid resolution
int grid_resolution = 100;
scalar grid_width = 100.0;

// We fix the refresh rate, and find the maximal available time step
scalar cfl_number = 3.0;
scalar refresh_rate = 60.0;
scalar frame_time = 1.0 / refresh_rate;
scalar step_limit = 0.01;

FluidSim sim;

#ifdef RENDER
#else
// Gluvi stuff
//-------------
Gluvi::PanZoom2D cam(0.0, 0.0, 1.2);
void display();
void timer(int junk);

#endif
using Clock = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;

TimePoint old_mouse_time;
Vector2s old_mouse = Vector2s::Zero();

// Boundary definition - several circles in a circular domain.

Vector2s c0(50, 50), c1(70, 50), c2(30, 35), c3(50, 70);
Vector2s s0(10, 5);
scalar rad0 = 40;
Vector2s o0(0.0, 0.0);
scalar brush_radius = 2.5;

#ifdef RENDER
#else

float GetDPIScale() {
  HDC screen = GetDC(0);
  float dpi = GetDeviceCaps(screen, LOGPIXELSX);
  ReleaseDC(0, screen);
  return dpi / 96.0f;  // 96为100%缩放
}

void MouseFunc(int button, int state, int x, int y) {
  cam.transform_mouse(x, y, old_mouse.data());
  old_mouse_time = Clock::now();
}

void DragFunc(int x, int y) {
  Vector2s new_mouse;
  cam.transform_mouse(x, y, new_mouse.data());
  TimePoint new_mouse_time = Clock::now();
  scalar duration = static_cast<scalar>(std::chrono::duration_cast<std::chrono::nanoseconds>(new_mouse_time - old_mouse_time).count()) * 1e-9;
  if (duration < 0.001) return;
  Vector2s vel = (new_mouse - old_mouse) * grid_width / duration;
  Vector2s center = new_mouse * grid_width + sim.get_origin();
  sim.paint_velocity(center, brush_radius, vel);
  old_mouse = new_mouse;
  old_mouse_time = new_mouse_time;
}

// Forward declarations of helper functions
void MainLoopStep();

static ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

void display_new() {
  sim.render();
  // Start the Dear ImGui frame
  ImGui_ImplOpenGL2_NewFrame();
  ImGui_ImplGLUT_NewFrame();
  ImGui::NewFrame();
  ImGuiIO& io = ImGui::GetIO();

  sim.renderImGuiStatusBar();
  sim.renderImGuiSidebar();
  // Rendering
  ImGui::Render();
  glViewport(0, 0, (GLsizei)io.DisplaySize.x, (GLsizei)io.DisplaySize.y);
  // glUseProgram(0); // You may want this if using this code in an OpenGL 3+ context where shaders may be bound, but prefer using the GL3+ code.
  ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
   glutSwapBuffers();
   glutPostRedisplay();
}

void MainLoopStep() {
  // Start the Dear ImGui frame
  ImGui_ImplOpenGL2_NewFrame();
  ImGui_ImplGLUT_NewFrame();
  ImGui::NewFrame();
  ImGuiIO& io = ImGui::GetIO();
  // Rendering
  ImGui::Render();
  glViewport(0, 0, (GLsizei)io.DisplaySize.x, (GLsizei)io.DisplaySize.y);
  glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
  glClear(GL_COLOR_BUFFER_BIT);
  // glUseProgram(0); // You may want this if using this code in an OpenGL 3+ context where shaders may be bound, but prefer using the GL3+ code.
  ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

  /*glutSwapBuffers();
  glutPostRedisplay();*/
}
#endif
// Main testing code
//-------------
int main(int argc, char **argv) {
#ifdef RENDER
  RenderWidget renderer;
  renderer.Init();


  // Set up the simulation
  sim.initialize(o0, grid_width, grid_resolution, grid_resolution, 1.0);
  sim.set_root_boundary(std::move(FluidSim::Boundary(c0, Vector2s(rad0, 0.0), FluidSim::BT_CIRCLE, true)));
  sim.update_boundary();
  sim.init_random_particles_2();

  while (!renderer.ShouldClose()) {
    renderer.PollEvents(); // 处理输入事件（如键盘、鼠标）

    // 开始 ImGui 帧
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    renderer.ProcessInput();  // 处理输入事件

    scalar max_timestep = std::min(step_limit, sim.compute_cfl() * cfl_number);
    int num_substeps = static_cast<int>(std::ceil(frame_time / max_timestep));
    scalar timestep = frame_time / static_cast<scalar>(num_substeps);
    for (int i = 0; i < num_substeps; i++) {  // 求解
      sim.advance(timestep);
    }

    // 提取粒子位置和密度信息并转换为LoadVertexes_new所需的参数类型  
    std::vector<glm::vec2> particlePositions;  
    std::vector<float> particleDensities;  
    // 获取粒子列表  
    const auto& particles = sim.get_particles();
    // 遍历粒子并提取位置和密度  
    for (const auto& particle : particles) {  
        if (particle.type_ == Particle::PT_AIR) continue;  // 跳过空气粒子的渲染
        // glm::vec2 position(particle.x_[0], particle.x_[1]);  
        glm::vec2 position(particle.x_[0] * 0.02f - 1.0f, particle.x_[1] * 0.02f - 1.0f);  
        float density = static_cast<float>(particle.dens_);  

        particlePositions.push_back(position);  
        particleDensities.push_back(density);
    }
    // 调用LoadVertexes_new方法  
    renderer.LoadVertexes_new(particlePositions, particleDensities);

    renderer.Update();
    // 4. 渲染 ImGui
    renderer.DrawImGuiSidebar();  // 仅记录 ImGui 绘制命令

  }
#else
  // Setup GLUT display function
  Gluvi::init("Basic Fluid Solver with Static Variational Boundaries", &argc, argv);
#ifdef __FREEGLUT_EXT_H__
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
#endif
  // Setup Dear ImGui context 
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  (void)io;
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls
  
  // Setup Dear ImGui style
  ImGui::StyleColorsDark();
  // ImGui::StyleColorsLight();
  // 
  // Setup Platform/Renderer backends
  ImGui_ImplGLUT_Init();
  ImGui_ImplOpenGL2_Init();
  // Install GLUT handlers (glutReshapeFunc(), glutMotionFunc(), glutPassiveMotionFunc(), glutMouseFunc(), glutKeyboardFunc() etc.)
  // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
  // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
  // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
  // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
  ImGui_ImplGLUT_InstallFuncs();
  // Load Fonts
  // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
  // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
  // - If the file cannot be loaded, the function will return a nullptr. Please handle those errors in your application (e.g. use an assertion, or display an
  // error and quit).
  // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which
  // ImGui_ImplXXXX_NewFrame below will call.
  // - Use '#define IMGUI_ENABLE_FREETYPE' in your imconfig file to use Freetype for higher quality font rendering.
  // - Read 'docs/FONTS.md' for more instructions and details.
  // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
  io.Fonts->Clear();
  float dpi_scale = GetDPIScale();
  ImFont* font = io.Fonts->AddFontFromFileTTF("C:\\Windows\\Fonts\\msyh.ttc",  // 字体路径
                                              18.0f,                           // 字号
                                              nullptr,                         // 字体配置
                                              io.Fonts->GetGlyphRangesChineseFull()  // 字形范围
  );
  io.FontGlobalScale = 1.0f;  // 字体已放大，不再需要全局缩放
  ImGui::GetStyle().ScaleAllSizes(dpi_scale);
  IM_ASSERT(font != nullptr);
  // Setup viewer stuff
  Gluvi::camera = &cam;
  Gluvi::userDisplayFunc = display_new;
  glClearColor(1, 1, 1, 1);
  scalar rho_liquid = 1.0f;
  //scalar rho_air = 0.001f;
  // Set up the simulation
  sim.initialize(o0, grid_width, grid_resolution, grid_resolution, rho_liquid);
  sim.set_root_boundary(std::move(FluidSim::Boundary(c0, Vector2s(rad0, 0.0), FluidSim::BT_CIRCLE, true)));
  sim.update_boundary();
  //sim.init_random_particles();
  sim.init_random_particles_2();  // 采用气液混合的初始粒子分布

  Gluvi::userMouseFunc = MouseFunc;
  Gluvi::userDragFunc = DragFunc;

  glutTimerFunc(1000, timer, 0);
  // Main loop
  Gluvi::run(); //glutMainLoop();
  //  Cleanup
  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGLUT_Shutdown();
  ImGui::DestroyContext();

#endif

  return 0;
}

void display(void) { sim.render(); }

#ifdef RENDER
#else
void specialKeys(int key, int x, int y) {
  if (key == GLUT_KEY_F1) {
    sim.OutputPointDataBgeo("D:/FluidSimulator/new_apic2d/point/pointData", sim.outframe_);
    sim.OutputGridXDataBgeo("D:/FluidSimulator/new_apic2d/gridX/gridXData", sim.outframe_);
    sim.OutputGridYDataBgeo("D:/FluidSimulator/new_apic2d/gridY/gridYData", sim.outframe_);
    sim.OutputGridDataBgeo("D:/FluidSimulator/new_apic2d/grid/gridData", sim.outframe_++);
  }
}

void timer(int junk) {
  glutPostRedisplay();
  glutTimerFunc(static_cast<int>(ceil(frame_time * 1000.0)), timer, 0);
  scalar max_timestep = std::min(step_limit, sim.compute_cfl() * cfl_number);
  int num_substeps = static_cast<int>(std::ceil(frame_time / max_timestep));
  scalar timestep = frame_time / static_cast<scalar>(num_substeps);
  // std::cout << "max_timestep = " << max_timestep << ", num_substeps = " << num_substeps << ", timestep = " << timestep << std::endl;
  for (int i = 0; i < num_substeps; ++i) sim.advance(timestep); //执行模拟num_substeps次后渲染一帧

  if (sim.export_enabled_ && sim.is_paused == false) {
    sim.OutputPointDataBgeo(sim.savePath_ + "/point", sim.outframe_);
    sim.OutputGridXDataBgeo(sim.savePath_ + "/gridX", sim.outframe_);
    sim.OutputGridYDataBgeo(sim.savePath_ + "/gridY", sim.outframe_);
    sim.OutputGridDataBgeo(sim.savePath_ + "/grid", sim.outframe_++);
  }
  glutSpecialFunc(specialKeys);
}
#endif
