#include <cfloat>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "array2_utils.h"
#include "fluidsim.h"

// #define RENDER

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
// Our state
static bool show_demo_window = true;
static bool show_another_window = false;
static ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

void display_new() {
  sim.render();
  // Start the Dear ImGui frame
  ImGui_ImplOpenGL2_NewFrame();
  ImGui_ImplGLUT_NewFrame();
  ImGui::NewFrame();
  ImGuiIO& io = ImGui::GetIO();
  // 3. Show another simple window.
  if (show_another_window) {
    ImGui::Begin("Another Window",
                 &show_another_window);  // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
    ImGui::Text("Hello from another window!");
    if (ImGui::Button("Close Me")) show_another_window = false;
    ImGui::End();
  }

  sim.renderImGuiSidebar();
  // Rendering
  ImGui::Render();
  glViewport(0, 0, (GLsizei)io.DisplaySize.x, (GLsizei)io.DisplaySize.y);

  /*glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
  glClear(GL_COLOR_BUFFER_BIT);*/
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

  // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
  if (show_demo_window) ImGui::ShowDemoWindow(&show_demo_window);

  // 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
  {
    static float f = 0.0f;
    static int counter = 0;

    ImGui::Begin("Hello, world!");  // Create a window called "Hello, world!" and append into it.

    ImGui::Text("This is some useful text.");           // Display some text (you can use a format strings too)
    ImGui::Checkbox("Demo Window", &show_demo_window);  // Edit bools storing our window open/close state
    ImGui::Checkbox("Another Window", &show_another_window);

    ImGui::SliderFloat("float", &f, 0.0f, 1.0f);             // Edit 1 float using a slider from 0.0f to 1.0f
    ImGui::ColorEdit3("clear color", (float*)&clear_color);  // Edit 3 floats representing a color

    if (ImGui::Button("Button"))  // Buttons return true when clicked (most widgets return true when edited/activated)
      counter++;
    ImGui::SameLine();
    ImGui::Text("counter = %d", counter);

    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
    ImGui::End();
  }

  // 3. Show another simple window.
  if (show_another_window) {
    ImGui::Begin("Another Window",
                 &show_another_window);  // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
    ImGui::Text("Hello from another window!");
    if (ImGui::Button("Close Me")) show_another_window = false;
    ImGui::End();
  }

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

  //Fluid2d::ParticalSystem ps;
  //ps.SetContainerSize(glm::vec2(-1.0, -1.0), glm::vec2(2.0, 2.0));
  //ps.AddFluidBlock(glm::vec2(-0.2, -0.2), glm::vec2(0.4, 0.4), glm::vec2(-2.0f, -10.0f), 0.01f * 0.7f);
  //std::cout << "partical num = " << ps.mPositions.size() << std::endl;

  // Set up the simulation
  sim.initialize(o0, grid_width, grid_resolution, grid_resolution, 1.0);
  sim.set_root_boundary(std::move(FluidSim::Boundary(c0, Vector2s(rad0, 0.0), FluidSim::BT_CIRCLE, true)));
  sim.update_boundary();
  sim.init_random_particles_2();

  while (!renderer.ShouldClose()) {
    renderer.PollEvents(); // ���������¼�������̡���꣩

    // ��ʼ ImGui ֡
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    renderer.ProcessInput();  // ���������¼�

    scalar max_timestep = std::min(step_limit, sim.compute_cfl() * cfl_number);
    int num_substeps = static_cast<int>(std::ceil(frame_time / max_timestep));
    scalar timestep = frame_time / static_cast<scalar>(num_substeps);
    for (int i = 0; i < num_substeps; i++) {  // ���
      sim.advance(timestep);
    }

    // ��ȡ����λ�ú��ܶ���Ϣ��ת��ΪLoadVertexes_new����Ĳ�������  
    std::vector<glm::vec2> particlePositions;  
    std::vector<float> particleDensities;  
    // ��ȡ�����б�  
    const auto& particles = sim.get_particles();
    // �������Ӳ���ȡλ�ú��ܶ�  
    for (const auto& particle : particles) {  
        if (particle.type_ == Particle::PT_AIR) continue;  // �����������ӵ���Ⱦ
        // glm::vec2 position(particle.x_[0], particle.x_[1]);  
        glm::vec2 position(particle.x_[0] * 0.02f - 1.0f, particle.x_[1] * 0.02f - 1.0f);  
        float density = static_cast<float>(particle.dens_);  

        particlePositions.push_back(position);  
        particleDensities.push_back(density);
        /*if (particle.x_[0] != 0) {
          std::cout << "particle x = " << position[0] << ", y = " << position[1] << std::endl;
        }*/
    }
    // ����LoadVertexes_new����  
    renderer.LoadVertexes_new(particlePositions, particleDensities);

    renderer.Update();
    // 4. ��Ⱦ ImGui
    renderer.DrawImGuiSidebar();  // ����¼ ImGui ��������

  }
#else
  Gluvi::init("Basic Fluid Solver with Static Variational Boundaries", &argc, argv);
  // Create GLUT window
  //glutInit(&argc, argv);
#ifdef __FREEGLUT_EXT_H__
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
#endif
  //glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_MULTISAMPLE);
  //glutInitWindowSize(1280, 720);
  //glutCreateWindow("Dear ImGui GLUT+OpenGL2 Example");

  // Setup GLUT display function
  // We will also call ImGui_ImplGLUT_InstallFuncs() to get all the other functions installed for us,
  // otherwise it is possible to install our own functions and call the imgui_impl_glut.h functions ourselves.
  //glutDisplayFunc(MainLoopStep);
  //glutDisplayFunc(display_new);
  // Setup Dear ImGui context 
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  (void)io;
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls
  
  // Setup Dear ImGui style
  ImGui::StyleColorsDark();
  // ImGui::StyleColorsLight();

  // Setup Platform/Renderer backends
  ImGui_ImplGLUT_Init(); // �ƺ���ԭ��Ⱦ���������ͻ��
  ImGui_ImplOpenGL2_Init();

  // Install GLUT handlers (glutReshapeFunc(), glutMotionFunc(), glutPassiveMotionFunc(), glutMouseFunc(), glutKeyboardFunc() etc.)
  // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
  // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
  // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
  // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
  ImGui_ImplGLUT_InstallFuncs();
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
  sim.init_random_particles_2();  // ������Һ��ϵĳ�ʼ���ӷֲ�

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
  /*void display(void) {
  // ����
  // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  sim.render();
  // ��ʼ ImGui ֡
  //ImGui_ImplGLUT_NewFrame();  // ������ImGui::NewFrameǰ����
  //ImGui_ImplOpenGL2_NewFrame();
  ImGui_ImplOpenGL2_NewFrame();
  ImGui::NewFrame();
  sim.renderImGuiSidebar();
  ImGui::Render();
  ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
  // glutSwapBuffers();  // ֻ�����ｻ������
}*/ 

#ifdef RENDER
#else
void specialKeys(int key, int x, int y) {
  if (key == GLUT_KEY_F1) {
    //if (true) {
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

  for (int i = 0; i < num_substeps; ++i) sim.advance(timestep); //ִ��ģ��num_substeps�κ���Ⱦһ֡

  glutSpecialFunc(specialKeys);
}
#endif
