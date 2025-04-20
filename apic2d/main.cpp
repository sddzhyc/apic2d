// Copyright 2016 Raymond Yun Fei, Christopher Batty, Robert Bridson
//
// Licensed under the Apache License,
// Version 2.0(the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <cfloat>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "array2_utils.h"
#include "fluidsim.h"

//#define RENDER

#ifdef RENDER
#include "RenderWidget.h"
#else
    #include "gluvi.h"
    #include "openglutils.h"
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
  sim.init_random_particles();

  while (!renderer.ShouldClose()) {
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
        glm::vec2 position(particle.x_[0] * 0.02f - 1.0f, particle.x_[1] * 0.02f - 1.0f);  
        float density = static_cast<float>(particle.dens_);  

        particlePositions.push_back(position);  
        particleDensities.push_back(density);
        /*if (particle.x_[0] != 0) {
          std::cout << "particle x = " << position[0] << ", y = " << position[1] << std::endl;
        }*/
    }
    // 调用LoadVertexes_new方法  
    renderer.LoadVertexes_new(particlePositions, particleDensities);

    renderer.Update();
    renderer.PollEvents();
  }
#else
  // Setup viewer stuff
  Gluvi::init("Basic Fluid Solver with Static Variational Boundaries", &argc, argv);
  Gluvi::camera = &cam;
  Gluvi::userDisplayFunc = display;
  glClearColor(1, 1, 1, 1);

  // Set up the simulation
  sim.initialize(o0, grid_width, grid_resolution, grid_resolution, 1.0);
  sim.set_root_boundary(std::move(FluidSim::Boundary(c0, Vector2s(rad0, 0.0), FluidSim::BT_CIRCLE, true)));
  sim.update_boundary();
  sim.init_random_particles();

  Gluvi::userMouseFunc = MouseFunc;
  Gluvi::userDragFunc = DragFunc;

  glutTimerFunc(1000, timer, 0);
  Gluvi::run();
#endif

  return 0;
}

void display(void) { sim.render(); }

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

  for (int i = 0; i < num_substeps; ++i) sim.advance(timestep); //执行模拟num_substeps次后渲染一帧

  glutSpecialFunc(specialKeys);
}
#endif
