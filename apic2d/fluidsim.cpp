
#define _INIT_TEMP 370.0f

#ifdef WIN32
#define NOMINMAX
#endif

#include "fluidsim.h"

#include "array2_utils.h"
#include "kernel.h"
#include "sorter.h"
#include "sparse_matrix.h"

#ifdef __APPLE__
#include <GLUT/glut.h>  // why does Apple have to put glut.h here...
#else
#include <GL/glut.h>  // ...when everyone else puts it here?
#endif

#include "openglutils.h"

#include <fstream>
// Change here to try different integration scheme, options:
// IT_PIC: original particle-in-cell (PIC)
// IT_FLIP: original fluid-implicit-particle (FLIP)
// IT_RPIC: rotational particle-in-cell (RPIC)
// IT_APIC: affine particle-in-cell (APIC)
// IT_AFLIP: affine fluid-implicit-particle (AFLIP)
// IT_ASFLIP: affine separable fluid-implicit-particle (ASFLIP)
const FluidSim::INTEGRATOR_TYPE integration_scheme = FluidSim::IT_FLIP;

// Change here to try different order for velocity evaluation from grid,
// options:
// VO_EULER: first order evaluation
// VO_RA2: Ralston's second order
//         evaluation
// VO_RK3: Runge Kutta's 3rd-order method
// VO_RK4: Runge Kutta's
//         4rd-order method
const FluidSim::VELOCITY_ORDER velocity_order = FluidSim::VO_EULER;

// Change here to try different order for interpolation
// options:
// IO_LINEAR: linear interpolation
// IO_QUADRATIC: quadratic interpolation
//TODO: IO_QUADRATIC的气体粒子g2p插值功能尚未实现，因为get_air_velocity_quadratic待实现
const FluidSim::INTERPOLATION_ORDER interpolation_order = FluidSim::IO_LINEAR;

const scalar lagrangian_ratio = 0.97f;
// const scalar lagrangian_ratio = 1.0f;
const int particle_correction_step = 1;

// whether collision is resolved with another pass on each particle (instead of
// grid)
const bool use_particle_boundary_collisions = false;

scalar fraction_inside(scalar phi_left, scalar phi_right);
void extrapolate(Array2s& grid, Array2s& old_grid, const Array2s& grid_weight, const Array2s& grid_liquid_weight, Array2c& valid_, Array2c old_valid_,
                 const Vector2i& offset, int num_layers);

FluidSim::~FluidSim() { delete m_sorter_; }

void FluidSim::initialize(const Vector2s& origin, scalar width, int ni, int nj, scalar rho, bool draw_grid, bool draw_particles, bool draw_velocities,
                          bool draw_boundaries, bool print_timing) {
  rho_air_ = 0.001f;
  rho_ = rho;
  draw_grid_ = draw_grid;
  draw_particles_ = draw_particles;
  draw_velocities_ = draw_velocities;
  draw_boundaries_ = draw_boundaries;
  print_timing_ = print_timing;
  origin_ = origin;
  ni_ = ni;
  nj_ = nj;
  dx_ = width / static_cast<scalar>(ni_);
  outframe_ = 0;
  u_.resize(ni_ + 1, nj_);
  temp_u_.resize(ni_ + 1, nj_);
  u_weights_.resize(ni_ + 1, nj_);
  u_valid_.resize(ni_ + 1, nj_);
  saved_u_.resize(ni_ + 1, nj_);
  v_.resize(ni_, nj_ + 1);
  temp_v_.resize(ni_, nj_ + 1);
  v_weights_.resize(ni_, nj_ + 1);
  v_valid_.resize(ni_, nj_ + 1);
  saved_v_.resize(ni_, nj_ + 1);
  u_.set_zero();
  v_.set_zero();

  u_a_.resize(ni_ + 1, nj_);
  v_a_.resize(ni_, nj_ + 1);
  u_a_.set_zero();
  v_a_.set_zero();

  saved_u_a_.resize(ni_ + 1, nj_);
  saved_v_a_.resize(ni_, nj_ + 1);
  temp_u_a_.resize(ni_ + 1, nj_);
  temp_v_a_.resize(ni_, nj_ + 1);

  nodal_solid_phi_.resize(ni_ + 1, nj_ + 1);
  valid_.resize(ni_ + 1, nj_ + 1);
  old_valid_.resize(ni_ + 1, nj_ + 1);
  valid_air_.resize(ni_ + 1, nj_ + 1);
  old_valid_air_.resize(ni_ + 1, nj_ + 1);
  liquid_phi_.resize(ni_, nj_);
  air_phi_.resize(ni_, nj_);
  merged_phi_.resize(ni_, nj_);
  m_sorter_ = new sorter(ni_, nj_);

  // compressible fluid
  comp_rho_.resize(ni_, nj_);
  saved_comp_rho_.resize(ni_, nj_);
  comp_pressure_.resize(ni_, nj_); // 需要初始化其大小

  grid_temp_.resize(ni_, nj_);
  grid_rho_.resize(ni_, nj_);

  // debug arrays
  laplacianP_.resize(ni_,nj_);

  // hardcode parameters
  R = 461.5;
  b = 949.7;
  a = 2977.4;  // at room temperature
  T_ = _INIT_TEMP;  // 初始温度27度

  // 在0.9到1.1区间上对coef_B函数绘制曲线
  /*for (int i = 0; i < 0.2 / 0.00005; i++) {
    scalar x = 0.9 + i * 0.00005;
    std::cout << x << ", " << compute_coef_A(x) << ", " << compute_coef_B(x) << ", " << get_pressure(x) << std::endl;
  }*/
}

/*!
  \brief  Initialize the grid-based signed distance field that dictates the position of
          the solid boundary
*/
void FluidSim::update_boundary() {
  parallel_for(0, nj_ + 1, [this](int j) {
    for (int i = 0; i < ni_ + 1; ++i) {
      Vector2s pos(i * dx_, j * dx_);
      nodal_solid_phi_(i, j) = solid_distance(pos + origin_);
    }
  });
}

void FluidSim::paint_velocity(const Vector2s& brush_center, const scalar brush_radius, const Vector2s& vel) {
  velocity_brush_data_.emplace_back(brush_center, vel, brush_radius);
}

void FluidSim::relaxation(scalar dt) {
  int np = static_cast<int>(particles_.size());

  const scalar re = dx_ / sqrt(2.0) * 1.1;

  int offset = rand() % particle_correction_step;
  // Compute Pseudo Moved Point
  parallel_for(0, np, [&](int n) {
    if (n % particle_correction_step != offset) return;

    Particle& p = particles_[n];

    Vector2s spring = Vector2s::Zero();

    int ix = std::max(0, std::min(static_cast<int>((p.x_(0) - origin_(0)) / dx_), ni_));
    int iy = std::max(0, std::min(static_cast<int>((p.x_(1) - origin_(1)) / dx_), nj_));

    m_sorter_->getNeigboringParticles_cell(ix, iy, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
      for (const Particle* np : neighbors) {
        if (&p != np) {
          scalar dist = (p.x_ - np->x_).norm();
          scalar w = 50.0 * kernel::smooth_kernel(dist * dist, re);
          if (dist > 0.01 * re) {
            spring += w * (p.x_ - np->x_) / dist * re;
          } else {
            spring(0) += 0.01 * re / dt * (rand() & 0xFF) / 255.0;
            spring(1) += 0.01 * re / dt * (rand() & 0xFF) / 255.0;
          }
        }
      }
    });

    p.buf0_ = p.x_ + dt * spring;

    Vector2s pp = (p.buf0_ - origin_) / dx_;
    scalar phi_value = interpolate_value(pp, nodal_solid_phi_);

    if (phi_value < 0) {
      Vector2s normal;
      interpolate_gradient(normal, pp, nodal_solid_phi_);
      normal.normalize();
      p.buf0_ -= phi_value * normal;
    }
  });

  // Update
  parallel_for(0, np, [&](int n) {
    if (n % particle_correction_step != offset) return;
    Particle& p = particles_[n];

    p.x_ = p.buf0_;
  });
}


void FluidSim::advance(scalar dt) {
  // Passively advect particles_
  tick();
  m_sorter_->sort(particles_, origin_, dx_);
  tock("sort");

  tick();
#ifdef COMPRESSIBLE_FLUID
  map_p2g_compressible();
#else
  map_p2g();
#endif
  tock("p2g");
  save_velocity(); //更新前保存原网格速度

  tick();
  add_force(dt);
  tock("add force");  // 加重力

  tick();
  // compute_liquid_distance();
  // compute_air_distance();
  compute_merged_distance();  // 联合计算liquid_phi_, air_phi_, merged_phi_
  tock("compute phi");

  // Compute finite-volume type_ face area weight for each velocity sample.
  tick();
  compute_weights();
  tock("compute weights");  // 求level set

  // Set up and solve the variational pressure solve.
  tick();
#ifdef COMPRESSIBLE_FLUID
  solve_compressible_density_new(dt);  // 启用不可压缩求解，更新comp_rho_
#else 
  solve_pressure_with_air(dt);
#endif
  solve_temperature(dt);  // 计算温度场
  tock("solve pressure");

  // Pressure projection only produces valid velocities in faces with non-zero
  // associated face area. Because the advection step may interpolate from these
  // invalid faces, we must extrapolate velocities from the fluid domain into
  // these zero-area faces.
  tick();
  extrapolate(u_, temp_u_, u_weights_, liquid_phi_, valid_, old_valid_, Vector2i(-1, 0), 2);
  extrapolate(v_, temp_v_, v_weights_, liquid_phi_, valid_, old_valid_, Vector2i(0, -1), 2);

  extrapolate(u_a_, temp_u_a_, u_weights_, air_phi_, valid_air_, old_valid_air_, Vector2i(-1, 0), 2);
  extrapolate(v_a_, temp_v_a_, v_weights_, air_phi_, valid_air_, old_valid_air_, Vector2i(0, -1), 2);
  tock("extrapolate");

  // For extrapolated velocities, replace the normal component with
  // that of the object.
  tick();
  constrain_velocity();
  tock("constrain velocity");

  tick();
  relaxation(dt);
  tock("relaxation");

  tick();
  switch (integration_scheme) {
    case IT_PIC:
      // The original PIC, which is a specific case of a more general (Affine)
      // FLIP scheme
      map_g2p_flip_general(dt, 0.0, 0.0, 0.0, 0.0);  // TODO:此处调用了interpolation 函数来获取粒子位置的密度插值，但是使用的却是更新前的密度，why？
      break;

    case IT_FLIP:
      // FLIP scheme
      map_g2p_flip_general(dt, lagrangian_ratio, 0.0, 0.0, 0.0);
      break;

    case IT_RPIC:
      // RPIC is a rotation-augmented PIC scheme
      map_g2p_flip_general(dt, 0.0, 0.0, 0.0, 1.0);
      break;

    case IT_APIC:
      // APIC is an affine-augmented PIC scheme
      map_g2p_flip_general(dt, 0.0, 0.0, 1.0, 1.0);
      break;

    case IT_AFLIP:
      // AFLIP is an affine-augmented FLIP scheme
      map_g2p_flip_general(dt, lagrangian_ratio, 0.0, 1.0, 1.0);
      break;

    case IT_ASFLIP:
      // ASFLIP is an affine-augmented PIC scheme with easier particle
      // separation
      map_g2p_flip_general(dt, lagrangian_ratio, 1.0, 1.0, 1.0);
      break;

    default:
      std::cerr << "Unknown integrator type_!" << std::endl;
      break;
  }
  tock("g2p");

  if (use_particle_boundary_collisions) {
    tick();
    particle_boundary_collision(dt);
    tock("particle collision");
  }
}

void FluidSim::save_velocity() {
  if (lagrangian_ratio > 0.0) {
    saved_u_ = u_;
    saved_v_ = v_;
    saved_comp_rho_ = comp_rho_;

    saved_u_a_ = u_a_;
    saved_v_a_ = v_a_;
  }
}

void FluidSim::add_force(scalar dt) {
  // gravity
  for (int j = 0; j < nj_ + 1; ++j) {
    for (int i = 0; i < ni_; ++i) {
      v_(i, j) += -9.810 * dt; //减少重力加速度？
      // TODO:是否需要为气体粒子设置不同的重力加速度？
      v_a_(i, j) += -9.810 * dt;  
    }
  }

  // user-specified velocity
  for (const BrushFootprint& p : velocity_brush_data_) {
    scalar radius_in_grid = p.radius_ / dx_;
    Vector2s center_in_grid_u = (p.center_ - origin_) / dx_ - Vector2s(0.0, 0.5);
    Vector2i low_range_u = Vector2i(std::max(0, std::min(u_.ni - 1, static_cast<int>(floor(center_in_grid_u(0) - radius_in_grid)))),
                                    std::max(0, std::min(u_.nj - 1, static_cast<int>(floor(center_in_grid_u(1) - radius_in_grid)))));
    Vector2i high_range_u = Vector2i(std::max(0, std::min(u_.ni - 1, static_cast<int>(ceil(center_in_grid_u(0) + radius_in_grid)))),
                                     std::max(0, std::min(u_.nj - 1, static_cast<int>(ceil(center_in_grid_u(1) + radius_in_grid)))));
    for (int i = low_range_u(0); i <= high_range_u(0); ++i) {
      for (int j = low_range_u(1); j <= high_range_u(1); ++j) {
        Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
        if ((pos - p.center_).norm() <= p.radius_) {
          u_(i, j) = p.vel_(0);
        }
      }
    }

    Vector2s center_in_grid_v = (p.center_ - origin_) / dx_ - Vector2s(0.5, 0.0);
    Vector2i low_range_v = Vector2i(std::max(0, std::min(v_.ni - 1, static_cast<int>(floor(center_in_grid_v(0) - radius_in_grid)))),
                                    std::max(0, std::min(v_.nj - 1, static_cast<int>(floor(center_in_grid_v(1) - radius_in_grid)))));
    Vector2i high_range_v = Vector2i(std::max(0, std::min(v_.ni - 1, static_cast<int>(ceil(center_in_grid_v(0) + radius_in_grid)))),
                                     std::max(0, std::min(v_.nj - 1, static_cast<int>(ceil(center_in_grid_v(1) + radius_in_grid)))));
    for (int i = low_range_v(0); i <= high_range_v(0); ++i) {
      for (int j = low_range_v(1); j <= high_range_v(1); ++j) {
        Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
        if ((pos - p.center_).norm() <= p.radius_) {
          v_(i, j) = p.vel_(1);
        }
      }
    }
  }

  velocity_brush_data_.clear();
}

scalar FluidSim::compute_cfl() const {
  scalar max_vel_squared = 0.0;
  for (const Particle& p : particles_) {
    max_vel_squared = std::max(max_vel_squared, p.v_.squaredNorm());
  }
  return dx_ / sqrt(max_vel_squared);
}

/*!
  \brief  For extrapolated points, replace the normal component
          of velocity with the object velocity (in this case zero).
*/
void FluidSim::constrain_velocity() {
  // At lower grid resolutions, the normal estimate from the signed
  // distance function is poor, so it doesn't work quite as well.
  // An exact normal would do better.

  // constrain u
  parallel_for(0, std::max(u_.nj, v_.nj), [this](int j) {
    if (j < u_.nj) {
      for (int i = 0; i < u_.ni; ++i) {
        if (u_weights_(i, j) == 0) {
          // apply constraint
          Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
          Vector2s vel = get_velocity(pos);
          Vector2s normal(0, 0);
          interpolate_gradient(normal, Vector2s(i, j + 0.5), nodal_solid_phi_);
          normal.normalize();
          scalar perp_component = vel.dot(normal);
          temp_u_(i, j) = vel[0] - perp_component * normal[0];
        }
      }
    }

    if (j < v_.nj) {
      for (int i = 0; i < v_.ni; ++i) {
        if (v_weights_(i, j) == 0) {
          // apply constraint
          Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
          Vector2s vel = get_velocity(pos);
          Vector2s normal(0, 0);
          interpolate_gradient(normal, Vector2s(i + 0.5, j), nodal_solid_phi_);
          normal.normalize();
          scalar perp_component = vel.dot(normal);
          temp_v_(i, j) = vel[1] - perp_component * normal[1];
        }
      }
    }
  });

  parallel_for(0, std::max(u_.nj, v_.nj), [this](int j) {
    if (j < u_.nj) {
      for (int i = 0; i < u_.ni; ++i) {
        if (u_weights_(i, j) == 0) {
          u_(i, j) = temp_u_(i, j);
        }
      }
    }
    if (j < v_.nj) {
      for (int i = 0; i < v_.ni; ++i) {
        if (v_weights_(i, j) == 0) {
          v_(i, j) = temp_v_(i, j);
        }
      }    
    }
  });
  // constrain u_a_
  parallel_for(0, std::max(u_a_.nj, v_a_.nj), [this](int j) {
    if (j < u_a_.nj) {
      for (int i = 0; i < u_a_.ni; ++i) {
        if (u_weights_(i, j) == 0) {
          // apply constraint
          Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
          Vector2s vel = get_air_velocity(pos);
          Vector2s normal(0, 0);
          interpolate_gradient(normal, Vector2s(i, j + 0.5), nodal_solid_phi_);
          normal.normalize();
          scalar perp_component = vel.dot(normal);
          temp_u_a_(i, j) = vel[0] - perp_component * normal[0];
        }
      }
    }
    if (j < v_a_.nj) {
      for (int i = 0; i < v_a_.ni; ++i) {
        if (v_weights_(i, j) == 0) {
          // apply constraint
          Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
          Vector2s vel = get_air_velocity(pos);
          Vector2s normal(0, 0);
          interpolate_gradient(normal, Vector2s(i + 0.5, j), nodal_solid_phi_);
          normal.normalize();
          scalar perp_component = vel.dot(normal);
          temp_v_a_(i, j) = vel[1] - perp_component * normal[1];
        }
      }
    }
  });
  parallel_for(0, std::max(u_a_.nj, v_a_.nj), [this](int j) {
    if (j < u_a_.nj) {
      for (int i = 0; i < u_a_.ni; ++i) {
        if (u_weights_(i, j) == 0) {
          u_a_(i, j) = temp_u_a_(i, j);
        }
      }
    }
    if (j < v_a_.nj) {
      for (int i = 0; i < v_a_.ni; ++i) {
        if (v_weights_(i, j) == 0) {
          v_a_(i, j) = temp_v_a_(i, j);
        }
      }
    }
  });
}

void FluidSim::compute_liquid_distance() {
  const scalar min_radius = dx_ / sqrtf(2.0);
  parallel_for(0, static_cast<int>(nj_), [&](int j) {
    for (int i = 0; i < ni_; ++i) {
      Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_; //网格中心位置
      // Estimate from particles
      scalar min_liquid_phi = dx_;
      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
        for (const Particle* p : neighbors) {
          if (p->type_ == Particle::PT_AIR) continue;  // 只考虑液体粒子
          scalar phi_temp = (pos - p->x_).norm() - std::max(p->radii_, min_radius);
          min_liquid_phi = std::min(min_liquid_phi, phi_temp); //网格中心到最近邻粒子边缘的距离
        }
      });

      // "extrapolate" phi into solids if nearby
      scalar solid_phi_val = solid_distance(pos);
      liquid_phi_(i, j) = std::min(min_liquid_phi, solid_phi_val);
    }
  });
}

void FluidSim::compute_air_distance() {
  const scalar min_radius = dx_ / sqrtf(2.0);
  parallel_for(0, static_cast<int>(nj_), [&](int j) {
    for (int i = 0; i < ni_; ++i) {
      Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_;  // 网格中心位置
      // Estimate from particles
      scalar min_air_phi = dx_;  // 初始化air_phi的初始值，也是没有检索到粒子时的最大值
      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
        for (const Particle* p : neighbors) {
          if (p->type_ == Particle::PT_LIQUID) continue;
          scalar phi_temp = (pos - p->x_).norm() - std::max(p->radii_, min_radius);
          min_air_phi = std::min(min_air_phi, phi_temp);  // 网格中心到最近邻粒子边缘的距离
        }
      });
      // "extrapolate" phi into solids if nearby
      scalar solid_phi_val = solid_distance(pos);
      air_phi_(i, j) = std::min(min_air_phi, solid_phi_val);
    }
  });
}


void FluidSim::compute_merged_distance() {
  const scalar min_radius = dx_ / sqrtf(2.0);
  parallel_for(0, static_cast<int>(nj_), [&](int j) {
    for (int i = 0; i < ni_; ++i) {
      Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_;  // 网格中心位置
      // Estimate from particles 
      scalar min_liquid_phi = 4 * dx_;// 初始化phi的初始值，也就是没有检索到粒子时的最大值
      scalar min_air_phi = 4 * dx_;
      m_sorter_->getNeigboringParticles_cell(i, j, -4, 4, -4, 4, [&](const NeighborParticlesType& neighbors) {
          for (const Particle* p : neighbors) {
            scalar phi_temp = (pos - p->x_).norm() - std::max(p->radii_, min_radius);
            if (p->type_ == Particle::PT_LIQUID) {
              min_liquid_phi = std::min(min_liquid_phi, phi_temp); 
            } else {
            min_air_phi = std::min(min_air_phi, phi_temp);  // 网格中心到最近邻粒子边缘的距离
            }
          }
      });
      scalar min_merged_phi = 4 * dx_; 
      // if (min_liquid_phi != 2 * dx_ || min_merged_phi != 2 * dx_) { // 将真空区域的phi值设置为2dx后，固液界面似乎变得更不稳定？
        min_merged_phi = (min_liquid_phi - min_air_phi) * 0.5f;
      //}
      // "extrapolate" phi into solids if nearby
      scalar solid_phi_val = solid_distance(pos);
      merged_phi_(i, j) = std::min(min_merged_phi, solid_phi_val);
      liquid_phi_(i, j) = std::min(min_liquid_phi, solid_phi_val);
      air_phi_(i, j) = std::min(min_air_phi, solid_phi_val);
    }
  });
}

// 计算曲率, 使用 standard central difference for the Laplacian
scalar FluidSim::compute_curvature(int i, int j, int is_x, int is_y) { 
  // scalar nabla_phi = merged_phi_(i + 1, j) + merged_phi_(i - 1, j) + merged_phi_(i, j + 1) + merged_phi_(i, j - 1) - 4 * merged_phi_(i, j);
  // 2D模拟中的界面曲率计算
  scalar nabla_phi = is_x * (merged_phi_(i + 1, j) + merged_phi_(i - 1, j) - 2 * merged_phi_(i, j)) +
                     is_y * (merged_phi_(i, j + 1) + merged_phi_(i, j - 1) - 2 * merged_phi_(i, j));
  nabla_phi /= (dx_ * dx_);
  return nabla_phi;
}
scalar FluidSim::compute_face_fraction(int phi_0, int phi_1) {
  float d = sqrt(dx_ * dx_ - sqr(phi_1 - phi_0));
  float vfrac = 0.5f - 0.5f * (phi_1 + phi_0) / d;
  vfrac = std::max(0.0f, std::min(1.0f, vfrac));  // clamp to [0, 1]
  return vfrac;
}

scalar FluidSim::get_merged_phi(int i, int j) { return (liquid_phi_(i, j) - air_phi_(i, j)) * 0.5f; }
    /*!
  \brief  Add a tracer particle for visualization
*/
void FluidSim::add_particle(const Particle& p) { particles_.push_back(p); }

/*!
  \brief  Particles can still occasionally leave the domain due to truncation
          errors, interpolation error, or large timesteps, so we project them
          back in for good measure.
*/
void FluidSim::particle_boundary_collision(scalar dt) {
  parallel_for(0, static_cast<int>(particles_.size()), [&](int p) {
    Vector2s pp = (particles_[p].x_ - origin_) / dx_;

    // Try commenting this section out to see the degree of accumulated error.
    scalar phi_value = interpolate_value(pp, nodal_solid_phi_);
    if (phi_value < 0) {
      Vector2s normal;
      interpolate_gradient(normal, pp, nodal_solid_phi_);
      normal.normalize();
      particles_[p].x_ -= phi_value * normal;
    }
  });

  m_sorter_->sort(particles_, origin_, dx_);

  auto remover = [&](const Particle& p) {
    return p.x_(0) < origin_(0) - 0.5 * dx_ || p.x_(0) > origin_(0) + (static_cast<scalar>(ni_) + 1.5) * dx_ || p.x_(1) < origin_(1) - 0.5 * dx_ ||
           p.x_(1) > origin_(1) + (static_cast<scalar>(nj_) + 1.5) * dx_;
  };
  particles_.erase(std::remove_if(particles_.begin(), particles_.end(), remover), particles_.end());
}

/*!
  \brief 对温度、密度等存储于网格中心的标量场插值回粒子位置
*/
scalar FluidSim::get_temperature_quadratic(const Vector2s& position, const Array2s& grid_temp) {
  Vector2s p = (position - origin_) / dx_;
  //Vector2s p0 = p - Vector2s(0.5, 0.5); 为了直观，不转换数学坐标的位置，而是对网格坐标加0.5dx的偏移
  Vector2s p0 = p;

  scalar ret = 0.0f;
  Vector2i ip = Vector2i(static_cast<int>(p0(0)), static_cast<int>(p0(1)));
  for (int i = ip(0) - 1; i <= ip(0) + 1; ++i) {
    for (int j = ip(1) - 1; j <= ip(1) + 1; ++j) {
      if (i < 0 || i > ni_ || j < 0 || j >= nj_) {
        continue;
      }
      Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_;
      scalar w = kernel::quadratic_kernel(position - pos, dx_);
      // scalar w = (velocity_order == IO_LINEAR) ? kernel::linear_kernel(position - pos, dx_) : kernel::quadratic_kernel(position - pos, dx_);
      ret += grid_temp(i, j) * w;
    }
  }

  return ret;
}

Vector2s FluidSim::get_velocity_quadratic_impl(const Vector2s& position, const Array2s& uu, const Array2s& vv) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);

  Vector2s ret = Vector2s::Zero();
  Vector2i ip0 = Vector2i(static_cast<int>(p0(0)), static_cast<int>(p0(1)));
  for (int i = ip0(0) - 1; i <= ip0(0) + 2; ++i) {
    for (int j = ip0(1) - 1; j <= ip0(1) + 2; ++j) {
      if (i < 0 || i > ni_ || j < 0 || j >= nj_) {
        continue;
      }
      Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
      scalar w = kernel::quadratic_kernel(position - pos, dx_);
      ret(0) += uu(i, j) * w;
    }
  }

  Vector2i ip1 = Vector2i(static_cast<int>(p1(0)), static_cast<int>(p1(1)));
  for (int i = ip1(0) - 1; i <= ip1(0) + 2; ++i) {
    for (int j = ip1(1) - 1; j <= ip1(1) + 2; ++j) {
      if (i < 0 || i >= ni_ || j < 0 || j > nj_) {
        continue;
      }
      Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
      scalar w = kernel::quadratic_kernel(position - pos, dx_);
      ret(1) += vv(i, j) * w;
    }
  }

  return ret;
}

Matrix2s FluidSim::get_affine_matrix_quadratic_impl(const Vector2s& position, const Array2s& uu, const Array2s& vv) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);

  Matrix2s ret = Matrix2s::Zero();
  scalar invD = 4.0 / (dx_ * dx_);
  Vector2i ip0 = Vector2i(static_cast<int>(p0(0)), static_cast<int>(p0(1)));
  for (int i = ip0(0) - 1; i <= ip0(0) + 2; ++i) {
    for (int j = ip0(1) - 1; j <= ip0(1) + 2; ++j) {
      if (i < 0 || i > ni_ || j < 0 || j >= nj_) {
        continue;
      }
      Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
      scalar w = kernel::quadratic_kernel(position - pos, dx_);
      ret.block<2, 1>(0, 0) += u_(i, j) * w * (position - pos) * invD;
    }
  }

  Vector2i ip1 = Vector2i(static_cast<int>(p1(0)), static_cast<int>(p1(1)));
  for (int i = ip1(0) - 1; i <= ip1(0) + 2; ++i) {
    for (int j = ip1(1) - 1; j <= ip1(1) + 2; ++j) {
      if (i < 0 || i >= ni_ || j < 0 || j > nj_) {
        continue;
      }
      Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
      scalar w = kernel::quadratic_kernel(position - pos, dx_);
      ret.block<2, 1>(0, 1) += v_(i, j) * w * (position - pos) * invD;
    }
  }

  return ret;
}

Vector2s FluidSim::get_velocity_quadratic(const Vector2s& position) { return get_velocity_quadratic_impl(position, u_, v_); }

Matrix2s FluidSim::get_affine_matrix_quadratic(const Vector2s& position) { return get_affine_matrix_quadratic_impl(position, u_, v_); }

Vector2s FluidSim::get_saved_velocity_quadratic(const Vector2s& position) { return get_velocity_quadratic_impl(position, saved_u_, saved_v_); }

Matrix2s FluidSim::get_saved_affine_matrix_quadratic(const Vector2s& position) { return get_affine_matrix_quadratic_impl(position, saved_u_, saved_v_); }

Vector2s FluidSim::get_velocity_and_affine_matrix_with_order(const Vector2s& position, scalar dt, FluidSim::VELOCITY_ORDER v_order,
                                                             FluidSim::INTERPOLATION_ORDER i_order, Matrix2s* affine_matrix) {
  auto get_velocity_func = (i_order == IO_LINEAR) ? &FluidSim::get_velocity : &FluidSim::get_velocity_quadratic;
  auto get_affine_matrix_func = (i_order == IO_LINEAR) ? &FluidSim::get_affine_matrix : &FluidSim::get_affine_matrix_quadratic;

  switch (v_order) {
    case FluidSim::VO_EULER:
      if (affine_matrix) {
        *affine_matrix = (this->*get_affine_matrix_func)(position);
      }
      return (this->*get_velocity_func)(position);
    case FluidSim::VO_RA2: {
      Vector2s v0 = (this->*get_velocity_func)(position);
      Vector2s v1 = (this->*get_velocity_func)(position + v0 * 2.0 / 3.0 * dt);
      if (affine_matrix) {
        Matrix2s a0 = (this->*get_affine_matrix_func)(position);
        Matrix2s a1 = (this->*get_affine_matrix_func)(position + v0 * 2.0 / 3.0 * dt);
        *affine_matrix = a0 * 0.25 + a1 * 0.75;
      }
      return v0 * 0.25 + v1 * 0.75;
    }
    case FluidSim::VO_RK3: {
      Vector2s v0 = (this->*get_velocity_func)(position);
      Vector2s v1 = (this->*get_velocity_func)(position + v0 * 0.5 * dt);
      Vector2s v2 = (this->*get_velocity_func)(position + (v1 * 2.0 - v0) * dt);
      if (affine_matrix) {
        Matrix2s a0 = (this->*get_affine_matrix_func)(position);
        Matrix2s a1 = (this->*get_affine_matrix_func)(position + v0 * 0.5 * dt);
        Matrix2s a2 = (this->*get_affine_matrix_func)(position + (v1 * 2.0 - v0) * dt);
        *affine_matrix = a0 / 6.0 + a1 * (2.0 / 3.0) + a2 / 6.0;
      }
      return v0 / 6.0 + v1 * (2.0 / 3.0) + v2 / 6.0;
    }
    case FluidSim::VO_RK4: {
      Vector2s v0 = (this->*get_velocity_func)(position);
      Vector2s v1 = (this->*get_velocity_func)(position + v0 * 0.5 * dt);
      Vector2s v2 = (this->*get_velocity_func)(position + v1 * 0.5 * dt);
      Vector2s v3 = (this->*get_velocity_func)(position + v2 * dt);
      if (affine_matrix) {
        Matrix2s a0 = (this->*get_affine_matrix_func)(position);
        Matrix2s a1 = (this->*get_affine_matrix_func)(position + v0 * 0.5 * dt);
        Matrix2s a2 = (this->*get_affine_matrix_func)(position + v1 * 0.5 * dt);
        Matrix2s a3 = (this->*get_affine_matrix_func)(position + v2 * dt);
        *affine_matrix = a0 / 6.0 + a1 / 3.0 + a2 / 3.0 + a3 / 6.0;
      }
      return v0 / 6.0 + v1 / 3.0 + v2 / 3.0 + v3 / 6.0;
    }
    default:
      return Vector2s::Zero();
  }
}


Vector2s FluidSim::get_air_velocity_and_affine_matrix_with_order(const Vector2s& position, scalar dt, FluidSim::VELOCITY_ORDER v_order,
                                                             FluidSim::INTERPOLATION_ORDER i_order, Matrix2s* affine_matrix) {
  // auto get_velocity_func = (i_order == IO_LINEAR) ? &FluidSim::get_air_velocity : &FluidSim::get_velocity_quadratic; //TODO:添加get_velocity_quadratic的气体实现
  // auto get_affine_matrix_func = (i_order == IO_LINEAR) ? &FluidSim::get_air_affine_matrix : &FluidSim::get_affine_matrix_quadratic;
  auto get_velocity_func = &FluidSim::get_air_velocity;  // TODO: 添加get_velocity_quadratic的气体实现
  auto get_affine_matrix_func = &FluidSim::get_air_affine_matrix;

  switch (v_order) {
    case FluidSim::VO_EULER:
      if (affine_matrix) {
        *affine_matrix = (this->*get_affine_matrix_func)(position);
      }
      return (this->*get_velocity_func)(position);
    case FluidSim::VO_RA2: {
      Vector2s v0 = (this->*get_velocity_func)(position);
      Vector2s v1 = (this->*get_velocity_func)(position + v0 * 2.0 / 3.0 * dt);
      if (affine_matrix) {
        Matrix2s a0 = (this->*get_affine_matrix_func)(position);
        Matrix2s a1 = (this->*get_affine_matrix_func)(position + v0 * 2.0 / 3.0 * dt);
        *affine_matrix = a0 * 0.25 + a1 * 0.75;
      }
      return v0 * 0.25 + v1 * 0.75;
    }
    case FluidSim::VO_RK3: {
      Vector2s v0 = (this->*get_velocity_func)(position);
      Vector2s v1 = (this->*get_velocity_func)(position + v0 * 0.5 * dt);
      Vector2s v2 = (this->*get_velocity_func)(position + (v1 * 2.0 - v0) * dt);
      if (affine_matrix) {
        Matrix2s a0 = (this->*get_affine_matrix_func)(position);
        Matrix2s a1 = (this->*get_affine_matrix_func)(position + v0 * 0.5 * dt);
        Matrix2s a2 = (this->*get_affine_matrix_func)(position + (v1 * 2.0 - v0) * dt);
        *affine_matrix = a0 / 6.0 + a1 * (2.0 / 3.0) + a2 / 6.0;
      }
      return v0 / 6.0 + v1 * (2.0 / 3.0) + v2 / 6.0;
    }
    case FluidSim::VO_RK4: {
      Vector2s v0 = (this->*get_velocity_func)(position);
      Vector2s v1 = (this->*get_velocity_func)(position + v0 * 0.5 * dt);
      Vector2s v2 = (this->*get_velocity_func)(position + v1 * 0.5 * dt);
      Vector2s v3 = (this->*get_velocity_func)(position + v2 * dt);
      if (affine_matrix) {
        Matrix2s a0 = (this->*get_affine_matrix_func)(position);
        Matrix2s a1 = (this->*get_affine_matrix_func)(position + v0 * 0.5 * dt);
        Matrix2s a2 = (this->*get_affine_matrix_func)(position + v1 * 0.5 * dt);
        Matrix2s a3 = (this->*get_affine_matrix_func)(position + v2 * dt);
        *affine_matrix = a0 / 6.0 + a1 / 3.0 + a2 / 3.0 + a3 / 6.0;
      }
      return v0 / 6.0 + v1 / 3.0 + v2 / 3.0 + v3 / 6.0;
    }
    default:
      return Vector2s::Zero();
  }
}
/*!
  \brief  Interpolate velocity from the MAC grid.
*/
Vector2s FluidSim::get_velocity(const Vector2s& position) {
  // Interpolate the velocity from the u_ and v_ grids
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);
  scalar u_value = interpolate_value(p0, u_);
  scalar v_value = interpolate_value(p1, v_);

  return Vector2s(u_value, v_value);
}

Vector2s FluidSim::get_air_velocity(const Vector2s& position) {
  // Interpolate the velocity from the u_a_ and v_a_ grids
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);
  scalar u_value = interpolate_value(p0, u_a_);
  scalar v_value = interpolate_value(p1, v_a_);

  return Vector2s(u_value, v_value);
}

Matrix2s FluidSim::get_affine_matrix(const Vector2s& position) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);

  Matrix2s c_;
  c_.col(0) = affine_interpolate_value(p0, u_) / dx_;
  c_.col(1) = affine_interpolate_value(p1, v_) / dx_;

  return c_;
}

Matrix2s FluidSim::get_air_affine_matrix(const Vector2s& position) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);

  Matrix2s c_;
  c_.col(0) = affine_interpolate_value(p0, u_a_) / dx_;
  c_.col(1) = affine_interpolate_value(p1, v_a_) / dx_;

  return c_;
}

Matrix2s FluidSim::get_saved_affine_matrix(const Vector2s& position) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);

  Matrix2s c_;
  c_.col(0) = affine_interpolate_value(p0, saved_u_);
  c_.col(1) = affine_interpolate_value(p1, saved_v_);

  return c_;
}

Vector2s FluidSim::get_saved_velocity(const Vector2s& position) {
  // Interpolate the velocity from the u_ and v_ grids
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);
  scalar u_value = interpolate_value(p0, saved_u_);
  scalar v_value = interpolate_value(p1, saved_v_);

  return Vector2s(u_value, v_value);
}

Vector2s FluidSim::get_saved_air_velocity(const Vector2s& position) {
  // Interpolate the velocity from the u_ and v_ grids
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);
  scalar u_value = interpolate_value(p0, saved_u_a_);
  scalar v_value = interpolate_value(p1, saved_v_a_);

  return Vector2s(u_value, v_value);
}

Vector2s FluidSim::get_saved_velocity_with_order(const Vector2s& position, FluidSim::INTERPOLATION_ORDER i_order) {
  return (i_order == IO_LINEAR) ? get_saved_velocity(position) : get_saved_velocity_quadratic(position);
}

Vector2s FluidSim::get_saved_air_velocity_with_order(const Vector2s& position, FluidSim::INTERPOLATION_ORDER i_order) {
  // return (i_order == IO_LINEAR) ? get_saved_air_velocity(position) : get_saved_velocity_quadratic(position);
  return get_saved_air_velocity(position); //TODO:待实现二次插值的情况
}


/*! compressible fluids settings */
scalar FluidSim::get_density(const Vector2s& position) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0.5, 0.5);
  scalar dens_value = interpolate_value(p0, comp_rho_);

  return dens_value;
}

scalar FluidSim::get_temperature(const Vector2s& position) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0.5, 0.5);
  scalar temp_value = interpolate_value(p0, grid_temp_);
  return temp_value;
}

scalar FluidSim::get_saved_density(const Vector2s& position) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0.5, 0.5);  // 交错网格，密度网格中的(i, j)是压强网格坐标中的 (i - 1/2, j - 1/2)
  scalar dens_value = interpolate_value(p0, saved_comp_rho_);

  return dens_value;
}

/*!
  \brief  Given two signed distance values, determine what fraction of a connecting
          segment is "inside"
          插值计算两个异号距离场函数采样点间，流体边界处距离场函数=0的位置
          返回一个theta值，等于液体一侧（phi < 0）到0值间的距离与采样距离的比值
          两点均小于0时，即两点间在液体内部，返回1；
          两点均大于0时，即两点均在气体中，返回0；
          返回值在[0, 1]之间
*/
scalar fraction_inside(scalar phi_left, scalar phi_right) {
  if (phi_left < 0 && phi_right < 0) return 1;
  if (phi_left < 0 && phi_right >= 0) return phi_left / (phi_left - phi_right);
  if (phi_left >= 0 && phi_right < 0)
    return phi_right / (phi_right - phi_left);
  else
    return 0;
}

/*!
  \brief  Compute finite-volume style face-weights for fluid from nodal signed
          distances
*/
void FluidSim::compute_weights() {
  parallel_for(0, std::max(u_weights_.nj, v_weights_.nj), [this](int j) {
    if (j < u_weights_.nj) {
      for (int i = 0; i < u_weights_.ni; ++i) {
        u_weights_(i, j) = 1 - fraction_inside(nodal_solid_phi_(i, j + 1), nodal_solid_phi_(i, j));
        u_weights_(i, j) = clamp(u_weights_(i, j), 0.0f, 1.0f);
      }
    }
    if (j < v_weights_.nj) {
      for (int i = 0; i < v_weights_.ni; ++i) {
        v_weights_(i, j) = 1 - fraction_inside(nodal_solid_phi_(i + 1, j), nodal_solid_phi_(i, j));
        v_weights_(i, j) = clamp(v_weights_(i, j), 0.0f, 1.0f);
      }
    }
  });
}

/*!
  \brief  An implementation of the variational pressure projection solve for static
          geometry
*/
void FluidSim::solve_pressure(scalar dt) {
  // This linear system could be simplified, but I've left it as is for clarity
  // and consistency with the standard naive discretization
  int system_size = ni_ * nj_;
  if (rhs_.size() != system_size) {
    rhs_.resize(system_size);
    pressure_.resize(system_size);
    matrix_.resize(system_size);
  }
  matrix_.zero();

  // Build the linear system for pressure
  parallel_for(1, nj_ - 1, [&](int j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      int index = i + ni_ * j;
      rhs_[index] = 0;
      pressure_[index] = 0;
      float centre_phi = liquid_phi_(i, j);
      if (centre_phi < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
        // right neighbour
        float term = u_weights_(i + 1, j) * dt / sqr(dx_); // TODO: 与原差分方程比，少乘了一个密度，对求解结果有影响吗？
        float right_phi = liquid_phi_(i + 1, j);
        if (right_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index + 1, -term);
        } else {  // liquid_phi_ >=0 时：邻居单元格为空气或固体边界时，其对应的参数(index, index + 1)的系数变化量为0
          float theta = fraction_inside(centre_phi, right_phi);  //theta 是什么？
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);  // TODO:theta是什么？
                                                               // theta 无穷大时，系数变化量为0，对应邻居单元格为固体边界
        }                                                      // theta = 1时，系数变化量为原值，对应邻居单元格为空气
        rhs_[index] -= u_weights_(i + 1, j) * u_(i + 1, j) / dx_; // 交错网格，(i, j)就是 (i - 1/2, j)
                                                                  // u_weights_(i, j) = 0时，左邻单元格是固体，该邻居单元格边上的速度项为0
        // left neighbour
        term = u_weights_(i, j) * dt / sqr(dx_);
        float left_phi = liquid_phi_(i - 1, j);
        if (left_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index - 1, -term);
        } else {
          float theta = fraction_inside(centre_phi, left_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);
        }
        rhs_[index] += u_weights_(i, j) * u_(i, j) / dx_;

        // top neighbour
        term = v_weights_(i, j + 1) * dt / sqr(dx_);
        float top_phi = liquid_phi_(i, j + 1);
        if (top_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index + ni_, -term);
        } else {
          float theta = fraction_inside(centre_phi, top_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);
        }
        rhs_[index] -= v_weights_(i, j + 1) * v_(i, j + 1) / dx_;

        // bottom neighbour
        term = v_weights_(i, j) * dt / sqr(dx_);
        float bot_phi = liquid_phi_(i, j - 1);
        if (bot_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index - ni_, -term);
        } else {
          float theta = fraction_inside(centre_phi, bot_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);
        }
        rhs_[index] += v_weights_(i, j) * v_(i, j) / dx_;
      }
    }
  });

  /*if (outframe_ == 10 ) {
    output_matrix_and_rhs_to_csv("D:/FluidSimulator/new_apic2d/matrix_" + std::to_string(outframe_) + ".csv",
                                 "D:/FluidSimulator/new_apic2d/rhs_" + std::to_string(outframe_) + ".csv");
  }*/

  // Solve the system using Robert Bridson's incomplete Cholesky PCG solver
  scalar residual;
  int iterations;
  bool success = solver_.solve(matrix_, rhs_, pressure_, residual, iterations);
  if (!success) {
    std::cout << "WARNING: Pressure solve failed! residual = " << residual << ", iters = " << iterations << std::endl;
  }

  // Apply the velocity update
  parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
    if (j < u_.nj) {
      for (int i = 1; i < u_.ni - 1; ++i) {
        int index = i + j * ni_;
        if (u_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i - 1, j) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i - 1, j) >= 0) theta = fraction_inside(liquid_phi_(i - 1, j), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;
            u_(i, j) -= dt * (pressure_[index] - pressure_[index - 1]) / dx_ / theta;
            u_valid_(i, j) = 1;
          } else {
            u_valid_(i, j) = 0;
          }
        } else {
          u_(i, j) = 0;
          u_valid_(i, j) = 0;
        }
      }    
    }
    if (j >= 1 && j < v_.nj - 1) {
      for (int i = 0; i < v_.ni; ++i) {
        int index = i + j * ni_;
        if (v_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i, j - 1) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i, j - 1) >= 0) theta = fraction_inside(liquid_phi_(i, j - 1), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;
            v_(i, j) -= dt * (pressure_[index] - pressure_[index - ni_]) / dx_ / theta;
            v_valid_(i, j) = 1;
          } else {
            v_valid_(i, j) = 0;
          }
        } else {
          v_(i, j) = 0;
          v_valid_(i, j) = 0;
        }
      }
    }
  });
}


void FluidSim::solve_pressure_with_air(scalar dt) {
  // This linear system could be simplified, but I've left it as is for clarity
  // and consistency with the standard naive discretization
  int system_size = ni_ * nj_;
  if (rhs_.size() != system_size) {
    rhs_.resize(system_size);
    pressure_.resize(system_size);
    matrix_.resize(system_size);
  }
  matrix_.zero();

  // Build the linear system for pressure
  parallel_for(1, nj_ - 1, [&](int j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      int index = i + ni_ * j;
      rhs_[index] = 0;
      pressure_[index] = 0;

      // scalar gema = .0f;
      scalar gema = 7.3f;  // γ for water and air at normal conditions is approximately 0.073J/m^2
      
      // float centre_phi = liquid_phi_(i, j);
      float centre_phi = merged_phi_(i, j);
      float right_phi = merged_phi_(i + 1, j);
      float left_phi = merged_phi_(i - 1, j);
      float top_phi = merged_phi_(i, j + 1);
      float bot_phi = merged_phi_(i, j - 1);

      if ((liquid_phi_(i, j) < 0 || air_phi_(i, j) < 0) &&
          (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
        // right neighbour
        float rho_inter = rho_;
        float theta = fraction_inside(centre_phi, right_phi);
        rho_inter = rho_ * theta + (1 - theta) * rho_air_;
        float term = u_weights_(i + 1, j) * dt / sqr(dx_) / rho_inter;

        if (liquid_phi_(i + 1, j) < 0 || air_phi_(i + 1, j) < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index + 1, -term);
          // 添加表面张力压强 
          // TODO: 1. 检查求表面曲率的维度是否正确 2. 注意center_phi为液体和气体两种情况时的符号是否需要改变
          if (theta > 0 && theta < 1) {
            int sign = (right_phi > 0) - (right_phi < 0);
            rhs_[index] += sign * gema * compute_curvature(i + 1, j, 0, 1) * dt / sqr(dx_) / rho_inter;  // sign为正时，对应right_phi > 0的情况
          }
        } else {  // 既无气体也无液体的真空单元格单独处理
          matrix_.add_to_element(index, index, term);
          //float term = u_weights_(i + 1, j) * dt / sqr(dx_);
          //// float theta = fraction_inside(centre_phi, right_phi);                                                          // theta 是什么？
          //if (theta < 0.01) theta = 0.01;
          //matrix_.add_to_element(index, index, term / theta);  
          // rhs_[index] -= u_weights_(i + 1, j) * u_(i + 1, j) / dx_;  
        }
        float face_frac = compute_face_fraction(centre_phi, right_phi);
        rhs_[index] -= u_weights_(i + 1, j) * (face_frac * u_(i + 1, j) + (1.0f - face_frac) * u_a_(i + 1, j)) / dx_;  // 交错网格，(i, j)就是 (i - 1/2, j)
        
        // left neighbour
        theta = fraction_inside(centre_phi, left_phi);
        rho_inter = rho_ * theta + (1 - theta) * rho_air_;
        term = u_weights_(i, j) * dt / sqr(dx_) / rho_inter;
        
        if (liquid_phi_(i - 1, j) < 0 || air_phi_(i - 1, j) < 0) {
          matrix_.add_to_element(index, index, term);  // TODO:注意差分方程的符号：p_i_j项为正
          matrix_.add_to_element(index, index - 1, -term);

        if (theta > 0 && theta < 1) {
          int sign = (left_phi > 0) - (left_phi < 0);
          rhs_[index] += sign * gema * compute_curvature(i - 1, j, 0, 1) * dt / sqr(dx_) / rho_inter;
          }
        } else {
          matrix_.add_to_element(index, index, term);
          /*float term = u_weights_(i + 1, j) * dt / sqr(dx_);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);*/  
        }

        face_frac = compute_face_fraction(centre_phi, left_phi);
        rhs_[index] += u_weights_(i, j) * (face_frac * u_(i, j) / dx_ + (1.0f - face_frac) * u_a_(i, j) / dx_);

        // top neighbour
        theta = fraction_inside(centre_phi, top_phi);
        rho_inter = rho_ * theta + (1 - theta) * rho_air_;
        term = v_weights_(i, j + 1) * dt / sqr(dx_) / rho_inter;

        if (liquid_phi_(i, j + 1) < 0 || air_phi_(i, j + 1) < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index + ni_, -term);

          if (theta > 0 && theta < 1) {
            int sign = (top_phi > 0) - (top_phi < 0);
            rhs_[index] += sign * gema * compute_curvature(i, j + 1, 1, 0) * dt / sqr(dx_) / rho_inter;
            // std::cout << "tension_term: " << gema * compute_curvature(i, j + 1) * dt / sqr(dx_) / rho_inter << std::endl;
          }
        } else {
          matrix_.add_to_element(index, index, term);
        }
        /*if (i == 41 && j == 10) {
          std::cout << "i: " << i << ", j: " << j << ", term: " << term << ", rho_inter: " << rho_inter << std::endl;
        }*/
        face_frac = compute_face_fraction(centre_phi, top_phi);
        rhs_[index] -= v_weights_(i, j + 1) * (face_frac * v_(i, j + 1) / dx_ + (1.0f - face_frac) * v_a_(i, j + 1) / dx_);
        
        // bottom neighbour
        theta = fraction_inside(centre_phi, bot_phi);
        rho_inter = rho_ * theta + (1 - theta) * rho_air_;
        term = v_weights_(i, j) * dt / sqr(dx_) / rho_inter;
        /*if (i == 41 && j == 11) {
          std::cout << "i: " << i << ", j: " << j << ", term: " << term << ", rho_inter: " << rho_inter << std::endl;
        }*/
        if (liquid_phi_(i, j - 1) < 0 || air_phi_(i, j - 1) < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index - ni_, -term);
          if (theta > 0 && theta < 1) {
            int sign = (bot_phi > 0) - (bot_phi < 0);
            rhs_[index] += sign * gema * compute_curvature(i, j - 1, 1, 0) * dt / sqr(dx_) / rho_inter;
            // std::cout << "tension_term: " << gema * compute_curvature(i, j - 1) * dt / sqr(dx_) / rho_inter << std::endl;
          }
        } else {
          matrix_.add_to_element(index, index, term);
          /*float term = v_weights_(i, j) * dt / sqr(dx_);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);*/
        }

        face_frac = compute_face_fraction(centre_phi, bot_phi);
        rhs_[index] += v_weights_(i, j) * (face_frac * v_(i, j) / dx_ + (1.0f - face_frac) * v_a_(i, j) / dx_);
      }
    }
  });

  /* //检查矩阵是否正定
  if (!is_symmetric(matrix_)) {
    // std::cout << std::endl << "Matrix is not symmetric." << std::endl;
  } else {
    std::cout << "Matrix is symmetric." << std::endl;
  }
  */

  // Solve the system using Robert Bridson's incomplete Cholesky PCG solver
  scalar residual;
  int iterations;
  bool success = solver_.solve(matrix_, rhs_, pressure_, residual, iterations);
  if (!success) {
    std::cout << "WARNING: Pressure solve failed! residual = " << residual << ", iters = " << iterations << std::endl;
  } else {
    // std::cout << "INFO: residual = " << residual << ", iters = " << iterations << std::endl;
  }

  // Apply the velocity update
  parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
    if (j < u_.nj) {
      for (int i = 1; i < u_.ni - 1; ++i) {
        int index = i + j * ni_;
        if (u_weights_(i, j) > 0) {
          float theta = fraction_inside(merged_phi_(i - 1, j), merged_phi_(i, j));
          float face_frac = compute_face_fraction(merged_phi_(i - 1, j), merged_phi_(i, j));
          float rho_inter = rho_ * theta + (1.0f - theta) * rho_air_;
          if (face_frac > 0) {
            u_(i, j) -= dt * (pressure_[index] - pressure_[index - 1]) / dx_ / rho_inter;
            u_valid_(i, j) = 1;
          }
          if (face_frac < 1) {
            u_a_(i, j) -= dt * (pressure_[index] - pressure_[index - 1]) / dx_ / rho_inter;
            // 需要增加u_a_valid_吗？
          }
        } else {
          u_(i, j) = 0;
          u_valid_(i, j) = 0;
        }
      }
    }
    if (j >= 1 && j < v_.nj - 1) {
      for (int i = 0; i < v_.ni; ++i) {
        int index = i + j * ni_;
        if (v_weights_(i, j) > 0) {
          /*if (liquid_phi_(i, j) < 0 || liquid_phi_(i, j - 1) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i, j - 1) >= 0) theta = fraction_inside(liquid_phi_(i, j - 1), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;

            float rho_inter = rho_ * theta + (1.0f - theta) * rho_air_;
            v_(i, j) -= dt * (pressure_[index] - pressure_[index - ni_]) / dx_ / rho_inter;
            // v_(i, j) -= dt * (pressure_[index] - pressure_[index - ni_]) / dx_ / theta;
            v_valid_(i, j) = 1;
          } else {
            v_valid_(i, j) = 0;
          }*/ 
          float theta = fraction_inside(merged_phi_(i, j - 1), merged_phi_(i, j));
          float face_frac = compute_face_fraction(merged_phi_(i, j - 1), merged_phi_(i, j));
          float rho_inter = rho_ * theta + (1.0f - theta) * rho_air_;
          if (face_frac > 0) {
            v_(i, j) -= dt * (pressure_[index] - pressure_[index - ni_]) / dx_ / rho_inter;
            v_valid_(i, j) = 1;
          }
          if (face_frac < 1) {
            v_a_(i, j) -= dt * (pressure_[index] - pressure_[index - ni_]) / dx_ / rho_inter;
          }
        } else {
          v_(i, j) = 0;
          v_valid_(i, j) = 0;
        }
      }
    }
  });
}

void FluidSim::solve_pressure_with_rho(scalar dt) {
  // This linear system could be simplified, but I've left it as is for clarity
  // and consistency with the standard naive discretization
  int system_size = ni_ * nj_;
  if (rhs_.size() != system_size) {
    rhs_.resize(system_size);
    pressure_.resize(system_size);
    matrix_.resize(system_size);
  }
  matrix_.zero();

  parallel_for(1, nj_ - 1, [&](int j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      int index = i + ni_ * j;
      rhs_[index] = 0;
      pressure_[index] = 0;
      float centre_phi = liquid_phi_(i, j);
      scalar rho = grid_rho_(i, j); // 密度
      if (centre_phi < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
        // right neighbour
        float term = u_weights_(i + 1, j) * dt / sqr(dx_) / rho;  // TODO: 与原差分方程比，少乘了一个密度，对求解结果有影响吗？
        float right_phi = liquid_phi_(i + 1, j);
        if (right_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index + 1, -term);
        } else {  // liquid_phi_ >=0 时：邻居单元格为空气或固体边界时，其对应的参数(index, index + 1)的系数变化量为0
          float theta = fraction_inside(centre_phi, right_phi);  // theta 是什么？
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);  // TODO:theta是什么？
                                                               // theta 无穷大时，系数变化量为0，对应邻居单元格为固体边界
        }  // theta = 1时，系数变化量为原值，对应邻居单元格为空气
        rhs_[index] -= u_weights_(i + 1, j) * u_(i + 1, j) / dx_;  // 交错网格，(i, j)就是 (i - 1/2, j)
                                                                   // u_weights_(i, j) = 0时，左邻单元格是固体，该邻居单元格边上的速度项为0
        // left neighbour
        term = u_weights_(i, j) * dt / sqr(dx_) / rho;
        float left_phi = liquid_phi_(i - 1, j);
        if (left_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index - 1, -term);
        } else {
          float theta = fraction_inside(centre_phi, left_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);
        }
        rhs_[index] += u_weights_(i, j) * u_(i, j) / dx_;

        // top neighbour
        term = v_weights_(i, j + 1) * dt / sqr(dx_) / rho;
        float top_phi = liquid_phi_(i, j + 1);
        if (top_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index + ni_, -term);
        } else {
          float theta = fraction_inside(centre_phi, top_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);
        }
        rhs_[index] -= v_weights_(i, j + 1) * v_(i, j + 1) / dx_;

        // bottom neighbour
        term = v_weights_(i, j) * dt / sqr(dx_) / rho;
        float bot_phi = liquid_phi_(i, j - 1);
        if (bot_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index - ni_, -term);
        } else {
          float theta = fraction_inside(centre_phi, bot_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);
        }
        rhs_[index] += v_weights_(i, j) * v_(i, j) / dx_;
      }
    }
  });

  /*if (outframe_ == 10 ) {
    output_matrix_and_rhs_to_csv("D:/FluidSimulator/new_apic2d/matrix_" + std::to_string(outframe_) + ".csv",
                                 "D:/FluidSimulator/new_apic2d/rhs_" + std::to_string(outframe_) + ".csv");
  }*/

  // Solve the system using Robert Bridson's incomplete Cholesky PCG solver
  scalar residual;
  int iterations;
  bool success = solver_.solve(matrix_, rhs_, pressure_, residual, iterations);
  if (!success) {
    std::cout << "WARNING: Pressure solve failed! residual = " << residual << ", iters = " << iterations << std::endl;
  }

  // Apply the velocity update
  parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
    if (j < u_.nj) {
      for (int i = 1; i < u_.ni - 1; ++i) {
        int index = i + j * ni_;
        if (u_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i - 1, j) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i - 1, j) >= 0) theta = fraction_inside(liquid_phi_(i - 1, j), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;
            u_(i, j) -= dt * (pressure_[index] - pressure_[index - 1]) / dx_ / theta / grid_rho_(i, j);
            u_valid_(i, j) = 1;
          } else {
            u_valid_(i, j) = 0;
          }
        } else {
          u_(i, j) = 0;
          u_valid_(i, j) = 0;
        }
      }
    }
    if (j >= 1 && j < v_.nj - 1) {
      for (int i = 0; i < v_.ni; ++i) {
        int index = i + j * ni_;
        if (v_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i, j - 1) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i, j - 1) >= 0) theta = fraction_inside(liquid_phi_(i, j - 1), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;
            v_(i, j) -= dt * (pressure_[index] - pressure_[index - ni_]) / dx_ / theta;
            v_valid_(i, j) = 1;
          } else {
            v_valid_(i, j) = 0;
          }
        } else {
          v_(i, j) = 0;
          v_valid_(i, j) = 0;
        }
      }
    }
  });
}

scalar FluidSim::solid_distance(const Vector2s& pos, const Boundary& b) const {
  switch (b.type_) {
    case BT_BOX:
      return b.sign_ * box_distance(pos, b.center_, b.parameter_);
    case BT_CIRCLE:
      return b.sign_ * circle_distance(pos, b.center_, b.parameter_(0));
    case BT_TORUS:
      return b.sign_ * torus_distance(pos, b.center_, b.parameter_(0), b.parameter_(1));
    case BT_TRIANGLE:
      return b.sign_ * triangle_distance(pos, b.center_, b.parameter_(0));
    case BT_HEXAGON:
      return b.sign_ * hexagon_distance(pos, b.center_, b.parameter_(0));
    case BT_CYLINDER:
      return b.sign_ * cylinder_distance(pos, b.center_, b.parameter_(0), b.parameter_(1));
    case BT_UNION:
      return union_distance(solid_distance(pos, *b.op0_), solid_distance(pos, *b.op1_));
    case BT_INTERSECTION:
      return intersection_distance(solid_distance(pos, *b.op0_), solid_distance(pos, *b.op1_));
    default:
      return 1e+20;
  }
}

scalar FluidSim::solid_distance(const Vector2s& pos) const { return solid_distance(pos, *root_boundary_); }

void FluidSim::init_random_particles() {
  int num_particle = ni_ * nj_;
  for (int i = 0; i < ni_; ++i) {
    for (int j = 0; j < nj_; ++j) {
      for (int k = 0; k < 2; ++k) {
        scalar x_ = (static_cast<scalar>(i) + 0.5 + ((static_cast<scalar>(rand()) / static_cast<scalar>(RAND_MAX)) * 2.0 - 1.0)) * dx_;
        scalar y = (static_cast<scalar>(j) + 0.5 + ((static_cast<scalar>(rand()) / static_cast<scalar>(RAND_MAX)) * 2.0 - 1.0)) * dx_;
        Vector2s pt = Vector2s(x_, y) + origin_;

        scalar phi = solid_distance(pt);
        // 用emplace_back调用了Particle的构造函数
        if (phi > dx_ * ni_ / 5) particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), rho_, T_, Particle::PT_LIQUID); 
        //if (phi > dx_ * ni_ / 5) particles_.emplace_back(pt, Vector2s::Zero(), dx_ / sqrt(2.0), rho_, T_); 
      }
    }
  }
}

void FluidSim::init_random_particles_2() {
  T_ = _INIT_TEMP;
  int num_particle = ni_ * nj_;
  for (int i = 0; i < ni_; ++i) {
    for (int j = 0; j < nj_; ++j) {
      for (int k = 0; k < 2; ++k) {
        scalar x_ = (static_cast<scalar>(i) + 0.5 + ((static_cast<scalar>(rand()) / static_cast<scalar>(RAND_MAX)) * 2.0 - 1.0)) * dx_;
        scalar y = (static_cast<scalar>(j) + 0.5 + ((static_cast<scalar>(rand()) / static_cast<scalar>(RAND_MAX)) * 2.0 - 1.0)) * dx_;
        Vector2s pt = Vector2s(x_, y) + origin_;

        scalar phi = solid_distance(pt);

       auto is_inside_circle = [](const scalar x, const scalar y, const scalar circle_center_x, const scalar circle_center_y, const scalar radius) {  
          scalar dx = x - circle_center_x;  
          scalar dy = y - circle_center_y;  
          return (dx * dx + dy * dy) < (radius * radius);  
       };
        //if (phi > dx_ * ni_ / 5) particles_.emplace_back(pt, Vector2s::Zero(), dx_ / sqrt(2.0), rho_, T_);
         // 中心生成气体被液体包裹的情况
        /*if (phi > dx_ * ni_ * 0.2 && phi <= dx_ * ni_ * 0.3) particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), rho_, T_, Particle::PT_LIQUID);
        if (phi > dx_ * ni_ * 0.3) particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), 0.01, 373.0f, Particle::PT_AIR);
        */
        // 场景：全液体粒子情况
        // if (phi > dx_ * ni_ * 0.2f) particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), rho_, T_, Particle::PT_LIQUID);
        // 场景：全气体粒子情况 
        // if (phi > dx_ * ni_ * 0.2f) particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), 0.01, 373.0f, Particle::PT_AIR);
        
         // 场景：流体静止，上层液体下层气体
         /*if (phi > 0 && y < nj_/ 2 && y > nj_ / 4) particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), rho_, T_, Particle::PT_LIQUID);
         if (phi > 0 && y <= nj_ / 4) particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), 0.01, 373.0f, Particle::PT_AIR);
        */
        // 场景：流体静止，上层气体下层液体
        /*if (phi > 0 && y < nj_ / 2 && y > nj_ / 2.8f) particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), 0.01, 373.0f, Particle::PT_AIR);
        if (phi > 0 && y <= nj_ / 3) particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), rho_, T_, Particle::PT_LIQUID);*/

        if (phi > 0 && y < nj_ * 0.5f && !is_inside_circle(x_, y, ni_ * 0.5f, nj_ * 0.2f, 2.0f))
         particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), rho_, T_, Particle::PT_LIQUID);
        // 在以 ni_2 , nj_ / 3 为圆心，半径为nj_ / 3内生成气体粒子
        if (is_inside_circle(x_, y, ni_ * 0.5f, nj_ * 0.2f, 2.0f))
          particles_.emplace_back(pt, Vector2s(.0f, .0f), dx_ / sqrt(2.0), rho_air_, 373.0f, Particle::PT_AIR);
      }
    }
  }
} 

void FluidSim::map_p2g() {
  if (interpolation_order == IO_LINEAR)
    map_p2g_linear();
  else
    map_p2g_quadratic();
}

void FluidSim::map_p2g_linear() {
  // u-component of velocity
  parallel_for(0, nj_ + 1, [this](int j) {
    if (j < nj_) {
      for (int i = 0; i < ni_ + 1; ++i) {
        Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
        scalar sumw = 0.0;
        scalar sumu = 0.0;

        scalar sumw_air = 0.0;
        scalar sumu_air = 0.0;
        m_sorter_->getNeigboringParticles_cell(i, j, -1, 0, -1, 1, [&](const NeighborParticlesType& neighbors) {
          for (const Particle* p : neighbors) {
            scalar w = p->mass_ * kernel::linear_kernel(p->x_ - pos, dx_);
            if (p->type_ == Particle::PT_AIR) {
              sumu_air += w * (p->v_(0) + p->c_.col(0).dot(pos - p->x_));
              sumw_air += w;
            } else {
            sumu += w * (p->v_(0) + p->c_.col(0).dot(pos - p->x_));  // p->c ？？
            sumw += w;
            }
          }
        });

        u_(i, j) = sumw ? sumu / sumw : 0.0;
        u_a_(i, j) = sumw_air ? sumu_air / sumw_air : 0.0;
      }
    }

    for (int i = 0; i < ni_; ++i) {
      Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
      scalar sumw = 0.0;
      scalar sumu = 0.0;

      scalar sumw_air = 0.0;
      scalar sumu_air = 0.0;
      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 0, [&](const NeighborParticlesType& neighbors) {
        for (const Particle* p : neighbors) {
          scalar w = p->mass_ * kernel::linear_kernel(p->x_ - pos, dx_);
          if (p->type_ == Particle::PT_AIR) {
            sumu_air += w * (p->v_(1) + p->c_.col(1).dot(pos - p->x_));
            sumw_air += w;
          } else {
            sumu += w * (p->v_(1) + p->c_.col(1).dot(pos - p->x_));
            sumw += w;
          }
        }
      });

      v_(i, j) = sumw ? sumu / sumw : 0.0;
      v_a_(i, j) = sumw_air ? sumu_air / sumw_air : 0.0;
    }
  
    // 温度传递
    for (int i = 0; i < ni_; ++i) {
      Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_;
      scalar grid_mass = 0.0;
      scalar sum_T_x_mass = 0.0;
      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
        for (const Particle* p : neighbors) {
          scalar w = kernel::quadratic_kernel(p->x_ - pos, dx_);  // TODO:温度的p2g和g2p目前全部使用的是quadratic_kernel
          grid_mass += w * p->mass_;
          sum_T_x_mass += p->temp_ * p->mass_ * w;
        }
      });
      grid_temp_(i, j) = grid_mass ? sum_T_x_mass / grid_mass : 0.0;
      /*if (grid_temp_(i, j) > 0.1 && grid_temp_(i, j) < 299.99f)
        std::cout << "grid_temp_:" << grid_temp_(i, j) << std::endl;*/
    }
  });

  //std::cout <<"一轮模拟结束！" << std::endl;
}

void FluidSim::map_p2g_quadratic() {
  // u-component of velocity
  parallel_for(0, nj_ + 1, [this](int j) {
    if (j < nj_) {
      for (int i = 0; i < ni_ + 1; ++i) {
        Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
        scalar sumw = 0.0;
        scalar sumu = 0.0;

        scalar sumw_air = 0.0;
        scalar sumu_air = 0.0;
        m_sorter_->getNeigboringParticles_cell(i, j, -2, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
          for (const Particle* p : neighbors) {
            scalar w = p->mass_ * kernel::quadratic_kernel(p->x_ - pos, dx_);  // w*m
            if (p->type_ == Particle::PT_AIR) {
              sumu_air += w * (p->v_(0) + p->c_.col(0).dot(pos - p->x_));
              sumw_air += w;
            } else {
              sumu += w * (p->v_(0) + p->c_.col(0).dot(pos - p->x_));
              sumw += w;
            }
          }
        });

        u_(i, j) = sumw > 0.0 ? sumu / sumw : 0.0;
        u_a_(i, j) = sumw_air > 0.0 ? sumu_air / sumw_air : 0.0;
      }
    }

    for (int i = 0; i < ni_; ++i) {
      Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
      scalar sumw = 0.0;
      scalar sumu = 0.0;

      scalar sumw_air = 0.0;
      scalar sumu_air = 0.0;
      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -2, 1, [&](const NeighborParticlesType& neighbors) {
        for (const Particle* p : neighbors) {
          scalar w = p->mass_ * kernel::quadratic_kernel(p->x_ - pos, dx_);
          if (p->type_ == Particle::PT_AIR) {
            sumu_air += w * (p->v_(1) + p->c_.col(1).dot(pos - p->x_));
            sumw_air += w;
          } else {
            sumu += w * (p->v_(1) + p->c_.col(1).dot(pos - p->x_));
            sumw += w;
          }
        }
      });

      v_(i, j) = sumw > 0.0 ? sumu / sumw : 0.0;
      v_a_(i, j) = sumw_air > 0.0 ? sumu_air / sumw_air : 0.0;
    }
    // 温度传递
    for (int i = 0; i < ni_; ++i) {
      Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_;
      scalar grid_mass = 0.0;
      scalar sum_T_x_mass = 0.0;
      // m_sorter_->getNeigboringParticles_cell(i, j, -1, 0, -1, 1, [&](const NeighborParticlesType& neighbors) { // 原来查找的是2x3邻域网格中的粒子，不知为何要这样写
      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
      for (const Particle* p : neighbors) {
          scalar w = kernel::quadratic_kernel(p->x_ - pos, dx_);
          grid_mass += w * p->mass_;
          sum_T_x_mass += p->temp_ * p->mass_ * w;
        }
      });
      grid_temp_(i, j) = grid_mass ? sum_T_x_mass / grid_mass : grid_temp_(i, j);
    }
    // 密度传递
    for (int i = 0; i < ni_; ++i) {
        Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_;
        scalar sumw = 0.0; // 问题：对于密度这样和mass_之间挂钩的量，需要gather后除以总质量吗？
        scalar grid_rho = 0.0; 
        m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
          for (const Particle* p : neighbors) {
            scalar w = kernel::quadratic_kernel(p->x_ - pos, dx_);
            grid_rho += p->dens_ * w;
          }
        });
        grid_rho_(i, j) = grid_rho ? grid_rho : grid_rho_(i, j);
        // grid_rho_(i, j) = sumw ? grid_rho / sumw : grid_rho_(i, j);
      }
 });
}

void FluidSim::map_p2g_compressible() {
    // u-component of momentum
    parallel_for(0, nj_ + 1, [this](int j) {
        if (j < nj_) {
            for (int i = 0; i < ni_ + 1; ++i) {
                Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
                scalar sumw = 0.0;
                scalar sumu = 0.0;
                m_sorter_->getNeigboringParticles_cell(i, j, -1, 0, -1, 1, [&](const NeighborParticlesType& neighbors) {
                    for (const Particle* p : neighbors) {
                        scalar w = p->mass_ * kernel::linear_kernel(p->x_ - pos, dx_);
                        sumu += w * (p->v_(0) + p->c_.col(0).dot(pos - p->x_));
                        sumw += w;
                    }
                });
                u_(i, j) = sumw ? sumu / sumw : 0.0;
            }
        }

        for (int i = 0; i < ni_; ++i) {
            Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
            scalar sumw = 0.0;
            scalar sumu = 0.0;
            m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 0, [&](const NeighborParticlesType& neighbors) {
                    for (const Particle* p : neighbors) {
                        scalar w = p->mass_ * kernel::linear_kernel(p->x_ - pos, dx_);
                        sumu += w * (p->v_(1) + p->c_.col(1).dot(pos - p->x_));
                        sumw += w;
                    }
            });
            v_(i, j) = sumw ? sumu / sumw : 0.0;
        }

        // rho  //注意：这里对rho的插值完全没有使用粒子的质量
        if (j < nj_) {
            for (int i = 0; i < ni_; ++i) {
            Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_;  // 交错网格，密度在压强单元格中心
                scalar sumw = 0.0;
                scalar sumrho = 0.0;
                m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
                    for (const Particle* p : neighbors) {
                        scalar w = kernel::linear_kernel(p->x_ - pos, dx_);
                        sumw += w;
                        sumrho += p->dens_* w;
                    }
                });
                comp_rho_(i,j) = sumw ? sumrho / sumw : 0.0;// TODO: 添加体积元不同带来的密度变化？
            }
        }
    });

    
}

/*!
  \brief  A general affine FLIP scheme that unifies all the other FLIP schemes
          used in this code
*/
void FluidSim::map_g2p_flip_general(float dt, const scalar lagrangian_ratio, const scalar lagrangian_symplecticity, const scalar affine_stretching_ratio,
                                    const scalar affine_rotational_ratio) {
    bool use_affine = affine_stretching_ratio > 0. || affine_rotational_ratio > 0.;
    scalar sum_rho = 0.0;
    scalar sum_temp = 0.0;
    parallel_for(0, static_cast<int>(particles_.size()), [&](int i) {
        auto& p = particles_[i];

        Matrix2s C = Matrix2s::Zero();  // APIC使用
        Vector2s next_grid_velocity = (p.type_ == Particle::PT_LIQUID ? get_velocity_and_affine_matrix_with_order(p.x_, dt, velocity_order, interpolation_order, use_affine ? (&C) : nullptr)
                        : get_air_velocity_and_affine_matrix_with_order(p.x_, dt, velocity_order, interpolation_order, use_affine ? (&C) : nullptr));
        Vector2s lagrangian_velocity = p.v_;
        Vector2s original_grid_velocity;

        if (lagrangian_ratio > 0.0) {
          original_grid_velocity = (p.type_ == Particle::PT_LIQUID ? get_saved_velocity_with_order(p.x_, interpolation_order)
                                                                   : get_saved_air_velocity_with_order(p.x_, interpolation_order));
            p.v_ = next_grid_velocity + (lagrangian_velocity - original_grid_velocity) * lagrangian_ratio;  // flip：只转换增量部分
        } else {
            p.v_ = next_grid_velocity;
        }
        // 温度传递回粒子
        scalar lagrangian_temp = p.temp_;
        scalar next_grid_temp = get_temperature_quadratic(p.x_, grid_temp_);
        if (lagrangian_ratio > 0.0) {
          scalar original_grid_temp = next_grid_temp;  // TODO:实现flip需要的saved_temp存储，目前暂时用next_grid_temp代替
          p.temp_ = next_grid_temp + (lagrangian_temp - original_grid_temp) * lagrangian_ratio;
        } else {
          p.temp_ = next_grid_temp;
        }
        sum_temp += p.temp_;

        // 根据温度发生相变
        if (p.temp_ > 373.0f && p.type_ == Particle::PT_LIQUID) {
          p.type_ = Particle::PT_AIR;
          p.dens_ = rho_air_;
        } else if (p.temp_ < 373.0f && p.type_ == Particle::PT_AIR) {
          p.type_ = Particle::PT_LIQUID;
          p.dens_ = rho_;
        }

#ifdef COMPRESSIBLE_FLUID
        // 密度传递回粒子
        scalar lagrangian_rho = p.dens_;
        scalar next_grid_rho = get_temperature_quadratic(p.x_,comp_rho_);
        if (lagrangian_ratio > 0.0) {
          scalar original_grid_rho = get_temperature_quadratic(p.x_, saved_comp_rho_);
            p.dens_ = next_grid_rho + (lagrangian_rho - original_grid_rho) * lagrangian_ratio;
            /*if (next_grid_rho - original_grid_rho > 0) {
              std::cout << "网格密度更新值：" << next_grid_rho - original_grid_rho << std::endl;
            }*/
        } else {
            p.dens_ = next_grid_rho;
        }
        
        sum_rho += p.dens_;
        //TODO:发现使用任何公式更新质量后，密度都会迅速衰减为0（具体表现：在流体自由落体到边界前，程序就崩了）
        //p.mass_ = p.dens_;
        //p.mass_ = M_PI * p.radii_ * p.radii_ * p.dens_;
        //p.mass_ = 4.0 / 3.0 * M_PI * p.radii_ * p.radii_ * p.dens_;
#endif

        p.c_ = C * (affine_stretching_ratio + affine_rotational_ratio) * 0.5 + C.transpose() * (affine_stretching_ratio - affine_rotational_ratio) * 0.5;

        if (lagrangian_symplecticity > 0.0 && lagrangian_ratio > 0.0) {
            scalar logdJ = C.trace() * dt;
            p.logJ_ += logdJ;
            scalar pos_adj = lagrangian_symplecticity;
            Vector2s pp = (p.x_ - origin_) / dx_;
            scalar phi_value = interpolate_value(pp, nodal_solid_phi_);
            if (phi_value < 0) {
                Vector2s normal;
                interpolate_gradient(normal, pp, nodal_solid_phi_);
                normal.normalize();
                float dotnv = -normal.dot(lagrangian_velocity);
                if (dotnv > 0.0) {
                    pos_adj = 0.0;
                }
            }
            if (p.logJ_ < -0.001) {
                pos_adj = 0.0;
            }
            p.x_ += (next_grid_velocity + (lagrangian_velocity - original_grid_velocity) * lagrangian_ratio * lagrangian_symplecticity) * dt;
        } 
        else {
            p.x_ += next_grid_velocity * dt;
        }
    });
    // std::cout << "sum_temp: " << sum_temp << ", average: " << sum_temp / particles_.size() << std::endl;
#ifdef COMPRESSIBLE_FLUID
    std::cout << "sum_rho in particles: " << sum_rho << std::endl;
    ////奇怪现象：加了个循环语句后，求解反而稳定了一些？
    //for (int i = 0; i < particles_.size(); ++i) {
    //  // 质量守恒
    // particles_[i].dens_ = particles_[i].dens_ / sum_rho * (particles_.size() * 1.0); //原密度为1
    // particles_[i].mass_ = M_PI * particles_[i].radii_ * particles_[i].radii_ * particles_[i].dens_;
    //}
#endif
}
void FluidSim::render_boundaries(const Boundary& b) {
  switch (b.type_) {
    case BT_CIRCLE:
      draw_circle2d(b.center_, b.parameter_(0), 64);
      break;
    case BT_BOX:
      draw_box2d(b.center_, b.parameter_(0), b.parameter_(1));
      break;
    case BT_TORUS:
      draw_circle2d(b.center_, b.parameter_(0), 64);
      draw_circle2d(b.center_, b.parameter_(1), 64);
      break;
    case BT_HEXAGON:
      draw_circle2d(b.center_, b.parameter_(0), 6);
      break;
    case BT_TRIANGLE:
      draw_circle2d(b.center_, b.parameter_(0), 3);
      break;
    case BT_UNION:
    case BT_INTERSECTION:
      render_boundaries(*b.op0_);
      render_boundaries(*b.op1_);
      break;
    default:
      break;
  }
}

void FluidSim::render() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // 清除颜色和深度缓冲区
  glPushMatrix();
  glScaled(1.0 / (dx_ * ni_), 1.0 / (dx_ * ni_), 1.0 / (dx_ * ni_));

  if (draw_grid_) {
    glColor3f(0.75, 0.75, 0.75);
    glLineWidth(1);
    draw_grid2d(origin_, dx_, ni_, nj_);
  }

  if (draw_boundaries_) {
    render_boundaries(*root_boundary_);
  }

  if (draw_velocities_) {
    glColor3f(1, 0, 0);
    scalar crit = dx_ * dx_ * 100.0f;
    glBegin(GL_LINES);
    for (int j = 0; j < nj_; ++j) {
      for (int i = 0; i < ni_; ++i) {
        Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_;
        Vector2s vel = get_velocity(pos);
        if (vel.squaredNorm() > crit) {
          Vector2s pos_v = pos + 0.01 * vel;
          glVertex2fv(pos.data());
          glVertex2fv(pos_v.data());
        }
      }
    }
    glEnd();
  }

  if (draw_particles_) {
    glColor3f(0, 0, 1);
    glPointSize(5);
    glBegin(GL_POINTS);
    // 线性插值计算颜色
    bool draw_particles_with_temp_color = false;
    bool draw_particles_by_density = true;
    if (draw_particles_with_temp_color) {
      glBegin(GL_POINTS);
      for (unsigned int i = 0; i < particles_.size(); ++i) {
        // 假设温度范围为 300K（蓝色）到 373K（红色）
        float t = (particles_[i].temp_ - 273.0f) / (373.0f - 273.0f);  // 归一化温度到[0, 1]
        t = std::max(std::min(t, 1.0f), 0.0f);                         // 确保t在[0, 1]范围内
        float r, g, b;
        if (t < 0.5f) {
          // 从蓝色到绿色
          r = 0.0f;
          g = t * 2.0f;         // 绿色分量从 0 增加到 1
          b = 1.0f - t * 2.0f;  // 蓝色分量从 1 减少到 0
        } else {
          // 从绿色到红色
          r = (t - 0.5f) * 2.0f;         // 红色分量从 0 增加到 1
          g = 1.0f - (t - 0.5f) * 2.0f;  // 绿色分量从 1 减少到 0
          b = 0.0f;
        }
        glColor3f(r, g, b);               // 设置颜色
        glVertex2fv(particles_[i].x_.data());  // 绘制粒子
      }
    }
    else if (draw_particles_by_density) {
      if (draw_particles_by_density) {
        glBegin(GL_POINTS);
        for (unsigned int i = 0; i < particles_.size(); ++i) {
          float t = (particles_[i].dens_) / (1.2f - .0f);  // 归一化密度到[0, 1]
          t = std::max(std::min(t, 1.0f), 0.0f);                    // 确保t在[0, 1]范围内
          float r = 0, g = 0, b = 0;
          b = t;                                 // 蓝色分量从 0 增加到 1
          glColor3f(r, g, b);                    
          glVertex2fv(particles_[i].x_.data()); 
        }
      }
    } else {
        for (unsigned int i = 0; i < particles_.size(); ++i) {
          glVertex2fv(particles_[i].x_.data());
        }
    }
    glEnd();
  }
  glPopMatrix();
  //glFlush();
}
void FluidSim::render2() {
  }

FluidSim::Boundary::Boundary(const Vector2s& center, const Vector2s& parameter, BOUNDARY_TYPE type, bool inside)
    : center_(center), parameter_(parameter), type_(type), sign_(inside ? -1.0 : 1.0) {}

FluidSim::Boundary::Boundary(Boundary* op0, Boundary* op1, BOUNDARY_TYPE type) : op0_(op0), op1_(op1), type_(type), sign_(op0 ? op0->sign_ : false) {}

Particle::Particle(const Vector2s& x, const Vector2s& v, const scalar& radii, const scalar& density, const scalar& tempreture,
                   const PARTICLE_TYPE type)
    : x_(x), v_(v), radii_(radii), mass_(4.0 / 3.0 * M_PI * radii * radii * radii * density), logJ_(0), dens_(density), temp_(tempreture), type_(type) {
  c_.setZero();
  buf0_.setZero();
}

Particle::Particle() : x_(Vector2s::Zero()), v_(Vector2s::Zero()), radii_(0.0), mass_(0.0), dens_(0.0), logJ_(0.0), temp_(0), type_(PT_LIQUID) {
  c_.setZero();
  buf0_.setZero();
}

Particle::Particle(const Particle& p) : x_(p.x_), v_(p.v_), radii_(p.radii_), mass_(p.mass_), dens_(p.dens_), logJ_(0.0), temp_(p.temp_), type_(p.type_) {
  c_.setZero();
  buf0_.setZero();
}

/*!
  \brief Apply several iterations of a very simple "Jacobi"-style propagation of
         valid velocity data in all directions
*/
void extrapolate(Array2s& grid, Array2s& old_grid, const Array2s& grid_weight, const Array2s& grid_liquid_weight, Array2c& valid, Array2c old_valid_,
                 const Vector2i& offset, int num_layers) {
  if (!num_layers) return;

  // Initialize the list of valid_ cells
  for (int j = 0; j < valid.nj; ++j) valid(0, j) = valid(valid.ni - 1, j) = 0;
  for (int i = 0; i < valid.ni; ++i) valid(i, 0) = valid(i, valid.nj - 1) = 0;

  Array2s* pgrids[] = {&grid, &old_grid};
  Array2c* pvalids[] = {&valid, &old_valid_};

  // 每次使用valid_前会重新赋值？
  for (int j = 1; j < grid.nj - 1; ++j) {
    for (int i = 1; i < grid.ni - 1; ++i) {
      valid(i, j) = grid_weight(i, j) > 0 && (grid_liquid_weight(i, j) < 0 || grid_liquid_weight(i + offset(0), j + offset(1)) < 0);
    }
  }

  old_valid_ = valid;
  old_grid = grid;

  for (int layers = 0; layers < num_layers; ++layers) {
    Array2s* pgrid_source = pgrids[layers & 1];
    Array2s* pgrid_target = pgrids[!(layers & 1)];

    Array2c* pvalid_source = pvalids[layers & 1];
    Array2c* pvalid_target = pvalids[!(layers & 1)];

    parallel_for(1, grid.nj - 1, [&](int j) {
      for (int i = 1; i < grid.ni - 1; ++i) {
        scalar sum = 0;
        int count = 0;

        if (!(*pvalid_source)(i, j)) {
          if ((*pvalid_source)(i + 1, j)) {
            sum += (*pgrid_source)(i + 1, j);
            ++count;
          }
          if ((*pvalid_source)(i - 1, j)) {
            sum += (*pgrid_source)(i - 1, j);
            ++count;
          }
          if ((*pvalid_source)(i, j + 1)) {
            sum += (*pgrid_source)(i, j + 1);
            ++count;
          }
          if ((*pvalid_source)(i, j - 1)) {
            sum += (*pgrid_source)(i, j - 1);
            ++count;
          }

          // If any of neighbour cells were valid_,
          // assign the cell their average value and tag it as valid_
          if (count > 0) {
            (*pgrid_target)(i, j) = sum / static_cast<scalar>(count);
            (*pvalid_target)(i, j) = 1;
          }
        }
      }
    });

    if (layers != num_layers - 1) {
      for (int j = 1; j < grid.nj - 1; ++j) {
        for (int i = 1; i < grid.ni - 1; ++i) {
          if (!(*pvalid_source)(i, j)) {
            (*pgrid_source)(i, j) = (*pgrid_target)(i, j);
            (*pvalid_source)(i, j) = (*pvalid_target)(i, j);
          }
        }
      }
    }
  }
}


scalar FluidSim::compute_coef_A(const scalar& rho) { 
	scalar background_T = 300.0f;
	scalar denominator = (b * b + 2.0f * rho * b - rho * rho);
	return R * background_T * b * b / ((b - rho) * (b - rho)) - 2.0f * a * b * b * b * rho * (rho + b) / (denominator * denominator); 
}

scalar FluidSim::compute_coef_B(const scalar& rho) {
	scalar background_T = 300.0f;
	scalar numerator = b * b * b + 3.0f * rho * rho * b + 2.0f * rho * rho * rho;
	scalar denominator = (b * b + 2.0f * rho * b - rho * rho);
    return 2.0f * R * background_T * b * b / ((b - rho) * (b - rho) * (b - rho)) - 2.0f * a * b * b * b * numerator / (denominator * denominator * denominator);
}

scalar FluidSim::compute_coef_A2(const scalar& rho) { 
    scalar A_l = 0, rho_l = 1.0, xi = 0.01, eta = 10.0;
    A_l = xi * (1 - eta) / eta * compute_coef_A(rho_l);
    scalar numerator = A_l * rho_l * rho_l * xi * (1 - xi);
    scalar denominator = (1 - xi) * rho_l - rho;
    return numerator / denominator * denominator + A_l;
}

scalar FluidSim::compute_coef_B2(const scalar& rho) { 
    scalar A_l = 0, rho_l = 1.0, xi = 0.01, eta = 10.0; 
    A_l = xi * (1 - eta) / eta * compute_coef_A(rho_l);
    scalar numerator = A_l * rho_l * rho_l * xi * (1 - xi);
    scalar denominator = (1 - xi) * rho_l - rho;
    return -2 * numerator / denominator * denominator * denominator;
}

/*! Output Data Bgeo */
void FluidSim::OutputPointDataBgeo(const std::string& s, const int frame) {
    std::string file = s + std::to_string(frame) + ".bgeo";
    std::cout << "Writing to: " << file << std::endl;

    Partio::ParticlesDataMutable* parts = Partio::create();
  Partio::ParticleAttribute pos, vel, rho, mass, temperature;
    pos = parts->addAttribute("position", Partio::VECTOR, 3);
    vel = parts->addAttribute("velocity", Partio::VECTOR, 3);
    rho = parts->addAttribute("rho", Partio::FLOAT, 1);
    mass = parts->addAttribute("mass", Partio::FLOAT, 1);
    temperature = parts->addAttribute("temperature", Partio::FLOAT, 1);

    for (int i = 0; i < particles_.size(); i++) {
        auto& p = particles_[i];
        Vector2s p_x = p.x_;
        Vector2s p_v = p.v_;
        scalar p_rho = p.dens_;
        scalar p_temp = p.temp_;
        scalar p_mass = p.mass_;

		int idx = parts->addParticle();
		float* x = parts->dataWrite<float>(pos, idx);
		float* v = parts->dataWrite<float>(vel, idx);
        float* dens = parts->dataWrite<float>(rho, idx);
        float* m = parts->dataWrite<float>(mass, idx);
        float* temp = parts->dataWrite<float>(temperature, idx);

		x[0] = p_x[0];
		x[1] = p_x[1];
		x[2] = 0.0f;
		v[0] = p_v[0];
		v[1] = p_v[1];
		v[2] = 0.0f;
        *dens = p_rho;
        *m = p_mass;
        *temp = p_temp;
	}

  Partio::write(file.c_str(), *parts);
  parts->release();
}

void FluidSim::OutputGridDataBgeo(const std::string& s, const int frame) {
  std::string file = s + std::to_string(frame) + ".bgeo";
  std::cout << "Writing to: " << file << std::endl;

  Partio::ParticlesDataMutable* parts = Partio::create();
  Partio::ParticleAttribute pos, rho, press, lapP, liquidPhi;
  pos = parts->addAttribute("position", Partio::VECTOR, 3);
  rho = parts->addAttribute("rho", Partio::FLOAT, 1);
  press = parts->addAttribute("pressure", Partio::FLOAT, 1);
  lapP = parts->addAttribute("laplacianP", Partio::FLOAT, 1);
  liquidPhi = parts->addAttribute("liquid_phi", Partio::FLOAT, 1);
  Partio::ParticleAttribute airPhi = parts->addAttribute("air_phi", Partio::FLOAT, 1);
  Partio::ParticleAttribute mergedPhi = parts->addAttribute("merged_phi", Partio::FLOAT, 1);
  
  for (int j = 0; j < nj_; j++)
    for (int i = 0; i < ni_; i++) {
      Vector2s p_x = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_);
      scalar p_rho = grid_rho_(i, j);
      scalar p_press = pressure_[i + j * ni_];
	  scalar p_lapP = laplacianP_(i, j);
      scalar p_liquid_phi = liquid_phi_(i, j);
      scalar p_air_phi = air_phi_(i, j);
      scalar p_merged_phi = merged_phi_(i, j);

      int idx = parts->addParticle();
      float* x = parts->dataWrite<float>(pos, idx);
      float* dens = parts->dataWrite<float>(rho, idx);
      float* pressure = parts->dataWrite<float>(press, idx);
	  float* laplacianP = parts->dataWrite<float>(lapP, idx);
      float* liquid_phi = parts->dataWrite<float>(liquidPhi, idx);
      float* air_phi = parts->dataWrite<float>(airPhi, idx);
      float* merged_phi = parts->dataWrite<float>(mergedPhi, idx);

      x[0] = p_x[0];
      x[1] = p_x[1];
      x[2] = 0.0f;
      *dens = p_rho;
      *pressure = p_press;
	  *laplacianP = p_lapP;
      *liquid_phi = p_liquid_phi;
      *air_phi = p_air_phi;
      *merged_phi = p_merged_phi;
    }

  Partio::write(file.c_str(), *parts);
  parts->release();
}

void FluidSim::OutputGridXDataBgeo(const std::string& s, const int frame) {
  std::string file = s + std::to_string(frame) + ".bgeo";
  std::cout << "Writing to: " << file << std::endl;

  Partio::ParticlesDataMutable* parts = Partio::create();
  Partio::ParticleAttribute pos, uf, velx, velxa, u_weight_;
  pos = parts->addAttribute("position", Partio::VECTOR, 3);
  velx = parts->addAttribute("u", Partio::FLOAT, 1);
  u_weight_ = parts->addAttribute("u_weight_", Partio::FLOAT, 1);
  Partio::ParticleAttribute air_velx = parts->addAttribute("u_a_", Partio::FLOAT, 1);

  for (int j = 0; j < nj_; j++)
    for (int i = 0; i < ni_ + 1; i++) {
      Vector2s p_x = Vector2s((i)*dx_, (j + 0.5) * dx_);
      scalar p_u = u_(i, j);
      scalar p_u_a = u_a_(i, j);

      int idx = parts->addParticle();
      float* x = parts->dataWrite<float>(pos, idx);
      float* u = parts->dataWrite<float>(velx, idx);
      float* u_weight = parts->dataWrite<float>(u_weight_, idx);
      float* u_a = parts->dataWrite<float>(air_velx, idx);

      x[0] = p_x[0];
      x[1] = p_x[1];
      x[2] = 0.0f;
      *u = p_u;
      *u_weight = u_weights_(i, j);
      *u_a = p_u_a;
    }

  Partio::write(file.c_str(), *parts);
  parts->release();
}

void FluidSim::OutputGridYDataBgeo(const std::string& s, const int frame)  {
  std::string file = s + std::to_string(frame) + ".bgeo";
  std::cout << "Writing to: " << file << std::endl;

  Partio::ParticlesDataMutable* parts = Partio::create();
  Partio::ParticleAttribute pos, vf, vely, velya, v_weight_;
  pos = parts->addAttribute("position", Partio::VECTOR, 3);
  vely = parts->addAttribute("v", Partio::FLOAT, 1);
  v_weight_ = parts->addAttribute("v_weight_", Partio::FLOAT, 1);
  Partio::ParticleAttribute air_vely = parts->addAttribute("v_a_", Partio::FLOAT, 1);

  for (int j = 0; j < nj_ + 1; j++)
    for (int i = 0; i < ni_; i++) {
      Vector2s p_x = Vector2s((i + 0.5) * dx_, (j)*dx_);
      scalar p_v = v_(i, j);
      scalar p_v_a = v_a_(i, j);

      int idx = parts->addParticle();
      float* x = parts->dataWrite<float>(pos, idx);
      float* v = parts->dataWrite<float>(vely, idx);
      float* v_weight = parts->dataWrite<float>(v_weight_, idx);
      float* v_a = parts->dataWrite<float>(air_vely, idx);
      x[0] = p_x[0];
      x[1] = p_x[1];
      x[2] = 0.0f;
      *v = p_v;
      *v_weight = v_weights_(i, j);
      *v_a = p_v_a;
    }

  Partio::write(file.c_str(), *parts);
  parts->release();
}

scalar FluidSim::get_pressure(const scalar& rho)
{
	scalar background_T = 300.0f;
	scalar p_R = R, p_a = a, p_b = b;
	// total_pressure = dynamic_pressure + static_pressure
    auto total_pressure = [p_R, p_a, p_b](scalar T, scalar dens) {
		return p_R * p_b * T * dens / (p_b - dens) - p_a * p_b * p_b * dens * dens / (p_b * p_b + 2.0f * dens * p_b - dens * dens);
    };

	// scalar static_pressure = total_pressure(background_T, rho_);
    // return max(0.0f, total_pressure(background_T, rho) - static_pressure);
    return max(0.0f, (rho - rho_) * 100.0f);
}

//void FluidSim::build_density_equation(int i, int j, ) { 
//    int index = i + ni_ * j;
//}

// this is incomplete version
void FluidSim::solve_compressible_density(scalar dt) {
	int system_size = ni_ * nj_;
	if (rhs_.size() != system_size) {
		rhs_.resize(system_size);
		comp_rho_solution_.resize(system_size);
		matrix_.resize(system_size);
	}
	matrix_.zero();
     
	// Build the linear system for density 
	parallel_for(1, nj_ - 1, [&](int j) {
          for (int i = 1; i < ni_ - 1; ++i) {
				laplacianP_(i, j) = 0.0f;
				int index = i + ni_ * j;
				rhs_[index] = 0;
				comp_rho_solution_[index] = 0;
				float centre_phi = liquid_phi_(i, j);
				float t2dx2 = dt * dt / (dx_ * dx_);
				float coef_A = t2dx2 * (compute_coef_A(comp_rho_(i, j)) + compute_coef_A2(comp_rho_(i, j)));
                float coef_B = (t2dx2) * (compute_coef_B(comp_rho_(i, j)) + compute_coef_B2(comp_rho_(i, j)));
				if (centre_phi < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
					// right neighbour
                    coef_A = t2dx2 * compute_coef_A((comp_rho_(i, j) + comp_rho_(i + 1, j))/2);
                    coef_B = (t2dx2 ) * compute_coef_B((comp_rho_(i, j) + comp_rho_(i + 1, j)) / 2);
					float centre_term = u_weights_(i + 1, j) * coef_A;
					float vel_term = 0.0f;  // u_weights_(i + 1, j) * (t2dx2 * u_(i + 1, j) / (2.0f * dx_));  有什么用？
                    // u_weights_(i + 1, j) * u_(i + 1, j) * dt / dx_是速度散度项
					//float local_term = centre_term + vel_term + u_weights_(i + 1, j) * u_(i + 1, j) * dt / dx_;  
					float local_term = centre_term + vel_term + u_weights_(i + 1, j) * u_(i + 1, j) * dt / dx_;  

					float term = u_weights_(i + 1, j) * (coef_A + coef_B * (comp_rho_(i + 1, j) - comp_rho_(i-1, j))) - vel_term;
					float right_phi = liquid_phi_(i + 1, j);
                    // liquid_phi_小于0时，是液体单元格
                    if (right_phi < 0) {
						matrix_.add_to_element(index, index, local_term);
						matrix_.add_to_element(index, index + 1, -term);
					} else { 
                        //邻居单元格是气体时
						float theta = fraction_inside(centre_phi, right_phi);
						if (theta < 0.01) theta = 0.01;
						matrix_.add_to_element(index, index, centre_term / theta);  // doesm't consider divergence term near interface
					}

					laplacianP_(i,j) += u_weights_(i + 1, j) * ((coef_A + coef_B * (comp_rho_(i + 1, j) - comp_rho_(i, j)))*comp_rho_(i+1,j)-coef_A * comp_rho_(i,j));

					// left neighbour
                    coef_A = t2dx2 * compute_coef_A((comp_rho_(i, j) + comp_rho_(i - 1, j)) * 0.5);
                    coef_B = (t2dx2 ) * compute_coef_B((comp_rho_(i, j) + comp_rho_(i - 1, j)) * 0.5);
					vel_term = 0.0f;//u_weights_(i, j) * (t2dx2 * u_(i, j) / (2.0f * dx_));
					term = u_weights_(i + 1, j) * (coef_A -  coef_B * (comp_rho_(i+1, j) - comp_rho_(i - 1, j))) + vel_term;
					//local_term = centre_term - vel_term -u_weights_(i,j) * u_(i,j) * dt / dx_;
					local_term = centre_term - vel_term -u_weights_(i,j) * u_(i,j) * dt / dx_;
					float left_phi = liquid_phi_(i - 1, j);
					if (left_phi < 0) {
					    matrix_.add_to_element(index, index, local_term);
						matrix_.add_to_element(index, index - 1, -term);
					} else {
						float theta = fraction_inside(centre_phi, left_phi);
						if (theta < 0.01) theta = 0.01;
						matrix_.add_to_element(index, index, centre_term / theta);
					}
					laplacianP_(i,j) += u_weights_(i, j) * ((coef_A + coef_B * (comp_rho_(i, j) - comp_rho_(i-1, j)))*comp_rho_(i-1,j)-coef_A * comp_rho_(i,j));

					// top neighbour
                    coef_A = t2dx2 * compute_coef_A((comp_rho_(i, j) + comp_rho_(i, j+1)) / 2);
                    coef_B = (t2dx2 ) * compute_coef_B((comp_rho_(i, j) + comp_rho_(i, j+1)) / 2);
					vel_term = 0.0f;//v_weights_(i, j + 1) * (t2dx2 * v_(i, j + 1) / (2.0f * dx_));
                    term = v_weights_(i, j + 1) * (coef_A + coef_B * (comp_rho_(i, j + 1) - comp_rho_(i, j-1))) - vel_term;
					// local_term = centre_term + vel_term + v_weights_(i,j+1) * v_(i,j+1) * dt / dx_;
					local_term = centre_term + vel_term + v_weights_(i,j+1) * v_(i,j+1) * dt / dx_;
					
					float top_phi = liquid_phi_(i, j + 1);
					if (top_phi < 0) {
						matrix_.add_to_element(index, index, local_term);
						matrix_.add_to_element(index, index + ni_, -term);
					} else {
						float theta = fraction_inside(centre_phi, top_phi);
						if (theta < 0.01) theta = 0.01;
						matrix_.add_to_element(index, index, centre_term / theta);
					}
					laplacianP_(i,j) += v_weights_(i, j+1) * ((coef_A + coef_B * (comp_rho_(i, j+1) - comp_rho_(i, j)))*comp_rho_(i,j+1)-coef_A * comp_rho_(i,j));

					// bottom neighbour
                    coef_A = t2dx2 * compute_coef_A((comp_rho_(i, j) + comp_rho_(i, j - 1)) / 2);
                    coef_B = (t2dx2 ) * compute_coef_B((comp_rho_(i, j) + comp_rho_(i, j - 1)) / 2);
					vel_term = 0.0f;//v_weights_(i, j - 1) * (t2dx2 * v_(i, j - 1) / (2.0f * dx_));
                    term = v_weights_(i, j - 1) * (coef_A - coef_B * (comp_rho_(i, j+1) - comp_rho_(i, j - 1))) + vel_term;
					local_term = centre_term - vel_term - v_weights_(i,j) * v_(i,j) * dt / dx_;
					float bot_phi = liquid_phi_(i, j - 1);
					if (bot_phi < 0) {
						matrix_.add_to_element(index, index, local_term);
						matrix_.add_to_element(index, index - ni_, -term);
					} else {
						float theta = fraction_inside(centre_phi, bot_phi);
						if (theta < 0.01) theta = 0.01;
						matrix_.add_to_element(index, index, centre_term / theta);
					}
					laplacianP_(i,j) += v_weights_(i, j) * ((coef_A + coef_B * (comp_rho_(i, j) - comp_rho_(i, j-1)))*comp_rho_(i,j-1)-coef_A * comp_rho_(i,j));

					// time dependent term
					matrix_.add_to_element(index, index, 1.0f);
					rhs_[index] += saved_comp_rho_(i,j); // 移项1/t，方程左边乘t^2
				}
		}
	});
    
    // 检查矩阵是否正定
    if (!is_symmetric(matrix_)) {
        std::cout << "Matrix is not symmetric." << std::endl;
    } else {
        std::cout << "Matrix is symmetric." << std::endl;
    }
	// Solve the system using Robert Bridson's incomplete Cholesky PCG solver
	scalar residual;
	int iterations;
	bool success = solver_.solve(matrix_, rhs_, comp_rho_solution_, residual, iterations);
	if (!success) {
		std::cout << "WARNING: Density solve failed! residual = " << residual << ", iters = " << iterations << std::endl;
       } else {
          //std::cout << "Density solve succeeded! residual = " << residual << ", iters = " << iterations << std::endl;
        }
    //std::cout << "初始密度下的压强："<< get_pressure(1.0) << std::endl; //结果为0
	comp_pressure_.set_zero();
	// Calculate pressure field using DVS equation
    parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
		if (j < u_.nj) {
			for (int i = 1; i < u_.ni - 1; ++i) {
				int index = i + j * ni_;
				comp_rho_(i,j) = comp_rho_solution_[index];
				comp_pressure_(i,j) = get_pressure(comp_rho_solution_[index]);
			}
		}
    });

    // 调试用：打印comp_pressure_的最大值
    scalar max_comp_pressure = 0.0;
    for (int j = 0; j < nj_; ++j) {
        for (int i = 0; i < ni_; ++i) {
        max_comp_pressure = std::max(max_comp_pressure, comp_pressure_(i, j));
        }
    }
    std::cout << "max_comp_pressure: " << max_comp_pressure << std::endl;

	// Apply the velocity update
	parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
		if (j < u_.nj) {
			for (int i = 1; i < u_.ni - 1; ++i) {
				int index = i + j * ni_;
				if (u_weights_(i, j) > 0) {
					if (liquid_phi_(i, j) < 0 || liquid_phi_(i - 1, j) < 0) {
						float theta = 1;
						if (liquid_phi_(i, j) >= 0 || liquid_phi_(i - 1, j) >= 0) theta = fraction_inside(liquid_phi_(i - 1, j), liquid_phi_(i, j));
						if (theta < 0.01) theta = 0.01;
                        u_(i, j) -= dt * (comp_pressure_(i, j) - comp_pressure_(i - 1, j)) / dx_ / theta;
                        // u_(i, j) -= dt * (pressure_[index] - pressure_[index - 1]) / dx_ / theta; //TODO: 是否应该再除以密度？
						u_valid_(i, j) = 1;
					} else {
						u_valid_(i, j) = 0;
					}
				} else {
					u_(i, j) = 0;
					u_valid_(i, j) = 0;
				}
			}
		}
		if (j >= 1 && j < v_.nj - 1) {
			for (int i = 0; i < v_.ni; ++i) {
				int index = i + j * ni_;
				if (v_weights_(i, j) > 0) {
					if (liquid_phi_(i, j) < 0 || liquid_phi_(i, j - 1) < 0) {
						float theta = 1;
						if (liquid_phi_(i, j) >= 0 || liquid_phi_(i, j - 1) >= 0) theta = fraction_inside(liquid_phi_(i, j - 1), liquid_phi_(i, j));
						if (theta < 0.01) theta = 0.01;
						v_(i, j) -= dt * (comp_pressure_(i, j) - comp_pressure_(i, j - 1)) / dx_ / theta;
						v_valid_(i, j) = 1;
					} else {
						v_valid_(i, j) = 0;
					}
				} else {
					v_(i, j) = 0;
					v_valid_(i, j) = 0;
				}
			}
		}
	});
}
//不使用||lapcian rho||的n+1时刻的近似
void FluidSim::solve_compressible_density_new(scalar dt) {
  int system_size = ni_ * nj_;
  if (rhs_.size() != system_size) {
    rhs_.resize(system_size);
    comp_rho_solution_.resize(system_size);
    matrix_.resize(system_size);
  }
  matrix_.zero();

  // Build the linear system for density
  parallel_for(1, nj_ - 1, [&](int j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      laplacianP_(i, j) = 0.0f;
      int index = i + ni_ * j;
      rhs_[index] = 0;
      comp_rho_solution_[index] = 0;
      float centre_phi = liquid_phi_(i, j);
      float t2dx2 = dt * dt / (dx_ * dx_);
      float coef_A = t2dx2 * (compute_coef_A(comp_rho_(i, j)) + compute_coef_A2(comp_rho_(i, j)));
      float coef_B = (t2dx2) * (compute_coef_B(comp_rho_(i, j)) + compute_coef_B2(comp_rho_(i, j)));
      float ave_rho = 0.0f;
      if (centre_phi < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
        // right neighbour
        ave_rho = (comp_rho_(i, j) + comp_rho_(i + 1, j)) * 0.5;  // 位于两个相邻单元格边界的平均密度
        coef_A = t2dx2 * (compute_coef_A(ave_rho) + compute_coef_A2(ave_rho));
        float centre_term = u_weights_(i + 1, j) * coef_A;
        // float local_term = centre_term + vel_term + u_weights_(i + 1, j) * u_(i + 1, j) * dt / dx_;
        float local_term = centre_term + u_weights_(i + 1, j) * u_(i + 1, j) * dt / dx_;

        float term = u_weights_(i + 1, j) * (coef_A /*+ coef_B * (comp_rho_(i + 1, j) - comp_rho_(i, j))*/);
        float right_phi = liquid_phi_(i + 1, j);
        // liquid_phi_小于0时，是液体单元格
        if (right_phi < 0) {
          matrix_.add_to_element(index, index, local_term);
          matrix_.add_to_element(index, index + 1, -term);
        } else {
          // 邻居单元格是气体时
          float theta = fraction_inside(centre_phi, right_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, centre_term / theta);  // doesm't consider divergence term near interface
        }

        laplacianP_(i, j) +=
            u_weights_(i + 1, j) * ((coef_A /*+ coef_B * (comp_rho_(i + 1, j) - comp_rho_(i, j))*/) * comp_rho_(i + 1, j) - coef_A * comp_rho_(i, j));

        // left neighbour
        ave_rho = (comp_rho_(i, j) + comp_rho_(i - 1, j)) * 0.5;  // 位于两个相邻单元格边界的平均密度
        coef_A = t2dx2 * (compute_coef_A(ave_rho) + compute_coef_A2(ave_rho));
     
        /*term = u_weights_(i + 1, j) * (coef_A - coef_B * (comp_rho_(i, j) - comp_rho_(i - 1, j)));*/ //TODO:原始写法的索引似乎本来就错了？
        term = u_weights_(i , j) * coef_A; 
        // local_term = centre_term - vel_term -u_weights_(i,j) * u_(i,j) * dt / dx_;
        local_term = centre_term - u_weights_(i, j) * u_(i, j) * dt / dx_;
        float left_phi = liquid_phi_(i - 1, j);
        if (left_phi < 0) {
          matrix_.add_to_element(index, index, local_term);
          matrix_.add_to_element(index, index - 1, -term);
        } else {
          float theta = fraction_inside(centre_phi, left_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, centre_term / theta);
        }
        laplacianP_(i, j) +=
            u_weights_(i, j) * ((coef_A /*+ coef_B * (comp_rho_(i, j) - comp_rho_(i-1, j))*/) * comp_rho_(i - 1, j) - coef_A * comp_rho_(i, j));

        // top neighbour
        ave_rho = (comp_rho_(i, j) + comp_rho_(i, j + 1)) * 0.5;  // 位于两个相邻单元格边界的平均密度
        coef_A = t2dx2 * (compute_coef_A(ave_rho) + compute_coef_A2(ave_rho));

        term = v_weights_(i, j + 1) * (coef_A /*+ coef_B * (comp_rho_(i, j + 1) - comp_rho_(i, j))*/);
        // local_term = centre_term + vel_term + v_weights_(i,j+1) * v_(i,j+1) * dt / dx_;
        local_term = centre_term + v_weights_(i, j + 1) * v_(i, j + 1) * dt / dx_;

        float top_phi = liquid_phi_(i, j + 1);
        if (top_phi < 0) {
          matrix_.add_to_element(index, index, local_term);
          matrix_.add_to_element(index, index + ni_, -term);
        } else {
          float theta = fraction_inside(centre_phi, top_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, centre_term / theta);
        }
        laplacianP_(i, j) +=
            v_weights_(i, j + 1) * ((coef_A /*+ coef_B * (comp_rho_(i, j+1) - comp_rho_(i, j))*/) * comp_rho_(i, j + 1) - coef_A * comp_rho_(i, j));

        // bottom neighbour
        ave_rho = (comp_rho_(i, j) + comp_rho_(i, j - 1)) * 0.5;
        coef_A = t2dx2 * (compute_coef_A(ave_rho) + compute_coef_A2(ave_rho));

        term = v_weights_(i, j) *coef_A;
        local_term = centre_term - v_weights_(i, j) * v_(i, j) * dt / dx_;
        float bot_phi = liquid_phi_(i, j - 1);
        if (bot_phi < 0) {
          matrix_.add_to_element(index, index, local_term);
          matrix_.add_to_element(index, index - ni_, -term);
        } else {
          float theta = fraction_inside(centre_phi, bot_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, centre_term / theta);
        }
        laplacianP_(i, j) +=
            v_weights_(i, j) * ((coef_A + coef_B * (comp_rho_(i, j) - comp_rho_(i, j-1))) * comp_rho_(i, j - 1) - coef_A * comp_rho_(i, j));

        // time dependent term
        matrix_.add_to_element(index, index, 1.0f);
        rhs_[index] += -saved_comp_rho_(i, j);  // 移项1/t，方程左边乘t^2

        laplacianP_(i, j) += (u_weights_(i + 1, j) * comp_pressure_(i + 1, j) + u_weights_(i, j) * comp_pressure_(i - 1, j) +
                              v_weights_(i, j + 1) * comp_pressure_(i, j + 1) + v_weights_(i, j) * comp_pressure_(i, j - 1) - 4 * comp_pressure_(i, j)) /
                             (dx_ * dx_);
        rhs_[index] += -coef_B *
                       ((u_weights_(i + 1, j) * comp_rho_(i + 1, j) - u_weights_(i, j) * comp_rho_(i - 1, j)) * 
                        (u_weights_(i + 1, j) * comp_rho_(i + 1, j) - u_weights_(i, j) * comp_rho_(i - 1, j)) +
                        (v_weights_(i, j + 1) * comp_rho_(i, j + 1) - v_weights_(i, j) * comp_rho_(i, j - 1)) *
                        (v_weights_(i, j + 1) * comp_rho_(i, j + 1) - v_weights_(i, j) * comp_rho_(i, j - 1)));
        /*rhs_[index] += t2dx2 * compute_coef_A(comp_rho_(i, j)) *
                       (comp_rho_(i + 1, j) + comp_rho_(i - 1, j) + comp_rho_(i, j + 1) + comp_rho_(i, j - 1) - 4 * comp_rho_(i, j));*/
        /*rhs_[index] +=
            t2dx2 * (comp_pressure_(i + 1, j) + comp_pressure_(i - 1, j) + comp_pressure_(i, j + 1) + comp_pressure_(i, j - 1) - 4 * comp_pressure_(i, j));*/
        

      }
    }
  });

  for (int j = 1; j < nj_ - 1; ++j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      float centre_phi = liquid_phi_(i, j);
      if (centre_phi < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
        int index = i + j * ni_;
      }
    }
  }

  // 检查矩阵是否正定
  /*if (!is_symmetric(matrix_)) {
    std::cout << "Matrix is not symmetric." << std::endl;
  } else {
    std::cout << "Matrix is symmetric." << std::endl;
  }*/

  // Solve the system using Robert Bridson's incomplete Cholesky PCG solver
  scalar residual;
  int iterations;
  bool success = solver_.solve(matrix_, rhs_, comp_rho_solution_, residual, iterations);
  if (!success) {
    std::cout << "WARNING: Density solve failed! residual = " << residual << ", iters = " << iterations << std::endl;
  } else {
    std::cout << "residual = " << residual << ", iters = " << iterations << std::endl;
  }

  comp_pressure_.set_zero();
  // Calculate pressure field using DVS equation
  scalar sum_rho = 0.0;
  parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
    if (j < u_.nj) {
      for (int i = 1; i < u_.ni - 1; ++i) {
        int index = i + j * ni_;
        //comp_rho_(i, j) = comp_rho_solution_[index];
        sum_rho += comp_rho_(i, j);
        comp_pressure_(i, j) = get_pressure(comp_rho_solution_[index]);
      }
    }
  });

  std::cout << "sum_rho in grid: " << sum_rho << std::endl;

  // 调试用：打印comp_pressure_的最大值
  scalar max_comp_pressure = 0.0;
  for (int j = 0; j < nj_; ++j) {
    for (int i = 0; i < ni_; ++i) {
      max_comp_pressure = std::max(max_comp_pressure, comp_pressure_(i, j));
    }
  }
  std::cout << "max_comp_pressure: " << max_comp_pressure << std::endl;

  // Apply the velocity update
  parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
    if (j < u_.nj) {
      for (int i = 1; i < u_.ni - 1; ++i) {
        int index = i + j * ni_;
        if (u_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i - 1, j) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i - 1, j) >= 0) theta = fraction_inside(liquid_phi_(i - 1, j), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;
            u_(i, j) -= dt * (comp_pressure_(i, j) - comp_pressure_(i - 1, j)) / dx_ / theta;
            // u_(i, j) -= dt * (pressure_[index] - pressure_[index - 1]) / dx_ / theta; //TODO: 是否应该再除以密度？
            u_valid_(i, j) = 1;
          } else {
            u_valid_(i, j) = 0;
          }
        } else {
          u_(i, j) = 0;
          u_valid_(i, j) = 0;
        }
      }
    }
    if (j >= 1 && j < v_.nj - 1) {
      for (int i = 0; i < v_.ni; ++i) {
        int index = i + j * ni_;
        if (v_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i, j - 1) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i, j - 1) >= 0) theta = fraction_inside(liquid_phi_(i, j - 1), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;
            v_(i, j) -= dt * (comp_pressure_(i, j) - comp_pressure_(i, j - 1)) / dx_ / theta;
            v_valid_(i, j) = 1;
          } else {
            v_valid_(i, j) = 0;
          }
        } else {
          v_(i, j) = 0;
          v_valid_(i, j) = 0;
        }
      }
    }
  });
}

void FluidSim::solve_compressible_density_new2(scalar dt) {
  int system_size = ni_ * nj_;
  if (rhs_.size() != system_size) {
    rhs_.resize(system_size);
    comp_rho_solution_.resize(system_size);
    matrix_.resize(system_size);
  }
  matrix_.zero();

  // Build the linear system for density
  parallel_for(1, nj_ - 1, [&](int j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      laplacianP_(i, j) = 0.0f;
      int index = i + ni_ * j;
      rhs_[index] = 0;
      comp_rho_solution_[index] = 0;
      float centre_phi = liquid_phi_(i, j);
      if (centre_phi < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
        laplacianP_(i, j) += (u_weights_(i + 1, j) * comp_pressure_(i + 1, j) + 
                              u_weights_(i, j    ) * comp_pressure_(i - 1, j) +
                              v_weights_(i, j + 1) * comp_pressure_(i, j + 1) +
                              v_weights_(i, j    ) * comp_pressure_(i, j - 1) - 4 * comp_pressure_(i, j)) /
                             (dx_*dx_);
        scalar div_vel =
            (u_weights_(i + 1, j) * u_(i + 1, j) - u_weights_(i, j) * u_(i, j) + v_weights_(i, j + 1) * v_(i, j + 1) - v_weights_(i, j) * v_(i, j)) / dx_;
        scalar rho_new = (comp_rho_(i, j) + dt * dt * laplacianP_(i, j)) / (1 + dt * div_vel);
        if (rho_new >= 1.1) {
          // std::cout << i << ", " << j << " rho_new: " << rho_new << " , diff: " << rho_new - comp_rho_(i, j) << ", laplacianP_: " << laplacianP_(i, j)
          //           << " , div_vel: " << div_vel << std::endl;
          /*std::cout << " rho_new: " << rho_new << " , diff: " << rho_new - comp_rho_(i, j) << ", laplacianP_: " << laplacianP_(i, j)
                    << " , div_vel: " << div_vel << std::endl;*/
        }
        comp_rho_(i, j) = rho_new;
      }
    }
  });

  for (int j = 1; j < nj_ - 1; ++j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      float centre_phi = liquid_phi_(i, j);
      if (centre_phi < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
        int index = i + j * ni_;
      }
    }
  }

  // Solve the system using Robert Bridson's incomplete Cholesky PCG solver
  /*scalar residual;
  int iterations;
  bool success = solver_.solve(matrix_, rhs_, comp_rho_solution_, residual, iterations);
  if (!success) {
    std::cout << "WARNING: Density solve failed! residual = " << residual << ", iters = " << iterations << std::endl;
  } else {
    std::cout << "residual = " << residual << ", iters = " << iterations << std::endl;
  }*/

  comp_pressure_.set_zero();
  // Calculate pressure field using DVS equation
  scalar sum_rho = 0.0;
  parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
    if (j < u_.nj) {
      for (int i = 1; i < u_.ni - 1; ++i) {
        sum_rho += comp_rho_(i, j);
        comp_pressure_(i, j) = get_pressure(comp_rho_(i, j)) * 0.0005; // get_pressure对于极小的密度增量都很敏感
      }
    }
  });

  std::cout << "sum_rho in grid: " << sum_rho << std::endl;

  // 调试用：打印comp_pressure_的最大值
  scalar max_comp_pressure = 0.0;
  for (int j = 0; j < nj_; ++j) {
    for (int i = 0; i < ni_; ++i) {
      max_comp_pressure = std::max(max_comp_pressure, comp_pressure_(i, j));
    }
  }
  std::cout << "max_comp_pressure: " << max_comp_pressure << std::endl;

  parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
    if (j < u_.nj) {
      for (int i = 1; i < u_.ni - 1; ++i) {
        int index = i + j * ni_;
        if (u_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i - 1, j) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i - 1, j) >= 0) theta = fraction_inside(liquid_phi_(i - 1, j), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;
            u_(i, j) -= dt * (comp_pressure_(i, j) - comp_pressure_(i - 1, j)) / dx_ / theta;
            // u_(i, j) -= dt * (pressure_[index] - pressure_[index - 1]) / dx_ / theta; //TODO: 是否应该再除以密度？
            u_valid_(i, j) = 1;
          } else {
            u_valid_(i, j) = 0;
          }
        } else {
          u_(i, j) = 0;
          u_valid_(i, j) = 0;
        }
      }
    }
    if (j >= 1 && j < v_.nj - 1) {
      for (int i = 0; i < v_.ni; ++i) {
        int index = i + j * ni_;
        if (v_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i, j - 1) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i, j - 1) >= 0) theta = fraction_inside(liquid_phi_(i, j - 1), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;
            v_(i, j) -= dt * (comp_pressure_(i, j) - comp_pressure_(i, j - 1)) / dx_ / theta;
            v_valid_(i, j) = 1;
          } else {
            v_valid_(i, j) = 0;
          }
        } else {
          v_(i, j) = 0;
          v_valid_(i, j) = 0;
        }
      }
    }
  });
}


void FluidSim::solve_temperature(scalar dt) {
  parallel_for(1, nj_ - 1, [&](int j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      scalar laplancianT = 0.0f;
      int index = i + ni_ * j;
      float centre_phi = liquid_phi_(i, j);

      scalar T0 = 473.0f;  // 底部边界加热温度
      scalar border_y = nj_ * 0.15;  // 加热边界
      // scalar border_y = 0; // 加热边界

      if (centre_phi < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {

        /*laplancianT = (u_weights_(i + 1, j) * grid_temp_(i + 1, j) + u_weights_(i, j) * grid_temp_(i - 1, j) + v_weights_(i, j + 1) * grid_temp_(i, j + 1) +
                        v_weights_(i, j) * grid_temp_(i, j - 1) - 4 * grid_temp_(i, j)) /
                       (dx_ * dx_);*/
        //定义一个匿名函数，若单元格为固体，则设置对应温度为T0
        auto set_temp = [&](scalar grid_temp, scalar border_weights) { 
            return j < border_y ? border_weights * grid_temp + (1 - border_weights) * T0 : grid_temp;
        };
        laplancianT = (set_temp(grid_temp_(i + 1, j), u_weights_(i + 1, j)) + set_temp(grid_temp_(i - 1, j), u_weights_(i, j)) +
                       set_temp(grid_temp_(i, j + 1), v_weights_(i, j + 1)) + set_temp(grid_temp_(i, j - 1), v_weights_(i, j)) - 4 * grid_temp_(i, j)) /
                      (dx_ * dx_);

        scalar D = 0.5;  // 热扩散率
        //scalar D = 0.000000145; //热扩散率
        /*if (laplancianT < 0) { //TODO: 发现相邻空气的流体并非完全绝热，laplancianT略为负数
          std::cout << "laplancianT: " << laplancianT << std::endl;
        }*/
        scalar rho_new = grid_temp_(i, j) + dt * laplancianT * D;
        grid_temp_(i, j) = rho_new;
      }
    }
  });


  for (int j = 1; j < nj_ - 1; ++j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      float centre_phi = liquid_phi_(i, j);
      if (centre_phi < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
        int index = i + j * ni_;
      }
    }
  }
}

void FluidSim::output_matrix_and_rhs_to_csv(const std::string& matrix_file, const std::string& rhs_file) {
  // Output matrix_ to CSV file
  std::ofstream matrix_out(matrix_file);
  if (matrix_out.is_open()) {
    for (int i = 0; i < matrix_.n; ++i) {
      for (int j = 0; j < matrix_.n; ++j) {
        matrix_out << matrix_(i, j);
        if (j < matrix_.n - 1) {
          matrix_out << ",";
        }
      }
      matrix_out << "\n";
    }
    matrix_out.close();
  } else {
    std::cerr << "Unable to open file " << matrix_file << std::endl;
  }

  // Output rhs_ to CSV file
  std::ofstream rhs_out(rhs_file);
  if (rhs_out.is_open()) {
    for (int i = 0; i < rhs_.size(); ++i) {
      rhs_out << rhs_[i] << "\n";
    }
    rhs_out.close();
  } else {
    std::cerr << "Unable to open file " << rhs_file << std::endl;
  }
}


bool FluidSim::is_symmetric(const robertbridson::SparseMatrix<scalar>& matrix) {
  int n = matrix.n;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (matrix(i, j) != matrix(j, i)) {
        // 计算其对应压强位置的坐标
        int x1 = i % ni_;
        int y1 = i / ni_;
        // 打印该位置的几何信息，如u_weights_和v_weights_, liquid_phi_
        std::cout << "Matrix is not symmetric with values: " << matrix(i, j) << " and " << matrix(j, i) << std::endl;
        std::cout << "(" << x1 << ", " << y1 << "): " << " ---u_weights_ : " << u_weights_(x1, y1) << ", v_weights_: " << v_weights_(x1, y1)
                  << ", liquid_phi_: " << liquid_phi_(x1, y1) << std::endl;
        int x2 = j % ni_;
        int y2 = j / ni_;
        std::cout << "(" << x2 << ", " << y2 << "): " << "--- u_weights_: " << u_weights_(x2, y2) << ", v_weights_: " << v_weights_(x2, y2)
                  << ", liquid_phi_: " << liquid_phi_(x2, y2)
                  << std::endl << std::endl;
        return false;
      }
    }
  }
  return true;
}

bool FluidSim::is_positive_definite(const robertbridson::SparseMatrix<scalar>& matrix) {
  int n = matrix.n;
  std::vector<scalar> L(n * n, 0.0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      scalar sum = 0.0;
      for (int k = 0; k < j; ++k) {
        sum += L[i * n + k] * L[j * n + k];
      }
      if (i == j) {
        L[i * n + j] = std::sqrt(matrix(i, i) - sum);
        if (L[i * n + j] <= 0) {
          return false;
        }
      } else {
        L[i * n + j] = (matrix(i, j) - sum) / L[j * n + j];
      }
    }
  }
  return true;
}