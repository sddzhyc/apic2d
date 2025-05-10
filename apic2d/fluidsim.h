
#ifndef APIC2D_FLUIDSIM_H_
#define APIC2D_FLUIDSIM_H_

// #define COMPRESSIBLE_FLUID

#define OUT_PUT

#include <chrono>
#include <memory>
#include <vector>

#include "array2.h"
#include "math_defs.h"
#include "sparse_matrix.h"
#include "pcg_solver.h"
#include "Partio.h"

class sorter;

struct Particle {

  enum PARTICLE_TYPE {
    PT_LIQUID,
    PT_AIR,
  };

  Particle(const Vector2s& x, const Vector2s& v, const scalar& radii, const scalar& density, const scalar& tempreture, const PARTICLE_TYPE type);
  Particle();
  Particle(const Particle&);

  Vector2s x_;
  Vector2s v_;
  Matrix2s c_;

  Vector2s buf0_;

  scalar dens_;
  scalar radii_;
  scalar mass_;
  scalar logJ_;
  scalar temp_;
  PARTICLE_TYPE type_;  // 粒子类型，PT_LIQUID表示流体粒子，PT_AIR表示空气粒子
};

class FluidSim {
 public:
  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = std::chrono::time_point<Clock>;
  using NeighborParticlesType = std::vector<const Particle*>;

  virtual ~FluidSim();

  scalar rho_;
  scalar rho_air_; // 气体初始密度
  scalar T_;
  int outframe_;

  enum INTEGRATOR_TYPE {
    IT_PIC,
    IT_FLIP,
    IT_RPIC,
    IT_APIC,
    IT_AFLIP,
    IT_ASFLIP,

    IT_COUNT
  };

  enum VELOCITY_ORDER {
    VO_EULER,
    VO_RA2,
    VO_RK3,
    VO_RK4,

    VO_COUNT
  };

  enum INTERPOLATION_ORDER {
    IO_LINEAR,
    IO_QUADRATIC,

    IO_COUNT
  };

  enum BOUNDARY_TYPE {
    BT_CIRCLE,
    BT_BOX,
    BT_HEXAGON,
    BT_TRIANGLE,
    BT_TORUS,
    BT_CYLINDER,

    BT_INTERSECTION,
    BT_UNION,

    BT_COUNT
  };

  struct Boundary {
    Boundary(const Vector2s& center, const Vector2s& parameter, BOUNDARY_TYPE type, bool inside);

    Boundary(Boundary* op0_, Boundary* op1_, BOUNDARY_TYPE type_);

    Vector2s center_;
    Vector2s parameter_;

    Boundary* op0_;
    Boundary* op1_;

    BOUNDARY_TYPE type_;
    scalar sign_;
  };

  struct BrushFootprint {
    BrushFootprint(const Vector2s& center, const Vector2s& vel, scalar radius) : center_(center), vel_(vel), radius_(radius) {}

    Vector2s center_;
    Vector2s vel_;
    scalar radius_;
  };

  void initialize(const Vector2s& origin, scalar width, int ni, int nj, scalar rho, bool draw_grid = true, bool draw_particles = true,
                  bool draw_velocities = true, bool draw_boundaries = true, bool print_timing = false);
  void advance(scalar dt);
  void update_boundary();
  void init_random_particles();
  void init_random_particles_2();
  void render();
  void render_boundaries(const Boundary& b);
  void render2();
  scalar compute_cfl() const;
  scalar solid_distance(const Vector2s& pos) const;
  void solve_pressure_with_rho(scalar dt);
  scalar solid_distance(const Vector2s& pos, const Boundary& b) const;

  Vector2s get_velocity_and_affine_matrix_with_order(const Vector2s& position, scalar dt, FluidSim::VELOCITY_ORDER v_order,
                                                     FluidSim::INTERPOLATION_ORDER i_order, Matrix2s* affine_matrix);
  Vector2s get_air_velocity_and_affine_matrix_with_order(const Vector2s& position, scalar dt, FluidSim::VELOCITY_ORDER v_order,
                                                         FluidSim::INTERPOLATION_ORDER i_order, Matrix2s* affine_matrix);
  Vector2s get_saved_velocity_with_order(const Vector2s& position, FluidSim::INTERPOLATION_ORDER i_order);

  Vector2s get_saved_air_velocity_with_order(const Vector2s& position, FluidSim::INTERPOLATION_ORDER i_order);

  /*! Quadratic interpolation kernels */
  Vector2s get_velocity_quadratic_impl(const Vector2s& position, const Array2s& uu, const Array2s& vv);
  Matrix2s get_affine_matrix_quadratic_impl(const Vector2s& position, const Array2s& uu, const Array2s& vv);
  Vector2s get_velocity_quadratic(const Vector2s& position);
  Matrix2s get_affine_matrix_quadratic(const Vector2s& position);
  Vector2s get_saved_velocity_quadratic(const Vector2s& position);
  Matrix2s get_saved_affine_matrix_quadratic(const Vector2s& position);

  /*! Linear interpolation kernels */
  Vector2s get_velocity(const Vector2s& position);
  Vector2s get_air_velocity(const Vector2s& position);
  Matrix2s get_affine_matrix(const Vector2s& position);
  Matrix2s get_air_affine_matrix(const Vector2s& position);
  Vector2s get_saved_velocity(const Vector2s& position);
  Vector2s get_saved_air_velocity(const Vector2s& position);
  Matrix2s get_saved_affine_matrix(const Vector2s& position);

  /*! Add particle to the system */
  void add_particle(const Particle& position);

  /*! P2G scheme */
  void map_p2g();
  void map_p2g_linear();
  void map_p2g_quadratic();

  /*! FLIP schemes */
  void map_g2p_flip_general(float dt, const scalar lagrangian_ratio, const scalar lagrangian_symplecticity, const scalar affine_stretching_ratio,
                            const scalar affine_rotational_ratio);

  void save_velocity();

  void relaxation(scalar dt);
  void set_root_boundary(const Boundary& b) { root_boundary_ = std::make_unique<Boundary>(b); }
  const std::vector<Particle>& get_particles() const { return particles_; }

  void paint_velocity(const Vector2s& brush_center, const scalar brush_radius, const Vector2s& vel);
  const Vector2s& get_origin() const { return origin_; }

  /*! Compressible fluid operations */
  void map_p2g_compressible();
  scalar get_density(const Vector2s& position);
  scalar get_temperature(const Vector2s& position);
  scalar get_saved_density(const Vector2s& position);
  


  /*! Output data bgeo */
  void OutputPointDataBgeo(const std::string& s, const int frame);
  void OutputGridDataBgeo(const std::string& s, const int frame);
  void OutputGridXDataBgeo(const std::string& s, const int frame);
  void OutputGridYDataBgeo(const std::string& s, const int frame);

  void output_matrix_and_rhs_to_csv(const std::string& matrix_file, const std::string& rhs_file);

  bool is_symmetric(const robertbridson::SparseMatrix<scalar>& matrix);

  bool is_positive_definite(const robertbridson::SparseMatrix<scalar>& matrix);

 protected:
  inline scalar circle_distance(const Vector2s& position, const Vector2s& centre, scalar radius) const { return ((position - centre).norm() - radius); }

  inline scalar box_distance(const Vector2s& position, const Vector2s& centre, const Vector2s& expand) const {
    scalar dx = fabs(position[0] - centre[0]) - expand[0];
    scalar dy = fabs(position[1] - centre[1]) - expand[1];
    scalar dax = max(dx, 0.0f);
    scalar day = max(dy, 0.0f);
    return min(max(dx, dy), 0.0f) + sqrt(dax * dax + day * day);
  }

  inline scalar hexagon_distance(const Vector2s& position, const Vector2s& centre, scalar radius) const {
    scalar dx = fabs(position[0] - centre[0]);
    scalar dy = fabs(position[1] - centre[1]);
    return max((dx * 0.866025f + dy * 0.5f), dy) - radius;
  }

  inline scalar triangle_distance(const Vector2s& position, const Vector2s& centre, scalar radius) const {
    scalar px = position[0] - centre[0];
    scalar py = position[1] - centre[1];
    scalar dx = fabs(px);
    return max(dx * 0.866025f + py * 0.5f, -py) - radius * 0.5f;
  }

  inline scalar cylinder_distance(const Vector2s& position, const Vector2s& centre, scalar theta, scalar radius) const {
    Vector2s nhat = Vector2s(cos(theta), sin(theta));
    Vector2s dx = position - centre;
    return sqrt(dx.transpose() * (Matrix2s::Identity() - nhat * nhat.transpose()) * dx) - radius;
  }

  inline scalar union_distance(const scalar& d1, const scalar& d2) const { return min(d1, d2); }

  inline scalar intersection_distance(const scalar& d1, const scalar& d2) const { return max(d1, d2); }

  inline scalar substraction_phi(const scalar& d1, const scalar& d2) const { return max(-d1, d2); }

  inline scalar torus_distance(const Vector2s& position, const Vector2s& centre, scalar radius0, scalar radius1) const {
    return max(-circle_distance(position, centre, radius0), circle_distance(position, centre, radius1));
  }

  inline void tick() {
    if (print_timing_) last_time_point_ = Clock::now();
  }

  inline void tock(const std::string& info) {
    if (print_timing_)
      std::cout << "[" << info << "] "
                << static_cast<scalar>(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - last_time_point_).count()) * 1e-6 << " ms"
                << std::endl;
  }

  Vector2s trace_rk2(const Vector2s& position, scalar dt);

  // tracer particle operations

  void particle_boundary_collision(scalar dt);

  scalar get_temperature_quadratic(const Vector2s& position, const Array2s& grid_temp);

  // fluid velocity operations
  void add_force(scalar dt);

  void compute_weights();
  void solve_pressure(scalar dt);
  void solve_pressure_with_air(scalar dt);
  void compute_liquid_distance();

  void compute_air_distance();

  void compute_merged_distance();

  scalar compute_curvature(int i, int j);

  scalar compute_face_fraction(int phi_0, int phi_1);

  scalar get_merged_phi(int i, int j);

  void constrain_velocity();

  // compressible fluid operations
  scalar compute_coef_A(const scalar& rho);
  scalar compute_coef_B(const scalar& rho);

  scalar compute_coef_A2(const scalar& rho);
  scalar compute_coef_B2(const scalar& rho);
  scalar get_pressure(const scalar& rho);
  void solve_compressible_density(scalar dt);
  void solve_compressible_density_new(scalar dt);

  void solve_compressible_density_new2(scalar dt);

  void solve_temperature(scalar dt);

 private:
  /*! Boundaries */
  std::unique_ptr<Boundary> root_boundary_;

  /*! Grid Origin */
  Vector2s origin_;

  /*! Grid dimensions */
  int ni_, nj_;
  scalar dx_;

  /*! Fluid velocity */
  Array2s u_, v_;
  Array2s temp_u_, temp_v_;
  Array2s saved_u_, saved_v_;

  /*! Air velocity */
  Array2s u_a_, v_a_;
  Array2s temp_u_a_, temp_v_a_;
  Array2s saved_u_a_, saved_v_a_;
  /*! Tracer particles */
  std::vector<Particle> particles_;

  /*! Static geometry representation */
  Array2s nodal_solid_phi_;
  Array2s liquid_phi_; // 空气与液体标识, 小于0大于-1时，是液体单元格，=1是空气单元格
  Array2s u_weights_, v_weights_;  // 0表示为固体边界，1表示流体可流动区域

  Array2s air_phi_; 
  Array2s merged_phi_;

  /*! Data arrays for extrapolation */
  Array2c valid_, old_valid_;
  Array2c u_valid_, v_valid_;

  /*! compressible fluid */
  Array2s comp_rho_;
  Array2s pre_comp_rho_;
  Array2s saved_comp_rho_;
  Array2s comp_pressure_;
  scalar a, b, R;

  Array2s grid_temp_;
  Array2s grid_rho_;

  /*! debug arrays*/
  Array2s laplacianP_;

  sorter* m_sorter_;

  /*! Solver data */
  robertbridson::PCGSolver<scalar> solver_;
  robertbridson::SparseMatrix<scalar> matrix_;
  std::vector<scalar> rhs_;
  std::vector<scalar> pressure_;
  std::vector<scalar> comp_rho_solution_;

  bool draw_grid_;
  bool draw_particles_;
  bool draw_velocities_;
  bool draw_boundaries_;
  bool print_timing_;

  TimePoint last_time_point_;

  std::vector<BrushFootprint> velocity_brush_data_;
};

#endif  // APIC2D_FLUIDSIM_H_
