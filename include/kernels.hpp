#ifndef HEXED_KERNELS_HPP_
#define HEXED_KERNELS_HPP_

#include "Basis.hpp"
#include "Sequence.hpp"
#include "Kernel_connection.hpp"
#include "Kernel_element.hpp"
#include "Refined_face.hpp"
#include "Stopwatch_tree.hpp"
#include "Transport_model.hpp"
#include "Face_permutation.hpp"

namespace hexed
{

struct Kernel_mesh
{
  int n_dim;
  int row_size;
  const Basis& basis;
  Sequence<Kernel_connection&>& car_cons;
  Sequence<Kernel_connection&>& def_cons;
  Sequence<Kernel_element&>& car_elems;
  Sequence<Kernel_element&>& def_elems;
  Sequence<Kernel_element&>& elems;
  Sequence<Refined_face&>& ref_faces;
};

struct Kernel_options
{
  Stopwatch_tree& sw_car;
  Stopwatch_tree& sw_def;
  Stopwatch_tree& sw_pr;
  double dt;
  int i_stage;
  bool compute_residual = false;
  bool use_filter = false;
};

void compute_euler(Kernel_mesh, Kernel_options);
void compute_advection(Kernel_mesh, Kernel_options, double advect_length);
void compute_navier_stokes(Kernel_mesh, Kernel_options, std::function<void()> flux_bc,
                           Transport_model visc, Transport_model therm_cond, bool laplacian_av);
void compute_smooth_av(Kernel_mesh, Kernel_options, std::function<void()> flux_bc, double diff_time, double chebyshev_step);
void compute_fix_therm_admis(Kernel_mesh, Kernel_options, std::function<void()> flux_bc);

double max_dt_euler(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety, bool local_time);
double max_dt_navier_stokes(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety, bool local_time,
                            Transport_model visc, Transport_model therm_cond, bool laplacian_av);
double max_dt_advection(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety, bool local_time, double advect_length);
double max_dt_smooth_av(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety, bool local_time);
double max_dt_fix_therm_admis(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety, bool local_time);

void compute_prolong(Kernel_mesh, bool scale = false, bool offset = false);
void compute_restrict(Kernel_mesh, bool scale = true, bool offset = false);
std::unique_ptr<Face_permutation_dynamic> face_permutation(int n_dim, int row_size, Connection_direction, double* data);
void compute_write_face(Kernel_mesh);
void compute_write_face_advection(Kernel_mesh);
void compute_write_face_smooth_av(Kernel_mesh);

}
#endif
