#define DEFINE_SLU

#include <time.h>
#include <stdio.h>

#include "solver.h"
#include "var_slu.h"
#include "fill_rhs.h"
#include "benchmark.hpp"

void init_solver(std::string file_name, int argc, char** argv) 
{

  // fill espilon and phi_bound
  set_solver_boundary(file_name, Eps); 
#if !NTRL_ONLY
  init_slu(file_name, Eps);
#endif
}

void calculate_potential( double* field_time, int r1, int r2, int z1, int z2) 
{
#if !NTRL_ONLY
  auto benchstart = Benchmark::start();
  int mpi_rank = get_rank(); 
  trans_t trans = NOTRANS;  // solving without tranposing
    fill_rhs(phi, qe);

  // solve step is just done by process 0
  if ( mpi_rank == 0 ) 
  {

    /* Solve the system, using L and U overwriting B with X. */
    dgstrs(trans, L, U, perm_c, perm_r, B, &stat_SLU, &info);
  }
  broadcast( phi, global_grid.mesh_dim );
  *field_time += Benchmark::stop(benchstart);

#endif
}

void destroy_solver()
{
#if !NTRL_ONLY
  destroy_slu();
#endif
}

