#define DEFINE_SLU 
#include <memory>
#include <vector>
#include <iostream>

#include "solver.h"
#include "domain_decomposition.h"
#include "fill_rhs.h"
#include "diag_flags.h"
#include "benchmark.hpp"

using namespace std;

void solve_SOR_reg(const double omega,
		   const double maxerror, double* phi /* result */,
		   const int NR, const int NZ, const std::vector<double>& densAll, /* input */ 
		   const int r1, const int r2, const int z1, const int z2 /* input */);

void init_SOR(int NG, int NR, int NZ, std::vector<std::vector<double>> Eps );


void init_solver(std::string file_name, int argc, char** argv)
{

  // fill espilon and phi_bound
  set_solver_boundary(file_name, Eps); 

  //init slu to obtain first guess for solution
  init_slu(file_name, Eps);

  init_SOR( global_grid.mesh_dim, global_grid.mesh_r_dim, global_grid.mesh_z_dim, Eps);
}

void calculate_potential( double* field_time,
		    int r1, int r2, int z1, int z2 )
{
#if !NTRL_ONLY
  auto begin_time = Benchmark::start();

  std::vector<double> densAll(global_grid.mesh_dim,0);

  fill_rhs(densAll.data(), qe);
  zero_array_outside_bounds_double(densAll.data());

  //on 3 or more cores omega>=1.975 apparantly leads to divergence
  //--> TODO automatic testing for best omega
  solve_SOR_reg( 1.20, 1e-9, phi, global_grid.mesh_r_dim, global_grid.mesh_z_dim, densAll, r1, r2, z1, z2);

  zero_array_outside_bounds_double(phi);
  allreduce( phi, global_grid.mesh_dim );

  *field_time += Benchmark::stop(begin_time);
#endif
}


struct CSM {
  int sz;
  int nnz;
  int maxNumInRow;
  int *ind = nullptr;
  int *ptr = nullptr;
  double *val = nullptr;

  CSM() : nnz(0), sz(0), maxNumInRow(0) {}
  CSM(int _sz, int _maxNumInRow) : CSM() { extend(_sz, _maxNumInRow); }

#if 0
  CSM(CSM&& src) noexcept : CSM() {
    (*this)=std::move(src);
  }
#endif

  ~CSM() {
    if (sz) {
      delete[] ind;
      delete[] ptr;
      delete[] val;
      ind = nullptr;
      ptr = nullptr;
      val = nullptr;
    }
  }

  CSM(CSM const &copy) {  // copy constructor
    sz = copy.sz;
    nnz = copy.nnz;
    maxNumInRow = copy.maxNumInRow;
    ind = new int[sz * maxNumInRow];
    ptr = new int[sz + 1];
    val = new double[sz * maxNumInRow];
    std::copy(&copy.ind[0], &copy.ind[copy.sz * copy.maxNumInRow], ind);
    std::copy(&copy.ptr[0], &copy.ptr[copy.sz + 1], ptr);
    std::copy(&copy.val[0], &copy.val[copy.sz * copy.maxNumInRow], val);
  }

  CSM &operator=(CSM rhs) {
    rhs.swap(*this);
    return *this;
  }

  void extend(int _sz, int _maxNumInRow) {
    ind = new int[_sz * _maxNumInRow];
    ptr = new int[_sz + 1];
    val = new double[_sz * _maxNumInRow];
    sz = _sz;
    maxNumInRow = _maxNumInRow;
  }

  void swap(CSM &s) {
    std::swap(this->sz, s.sz);
    std::swap(this->nnz, s.nnz);
    std::swap(this->maxNumInRow, s.maxNumInRow);
    std::swap(this->ind, s.ind);
    std::swap(this->ptr, s.ptr);
    std::swap(this->val, s.val);
  }
};

CSM crm;

void init_SOR(int NG, int NR, int NZ, std::vector<std::vector<double>> Eps) {

  double east, west, south, north, bingo;
  double e_iP1_jP1, e_iP1_jM1, e_iM1_jP1, e_iM1_jM1;
  int i, j, inz, i_row;

  double factorize_time = 0.0; 
  auto begin_time = Benchmark::start();

  // local crm 
  CSM lcrm(global_grid.mesh_dim,5);

  // count of non null entries
  inz = 0;

  for (i = 0; i < global_grid.mesh_r_dim; ++i)
    for (j = 0; j < global_grid.mesh_z_dim; ++j) {
      i_row = j + i * global_grid.mesh_z_dim;

      if (Eps[i][j] == -100.)  // Fixed potential on metal
      {
	lcrm.ptr[i_row] = inz;
	lcrm.ind[inz] = i_row;  // diagonal element
	lcrm.val[inz] = 1.;
	inz++;
      } else if (Eps[i][j] == -200.)  // Fixed gradient
      {
	lcrm.ptr[i_row] = inz;
	lcrm.ind[inz] = i_row - 1;  // c - west
	lcrm.val[inz] = -1.;

	lcrm.ind[inz + 1] = i_row;  // e - bingo
	lcrm.val[inz + 1] = 1.;
	inz += 2;
      }

      else if (i == 0)  // Axis
      {
	e_iP1_jM1 = Eps[i + 1][j - 1];
	if (e_iP1_jM1 <= 0.) e_iP1_jM1 = Eps[i][j];

	e_iP1_jP1 = Eps[i + 1][j + 1];
	if (e_iP1_jP1 <= 0.) e_iP1_jP1 = Eps[i][j];

	west = e_iP1_jM1;  // a,b,c,d,e    *= dz*dz

	east = e_iP1_jP1;

	north = 2. * (e_iP1_jM1 + e_iP1_jP1);

	bingo = -3. * (e_iP1_jM1 + e_iP1_jP1);

	lcrm.ptr[i_row] = inz;

	lcrm.ind[inz] = i_row - 1;  // c - west
	lcrm.val[inz] = west;

	lcrm.ind[inz + 1] = i_row;  // e - bingo
	lcrm.val[inz + 1] = bingo;

	lcrm.ind[inz + 2] = i_row + 1;  // d - east
	lcrm.val[inz + 2] = east;

	lcrm.ind[inz + 3] = i_row + global_grid.mesh_z_dim;  // b - north
	lcrm.val[inz + 3] = north;

	inz += 4;

      }

      else  // the rest: channell, dielectric, boundaries, etc
      {
	e_iP1_jM1 = Eps[i + 1][j - 1];
	if (e_iP1_jM1 <= 0.) e_iP1_jM1 = Eps[i][j];
	if (e_iP1_jM1 == 0) {
	  if (Eps[i + 1][j] > 0)
	    e_iP1_jM1 = Eps[i + 1][j];
	  else
	    e_iP1_jM1 = Eps[i][j - 1];
	}

	e_iP1_jP1 = Eps[i + 1][j + 1];
	if (e_iP1_jP1 <= 0.) e_iP1_jP1 = Eps[i][j];
	if (e_iP1_jP1 == 0) {
	  if (Eps[(i + 1)][j] > 0)
	    e_iP1_jP1 = Eps[(i + 1)][j];
	  else
	    e_iP1_jP1 = Eps[i][j + 1];
	}

	e_iM1_jM1 = Eps[(i - 1)][j - 1];
	if (e_iM1_jM1 <= 0.) e_iM1_jM1 = Eps[i][j];
	if (e_iM1_jM1 == 0) {
	  if (Eps[(i - 1)][j] > 0)
	    e_iM1_jM1 = Eps[(i - 1)][j];
	  else
	    e_iM1_jM1 = Eps[i][j - 1];
	}

	e_iM1_jP1 = Eps[(i - 1)][j + 1];
	if (e_iM1_jP1 <= 0.) e_iM1_jP1 = Eps[i][j];
	if (e_iM1_jP1 == 0) {
	  if (Eps[(i - 1)][j] > 0)
	    e_iM1_jP1 = Eps[(i - 1)][j];
	  else
	    e_iM1_jP1 = Eps[i][j + 1];
	}

	north =
	    0.25 * (e_iP1_jM1 + e_iP1_jP1) / i + 0.5 * (e_iP1_jM1 + e_iP1_jP1);

	south =
	    0.5 * (e_iM1_jM1 + e_iM1_jP1) - 0.25 * (e_iM1_jM1 + e_iM1_jP1) / i;

	bingo = 0.25 * (e_iM1_jM1 + e_iM1_jP1 - e_iP1_jM1 - e_iP1_jP1) / i -
		(e_iP1_jM1 + e_iP1_jP1 + e_iM1_jM1 + e_iM1_jP1);

	west = 0.5 * (e_iP1_jM1 + e_iM1_jM1);

	east = 0.5 * (e_iP1_jP1 + e_iM1_jP1);

	lcrm.ptr[i_row] = inz;

	lcrm.ind[inz] = i_row - global_grid.mesh_z_dim;  // a - south
	lcrm.val[inz] = south;

	lcrm.ind[inz + 1] = i_row - 1;  // c - west
	lcrm.val[inz + 1] = west;

	lcrm.ind[inz + 2] = i_row;  // e - bingo
	lcrm.val[inz + 2] = bingo;

	lcrm.ind[inz + 3] = i_row + 1;  // d - east
	lcrm.val[inz + 3] = east;

	lcrm.ind[inz + 4] = i_row + global_grid.mesh_z_dim;  // b - north
	lcrm.val[inz + 4] = north;

	inz += 5;
      }
    }

  lcrm.nnz = inz;
  lcrm.ptr[global_grid.mesh_dim] = lcrm.nnz;

  crm = lcrm;
#if DEBUG_MATRIX_CREATION
    //printf("sz %d\n",sz);
    printf("nnz %d \n", crm.nnz);
    printf("maxNumInRow %d\n",crm.maxNumInRow);
    for (int i = 0; i < global_grid.mesh_dim+1; ++i){
      printf("%d ", crm.ptr[i]);
    }
    printf("\n");

    for (int i = 0; i < global_grid.mesh_dim*5; ++i){
      printf("%d ", crm.ind[i]);
    }
    printf("\n");

    for (int i = 0; i < global_grid.mesh_dim*5; ++i){
      printf("%f ", crm.val[i]);
    }
    printf("\n");
#endif

  factorize_time += Benchmark::stop(begin_time);

}


void solve_SOR_reg(const double omega,
		   const double maxerror, double* phi /* result */,
		   const int NR, const int NZ, const std::vector<double>& densAll, /* input */ 
		   const int r1, const int r2, const int z1, const int z2 /* input */){

  double error = 0;
  int mpi_size = get_size();
  int sor_steps = 0;
  double diag_element = 0;

  do {
    error = 0;
    exchange_array_bounds_double(phi);

    auto doProcessing = [&](int i){ 
      double phi0 = densAll[i];
      for (auto j = crm.ptr[i]; j < crm.ptr[i + 1]; ++j) {
	if (crm.ind[j] != i){
	  phi0 -= crm.val[j] * phi[crm.ind[j]];
	}else{
	  diag_element = crm.val[j];
	}
      }
      phi0 /= diag_element;

      if (fabs(phi0 - phi[i]) > error) {
	error = fabs(phi0 - phi[i]);
      }

      phi[i] = (1.0 - omega) * phi[i] + omega * phi0;

      // constant bc
      if ((crm.ind[crm.ptr[i]] == i) &&
	  (crm.val[crm.ptr[i]] == 1) &&
	  (crm.val[crm.ptr[i] + 1] != -1)) {
	phi[i] = densAll[i];
      }
    };

    // loop over subdomain
    for (int r = r1; r < r2; ++r){
      for (int z = z1; z < z2; ++z){
	int i = r * global_grid.mesh_z_dim + z;
	doProcessing(i);
      }
    }

    ++sor_steps;

    //TODO better way for termination condition
    //testing purposes
    allreduce(&error,1);
    error/=mpi_size;

#if 0
    if (!(sor_steps % 50)){
      std::cout << std::endl << "error " << error << std::endl ;
    }
#endif
  }while (error > maxerror);
#if 1
  std::cout << std::endl << "error " << error << std::endl << "steps " << sor_steps << std::endl;
#endif
}


void destroy_solver()
{
#if !NTRL_ONLY
  destroy_slu();
#endif
}
