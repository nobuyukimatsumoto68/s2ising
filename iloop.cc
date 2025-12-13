#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>
#include <complex>

// #include <hdf5.h>
// #include "highfive/H5File.hpp"
// #include <highfive/H5File.hpp>
// #include <highfive/H5DataSet.hpp>
// #include <highfive/H5DataSpace.hpp>
#include <omp.h>

#include <algorithm>
#include <filesystem>


#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>


#include <cstdint>
// using Idx= std::uint_fast64_t;
using Idx= std::int_fast64_t;
using V3 = Eigen::Vector3d;
using V2 = Eigen::Vector2d;
using M3 = Eigen::Matrix3d;
using M2 = Eigen::Matrix2cd;
using Complex = std::complex<double>;



#include "sphere.h"
#include "icosahedron.h"
#include "lattices.h"
#include "groups.h"
// #include "dual_optimizer.h"

#include "connection.h"



// #include <omp.h>
#include "fermion.h"

#include <memory>
#include "loop.h"
#include "ising.h"

constexpr int L = 2;
constexpr Idx N = 10*L*L+2;

#ifndef _OPENMP
int omp_get_thread_num(){ return 0; }
#endif



int main(int argc, char* argv[]){
// int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  const int nparallel = 12;
  omp_set_num_threads(nparallel);

  // --------------------
  int k=0;
  if(argc>=2) k = atoi(argv[1]);
  int r=1;
  if(argc>=3) r = atoi(argv[2]);

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);
  Fermion D(spin);

  // const Idx nloops = std::pow(2,N);
  DualLoop<N> loop(D);

  Ising ising(loop);

  Idx q = std::pow(2, k) + r; // Binary: 00100101
  loop.set( q );
  std::cout << "config = " << loop.config() << std::endl;
  std::cout << ising.eval_loop( q ) << std::endl;
  std::cout << loop.eval() << std::endl;

  return 0;
}
