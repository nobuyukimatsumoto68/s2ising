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

#include <algorithm>
#include <filesystem>


#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>


#include <cstdint>
using Idx= std::uint_fast64_t;
using V3 = Eigen::Vector3d;
using V2 = Eigen::Vector2d;
using M3 = Eigen::Matrix3d;
using M2 = Eigen::Matrix2cd;
using Complex = std::complex<double>;


constexpr int nparallel = 1;

#include "sphere.h"
#include "icosahedron.h"
#include "lattices.h"
#include "groups.h"
// #include "dual_optimizer.h"

#include "connection.h"


#include "fermion.h"

// #include <omp.h>
#include "loop.h"

constexpr int L = 1;
constexpr Idx N = 10*L*L+2;





int main(int argc, char* argv[]){
// int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------
  int k=3;
  if(argc>=2) k = atoi(argv[1]);
  int r=0;
  if(argc>=3) r = atoi(argv[2]);

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);

  Fermion D(spin);

  Eigen::MatrixXcd Dmat = D.matrix();
  auto inverse=Dmat.inverse();
  std::cout << real( inverse(0,0)+inverse(1,1) ) << std::endl;


  return 0;
}
