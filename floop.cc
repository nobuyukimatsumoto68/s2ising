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



#include "sphere.h"
#include "icosahedron.h"
#include "lattices.h"
#include "groups.h"
// #include "dual_optimizer.h"

#include "connection.h"



// #include <omp.h>
#include "fermion.h"

#include "loop.h"

constexpr int L = 2;
constexpr Idx N = 10*L*L+2;





int main(int argc, char* argv[]){
// int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------
  int k=0;
  if(argc>=2) k = atoi(argv[1]);
  int r=1;
  if(argc>=3) r = atoi(argv[2]);

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  // FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  // Rotation rot;
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);

  Fermion D(spin);
  DualLoops<N> loops(D);
  // FermionVector phi(dual);

  Idx q = std::pow(2, k) + r; // Binary: 00100101
  loops.set( q );

  std::cout << "idx = " << loops() << std::endl;
  std::cout << "config = " << loops.config() << std::endl;
  std::cout << "loops = " << std::endl;
  std::cout << loops.printLoops() << std::endl;

  std::cout << "w_prod  = " << loops.eval() << std::endl;
  std::cout << "w_prod  = " << loops.eval_prod() << std::endl;


  return 0;
}
