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
#include "loop.h"
#include "fermion.h"

constexpr int L = 2;
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
  // FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  // Rotation rot;
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);

  DualLoops<N> loops(dual);

  Fermion D(spin);
  FermionVector phi(dual);

  phi( 0, 0 ) = 1.0;
  std::cout << "# phi = " << std::endl;
  std::cout << phi.print() << std::endl;

  Eigen::MatrixXcd Dmat = D.matrix();
  // std::cout << Dmat.real() << std::endl;
  // std::cout << Dmat.imag() << std::endl;
  std::cout << Dmat.determinant() << std::endl;

  // Eigen::VectorXcd v0(2*dual.NVertices());
  // v0[0] = 1.0;
  // v0 = Dmat.inverse() * v0;
  // std::cout << v0.array().abs() << std::endl;

  // FermionVector Dphi(dual);
  // D(Dphi, phi);

  // SpinMatrix mat = spin.sigma[0];
  // std::cout << "# mat = " << std::endl;
  // std::cout << mat << std::endl;
  // Dphi.insertMultBlock(1, mat, phi, 0);

  // std::cout << "# Dphi = " << std::endl;
  // std::cout << Dphi.print() << std::endl;

  return 0;
}
