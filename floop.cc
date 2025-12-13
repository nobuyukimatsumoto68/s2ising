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

constexpr int L = 1;
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
  // int k=0;
  // if(argc>=2) k = atoi(argv[1]);
  // int r=1;
  // if(argc>=3) r = atoi(argv[2]);

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  // FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  // Rotation rot;
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);

  // Fermion D(spin, true);
  Fermion D(spin);
  // FermionVector phi(dual);


  // Idx q = std::pow(2, k) + r; // Binary: 00100101
  // Idx q = 1025; // std::pow(2, N) + 0; // Binary: 00100101
  // std::unique_ptr<double> ws = new std::unique_ptr<vector_test::SimpleClass>();
  // ptr_to_unique_ptr->reset(new vector_test::SimpleClass(555));
  std::vector<double> ws( std::pow(2,N) );

  const Idx nloops = std::pow(2,N);
  // std::unique_ptr<double[]> ws = std::make_unique<double[]>(nloops);
  // const int length = numEntries;
  // int* arr = new int[length];

  // std::cout << "debug. pt1 " << std::endl;
  //   std::vector<double> tmp(nparallel, 0.0);
  // {
  //   DualLoop<N> loop(D);

  //   loop.set( 1 );
  //   std::cout << "idx = " << loop() << std::endl;
  //   std::cout << "config = " << loop.config() << std::endl;
  //   std::cout << "loop = " << std::endl;
  //   std::cout << loop.printLoops() << std::endl;
  //   std::cout << "w_prod  = " << loop.eval() << std::endl;

  //   loop.set( 2 );
  //   std::cout << "idx = " << loop() << std::endl;
  //   std::cout << "config = " << loop.config() << std::endl;
  //   std::cout << "loop = " << std::endl;
  //   std::cout << loop.printLoops() << std::endl;
  //   std::cout << "w_prod  = " << loop.eval() << std::endl;

  //   loop.set( 4 );
  //   std::cout << "idx = " << loop() << std::endl;
  //   std::cout << "config = " << loop.config() << std::endl;
  //   std::cout << "loop = " << std::endl;
  //   std::cout << loop.printLoops() << std::endl;
  //   std::cout << "w_prod  = " << loop.eval() << std::endl;

  //   loop.set( 1024 );
  //   std::cout << "idx = " << loop() << std::endl;
  //   std::cout << "config = " << loop.config() << std::endl;
  //   std::cout << "loop = " << std::endl;
  //   std::cout << loop.printLoops() << std::endl;
  //   std::cout << "w_prod  = " << loop.eval() << std::endl;


  //   loop.set( 2048 );
  //   std::cout << "idx = " << loop() << std::endl;
  //   std::cout << "config = " << loop.config() << std::endl;
  //   std::cout << "loop = " << std::endl;
  //   std::cout << loop.printLoops() << std::endl;
  //   std::cout << "w_prod  = " << loop.eval() << std::endl;

  // }
  // // return 1;

  // #ifdef _OPENMP
  // #pragma omp parallel for num_threads(nparallel) schedule(static)
  // #endif
  for(Idx q=0; q<nloops; q++){
    // for(Idx q=nloops-1; q>=0; q--){
    DualLoop<N> loop(D);
    loop.set( q );
    assert( std::abs(loop.eval()-loop.eval_prod())<1.0e-9 );
    // std::cout << "w_prod  = " << loop.eval() << std::endl;
    // std::cout << "w  = " << loop.eval_prod() << std::endl;
    ws[q] = loop.eval();
    if(ws[q]<0.0){
      std::cout << "idx = " << loop() << std::endl;
      std::cout << "config = " << loop.config() << std::endl;
      std::cout << "loop = " << std::endl;
      std::cout << loop.printLoops() << std::endl;
    }
    // tmp[omp_get_thread_num()] += loop.eval();
    // if(q%100==0) std::cout << "q = " << loop.config() << std::endl;
  }
  // for(auto elem : tmp) std::cout << elem << std::endl;
  double zloop = 0.0;
  for(Idx q=0; q<nloops; q++){
    zloop += ws[q];
    // std::cout << q << " " << ws[q] << std::endl;
  }

  double prod_mu = 1.0;
  for(auto mu : D.mus) prod_mu *= mu;

  zloop *= prod_mu / 2.0;
  std::cout << "Z (loop) = " << zloop << std::endl;

  auto matD = D.matrix();
  // std::cout << "det = " << matD.determinant() << std::endl;
  std::cout << "Z (Pf) = " << std::sqrt( matD.determinant().real() ) << std::endl;

  return 0;
}
