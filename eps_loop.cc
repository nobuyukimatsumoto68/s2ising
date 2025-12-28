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

#include <random>


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


const int nparallel = 1;

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

constexpr int L = 1; // 4
constexpr Idx N = 10*L*L+2;
constexpr Idx N2 = 20*L*L;

#ifndef _OPENMP
int omp_get_thread_num(){ return 0; }
#endif



int main(int argc, char* argv[]){
// int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  omp_set_num_threads(nparallel);

  int k=0;
  if(argc>=2) k = atoi(argv[1]);
  int r=0;
  if(argc>=3) r = atoi(argv[2]);


  // --------------------

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);
  // Fermion D(spin, true);
  Fermion D(spin);

  Idx q = std::pow(2, k) + r; // Binary: 00100101

  {
    DualLoop<N> loop(D);
    Ising ising(loop);

    double prod = 1.0;
    for(Idx il=0; il<2*dual.NLinks(); il++){
      const DirectedLink ell = dual.directedlinkidx2DirectedLink( il );
      const Idx if1 = ell.first;
      const int df = ell.second;
      FaceCoords f1 = dual.idx2FaceCoords( if1 );
      if(f1[3]!=XZ) continue;
      const double beta = ising.betas.at( ell );
      FaceCoords f2;
      dual.shift(f2, f1, df);
      const Idx if2 = dual.idx(f2);
      int sign = 1;
      if( if1==k || if2==k ) sign = -1;
      prod *= 1.0 + sign*std::tanh(beta);
    }
    std::cout << "prod = " << prod << std::endl;
  }
 
  { // this will only modify q loop
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    double prod = 0.0;
    double sum = 0.0;
    for(Idx q1=0; q1<std::pow(2,N2); q1++){
      s.set( q1 );
      prod += s.weight();
      sum += std::exp(s.S());
      // std::cout << ising.eval_loop( q ) << std::endl;
      // break;
    }

    double factor = 1.0;
    for(Idx il=0; il<2*dual.NLinks(); il++){
      const DirectedLink ell = dual.directedlinkidx2DirectedLink( il );
      const Idx if1 = ell.first;
      FaceCoords f = dual.idx2FaceCoords( ell.first );

      FaceCoords f2;
      dual.shift( f2, dual.idx2FaceCoords(ell.first), ell.second );
      const Idx if2 = dual.idx(f2);
      if(if1>=if2) continue;

      const double beta = ising.betas.at( ell );
      factor *= std::cosh(beta);
    }
    factor *= std::pow(2,N2);

    // std::cout << "debug. prod = " << prod << std::endl;
    // std::cout << "debug. sum = " << sum << std::endl;
    std::cout << "Z (ising, prod) = " << prod/std::pow(2,N2) << std::endl;
    std::cout << "Z (ising, sum) = " << sum/factor << std::endl;
  }
  // {
  //   double sum = 0.0;
  //   for(Idx q=0; q<std::pow(2,N); q++) sum += ising.eval_loop( q );
  //   std::cout << "Z (ising, loopsum) = " << sum/2.0 << std::endl;
  // }


  // dual loop expansion double counts
  {
    std::vector<double> ws( std::pow(2,N) );
    const Idx nloops = std::pow(2,N);
    double zloop = 0.0;
    double zeps = 0.0;

    for(Idx q=0; q<nloops; q++){
      DualLoop<N> loop(D);
      loop.set( q );
      double w = loop.eval();
      zloop += w;

      bool is_have0 = false;
      for(auto& loop : loop.loops){
        std::vector<Idx>::iterator it = std::find(loop.begin(), loop.end(), 0);
        if (it != loop.end()) {
          is_have0 = true;
          break;
        }
      }
      if(!is_have0) zeps += w;
    }
    std::cout << "Z (loop) = " << zloop/2.0 << std::endl;
    std::cout << "Z eps = " << zeps/2.0 << std::endl;
    std::cout << "<eps> = " << zeps/zloop << std::endl;
    std::cout << "mu = " << D.mus[0] << std::endl;

    // double prod_mu = 1.0;
    // for(auto mu : D.mus) prod_mu *= mu;

    // zloop *= prod_mu / 2.0;
    // std::cout << "Z (multed) = " << zloop << std::endl;
  }

  {
    auto matD = D.matrix();
    std::cout << "Z (Pf) = " << std::sqrt( matD.determinant().real() ) << std::endl;
  }

  // {
  //   DualLoop<N> loop(D);
  //   Ising ising(loop);
  //   Spin<N, N2> s(ising);
  //   s.random();
  //   std::cout << "s = " << s.print() << std::endl;
  //   s.ckpoint("check.dat");
  //   s.set0();
  //   std::cout << "s0 = " << s.print() << std::endl;
  //   s.read("check.dat");
  //   std::cout << "s = " << s.print() << std::endl;
  // }

  return 0;
}






  // {
  //   DualLoop<N> loops(D);
  //   loops.set( q );
  //   std::cout << loops.printLoops() << std::endl;

  //   for(const auto& loop : loops){
  //     for(Idx i=0; i<loop.size(); i++){
  //       const Idx if1 = loop[i];
  //       const Idx if2 = loop[(i+1)%loop.size()];
  //       const int df = dual.getDirection.at( Link{if1, if2} );
  //       D.kappas[ dual.directedlinkidx(if1, df) ] = 1.0;
  //     }
  //   }
  // }

  // {
  //   DualLoop<N> loop(D);
  //   Ising ising(loop);

  //   double prod = 1.0;
  //   for(Idx il=0; il<2*dual.NLinks(); il++){
  //     const DirectedLink ell = dual.directedlinkidx2DirectedLink( il );
  //     const Idx if1 = ell.first;
  //     const int df = ell.second;
  //     FaceCoords f1 = dual.idx2FaceCoords( if1 );
  //     if(f1[3]!=XZ) continue;
  //     const double beta = ising.betas.at( ell );
  //     FaceCoords f2;
  //     dual.shift(f2, f1, df);
  //     const Idx if2 = dual.idx(f2);
  //     int sign = 1;
  //     if( if1==k || if2==k ) sign = -1;
  //     prod *= 1.0 + sign*std::tanh(beta);
  //   }
  //   std::cout << "prod = " << prod << std::endl;
  // }
  // //

  // // {
  // //   DualLoop<N> loop(D);
  // //   Ising ising(loop);
  // //   Spin<N, N2> s(ising);

  // //   s.set_from_loop( 5 );
  // //   std::cout << "dual :" << s.ising.loops.config() << std::endl;
  // //   std::cout << "loops :" << s.ising.loops.printLoops() << std::endl;
  // //   std::cout << "spin :" << s.print() << std::endl;
  // // }

  // {
  //   DualLoop<N> loop(D);
  //   Ising ising(loop);
  //   Spin<N, N2> s(ising);

  //   // double sum = 0.0;
  //   double prod = 0.0;
  //   for(Idx q1=0; q1<std::pow(2,N2); q1++){
  //   // for(Idx q=0; q<std::pow(2,0); q++){
  //   // {
  //     // Idx q = std::pow(2,k);
  //     s.set( q1 );
  //     // std::cout << s.print() << std::endl;
  //     // std::cout << s.S() << std::endl;
  //     // sum += std::exp( s.S() );

  //     prod += s.weight();
  //     // prod += s[0]*s[1] * s.weight();

  //     // double w1 = s.weight();
  //     // s.flip();
  //     // double w2 = s.weight();
  //     // prod += w1+w2;
  //     // if( std::abs(w1-w2)<1.0e-14 ) counter++;
  //   }
  //   // std::cout << "debug. sum = " << sum << std::endl;
  //   std::cout << "debug. prod = " << prod << std::endl;
  //   std::cout << "Z (ising, prod) = " << prod/std::pow(2,N2) << std::endl;

  //   // double factor = 1.0;
  //   // for(Idx il=0; il<2*dual.NLinks(); il++){
  //   //   const DirectedLink ell = dual.linkidx2DirectedLink( il );
  //   //   FaceCoords f = dual.idx2FaceCoords( ell.first );
  //   //   if(f[3]!=XZ) continue;
  //   //   const double beta = ising.betas.at( ell );
  //   //   factor *= std::cosh(beta);
  //   // }
  //   // factor *= std::pow(2,N2);
  //   // std::cout << "Z (ising) = " << sum/factor << std::endl;
  // }




















  // // dual loop expansion double counts
  // {
  //   std::vector<double> ws( std::pow(2,N) );

  //   const Idx nloops = std::pow(2,N);
  //   // for(Idx q=0; q<nloops; q++){
  //   {
  //     // Idx q=0;
  //     DualLoop<N> loop(D);

  //     loop.set( q );
  //     if(std::abs(loop.eval()-loop.eval_prod())>1.0e-9){
  //       std::cout << loop.eval()
  //                 << " "
  //                 << loop.eval_()
  //                 << " "
  //                 << loop.eval_prod()
  //                 << std::endl;
  //     }
  //     // ws[q] = loop.eval();
  //     ws[q] = loop.eval_();
  //     // ws[q] = loop.eval_prod();
  //     // if(ws[q]<0.0){
  //     //   std::cout << "idx = " << loop() << std::endl;
  //     //   std::cout << "config = " << loop.config() << std::endl;
  //     //   std::cout << "loop = " << std::endl;
  //     //   std::cout << loop.printLoops() << std::endl;
  //     // }
  //   }

  //   double zloop = 0.0;
  //   for(Idx q1=0; q1<nloops; q1++){
  //     zloop += ws[q1];
  //     // std::cout << q << " " << ws[q] << std::endl;
  //   }
  //   std::cout << "Z (loop) = " << zloop/2.0 << std::endl;

  //   // double prod_mu = 1.0;
  //   // for(auto mu : D.mus) prod_mu *= mu;

  //   // zloop *= prod_mu / 2.0;
  //   // std::cout << "Z (multed) = " << zloop << std::endl;
  // }

  // {
  //   auto matD = D.matrix();
  //   std::cout << "Z (Pf) = " << std::sqrt( matD.determinant().real() ) << std::endl;
  // }
