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
    // std::cout << "prod = " << prod << std::endl;
  }
 
  { // this will only modify q loop
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    double prod = 0.0;
    double sum = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    for(Idx q1=0; q1<std::pow(2,N2); q1++){
      s.set( q1 );
      prod += s.weight();
      sum += std::exp(s.S());
      double eps = s.s.eps_hat(0, ising);
      sum2 += eps*std::exp(s.S());
      sum3 += s.s(0)*s.s(1) * std::exp(s.S());
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
    // std::cout << "Z (ising, prod) = " << prod/std::pow(2,N2) << std::endl;
    // std::cout << "Z (ising, sum) = " << sum/factor << std::endl;
    std::cout << "<eps> = " << sum2/sum << std::endl;
    // std::cout << "<E> = " << sum3/sum << std::endl;
    // std::cout << "dbeta_dmu = " << ising.dbeta_dmus[0] << std::endl;
    // std::cout << "tanh = " << std::tanh( ising.beta[0] ) << std::endl;
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
    // std::cout << "Z (loop) = " << zloop/2.0 << std::endl;
    // std::cout << "Z eps = " << zeps/2.0 << std::endl;
    // std::cout << "<eps> = " << zeps/zloop << std::endl;
    // std::cout << "mu = " << D.mus[0] << std::endl;
    std::cout << "<eps>/mu = " << -0.5*zeps/zloop/D.mus[0]*2.0 << std::endl;
  }

  // std::cout << "Pf. numerical deriv" << std::endl;
  {
    auto matD = D.matrix();
    double prod_mu = 1.0;
    for(auto mu : D.mus) prod_mu *= mu;

    // std::cout << "Z (Pf) = " << std::sqrt( matD.determinant().real() )/prod_mu << std::endl;

    auto inverse=matD.inverse();
    std::cout << "from inv = " << 0.5*real( inverse(0,0)+inverse(1,1) ) << std::endl;
  }
  const double eta = 1.0e-6;
  double zp, zm;
  {
    D.mus[0] += eta;
    auto matD = D.matrix();
    double prod_mu = 1.0;
    for(auto mu : D.mus) prod_mu *= mu;

    zp = std::sqrt( matD.determinant().real() )/prod_mu;
    // std::cout << "Z (Pf) = " << std::sqrt( matD.determinant().real() )/prod_mu << std::endl;
  }
  {
    D.mus[0] -= 2.0*eta;
    auto matD = D.matrix();
    double prod_mu = 1.0;
    for(auto mu : D.mus) prod_mu *= mu;

    zm = std::sqrt( matD.determinant().real() )/prod_mu;
    // std::cout << "Z (Pf) = " << std::sqrt( matD.determinant().real() )/prod_mu << std::endl;
  }
  // std::cout << "deriv = " << (std::log(zp)-std::log(zm))/(2.0*eta) << std::endl;
  std::cout << "1/mu+deriv = " << 1.0/D.mus[0]+(std::log(zp)-std::log(zm))/(2.0*eta) << std::endl;


  return 0;
}


