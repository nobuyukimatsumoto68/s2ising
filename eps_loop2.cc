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

  // dual loop expansion double counts
  const double eta = 1.0e-6;

  std::vector<double> ws( std::pow(2,N) );
  std::vector<double> wsp( std::pow(2,N) );
  std::vector<double> wsm( std::pow(2,N) );
  std::vector<bool> is_have0( std::pow(2,N), false );
  const Idx nloops = std::pow(2,N);
  for(Idx q=0; q<nloops; q++){
    DualLoop<N> loop(D);
    loop.set( q );
    ws[q] = loop.eval();

    for(auto& loop : loop.loops){
      std::vector<Idx>::iterator it = std::find(loop.begin(), loop.end(), 0);
      if (it != loop.end()) {
        is_have0[q] = true;
        break;
      }
    }
  }
  D.mus[0] += eta;
  for(Idx q=0; q<nloops; q++){
    DualLoop<N> loop(D);
    loop.set( q );
    wsp[q] = loop.eval();
  }
  D.mus[0] -= 2.0*eta;
  for(Idx q=0; q<nloops; q++){
    DualLoop<N> loop(D);
    loop.set( q );
    wsm[q] = loop.eval();
  }
  D.mus[0] += eta;


  double DeltaZ=0.0;
  for(Idx q=0; q<nloops; q++){
    double diff = wsp[q]-wsm[q];
    double iszero = D.mus[0]*diff/ws[q]/(2.0*eta) + 1.0;
    if(is_have0[q]){
      assert( std::abs(iszero)<eta );
      std::cout << iszero << std::endl;
    }
    else{
      assert( std::abs(diff)<eta );
    }
    // std::cout << q << ": " << is_have0[q] << " " << diff << " " <<  << std::endl;
    DeltaZ += diff;
  }


  double zwp=0.0, zwm=0.0;
  double zweps=0.0, zw=0.0;
  for(Idx q=0; q<nloops; q++){
    zwp += wsp[q];
    zwm += wsm[q];
    zw += ws[q];
    if(!is_have0[q]) zweps += ws[q];
  }
  std::cout << "mu[0] = " << D.mus[0] << std::endl;
  std::cout << "zweps = " << zweps << std::endl;
  std::cout << "zw = " << zw << std::endl;
  std::cout << "deriv. loop = " << (std::log(zwp)-std::log(zwm))/(2.0*eta) << std::endl;
  std::cout << "zweps/zw = " << zweps/zw << std::endl;
  std::cout << "zweps/zw/D.mus[0]*2.0 = " << zweps/zw/D.mus[0]*2.0 << std::endl;

  std::cout << "sum DeltaZ/eta/zw = " << DeltaZ/(2.0*eta)/zw << std::endl;
  std::cout << "sum DeltaZ/eta/zw/D.mus[0]*2.0 = " << DeltaZ/(2.0*eta)/zw/D.mus[0]*2.0 << std::endl;

  std::cout << "DeltaZ*mu = " << -DeltaZ/(2.0*eta) * D.mus[0] << std::endl;
  std::cout << "zweps = " <<  zweps << std::endl;
  std::cout << "check 1 = " << (-DeltaZ/(2.0*eta) * D.mus[0] + zweps)/zw << std::endl;



  std::cout << "Pf. numerical deriv" << std::endl;
  {
    auto matD = D.matrix();
    double prod_mu = 1.0;
    for(auto mu : D.mus) prod_mu *= mu;

    std::cout << "Z (Pf) = " << std::sqrt( matD.determinant().real() )/prod_mu << std::endl;

    auto inverse=matD.inverse();
    std::cout << "from inv = " << real( inverse(0,0)+inverse(1,1) ) << std::endl;
  }
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
  std::cout << "deriv = " << (std::log(zp)-std::log(zm))/(2.0*eta) << std::endl;
  std::cout << "1+deriv = " << 1+(std::log(zp)-std::log(zm))/(2.0*eta) << std::endl;


  return 0;
}


