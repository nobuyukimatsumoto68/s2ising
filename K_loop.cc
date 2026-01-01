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
  const Idx ell = 0;
  const DirectedLink link = dual.directedlinkidx2DirectedLink(ell);
  const Idx if1 = link.first;
  const int rho = link.second;
  const FaceCoords f1 = dual.idx2FaceCoords(if1);
  FaceCoords f2;
  dual.shift( f2, f1, rho );
  const Idx if2 = dual.idx(f2);
  const int sigma = dual.getDirection.at({if2,if1});
  const Idx em = dual.directedlinkidx(if2,sigma);



  std::vector<double> ws( std::pow(2,N) );
  std::vector<double> wsp( std::pow(2,N) );
  std::vector<double> wsm( std::pow(2,N) );
  std::vector<bool> is_have_ell( std::pow(2,N), false );
  const Idx nloops = std::pow(2,N);
  for(Idx q=0; q<nloops; q++){
    DualLoop<N> loop(D);
    loop.set( q );
    ws[q] = loop.eval();

    for(auto& loop : loop.loops){
      std::vector<Idx>::iterator it1 = std::find(loop.begin(), loop.end(), if1);
      std::vector<Idx>::iterator it2 = std::find(loop.begin(), loop.end(), if2);
      if ( it1 != loop.end() && it2 != loop.end() ) {
        is_have_ell[q] = true;
        break;
      }
    }
  }
  D.kappas[ell] += eta;
  D.kappas[em] += eta;
  for(Idx q=0; q<nloops; q++){
    DualLoop<N> loop(D);
    loop.set( q );
    wsp[q] = loop.eval();
  }
  D.kappas[ell] -= 2.0*eta;
  D.kappas[em] -= 2.0*eta;
  for(Idx q=0; q<nloops; q++){
    DualLoop<N> loop(D);
    loop.set( q );
    wsm[q] = loop.eval();
  }
  D.kappas[ell] += eta;
  D.kappas[em] += eta;


  // double DeltaZ=0.0;
  // for(Idx q=0; q<nloops; q++){
  //   double diff = wsp[q]-wsm[q];
  //   double iszero = D.mus[0]*diff/ws[q]/(2.0*eta) + 1.0;
  //   if(is_have_ell[q]){
  //     // assert( std::abs(iszero)<eta );
  //     std::cout << iszero << std::endl;
  //   }
  //   else{
  //     assert( std::abs(diff)<eta );
  //   }
  //   // std::cout << q << ": " << is_have0[q] << " " << diff << " " <<  << std::endl;
  //   DeltaZ += diff;
  // }


  double zwp=0.0, zwm=0.0;
  double zwK=0.0, zw=0.0;
  for(Idx q=0; q<nloops; q++){
    zwp += wsp[q];
    zwm += wsm[q];
    zw += ws[q];
    if(!is_have_ell[q]) zwK += ws[q];
  }
  // std::cout << "mu[0] = " << D.mus[0] << std::endl;
  std::cout << "zwK = " << zwK << std::endl;
  // std::cout << "zw = " << zw << std::endl;
  std::cout << "deriv. loop = " << (std::log(zwp)-std::log(zwm))/(2.0*eta) << std::endl;
  // std::cout << "zweps/zw = " << zweps/zw << std::endl;
  // std::cout << "zweps/zw/D.mus[0]*2.0 = " << zweps/zw/D.mus[0]*2.0 << std::endl;

  // std::cout << "sum DeltaZ/eta/zw = " << DeltaZ/(2.0*eta)/zw << std::endl;
  // std::cout << "sum DeltaZ/eta/zw/D.mus[0]*2.0 = " << DeltaZ/(2.0*eta)/zw/D.mus[0]*2.0 << std::endl;

  // std::cout << "DeltaZ*mu = " << -DeltaZ/(2.0*eta) * D.mus[0] << std::endl;
  // std::cout << "zweps = " <<  zweps << std::endl;
  // std::cout << "check 1 = " << (-DeltaZ/(2.0*eta) * D.mus[0] + zweps)/zw << std::endl;








  { // this will only modify q loop
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    double prod = 0.0;
    double sum = 0.0;
    double sum11 = 0.0;
    double sum22 = 0.0;
    double sum12 = 0.0;
    double sumT = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    for(Idx q1=0; q1<std::pow(2,N2); q1++){
      s.set( q1 );
      prod += s.weight();
      sum += std::exp(s.S());

      double K = s.s.K_hat(if1, rho, ising);
      double eps_hat = s.s.eps_hat(if1, ising);
      // double eps_hat2 = s.s.eps_hat(if2, ising);
      // std::cout << "K, eps " << K << " " << eps_hat << " " << eps_hat2 << std::endl;
      sum2 += K*std::exp(s.S());
      sum3 += eps_hat*std::exp(s.S());

      double T11 = s.s.T11(if1, ising);
      double T22 = s.s.T22(if1, ising);
      double T12 = s.s.T12(if1, ising);
      double trT = s.s.trT(if1, ising);
      sum11 += T11*std::exp(s.S());
      sum22 += T22*std::exp(s.S());
      sum12 += T12*std::exp(s.S());
      sumT += trT*std::exp(s.S());
      // sum3 += s.s(0)*s.s(1) * std::exp(s.S());
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

    std::cout << "<T11> = " << sum11/sum << std::endl;
    std::cout << "<T22> = " << sum22/sum << std::endl;
    std::cout << "<T12> = " << sum12/sum << std::endl;
    std::cout << "<trT> = " << sumT/sum << std::endl;
    std::cout << "<K> = " << sum2/sum << std::endl;
    std::cout << "<eps> = " << sum3/sum << std::endl;
  }



  std::cout << "Ising. numerical deriv" << std::endl;
  double isump=0.0, isumm=0.0;
  { // this will only modify q loop
    D.kappas[ell] += eta;
    D.kappas[em] += eta;

    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    for(Idx q1=0; q1<std::pow(2,N2); q1++){
      s.set( q1 );
      isump += std::exp(s.S());
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
    isump /= factor;
  }
  { // this will only modify q loop
    D.kappas[ell] -= 2.0*eta;
    D.kappas[em] -= 2.0*eta;

    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    for(Idx q1=0; q1<std::pow(2,N2); q1++){
      s.set( q1 );
      isumm += std::exp(s.S());
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
    isumm /= factor;
  }
  D.kappas[ell] += eta;
  D.kappas[em] += eta;
  std::cout << "deriv = " << (std::log(isump)-std::log(isumm))/(2.0*eta) << std::endl;




  std::cout << "Pf. numerical deriv" << std::endl;
  {
    auto matD = D.matrix();
    double prod_mu = 1.0;
    for(auto mu : D.mus) prod_mu *= mu;

    std::cout << "Z (Pf) = " << std::sqrt( matD.determinant().real() )/prod_mu << std::endl;

    auto inverse=matD.inverse();
    SpinMatrix P = spin.P(f1,rho);
    SpinMatrix Dinv;
    Dinv << inverse(2*if2, 2*if1)
      , inverse(2*if2, 2*if1+1)
      , inverse(2*if2+1, 2*if1)
      , inverse(2*if2+1, 2*if1+1);
    std::cout << "from inv = " << -(P*Dinv).trace().real() << std::endl;
    // std::cout << "from inv = " << real( inverse(2*if1, 2*if2+1) ) << std::endl;
    // std::cout << "from inv = " << real( inverse(2*if1, 2*if2) ) << std::endl;
    // std::cout << "from inv = " << real( inverse(2*if1, 2*if2) ) << std::endl;
  }
  double zp, zm;
  {
    D.kappas[ell] += eta;
    D.kappas[em] += eta;
    auto matD = D.matrix();
    double prod_mu = 1.0;
    for(auto mu : D.mus) prod_mu *= mu;

    zp = std::sqrt( matD.determinant().real() )/prod_mu;
    // std::cout << "Z (Pf) = " << std::sqrt( matD.determinant().real() )/prod_mu << std::endl;
  }
  {
    D.kappas[ell] -= 2.0*eta;
    D.kappas[em] -= 2.0*eta;
    auto matD = D.matrix();
    double prod_mu = 1.0;
    for(auto mu : D.mus) prod_mu *= mu;

    zm = std::sqrt( matD.determinant().real() )/prod_mu;
    // std::cout << "Z (Pf) = " << std::sqrt( matD.determinant().real() )/prod_mu << std::endl;
  }
  std::cout << "deriv = " << (std::log(zp)-std::log(zm))/(2.0*eta) << std::endl;
  // std::cout << "1+deriv = " << 1+(std::log(zp)-std::log(zm))/(2.0*eta) << std::endl;


  {
    DualLoop<N> loop(D);
    Ising ising(loop);

    std::vector<double> alpha(3); // [mu]; y+nu,y,y+rho; (mu,nu,rho)
    std::vector<double> dalpha(3); // [mu]; y+nu,y,y+rho; (mu,nu,rho)
    for(int mu=0; mu<3; mu++) {
      alpha[mu] = ising.spin.alphas[ising.dual.directedlinkidx(if1, mu)];
    }
    for(int mu=0; mu<3; mu++) {
      const double al_nu = ising.spin.alphas[ising.dual.directedlinkidx(if1, (mu+1)%3)];
      const double al_rho = ising.spin.alphas[ising.dual.directedlinkidx(if1, (mu+2)%3)];
      dalpha[mu] = M_PI+al_nu - al_rho;
    }

    Eigen::Matrix3d mm, mminv;
    mm << std::pow( std::cos(alpha[0]), 2 ), std::sin( 2.0*alpha[0] ), std::pow( std::sin(alpha[0]), 2 )
      , std::pow( std::cos(alpha[1]), 2 ), std::sin( 2.0*alpha[1] ), std::pow( std::sin(alpha[1]), 2 )
      , std::pow( std::cos(alpha[2]), 2 ), std::sin( 2.0*alpha[2] ), std::pow( std::sin(alpha[2]), 2 );

    mminv << -std::sin(alpha[1])*std::sin(alpha[2])/std::sin(dalpha[1])/std::sin(dalpha[2])
      , -std::sin(alpha[2])*std::sin(alpha[0])/std::sin(dalpha[2])/std::sin(dalpha[0])
      , -std::sin(alpha[0])*std::sin(alpha[1])/std::sin(dalpha[0])/std::sin(dalpha[1])
      , 0.5*std::sin(alpha[1]+alpha[2])/std::sin(dalpha[1])/std::sin(dalpha[2])
      , 0.5*std::sin(alpha[2]+alpha[0])/std::sin(dalpha[2])/std::sin(dalpha[0])
      , 0.5*std::sin(alpha[0]+alpha[1])/std::sin(dalpha[0])/std::sin(dalpha[1])
      , -std::cos(alpha[1])*std::cos(alpha[2])/std::sin(dalpha[1])/std::sin(dalpha[2])
      , -std::cos(alpha[2])*std::cos(alpha[0])/std::sin(dalpha[2])/std::sin(dalpha[0])
      , -std::cos(alpha[0])*std::cos(alpha[1])/std::sin(dalpha[0])/std::sin(dalpha[1]);

    std::cout << mm*mminv << std::endl;

    std::vector<double> coeffs(3); // [mu]; y+nu,y,y+rho; (mu,nu,rho)
    coeffs[0] = std::cos(dalpha[0])/std::sin(dalpha[1])/std::sin(dalpha[2]);
    coeffs[1] = std::cos(dalpha[1])/std::sin(dalpha[2])/std::sin(dalpha[0]);
    coeffs[2] = std::cos(dalpha[2])/std::sin(dalpha[0])/std::sin(dalpha[1]);

    std::cout << coeffs[0] << " " << std::cos(dalpha[0]) / std::sin(dalpha[1])/std::sin(dalpha[2]) << std::endl;
    std::cout << coeffs[1] << " " << -( std::sin(alpha[2])*std::sin(alpha[0])+std::cos(alpha[2])*std::cos(alpha[0]) ) / std::sin(dalpha[2])/std::sin(dalpha[0]) << std::endl;
    std::cout << coeffs[2] << " " << -( std::sin(alpha[0])*std::sin(alpha[1])+std::cos(alpha[0])*std::cos(alpha[1]) ) / std::sin(dalpha[0])/std::sin(dalpha[1]) << std::endl;
  }

  return 0;
}


