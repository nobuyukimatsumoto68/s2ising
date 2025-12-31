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

const int nparallel = 12;

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

constexpr int L = 4; // 4
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
  Idx iy=0;
  if(argc>=4) iy = atoi(argv[3]);
  int nu=0;
  if(argc>=5) nu = atoi(argv[4]);


  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);
  // Fermion D(spin, true);
  Fermion D(spin);


  Idx iyp;
  int nup;
  {
    FaceCoords f2;
    dual.shift( f2, dual.idx2FaceCoords(iy), nu );
    iyp = dual.idx(f2);
    nup = dual.getDirection.at({iyp,iy});
  }

  // --------------------

  double dbeta_dmus[3];
  double dbeta_dkappas[3];
  double betaps[3];
  double betams[3];
  double ts_augps[3];
  double ts_augms[3];
  double sum, sump, summ;

  double eps = 1.0e-6;
  Idx q = std::pow(2, k) + r; // Binary: 00100101
  // Idx iy=0;
  {
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);
    s.set( q );
    sum = std::exp(s.S());

    std::cout << "Z (ising, sum) = " << sum << std::endl;
  }
  {
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);
    s.set( q );

    double sum1 = 0.0;
    for(int mu=0; mu<3; mu++){
      double dbeta_dmu = -std::sinh( 2.0*ising.betas.at(DirectedLink{iy,mu}) ) /4.0/D.mus[iy];
      dbeta_dmus[mu] = dbeta_dmu;

      FaceCoords f2;
      dual.shift( f2, dual.idx2FaceCoords(iy), mu );
      const Idx if2 = dual.idx(f2);
      sum1 += dbeta_dmu * s.s(iy) * s.s(if2);
    }
    std::cout << "sum1 = " << sum1 << std::endl;
    std::cout << "ratio = " << sum1 / sum << std::endl;
  }
  {
    D.mus[iy] += eps;
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    for(int mu=0; mu<3; mu++){
      betaps[mu] = ising.betas.at(DirectedLink{iy,mu});
    }

    s.set( q );
    sump = std::exp(s.S());
  }
  {
    D.mus[iy] -= 2.0*eps;
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    for(int mu=0; mu<3; mu++){
      betams[mu] = ising.betas.at(DirectedLink{iy,mu});
    }

    s.set( q );
    summ = std::exp(s.S());
  }
  std::cout << "num deriv = " << (std::log(sump)-std::log(summ))/(2.0*eps) << std::endl;

  {
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    for(int mu=0; mu<3; mu++){
      // std::cout << dbeta_dmus[mu] << " " << (betaps[mu]-betams[mu])/(2.0*eps) << std::endl;
      std::cout << dbeta_dmus[mu] << " " << (betaps[mu]-betams[mu])/(2.0*eps) << " " << ising.dbeta_dmus[dual.directedlinkidx(iy,mu)] << std::endl;
    }
  }




  D.mus[iy] -= eps;
  {
    DualLoop<N> loop(D);
    Ising ising(loop);

    for(int mu=0; mu<3; mu++){
      double dbeta_dkappa = 0.5*std::sinh( 2.0*ising.betas.at(DirectedLink{iy,mu}) ) /D.kappas[dual.directedlinkidx(iy,mu)];
      dbeta_dkappas[mu] = dbeta_dkappa;
    }
  }
  {
    D.kappas[dual.directedlinkidx(iy,nu)] += eps;
    D.kappas[dual.directedlinkidx(iyp,nup)] += eps;
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    for(int mu=0; mu<3; mu++){
      betaps[mu] = ising.betas.at(DirectedLink{iy,mu});
      ts_augps[mu] = ising.ts_aug.at( DirectedLink{iy,mu} );
    }

    s.set( q );
    sump = std::exp(s.S());
  }
  {
    D.kappas[dual.directedlinkidx(iy,nu)] -= 2.0*eps;
    D.kappas[dual.directedlinkidx(iyp,nup)] -= 2.0*eps;
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    for(int mu=0; mu<3; mu++){
      betams[mu] = ising.betas.at(DirectedLink{iy,mu});
      ts_augms[mu] = ising.ts_aug.at( DirectedLink{iy,mu} );
    }

    s.set( q );
    summ = std::exp(s.S());
  }
  std::cout << "num deriv = " << (std::log(sump)-std::log(summ))/(2.0*eps) << std::endl;

  {
    DualLoop<N> loop(D);
    Ising ising(loop);
    Spin<N, N2> s(ising);

    for(int mu=0; mu<3; mu++){
      // std::cout << 0.5/D.kappas[dual.directedlinkidx(iy,mu)] << " " << (std::log(ts_augps[mu])-std::log(ts_augms[mu]))/(2.0*eps) << std::endl;
      std::cout << dbeta_dkappas[mu] << " " << (betaps[mu]-betams[mu])/(2.0*eps) << " " << ising.dbeta_dkappas[dual.directedlinkidx(iy,mu)] << std::endl;
    }
  }

  return 0;
}


