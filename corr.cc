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

  const int nparallel = 1;
  omp_set_num_threads(nparallel);

  int seed=0;
  if(argc>=2) seed = atoi(argv[1]);


  bool if_read = true;

  const std::string description = "L"+std::to_string(L);
  const std::string dir = "./data_"+description+"/";
  std::filesystem::create_directories( dir );

  // const int Nrepeat = 100;
  const int Nconf = 1e5;
  const int n_init = 1e2;

  // --------------------

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);
  Fermion D(spin);
  DualLoop<N> loop(D);
  Ising ising(loop);
  Spin<N, N2> s(ising,seed);




  FullIcosahedralGroup Ih( "multtablemathematica.dat",
                           3, 19, 60 );
  Rotation rot;

  Orbits orbits(dual.vertices,
                dual.basePoints,
                dual.baseTypes,
                lattice, Ih, rot);


  int n_max=0;
  {
    for(n_max=1; n_max<Nconf; n_max++ ){
      const std::string filepath = dir+std::to_string(n_max);
      const bool bool_lat = std::filesystem::exists(filepath+".lat");
      const bool bool_rng = std::filesystem::exists(filepath+".rng");
      if(!(bool_lat&&bool_rng)) break;
    }
    n_max -= 1;
  }


  std::vector<double> Nsq_mean(orbits.nbase(), 0.0);

  const Idx i0 = 0;

  std::vector<double> ss_mean(s.s.size(), 0.0);
  for(int n=n_init; n<=n_max; n++){
    const std::string filepath = dir+std::to_string(n);
    s.read(filepath);

    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const auto& b0_g = orbits.b0_g_pairs[if1];
      Nsq_mean[b0_g.first] += s[if1]*s[orbits.antipodal[if1]];
    }

    for(Idx g=0; g<NIh; g++){
      Spin<N,N2> ss(ising);
      ss.s = s.s;
      if(!ss.s[ orbits[i0][g] ]) ss.flip();
      for(Idx i=0; i<ss.s.size(); i++) ss_mean[ i ] += ss[ orbits[i][g] ];
    }
  }
  for(Idx b0=0; b0<orbits.nbase(); b0++) Nsq_mean[b0] /= n_max * orbits.npts[b0];
  for(Idx i=0; i<ss_mean.size(); i++) ss_mean[i] /= n_max * NIh;


  std::vector<double> ss_var(s.s.size(), 0.0);
  for(int n=n_init; n<=n_max; n++){
    const std::string filepath = dir+std::to_string(n);
    s.read(filepath);

    for(Idx g=0; g<NIh; g++){
      Spin<N,N2> ss(ising);
      ss.s = s.s;
      if(!ss.s[ orbits[i0][g] ]) ss.flip();
      for(Idx i=0; i<ss.s.size(); i++) ss_var[ i ] += std::pow(ss[ orbits[i][g] ]-ss_mean[i], 2);
    }
  }
  for(Idx i=0; i<ss_var.size(); i++) ss_var[i] /= n_max * NIh;



  std::cout << "# Nsq : " << std::endl << "# ";
  for(Idx b0=0; b0<orbits.nbase(); b0++) {
    std::cout << std::sqrt(Nsq_mean[b0]) << " ";
  }
  std::cout << std::endl;;


  std::cout << "# ss : " << std::endl;
  const V3 r0 = dual.vertices[i0];
  const Idx b0 = orbits.b0_g_pairs[i0].first;
  for(Idx i=0; i<ss_mean.size(); i++) {
    const Idx b1 = orbits.b0_g_pairs[i].first;
    const double norm = 1.0; // std::sqrt( Nsq_mean[b0]*Nsq_mean[b1] ); // @@@

    const V3 r1 = dual.vertices[i];
    double ell = arcLength( r0, r1 );
    if(isnan(ell)){
      std::cout << "# debug. ell = " << ell << std::endl;
      std::cout << "# debug. r0 = " << r0.transpose() << std::endl;
      std::cout << "# debug. r1 = " << r1.transpose() << std::endl;
      ell = M_PI;
    }
    std::cout << ell << " " << ss_mean[i]/norm << " " << std::sqrt(ss_var[i])/norm << std::endl;
  }
  std::cout << std::endl;

  return 0;
}



  // Fundamental = dual.basepoints
  // antipodal normalization
  // func find the corresp fund.

  // std::cout << "debug. orbits = " << std::endl;
  // std::cout << orbits.print() << std::endl;
  // std::cout << "debug. fund = " << std::endl;
  // for(const Idx elem : orbits.basePoints) std::cout << elem << " ";
  // std::cout << std::endl;




  // orbits.b0_g_pairs.clear();
  // orbits.b0_g_pairs.resize(dual.NVertices());
  // for(Idx if1=0; if1<dual.NVertices(); if1++){
  //   const auto orbit = orbits[if1];
  //   bool is_found = false;
  //   for(Idx g=0; g<NIh; g++){
  //     // for( const Idx b0 : orbits.basePoints ){
  //     for(Idx b0=0; b0<orbits.basePoints.size(); b0++){
  //       if(orbits.basePoints[b0]==orbit[g]) {
  //         orbits.b0_g_pairs[if1] = std::pair<Idx,Idx>{b0,g};
  //         is_found = true;
  //         break;
  //       }
  //     }
  //     if(is_found) break;
  //   }
  //   if(!is_found) {
  //     std::cout << "debug. orbit = " << std::endl;
  //     for(auto elem : orbit) std::cout << elem << " ";
  //     std::cout << std::endl;
  //     assert(false);
  //   }
  // } // end if1

  // // std::cout << "debug. pairs = " << std::endl;
  // // for(const auto& b0_g : orbits.b0_g_pairs) std::cout << b0_g.first << " " << b0_g.second << std::endl;
  // // std::cout << "debug. b0_g.size = " << orbits.b0_g_pairs.size() << std::endl;

  // orbits.npts.clear();
  // orbits.npts.resize(orbits.basePoints.size());
  // for(Idx& elem : orbits.npts) elem = 0;
  // Idx counter=0;
  // for(const auto& b0_g : orbits.b0_g_pairs) {
  //   orbits.npts[b0_g.first]+=1;
  //   counter++;
  // }
  // std::cout << "debug. counter = " << counter << std::endl;
  // std::cout << "debug. npts = " << std::endl;
  // for(Idx elem : orbits.npts) std::cout << elem << " ";
  // std::cout << std::endl;

  // Idx sum=0;
  // for(Idx elem : orbits.npts) sum+=elem;
  // std::cout << "debug. sum = " << sum << std::endl;
  // std::cout << "debug. NVertices() = " << dual.NVertices() << std::endl;
  // assert(sum==dual.NVertices());

  // return 1;
