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

constexpr int nparallel = 8;

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


// constexpr int nparallel = 1;

#include "obs.h"

constexpr int L = 64; // 4
constexpr Idx N = 10*L*L+2;
constexpr Idx N2 = 20*L*L;

#ifndef _OPENMP
int omp_get_thread_num(){ return 0; }
#endif



int main(int argc, char* argv[]){
// int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // omp_set_num_threads(nparallel);

  int seed=0;
  if(argc>=2) seed = atoi(argv[1]);

  bool is_correction=true;
  if(argc>=3) is_correction = atoi(argv[2]);


  bool if_read = true;

  int id=0;
  const std::string description = "L"+std::to_string(L)+"_"+std::to_string(id);
  const std::string dir = "./data_"+description+"/";
  std::filesystem::create_directories( dir );
  const std::string obsdir = "./obs_"+description+"/";
  std::filesystem::create_directories( obsdir );

  const int Nconf = 1e5;
  const int n_init = 1e3;

  // --------------------

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);
  Fermion D(spin);
  DualLoop<N> loop(D);
  FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
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

  using T1=Eigen::VectorXd;
  using T2=SpinField<N2>;
  // Jackknife<T1,T2> obs;
  Jackknife<T1,T2> obs(n_max-n_init+1);

  {
    T2 s;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(int n=n_init; n<=n_max; n++){
      const std::string filepath = dir+std::to_string(n);
      s.read(filepath);
      // obs.meas( s );
      obs.meas( n-n_init, s );
    }
  }

  std::cout << "# MEAS FINISHED. LEN = " << obs.size() << std::endl;

  const Idx i0 = 0;
  const T1 zero = T1::Zero( dual.NVertices() );

  auto f = [&](const std::vector<T2>& vs) {
    T1 mean = zero;
    std::vector<double> Nsq_mean(orbits.nbase(), 0.0);

    for(Idx k=0; k<vs.size(); k++) {
      const T2& s = vs[k];

      for(Idx if1=0; if1<dual.NVertices(); if1++){
        const auto& b0_g = orbits.b0_g_pairs[if1];
        Nsq_mean[b0_g.first] += s(if1) * s(orbits.antipodal[if1]);
      }

      for(Idx g=0; g<NIh; g++){
        T2 ss = s; // copy
        if( !ss[orbits[i0][g]] ) ss.flip();

        T1 tmp = zero;
        for(Idx i=0; i<dual.NVertices(); i++) tmp[ i ] = ss( orbits[i][g] );
        mean += tmp;
      }
    }
    mean /= n_max * NIh;

    if(is_correction){
      for(Idx b0=0; b0<orbits.nbase(); b0++) Nsq_mean[b0] /= n_max * orbits.npts[b0];
      const Idx b0 = orbits.b0_g_pairs[i0].first;
      for(Idx i=0; i<dual.NVertices(); i++) {
        const Idx b1 = orbits.b0_g_pairs[i].first;
        const double norm = std::sqrt( Nsq_mean[b0]*Nsq_mean[b1] ); // @@@
        mean[i] /= norm;
      }
    }

    return mean;
  };


  const int nbins = 10;
  const int binsize = obs.size()/nbins;
  obs.init( binsize );

  int ibin_min;
  for(ibin_min=0; ibin_min<obs.nbins; ibin_min++){
    // check ckpoints
    const std::string filepath = obsdir+"corr_"+std::to_string(ibin_min)+".dat";
    const bool bool_corr = std::filesystem::exists(filepath);
    if(!(bool_corr)) break;
  }

  std::cout << "# starting from ibin = " << ibin_min << std::endl;

#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
  for(int ibin=ibin_min; ibin<obs.nbins; ibin++){
    std::cout << "# debug. ibin = " << ibin << std::endl;
    const T1 jk_avg_corr = obs.jk_avg( ibin, f );
    // std::cout << "# debug. jk_corr" << std::endl;
    // for(double elem : jk_avg_corr) std::cout << elem << std::endl;
    const std::string filepath = obsdir+"corr_"+std::to_string(ibin)+".dat";
    obs.write( jk_avg_corr, filepath, jk_avg_corr.size() );
    // write ckpoints
  }

  // // in sep file
  // read ckpoints

  // obs.do_it( f, square, zero, binsize );

  // const T1 mean = obs.mean;
  // const T1 var = obs.var;

  // std::cout << "# ss : " << std::endl;
  // const V3 r0 = dual.vertices[i0];
  // const Idx b0 = orbits.b0_g_pairs[i0].first;
  // for(Idx i=0; i<dual.NVertices(); i++) {
  //   const Idx b1 = orbits.b0_g_pairs[i].first;
  //   const V3 r1 = dual.vertices[i];
  //   double ell = arcLength( r0, r1 );
  //   if(isnan(ell)){
  //     std::cout << "# debug. ell = " << ell << std::endl;
  //     std::cout << "# debug. r0 = " << r0.transpose() << std::endl;
  //     std::cout << "# debug. r1 = " << r1.transpose() << std::endl;
  //     ell = M_PI;
  //   }
  //   std::cout << ell << " " << mean[i] << " " << std::sqrt(var[i]) << std::endl;
  // }
  // std::cout << std::endl;

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
