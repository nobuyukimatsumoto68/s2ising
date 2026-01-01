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

#ifndef nparallel
#define nparallel 1
#endif

#ifndef LL
#define LL 2
#endif

#ifndef Nconf
#define Nconf 1e5
#endif


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



#include "obs.h"

constexpr int L = LL; // 4$
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

  // const int Nconf = 4e5;
  const int n_init = 1e3;

  // --------------------

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);
  Fermion D(spin);
  DualLoop<N> loop(D);
  Ising ising(loop);
  FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  Rotation rot;

  Orbits orbits(dual.vertices,
                dual.basePoints,
                dual.baseTypes,
                lattice, Ih, rot);


  int n_max=Nconf-1;
  // {
  //   for(n_max=1; n_max<Nconf; n_max++ ){
  //     const std::string filepath = dir+std::to_string(n_max);
  //     const bool bool_lat = std::filesystem::exists(filepath+".lat");
  //     const bool bool_rng = std::filesystem::exists(filepath+".rng");
  //     if(!(bool_lat&&bool_rng)) break;
  //   }
  //   n_max -= 1;
  // }

  std::cout << "# START MEAS" << std::endl;

  using T1=Eigen::VectorXd;
  using T2=SpinField<N2>;
  // Jackknife<T1,T2> obs;
  JackknifeSimp<T1,T2> obs(n_max-n_init+1);

  std::cout << "# READING" << std::endl;
  {
#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(int n=n_init; n<=n_max; n++){
      T2 s;
      const std::string filepath = dir+std::to_string(n);
      s.read(filepath);
      obs.meas( n-n_init, s );
    }
  }

  std::cout << "# MEAS FINISHED. LEN = " << obs.size() << std::endl;
  std::cout << "# ell_mean = " << dual.mean_ell << std::endl;

  const Idx i0 = 0;
  // const T1 zero = T1::Zero( dual.NVertices() );
  const T1 zero = T1::Zero( orbits.nbase() ); // T1::Zero( dual.NVertices() );

  auto mean12 = [&](const std::vector<T2>& vs) {
                  T1 T11_mean = zero;
                  for(Idx k=0; k<vs.size(); k++) {
                    const T2& s = vs[k];
                    for(Idx if1=0; if1<dual.NVertices(); if1++){
                      const auto& b0_g = orbits.b0_g_pairs[if1];
                      T11_mean[b0_g.first] += s.T11(if1, ising);
                    }
                  }
                  for(Idx b0=0; b0<orbits.nbase(); b0++) T11_mean[b0] /= 1.0 * vs.size() * orbits.npts[b0];

                  return T11_mean;
                };

  auto mean11 = [&](const std::vector<T1>& vs) {
                  T1 mean = zero;
                  for(Idx k=0; k<vs.size(); k++) mean += vs[k];
                  mean /= 1.0 * vs.size();
                  return mean;
                };

  auto square = [](const T1& x) { return x.array().square().matrix(); };


  const int nbins = 10;
  const int binsize = obs.size()/nbins;
  obs.init( binsize );

#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
  for(int ibin=0; ibin<obs.nbins; ibin++) {
    std::cout << "# ibin = " << ibin << std::endl;
    obs.bin_avg[ibin] = obs.get_bin_avg(ibin, mean12);
  }
  for(int ibin=0; ibin<obs.nbins; ibin++) obs.jack_avg[ibin] = obs.jk_avg(ibin, mean11);
  obs.finalize(square, zero);

  const T1 mean = obs.mean;
  const T1 var = obs.var;

  {
    const std::string filepath = obsdir+"T11_jk.dat";
    std::ofstream os( filepath, std::ios::out | std::ios::trunc );
    os << std::scientific << std::setprecision(25);
    if(!os) assert(false);
    os << "# ell_mean = " << dual.mean_ell << std::endl;
    std::cout << "# T11 : " << std::endl;
    for(Idx b0=0; b0<orbits.nbase(); b0++) {
      std::cout << mean[b0] << " " << std::sqrt(var[b0]) << std::endl;
      os << mean[b0] << " " << std::sqrt(var[b0]) << std::endl;
    }
    os.close();
  }

  {
    const std::string filepath = obsdir+"orbit_size.dat";
    std::ofstream os( filepath, std::ios::out | std::ios::trunc );
    os << std::scientific << std::setprecision(25);
    if(!os) assert(false);
    os << "# ell_mean = " << dual.mean_ell << std::endl;
    std::cout << "# orbit_size : " << std::endl;
    os << orbits.print_orbitsize();
    os.close();
  }


  return 0;
}

