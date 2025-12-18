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

  std::vector<double> ss_mean(s.s.size(), 0.0);
  for(int n=n_init; n<=n_max; n++){
    // std::cout << "read from n_max = " << n_max << std::endl;
    const std::string filepath = dir+std::to_string(n);
    s.read(filepath);
    std::cout << s.print() << std::endl;

    Spin<N,N2> ss(ising);
    ss.s = s.s;
    if(!ss[0]) ss.flip();
    for(Idx i=0; i<ss.s.size(); i++) ss_mean[i] += ss[i];
    //
    // ss corr -> s field
    // double object -> mean, var
    // icos symmetry shuffle; avg
  }
  for(Idx i=0; i<ss_mean.size(); i++) ss_mean[i] /= n_max;

  for(Idx i=0; i<ss_mean.size(); i++) std::cout << ss_mean[i] << " ";
  std::cout << std::endl;

  return 0;
}

