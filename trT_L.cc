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


const int nparallel = 8;

#include "sphere.h"
#include "icosahedron.h"
#include "lattices.h"
#include "groups.h"


#ifndef _OPENMP
int omp_get_thread_num(){ return 0; }
#endif


int main(int argc, char* argv[]){
// int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  omp_set_num_threads(nparallel);

  // --------------------

  FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  Rotation rot;

  for(int L : std::vector<int>{2,4,8,16,32,64}){
    if(L==1 || L%2==0){
      RefinedIcosahedron lattice(L);
      RefinedIcosahedronDual dual(lattice);

      Orbits orbits(dual.vertices,
                    dual.basePoints,
                    dual.baseTypes,
                    lattice, Ih, rot);
      const int norbits = orbits.nbase();

      // const Idx b0 = orbits.b0_g_pairs.back().first;
      // const Idx b0 = orbits.b0_g_pairs[norbits-1].first;
      Idx if1;
      for(if1=0; if1<dual.NVertices(); if1++){
        const auto& b0_g = orbits.b0_g_pairs[if1];
        if(b0_g.first==norbits-1) break;
      }
      std::cout << L << " " << dual.mean_ell << " " << std::pow(dual.mean_ell, 2) * 2.0*M_PI/(dual.site_volumes[if1]) << std::endl;
    }
  }



  return 0;
}


