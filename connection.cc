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

#include <algorithm>
#include <filesystem>


#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>


using Idx=std::size_t;
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
using Complex = std::complex<double>;

#include "icosahedron.h"
#include "lattices.h"
#include "groups.h"
// #include "dual_optimizer.h"

#include "connection.h"


// #include <omp.h>




int main(int argc, char* argv[]){
// int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------
  // const int L = 16;
  // const int L = 32;
  // const int L = 48;
  int L=1;
  if(argc==2) L = atoi(argv[1]);

  RefinedIcosahedron lattice(L);
  FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  Rotation rot;
  RefinedIcosahedronDual dual(lattice);

  // Orbits orbits(dual.vertices,
  //               dual.basePoints,
  //               dual.baseTypes,
  //               lattice, Ih, rot);

  DualSpinStructure spin(dual);

  V3 r1,r2;
  r1 << 0.3, 0.6, 0.5;
  r2 << 0.4, 0.2, 0.3;
  r1 = r1/r1.norm();
  r2 = r2/r2.norm();

  const double stilde = 0.2;
  std::cout << spin.ptilde(r1,r2,stilde).transpose() << std::endl;
  std::cout << spin.p(r1,r2,stilde).transpose() << std::endl;
  std::cout << spin.e(r1,r2,stilde).transpose() << std::endl;
  std::cout << spin.dtheta_dstilde(r1,r2,stilde) << std::endl
            << spin.dphi_dstilde(r1,r2,stilde) << std::endl;
  std::cout << spin.alpha(r1,r2,stilde) << std::endl;
  // std::cout << spin.p(r1,r2,stilde).transpose() << std::endl;

<<<<<<< HEAD
  std::cout << "faces = " << std::endl;
  for(auto& face : spin.faces){
    std::cout << "length = " << face.size() << std::endl;
    for(auto elem : face) std::cout << elem << " ";
    std::cout << std::endl;
  }

=======

  for(Idx if1=0; if1<dual.NVertices(); if1++){
    const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
    assert( if1==dual.idx(f1) );
    // const V3 rf1 = dual.vertices[if1];

    for(int df=0; df<3; df++){
      FaceCoords f2;
      dual.shift( f2, f1, df);
      const Idx if2 = dual.idx( f2 );
      // assert( if1==dual.idx(f1) );
      // const V3 rf2 = dual.vertices[if2];

      const double alpha1 = spin.alphas[ spin.idx(if1, df) ];
      const double alpha2 = spin.alphas[ spin.idx(if2, df) ];
      const double omega12 = spin.omegas[ spin.idx(if1, df) ];

      double diff = (alpha2 + M_PI + omega12) - alpha1;
      // while(diff>M_PI) diff -= 2.0*M_PI;
      // while(diff<=-M_PI) diff += 2.0*M_PI;
      // std::cout << "diff = " << diff << std::endl;
      assert( std::abs(diff)<1.0e-14 );
    }
  }










  std::cout << "faces = " << std::endl;
  for(auto& face : spin.faces){
    std::cout << "length = " << face.size() << std::endl;
    double omega_sum = 0.0;
    for(auto elem : face) {
      omega_sum += spin.omegas[ spin.idx( elem.first, elem.second ) ];
    }
    std::cout << omega_sum << std::endl;
  }

  {
    const double phi0 = M_PI/48.0;
    // const double phi0 = -M_PI/48.0;
    for(const auto& r : dual.vertices ) assert( std::abs( spin.phi(r) - phi0 ) >1.0e-14 );
    std::vector<DualSpinStructure::DirectedLink> section;

    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );
        const V3 rf2 = dual.vertices[if2];

        double phi1 = spin.phi(rf1);
        double phi2 = spin.phi(rf2);

        bool is_yes = false;

        const int nstep=20;
        for(double stilde=0.0; stilde<=1.0; stilde+=1.0/nstep){
          double phi = spin.phi( spin.p(rf1,rf2,stilde) );
          // std::cout << "phi = " << phi << std::endl;
          if( std::abs(phi-phi0)<2.0*M_PI/5.0/L/nstep ) is_yes=true;
        }

        // if( std::abs(phi1-M_PI)<1.0e-14 && phi2<0.0 ) phi1 -= 2.0*M_PI;
        // if( std::abs(phi1+M_PI)<1.0e-14 && phi2>0.0 ) phi1 += 2.0*M_PI;
        // if( std::abs(phi2-M_PI)<1.0e-14 && phi1<0.0 ) phi2 -= 2.0*M_PI;
        // if( std::abs(phi2+M_PI)<1.0e-14 && phi1>0.0 ) phi2 += 2.0*M_PI;

        // if(phi1<phi2) {
        //   if( std::abs(phi2-phi1-2.0*M_PI)<1.0e-14 ) continue;
        //   else if(phi1<phi0 && phi0<phi2) is_yes = true;
        // }
        // if(phi2<phi1) {
        //   if( std::abs(phi1-phi2-2.0*M_PI)<1.0e-14 ) continue;
        //   else if(phi2<phi0 && phi0<phi1) is_yes = true;
        // }
        // std::cout << "debug. " << std::endl
        //           << rf1.transpose() << std::endl
        //           << rf2.transpose() << std::endl;
        // std::cout << "phi1, phi0, phi2 = " << phi1 << " " << phi0 << " " << phi2 << " " << is_yes << std::endl;

        // double dphi1 = phi1-phi0;
        // double dphi2 = phi2-phi0;
        // double dphi = phi1-phi2;

        // if(dphi1 >= M_PI) dphi1 -= 2.0*M_PI;
        // else if(dphi1 < -M_PI) dphi1 += 2.0*M_PI;
        // if(dphi2 >= M_PI) dphi2 -= 2.0*M_PI;
        // else if(dphi2 < -M_PI) dphi2 += 2.0*M_PI;
        // if(dphi >= M_PI) dphi -= 2.0*M_PI;
        // else if(dphi < -M_PI) dphi += 2.0*M_PI;
        // std::cout << "dphi1, dphi2, dphi = " << dphi1 << " " << dphi2 << " " << dphi << std::endl;
        // if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < 1.0e-14 ){
        //   section.push_back( std::make_pair(if1, df) );
        //   std::cout << "yes" << std::endl;
        // }
        if( is_yes ){
          section.push_back( std::make_pair(if1, df) );
          std::cout << "yes" << std::endl;
        }
      }
    }
    for(const DualSpinStructure::DirectedLink& ell : section) {
      spin.omegas[spin.idx(ell)] -= 2.0*M_PI;
    }
    // for(auto elem : section) std::cout << elem << " ";
    // std::cout << std::endl;
    // std::cout << section.size() << std::endl;
    // std::cout << links.size() << std::endl;
  }


  std::cout << "faces = " << spin.faces.size() << std::endl;
  for(auto& face : spin.faces){
    std::cout << "length = " << face.size() << std::endl;
    double omega_sum = 0.0;
    for(auto elem : face) {
      // std::cout << elem.first << " "  << elem.second << std::endl;
      omega_sum += spin.omegas[ spin.idx( elem.first, elem.second ) ];
    }
    while(omega_sum>2.0*M_PI) omega_sum -= 4.0*M_PI;
    while(omega_sum<=-2.0*M_PI) omega_sum += 4.0*M_PI;

    std::cout << omega_sum << std::endl;
  }


>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41


  return 0;
}
