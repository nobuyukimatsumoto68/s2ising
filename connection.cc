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

  std::cout << "faces = " << std::endl;
  for(auto& face : spin.faces){
    std::cout << "length = " << face.size() << std::endl;
    for(auto elem : face) std::cout << elem << " ";
    std::cout << std::endl;
  }



  return 0;
}
