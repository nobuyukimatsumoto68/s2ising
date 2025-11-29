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
#include "dual_optimizer.h"








// int main(int argc, char* argv[]){
int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------
  const int L = 8;
  RefinedIcosahedron lattice(L);

  FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  Rotation rot;

  RefinedIcosahedronDual dual(lattice);
  Orbits orbits(dual.vertices,
                dual.basePoints,
                dual.baseTypes,
                lattice, Ih, rot);

  DualLatticeAngleCostEvaluator opt(dual, lattice, Ih, rot);
  std::cout << "cost = " << opt.cost() << std::endl;

  return 0;
}
