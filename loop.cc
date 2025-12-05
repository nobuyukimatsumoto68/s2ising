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


#include <cstdint>
using Idx= std::uint_fast64_t;
using V3 = Eigen::Vector3d;
using V2 = Eigen::Vector2d;
using M3 = Eigen::Matrix3d;
using Complex = std::complex<double>;

#include "icosahedron.h"
#include "lattices.h"
#include "groups.h"
// #include "dual_optimizer.h"

#include "connection.h"



// #include <omp.h>
#include "loop.h"

constexpr int L = 2;
constexpr Idx N = 10*L*L+2;





int main(int argc, char* argv[]){
// int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------
  int k=3;
  if(argc>=2) k = atoi(argv[1]);
  int r=0;
  if(argc>=3) r = atoi(argv[2]);

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  // FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  // Rotation rot;
  RefinedIcosahedronDual dual(lattice);

  DualSpinStructure spin(dual);

  DualLoops<N> loops(dual);
  {
    Idx q = std::pow(2, k) + r; // Binary: 00100101
    loops.set( q );

    std::cout << "idx = " << loops() << std::endl;

    std::cout << "config = " << loops.config() << std::endl;

    std::cout << "loops = " << std::endl;
    std::cout << loops.printLoops() << std::endl;



    {
      double eps = 0.02;

      // std::ofstream ofs("loops"+std::to_string(L)+"_"+std::to_string(q)+".dat");
      std::ofstream ofs("loops.dat");
      ofs << std::scientific << std::setprecision(25);

      for(const auto& loop : loops.loops){
        for(Idx i1=0; i1<loop.size(); i1++){
          const Idx i2 = (i1+1)%loop.size();

          const Idx if1 = loop[i1];
          const Idx if2 = loop[i2];
          std::cout << if1 << " " << if2 << std::endl;

          const V3 rf1 = dual.vertices[if1];
          const V3 rf2 = dual.vertices[if2];

          for(double t=0.0; t<1.0; t+=eps){
            V3 p = spin.p(rf1, rf2, t);
            ofs << p.transpose() << std::endl;
          }
        }}

    }
  }

  

  return 0;
}
