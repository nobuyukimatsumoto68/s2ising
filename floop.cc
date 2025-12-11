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
using M2 = Eigen::Matrix2cd;
using Complex = std::complex<double>;



#include "sphere.h"
#include "icosahedron.h"
#include "lattices.h"
#include "groups.h"
// #include "dual_optimizer.h"

#include "connection.h"



// #include <omp.h>
#include "loop.h"
#include "fermion.h"

constexpr int L = 2;
constexpr Idx N = 10*L*L+2;





int main(int argc, char* argv[]){
// int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------
  int k=0;
  if(argc>=2) k = atoi(argv[1]);
  int r=1;
  if(argc>=3) r = atoi(argv[2]);

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  // FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  // Rotation rot;
  RefinedIcosahedronDual dual(lattice);
  DualSpinStructure spin(dual);

  DualLoops<N> loops(dual);

  Fermion D(spin);
  // FermionVector phi(dual);

  Idx q = std::pow(2, k) + r; // Binary: 00100101
  loops.set( q );

  std::cout << "idx = " << loops() << std::endl;
  std::cout << "config = " << loops.config() << std::endl;
  std::cout << "loops = " << std::endl;
  std::cout << loops.printLoops() << std::endl;


  double w_prod = 1;
  for(const auto& loop : loops){
    double w = 1;
    double factor = 1;
    SpinMatrix mat = SpinMatrix::Identity();
    for(int i=0; i<loop.size(); i++){
      const Idx if0 = loop[(i+loop.size()-1)%loop.size()];
      const Idx if1 = loop[i];
      const Idx if2 = loop[(i+1)%loop.size()];

      const FaceCoords f1 = dual.idx2FaceCoords(if1);
      const int df2 = dual.getDirection.at(Link{if1,if2});

      std::cout << if1 << " " << if2 << " " << df2 << std::endl;
      factor *= D.kappas[dual.linkidx(if1,df2)];
      factor /= D.mus[if1];
      mat = mat * spin.P( f1, df2 );
      std::cout << D.kappas[dual.linkidx(if1,df2)] << " " << D.mus[if1] << std::endl;

      const int df0 = dual.getDirection.at(Link{if1,if0});

      double al10 = spin.alphas[dual.linkidx(if1, df0)];
      double al12 = spin.alphas[dual.linkidx(if1, df2)];
      double dalpha = (al10+M_PI)-al12;
      w *= D.kappas[dual.linkidx(if1,df2)];
      w /= D.mus[if1];
      while(dalpha>M_PI) dalpha -= 2.0*M_PI;
      while(dalpha<-M_PI) dalpha += 2.0*M_PI;
      w *= std::cos( dalpha/2 );
      std::cout << "debug. cos = " << std::cos( dalpha/2 ) << std::endl;
    }
    w *= -1.0;
    const Complex w_ = factor*mat.trace();
    std::cout << "w_ = " << w_ << std::endl;
    std::cout << "w  = " << w << std::endl;
    w_prod *= w;
  }
  std::cout << "w_prod  = " << w_prod << std::endl;


  return 0;
}
