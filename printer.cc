#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>
#include <complex>

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

#include "sphere.h"
#include "icosahedron.h"
#include "lattices.h"
#include "groups.h"
// #include "dual_optimizer.h"

#include "connection.h"



// #include <omp.h>


constexpr int L = 1;
constexpr Idx N = 10*L*L+2;
#include <bitset>
using SimpConfig = std::bitset<N>;
#include "loop.h"




int main(int argc, char* argv[]){
// int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------
  // const int L = 16;
  // const int L = 32;
  // const int L = 48;
  int k=3;
  if(argc>=2) k = atoi(argv[1]);
  int r=0;
  if(argc>=3) r = atoi(argv[2]);

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  Rotation rot;
  RefinedIcosahedronDual dual(lattice);

  {
    std::ofstream ofs("simp"+std::to_string(L)+".dat");
    ofs << std::scientific << std::setprecision(25) << lattice.print();
  }
  {
    std::ofstream ofs("dual"+std::to_string(L)+".dat");
    ofs << std::scientific << std::setprecision(25) << dual.print();
  }

  // int nsteps=20;
  double eps = 0.02;
  DualSpinStructure spin(dual);


  // simp
  {
    std::ofstream ofs("simplink"+std::to_string(L)+".dat");
    ofs << std::scientific << std::setprecision(25);

    for(Idx in1=0; in1<lattice.NVertices(); in1++){
      const Coords n1 = lattice.idx2Coords(in1); //{ n[0],n[1],n[2],type };
      const V3 rn1 = lattice.vertices[in1];

      for(int dn=0; dn<3; dn++){
        Coords n2;
        lattice.shift(n2, n1, dn);
        const V3 rn2 = lattice.vertices[lattice.idx(n2)];

        for(double t=0.0; t<1.0; t+=eps){
          V3 p = spin.p(rn1,rn2, t);
          ofs << p.transpose() << std::endl;
        }
      }
    }
  }

  // dual
  {
    std::ofstream ofs("duallink"+std::to_string(L)+".dat");
    ofs << std::scientific << std::setprecision(25);

    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      // if(f1[3]!=XZ) continue;
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );
        const V3 rf2 = dual.vertices[if2];

        for(double t=0.0; t<1.0; t+=eps){
          V3 p = spin.p(rf1,rf2, t);
          ofs << p.transpose() << std::endl;
        }
      }
    }
  }

  {
    // [5 +z ->7]
    const Idx in1 = 6; // in2 = 7;
    const int dn = lattice.nZ;
    const V3 rn1 = lattice.vertices[in1];

    Coords n1 = lattice.idx2Coords(in1), n2;
    lattice.shift(n2, n1, dn);
    const Idx in2 = lattice.idx(n2);
    const V3 rn2 = lattice.vertices[in2];

    auto duallink = dual.linkSimp2Dual.at( DirectedLink{in1, dn} );

    const Idx if1 = duallink.first;
    const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
    const V3 rf1 = dual.vertices[if1];

    FaceCoords f2;
    dual.shift( f2, f1, duallink.second);
    const V3 rf2 = dual.vertices[dual.idx(f2)];

    std::ofstream ofs("mix"+std::to_string(L)+".dat");
    ofs << std::scientific << std::setprecision(25);
    {
      for(double t=0.0; t<1.0; t+=eps){
        V3 p = spin.p(rf1,rn1, t);
        ofs << p.transpose() << std::endl;
      }
      for(double t=0.0; t<1.0; t+=eps){
        V3 p = spin.p(rf1,rn2, t);
        ofs << p.transpose() << std::endl;
      }
      for(double t=0.0; t<1.0; t+=eps){
        V3 p = spin.p(rf2,rn1, t);
        ofs << p.transpose() << std::endl;
      }
      for(double t=0.0; t<1.0; t+=eps){
        V3 p = spin.p(rf2,rn2, t);
        ofs << p.transpose() << std::endl;
      }
    }
  }

  return 0;
}
