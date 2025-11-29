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
  std::cout << "# points" << std::endl
            << lattice.print() << std::endl;

  {
    for(int s=0; s<NPatches; s++){
      for(int x=0; x<L; x++){
        for(int y=0; y<L; y++){
          const Coords n{x,y,s};
          const Idx in = lattice.idx(n);

          Coords nPX, nPY, nPZ;
          lattice.shiftPX( nPX, n );
          lattice.shiftPY( nPY, n );
          lattice.shiftPZ( nPZ, n );
          const Idx inPX = lattice.idx(nPX);
          const Idx inPY = lattice.idx(nPY);
          const Idx inPZ = lattice.idx(nPZ);

          const V3 vn = lattice.vertices[in];
          const V3 vnPX = lattice.vertices[inPX];
          const V3 vnPY = lattice.vertices[inPY];
          const V3 vnPZ = lattice.vertices[inPZ];

          {
            const double dx = (vn-vnPX).norm();
            const double dy = (vn-vnPY).norm();
            const double dz = (vn-vnPZ).norm();

            std::cout << "@@ in = " << in << std::endl;
            std::cout << "@@ x,y,s = " << x << " " << y << " " << s << std::endl;
            std::cout << "+x x,y,s = " << nPX[0] << " " << nPX[1] << " " << nPX[2] << std::endl;
            std::cout << "+y x,y,s = " << nPY[0] << " " << nPY[1] << " " << nPY[2] << std::endl;
            std::cout << "+z x,y,s = " << nPZ[0] << " " << nPZ[1] << " " << nPZ[2] << std::endl;
            std::cout << dx << std::endl;
            std::cout << dy << std::endl;
            std::cout << dz << std::endl;

            // // for L=2
            // assert( dx<0.7 );
            // assert( dy<0.7 );
            // assert( dz<0.7 );
            // // for L=4
            // assert( dx<0.4 );
            // assert( dy<0.4 );
            // assert( dz<0.4 );
          }

          Coords nMX, nMY, nMZ;
          lattice.shiftMX( nMX, n );
          lattice.shiftMY( nMY, n );
          lattice.shiftMZ( nMZ, n );

          if(nMX[2]!=PatchIdxVoid){
            const Idx inMX = lattice.idx(nMX);
            const V3 vnMX = lattice.vertices[inMX];
            const double dx = (vn-vnMX).norm();
            std::cout << "-x x,y,s = " << nMX[0] << " " << nMX[1] << " " << nMX[2] << std::endl;
            std::cout << dx << std::endl;
            // assert( dx<0.7 );
            // assert( dx<0.4 );
          }
          if(nMY[2]!=PatchIdxVoid){
            const Idx inMY = lattice.idx(nMY);
            const V3 vnMY = lattice.vertices[inMY];
            const double dy = (vn-vnMY).norm();
            std::cout << "-x x,y,s = " << nMY[0] << " " << nMY[1] << " " << nMY[2] << std::endl;
            std::cout << dy << std::endl;
            // assert( dy<0.7 );
            // assert( dy<0.4 );
          }
          if(nMZ[2]!=PatchIdxVoid){
            const Idx inMZ = lattice.idx(nMZ);
            const V3 vnMZ = lattice.vertices[inMZ];
            const double dz = (vn-vnMZ).norm();
            std::cout << "-x x,y,s = " << nMZ[0] << " " << nMZ[1] << " " << nMZ[2] << std::endl;
            std::cout << dz << std::endl;
            // assert( dz<0.7 );
            // assert( dz<0.4 );
          }
        }}}
  }

  FullIcosahedralGroup Ih( "multtablemathematica.dat",
                           3, 19, 60 );

  Rotation rot;

  FullIcosahedralGroup::Rep rep;
  std::cout << "# Ih elems" << std::endl;
  std::cout << Ih.print() << std::endl;

  {
    // BasePoints basePoints;
    // BaseTypes baseTypes;
    // lattice.GetBasePoints( basePoints, baseTypes );

    for(Idx i=0; i<lattice.basePoints.size(); i++){
      std::cout << "# basepoints" << std::endl;
      std::cout << "@@@ i = " << i << " type = " << lattice.baseTypes[i] << std::endl;
      const Coords n = lattice.idx2Coords( lattice.basePoints[i] );
      std::cout << "x,y,s = " << n[0] << " " << n[1] << " " << n[2] << std::endl;
    }
  }

  {
    // for simplicial
    Orbits orbits(lattice.vertices,
                  lattice.basePoints,
                  lattice.baseTypes,
                  lattice, Ih, rot);
    std::cout << "# orbits (simplicial)" << std::endl;
    std::cout << orbits.print() << std::endl;
  }



  RefinedIcosahedronDual dual(lattice);
  {
    for(int s=0; s<NPatches; s++){
      for(int x=0; x<L; x++){
        for(int y=0; y<L; y++){
          // for(const int type : std::vector<int>{XZ, ZY}){
          {
            const int type = XZ;
            const FaceCoords f{x,y,s,type};
            const Idx iff = dual.idx(f);

            FaceCoords fPA, fPB, fPC;
            dual.shiftPA( fPA, f );
            dual.shiftPB( fPB, f );
            dual.shiftPC( fPC, f );
            const Idx ifPA = dual.idx(fPA);
            const Idx ifPB = dual.idx(fPB);
            const Idx ifPC = dual.idx(fPC);

            const V3 vf = dual.vertices[iff];
            const V3 vfPA = dual.vertices[ifPA];
            const V3 vfPB = dual.vertices[ifPB];
            const V3 vfPC = dual.vertices[ifPC];

            const double dA = (vf-vfPA).norm();
            const double dB = (vf-vfPB).norm();
            const double dC = (vf-vfPC).norm();

            std::cout << "@@ iff = " << iff << std::endl;
            std::cout << "@@ x,y,s = " << x << " " << y << " " << s << " " << type << std::endl;
            std::cout << "+A x,y,s = " << fPA[0] << " " << fPA[1] << " " << fPA[2] << " " << fPA[3] << std::endl;
            std::cout << "+B x,y,s = " << fPB[0] << " " << fPB[1] << " " << fPB[2] << " " << fPB[3] << std::endl;
            std::cout << "+C x,y,s = " << fPC[0] << " " << fPC[1] << " " << fPC[2] << " " << fPC[3] << std::endl;
            std::cout << dA << std::endl;
            std::cout << dB << std::endl;
            std::cout << dC << std::endl;

            // // for L=4
            // assert( dA<0.2 );
            // assert( dB<0.2 );
            // assert( dC<0.2 );
          }

          {
            const int type = ZY;
            const FaceCoords f{x,y,s,type};
            const Idx iff = dual.idx(f);

            FaceCoords fMA, fMB, fMC;
            dual.shiftMA( fMA, f );
            dual.shiftMB( fMB, f );
            dual.shiftMC( fMC, f );

            const V3 vf = dual.vertices[iff];

            if(fMA[2]!=PatchIdxVoid){
              const Idx ifMA = dual.idx(fMA);
              const V3 vfMA = dual.vertices[ifMA];
              const double dx = (vf-vfMA).norm();
              std::cout << "-A x,y,s = " << fMA[0] << " " << fMA[1] << " " << fMA[2] << " " << fMA[3] << std::endl;
              std::cout << dx << std::endl;
              // assert( dx<0.7 );
              // assert( dx<0.2 );
            }
            if(fMB[2]!=PatchIdxVoid){
              const Idx ifMB = dual.idx(fMB);
              const V3 vfMB = dual.vertices[ifMB];
              const double dy = (vf-vfMB).norm();
              std::cout << "-B x,y,s = " << fMB[0] << " " << fMB[1] << " " << fMB[2] << " " << fMB[3] << std::endl;
              std::cout << dy << std::endl;
              // assert( dy<0.7 );
              // assert( dy<0.2 );
            }
            if(fMC[2]!=PatchIdxVoid){
              const Idx ifMC = dual.idx(fMC);
              const V3 vfMC = dual.vertices[ifMC];
              const double dz = (vf-vfMC).norm();
              std::cout << "-C x,y,s = " << fMC[0] << " " << fMC[1] << " " << fMC[2] << " " << fMC[3] << std::endl;
              std::cout << dz << std::endl;
              // assert( dz<0.7 );
              // assert( dz<0.2 );
            }
          }
        }}}
  }


  {
    // BasePoints basePoints;
    // BaseTypes baseTypes;
    // dual.GetBasePoints( basePoints, baseTypes );

    for(Idx i=0; i<dual.basePoints.size(); i++){
      std::cout << "# basepoints" << std::endl;
      std::cout << "@@@ i = " << i << " type = " << dual.baseTypes[i] << std::endl;
      const FaceCoords f = dual.idx2FaceCoords( dual.basePoints[i] );
      std::cout << "x,y,s,type = " << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << std::endl;
    }


    // for dual
    Orbits orbits(dual.vertices,
                  dual.basePoints,
                  dual.baseTypes,
                  lattice, Ih, rot);
    // std::cout << "# orbits (dual)" << std::endl;
    // std::cout << orbits.print() << std::endl;

    const Vertices fs_old = dual.vertices;
    orbits.RefreshOrbits();
    const Vertices fs_new = dual.vertices;
    const double norm2 = squaredNorm(fs_old-fs_new);
    std::cout << "norm = " << std::sqrt( norm2/fs_old.size() ) << std::endl;
    assert( std::sqrt( norm2/fs_old.size() )<1.0e-14 );

    {
      DualLatticeOptimizerBase opt(dual, lattice, Ih, rot);

      for(Idx i=0; i<opt.basePoints.size(); i++){
        const Idx ir = opt.basePoints[i];
        std::cout << "type: " << opt.baseTypes[i] << std::endl;

        // const Coords coords = .idx2Coords( ir );
        const V3 r = dual.vertices[ir];
        // const XiCoords xi = opt.XiProjector(r);
        // std::cout << "xy: " << xi[0] << " " << xi[1] << std::endl;

        const V3 check = opt.embedXi( opt.XiProjector(r) );
        assert( (check-r).norm() < 1.0e-14 );
        const V3 check2 = opt.embedZeta( opt.ZetaProjector(r) );
        assert( (check2-r).norm() < 1.0e-14 );

        std::cout << "types: " << opt.is_type1(r)
                  << " " << opt.is_type2(r) << std::endl;
      }
    }

  }


  return 0;
}




  // // only defined for regular points
  // void shiftMX( Coords& nP,
  //               const Coords& n ) const {
  //   const int x = n[0];
  //   const int y = n[1];
  //   const int s = n[2];
  //   assert(0<=x && x<L);
  //   assert(0<=y && y<L);
  //   assert(0<=s && s<NPatches);

  //   if(x==0) { // steps over the patch bdy
  //     if(s<NPatches/2) np = Coords{L-1, y, s+4}; // southern hemisphere
  //     else np = Coords{L-1-y, L-1, s-1}; // northern hemisphere
  //   }
  //   else np = Coords{x-1, y, s}; // within patch
  // }

  // // only defined for regular points
  // void shiftMY( Coords& np,
  //               const Coords& n ) const {
  //   const int x = n[0];
  //   const int y = n[1];
  //   const int s = n[2];
  //   assert(0<=x && x<L);
  //   assert(0<=y && y<L);
  //   assert(0<=s && s<NPatches);

  //   if(y==0){ // steps over the patch bdy
  //     if(s<NPatches/2) np = Coords{L-1, L-1-x, s-1}; // southern hemisphere
  //     else np = Coords{x, L-1, s-5}; // northern hemisphere
  //   }
  //   else np = Coords{x, y-1, s}; // within patch
  // }

  // // only defined for regular points
  // void shiftMZ( Coords& np,
  //               const Coords& n ) const {
  //   const int x = n[0];
  //   const int y = n[1];
  //   const int s = n[2];
  //   assert(0<=x && x<L);
  //   assert(0<=y && y<L);
  //   assert(0<=s && s<NPatches);

  //   if(s<NPatches/2){ // southern hemisphere
  //     if(x==0) {
  //       if(y==0) np = Coords{-1, -1, PatchIdxVoid};
  //       else np = Coords{L-1, y-1, s+4};
  //     }
  //     else if(y==0) np = Coords{L-1, L-x, s-1};
  //     else np = Coords{x-1, y-1, s};
  //   }
  //   else { // northern hemisphere
  //     if(x==0) {
  //       if(y==0) np = Coords{-1, -1, PatchIdxVoid};
  //       else np = Coords{L-y, L-1, s-1};
  //     }
  //     else if(y==0) np = Coords{x-1, L-1, s-5};
  //     else np = Coords{x-1, y-1, s};
  //   }
  // }
