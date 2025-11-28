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





class DualLatticeAngleCostEvaluator{
private:
  std::vector<double> coordsType1;
  std::vector<double> coordsType2;
  std::vector<std::array<double,2>> coordsType4;

public:
  RefinedIcosahedronDual& dual;
  Orbits orbits;

  BasePoints basePoints;
  BaseTypes baseTypes;

  int nType1;
  int nType2;
  int nType4;

  // pointers collecting all the coordinates
  std::vector<double*> p;

  DualLatticeAngleCostEvaluator
  (
   RefinedIcosahedronDual& dual_,
   const Icosahedron& icos_,
   const FullIcosahedralGroup& Ih_,
   const Rotation& rot_
   )
    : dual(dual_)
    , orbits( dual.vertices,
              dual.basePoints,
              dual.baseTypes,
              icos_, Ih_, rot_ )
  {
    nType1=0, nType2=0, nType4=0;
    for(const int type : dual.baseTypes){
      if(type==1) nType1++;
      else if(type==2) nType2++;
      else if(type==4) nType4++;
      else assert(false);
    }

    coordsType1.clear();
    coordsType1.resize(nType1);
    coordsType2.clear();
    coordsType2.resize(nType2);
    coordsType4.clear();
    coordsType4.resize(nType4);

    FillBaseCoordsConvenient();
    p.clear();
    for(Idx i=0; i<coordsType1.size(); i++) p.push_back( coordsType1.data()+i );
    for(Idx i=0; i<coordsType2.size(); i++) p.push_back( coordsType2.data()+i );
    for(Idx i=0; i<coordsType4.size(); i++) {
      p.push_back( coordsType4[i].data() );
      p.push_back( coordsType4[i].data()+1 );
    }

    FillBasePointsFromCoords();
    orbits.RefreshOrbits();
  }


  void FillBaseCoordsConvenient(){
    // vertices.clear();

    // Idx counter=0;
    // for(int s=0; s<NPatches; s++){
    //   for(int y=0; y<L; y++){
    //     for(int x=0; x<L; x++){
    //       const Coords n{x,y,s};
    //       const Idx in = simplicial.idx(n);

    //       Coords nPX, nPY, nPZ;
    //       simplicial.shiftPX( nPX, n );
    //       simplicial.shiftPY( nPY, n );
    //       simplicial.shiftPZ( nPZ, n );
    //       const Idx inPX = simplicial.idx(nPX);
    //       const Idx inPY = simplicial.idx(nPY);
    //       const Idx inPZ = simplicial.idx(nPZ);

    //       const V3 vn = simplicial.vertices[in];
    //       const V3 vnPX = simplicial.vertices[inPX];
    //       const V3 vnPY = simplicial.vertices[inPY];
    //       const V3 vnPZ = simplicial.vertices[inPZ];

    //       // xz face
    //       V3 r = project( circumcenter(vn, vnPX, vnPZ) );
    //       assert(vertices.size()==idx(FaceCoords{x,y,s,XZ}));
    //       vertices.push_back( r );

    //       // zy face
    //       r = project( circumcenter(vn, vnPZ, vnPY ) );
    //       assert(vertices.size()==idx(FaceCoords{x,y,s,ZY}));
    //       vertices.push_back( r );
    //     }}}

    // assert( vertices.size()==NVertices() );
  }

  void FillBasePointsFromCoords(){

  }


  double cost() const {
    double sqrd = 0.0;

    for(const Idx iff : dual.basePoints){

      const V3 rf = dual.vertices[iff];
      const FaceCoords nf = dual.idx2FaceCoords( iff );

      FaceCoords nfA, nfB, nfC;
      dual.shiftPA( nfA, nf );
      dual.shiftPB( nfB, nf );
      dual.shiftPC( nfC, nf );

      const V3 rfA = dual.vertices[dual.idx(nfA)];
      const V3 rfB = dual.vertices[dual.idx(nfB)];
      const V3 rfC = dual.vertices[dual.idx(nfC)];

      const double A0B = sphericalAngle( rfA, rf, rfB );
      const double B0C = sphericalAngle( rfB, rf, rfC );
      const double C0A = sphericalAngle( rfC, rf, rfA );

      assert( std::abs( A0B+B0C+C0A-M_PI ) < 1.0e-14 );

      const double diff = A0B - B0C;
      sqrd += diff*diff;
    }
    return std::sqrt( sqrd/dual.basePoints.size() );
  }

  

};




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
