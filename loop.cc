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
  int k=3;
  if(argc>=2) k = atoi(argv[1]);
  int r=0;
  if(argc>=3) r = atoi(argv[2]);

  RefinedIcosahedron lattice(L);
  assert(lattice.NVertices()==N);
  FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  Rotation rot;
  RefinedIcosahedronDual dual(lattice);

  // std::cout << "len = " << arcLength( lattice.vertices[0], lattice.vertices[1] ) << std::endl;
  // std::cout << "len = " << arcLength( dual.vertices[0], dual.vertices[1] ) << std::endl;
  // std::cout << "len = " << arcLength( dual.vertices[0], dual.vertices[2] ) << std::endl;
  // std::cout << "len = " << arcLength( dual.vertices[0], dual.vertices[3] ) << std::endl;

  // std::cout << "faces = " << std::endl;
  // for(auto& face : dual.faces){
  //   std::cout << "length = " << face.size() << std::endl;
  //   for(auto& elem : face){
  //     FaceCoords f = dual.idx2FaceCoords(elem.first);
  //     std::cout << "debug. f = " << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << std::endl;
  //     std::cout << "debug. df = " << elem.second << std::endl;
  //   }
  // }


  // Orbits orbits(dual.vertices,
  //               dual.basePoints,
  //               dual.baseTypes,
  //               lattice, Ih, rot);

  DualSpinStructure spin(dual);

  // DualLoopUtils loops;
  {
    Idx q = std::pow(2, k) + r; // Binary: 00100101

    SimpConfig c( q );

    std::cout << "c = " << c << std::endl;
    // for(Idx in=0; in<N; in++) std::cout << c[in] << std::endl;
    // std::cout << std::endl;

    // std::vector<Idx> ups;
    std::set<Idx> ups;
    // 1:up, 0:down
    for(Idx in=0; in<N; in++) {
      if(c[in]==1) ups.insert(in);
    }

    std::set<Link> links_in_loops;
    for(const Idx in : ups){
      const Coords n = lattice.idx2Coords(in);
      // std::cout << "debug. in = " << in << " n = " << n[0] << " " << n[1] << " " << n[2] << std::endl;

      for(int mu=0; mu<6; mu++){
        // std::cout << "debug. mu = " << mu << std::endl;
        Coords np;
        lattice.shift( np, n, mu );
        if(np[2]==PatchIdxVoid) continue;
        const Idx inp = lattice.idx(np);

        if( !ups.contains(inp) ){ // (n, np) : broken
          // std::cout << "debug. inp = " << inp << " np = " << np[0] << " " << np[1] << " " << np[2] << std::endl;

          const DirectedLink duallink = dual.linkSimp2Dual.at( DirectedLink{in,mu} );
          // const FaceCoords fA = dual.idx2FaceCoords( duallink.first );
          const FaceCoords fA = dual.idx2FaceCoords( duallink.first );
          const int df = duallink.second;
          FaceCoords fB;
          dual.shift( fB, fA, df );

          // dual.get_dual_link( fA, fB, n, mu );
          // std::cout << "debug. fA = " << fA[0] << " " << fA[1] << " " << fA[2] << " " << fA[3] << std::endl;
          // std::cout << "debug. fB = " << fB[0] << " " << fB[1] << " " << fB[2] << " " << fB[3] << std::endl;

          const Idx ifA = dual.idx(fA);
          const Idx ifB = dual.idx(fB);
          links_in_loops.insert( Link{ifA, ifB} );
          // if(ifA<ifB) ;
          // else links_in_loops.insert( Link{ifB, ifA} );
        }
      }
    }
    std::cout << "links_in_loops = " << std::endl;
    for(auto elem : links_in_loops) std::cout << elem.first << " " << elem.second << std::endl;
    std::cout << std::endl;

    // std::set<Idx> sites_in_loops;
    // for(auto elem : links_in_loops) {
    //   sites_in_loops.insert(elem.first);
    //   sites_in_loops.insert(elem.second);
    // }

    // group links
    using Loop = std::vector<Idx>; //  link;
    std::vector<Loop> loops;


    while(links_in_loops.size()!=0){
      // initialize a loop
      Loop loop;
      const Link ell = *links_in_loops.begin();
      const Idx if_end = ell.first, if_init = ell.second;

      // retrieve a loop starting from if_init; coming back to if_end
      Idx if_current = if_init;
      loop.push_back(if_current);

      while(*std::prev(loop.end())!=if_end){
        for(auto em=std::next(links_in_loops.begin()); em!=links_in_loops.end(); em++){
          // std::cout << "em = " << em->first << " " << em->second << std::endl; // always updated
          if(em->first==if_current){
            if_current = em->second;
            loop.push_back(if_current);
            // std::cout << "debug. insert second" << std::endl;
            links_in_loops.erase(em);
            break;
          }
          else if(em->second==if_current){
            if_current = em->first;
            loop.push_back(if_current);
            // std::cout << "debug. insert first" << std::endl;
            links_in_loops.erase(em);
            break;
          }
        } // end for
        //  std::cout << "debug. if_current = " << if_current << std::endl;
      }
      links_in_loops.erase(ell);
      loops.push_back(loop);
    }

    std::cout << "loops = " << std::endl;
    for(auto loop : loops){
      for(auto elem : loop) std::cout << elem << " ";
      std::cout << std::endl;
    }

  }

  return 0;
}
