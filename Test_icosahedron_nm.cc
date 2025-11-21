    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/debug/Test_icosahedron.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <Grid/Grid.h>
#include <Grid/stencil/IcosahedralStencil.h>


using Double = double;
using Idx = std::int32_t;
// using Idx = long int;
using Complex = std::complex<double>;

using Link = std::array<Idx,2>; // <int,int>;
using Face = std::vector<Idx>;

using MS=Eigen::Matrix2cd;
using VD=Eigen::Vector2d;
using VE=Eigen::Vector3d;
using VC=Eigen::VectorXcd;

static constexpr int NS = 2;
static constexpr int DIM = 2;
static constexpr Complex I = Complex(0.0, 1.0);


constexpr int FOUR=4;
constexpr int THREE=3;


// #define IS_DUAL
// #define IS_OVERLAP

// #define IsVerbose
// #define InfoForce
// #define InfoDelta

namespace Comp{
  constexpr bool is_compact=true;

  // d_DW.update() is always done independently
#ifdef IS_OVERLAP
  constexpr int NPARALLEL_DUPDATE=1;
  constexpr int NPARALLEL=12; // 12
  constexpr int NSTREAMS=4; // 4
#else
  constexpr int NPARALLEL_DUPDATE=12;
  constexpr int NPARALLEL=1; // 12
  constexpr int NSTREAMS=12; // for grad loop
#endif
  constexpr int NPARALLEL_GAUGE=12; // 12
  constexpr int NPARALLEL_SORT=12; // 12

  constexpr int N_REFINE=4;
  constexpr int NS=2;

  constexpr int Nt=16;
  // constexpr int Nt=24;
  // constexpr int Nt=48; // add 2
  // constexpr int Nt=64;
  // constexpr int Nt=96; // add 4
  // constexpr int Nt=120;
  // constexpr int Nt=144; // add 8
  // constexpr int Nt=160;

  // constexpr int Nt=24;
  // constexpr int Nt=192;
  // constexpr int Nt=1;
  // constexpr int Nt=16;

#ifdef IS_DUAL
  constexpr Idx N_SITES=20*N_REFINE*N_REFINE;
#else
  constexpr Idx N_SITES=10*N_REFINE*N_REFINE+2;
  constexpr int N_LINKS=30*N_REFINE*N_REFINE; // 30, 120, 480
#endif

  constexpr Idx Nx=NS*N_SITES; // matrix size of DW
  constexpr Idx N=Nx*Nt; // matrix size of DW

  const double TOL_INNER=1.0e-15;
  const double TOL_OUTER=1.0e-14;
}

const std::string dir = "/Users/nobuyukimatsumoto/grid4/qed3-1.0/geometry/data/";

// /Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header
#include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/s2n_simp.h"
// #include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/s2n_dual.h"
// #include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/rng.h"
#include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/valence.h"
#include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/gauge_ext.h"
#include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/action_ext.h"



// #include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/sparse_matrix.h"
class COOEntry;
#include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/geometry/geodesic.h"
#include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/dirac_base.h"
#include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/dirac_simp.h"
#include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/dirac_ext.h"
// #include <cuComplex.h>
// #include <cuda_runtime.h>
// #include <cublas_v2.h>
// #include <cublas_api.h>
// #include <cusolverDn.h>
// using CuC = cuDoubleComplex;
// #include "/Users/nobuyukimatsumoto/grid4/qed3-1.0/main/header/gpu_header.h"


using GridCoord = std::array<int,3>;

struct SphereUtils{
  std::vector<VE> IcosSites;
  std::vector<VE> RefinedIcosSites;

  const int L;

  // static constexpr auto comp = [](const GridCoord& a, const GridCoord& b){
  //   if(a[0]!=b[0]) return a[0]<b[0];
  //   else if(a[1]!=b[1]) return a[1]<b[1];
  //   else if(a[2]!=b[2]) return a[2]<b[2];
  //   else return a[3]<b[3];
  // };

  // ~~~~~~~~~~~~~~~~~~~
  // IcosIdx and Idx are arranged with the original convention of Evan's:
  //
  // IcosIdx: 1 dim idx for the icosahedron
  // pole exceptions for the Patch values: S = IcosahedralPatches = 10, N = 11;
  //
  // Idx: 1dim idx for the refined icosahedron
  // ~~~~~~~~~~~~~~~~~~~

  // Patch <-> IcosIdx
  std::map<const int, const int> Patch2IcosIdx;
  std::map<const int, const int> IcosIdx2Patch;

  // X, Y, D shifts for IcosIdx
  std::map<const int, const int> IcosIdxNbrForPlusX;
  std::map<const int, const int> IcosIdxNbrForPlusY;
  std::map<const int, const int> IcosIdxNbrForPlusDiagonal;

  // GridCoord -> R3Coord
  std::map<const GridCoord, const VE> GridCoord2R3Coord;
  // std::map<const GridCoord, const VE, decltype(comp)> GridCoord2R3Coord;

  // GridCoord <-> Idx
  std::map<const GridCoord, const Idx> GridCoord2Idx;
  std::map<const Idx, const GridCoord> Idx2GridCoord;

  Idx IdxN;
  Idx IdxS;

  SphereUtils( const int n_refine )
    : L(n_refine)
    // , GridCoord2R3Coord(comp)
    // , GridCoord2Idx(comp)
  {
    SetPatch2IcosIdx();
    SetIcosIdx2Patch();
    SetIcosIdxNbrForPlusX();
    SetIcosIdxNbrForPlusY();
    SetIcosIdxNbrForPlusDiagonal();
    ReadIcosCoords( IcosSites, 1 );
    ReadIcosCoords( RefinedIcosSites, n_refine );
    SetGridCoord2R3Coord();
    SetGridCoord2Idx();
    SetIdx2GridCoord();
    SetIdxPoles();
  }

  void SetPatch2IcosIdx(){
    // {Patch, IcosIdx}
    Patch2IcosIdx.insert( { 0, 7 } );
    Patch2IcosIdx.insert( { 1, 10 } );
    Patch2IcosIdx.insert( { 2, 1 } );
    Patch2IcosIdx.insert( { 3, 2 } );
    Patch2IcosIdx.insert( { 4, 4 } );
    Patch2IcosIdx.insert( { 5, 3 } );
    Patch2IcosIdx.insert( { 6, 5 } );
    Patch2IcosIdx.insert( { 7, 8 } );
    Patch2IcosIdx.insert( { 8, 9 } );
    Patch2IcosIdx.insert( { 9, 11 } );
    Patch2IcosIdx.insert( { 10, 6 } ); // S
    Patch2IcosIdx.insert( { 11, 0 } ); // N
  }

  void SetIcosIdx2Patch(){
    for(auto it = Patch2IcosIdx.begin(); it != Patch2IcosIdx.end(); it++){
      IcosIdx2Patch.insert( {it->second, it->first} );
    }
  }

  void SetIcosIdxNbrForPlusX(){
    // {IcosIdx, IcosIdx}
    IcosIdxNbrForPlusX.insert( { 3, 10 } );
    IcosIdxNbrForPlusX.insert( { 5, 1 } );
    IcosIdxNbrForPlusX.insert( { 8, 2 } );
    IcosIdxNbrForPlusX.insert( { 9, 4 } );
    IcosIdxNbrForPlusX.insert( { 11, 7 } );
    IcosIdxNbrForPlusX.insert( { 7, 6 } );
    IcosIdxNbrForPlusX.insert( { 10, 6 } );
    IcosIdxNbrForPlusX.insert( { 1, 6 } );
    IcosIdxNbrForPlusX.insert( { 2, 6 } );
    IcosIdxNbrForPlusX.insert( { 4, 6 } );
  }

  void SetIcosIdxNbrForPlusY(){
    // {IcosIdx, IcosIdx}
    IcosIdxNbrForPlusY.insert( { 3, 0 } );
    IcosIdxNbrForPlusY.insert( { 5, 0 } );
    IcosIdxNbrForPlusY.insert( { 8, 0 } );
    IcosIdxNbrForPlusY.insert( { 9, 0 } );
    IcosIdxNbrForPlusY.insert( { 11, 0 } );

    IcosIdxNbrForPlusY.insert( { 7, 3 } );
    IcosIdxNbrForPlusY.insert( { 10, 5 } );
    IcosIdxNbrForPlusY.insert( { 1, 8 } );
    IcosIdxNbrForPlusY.insert( { 2, 9 } );
    IcosIdxNbrForPlusY.insert( { 4, 11 } );
  }

  void SetIcosIdxNbrForPlusDiagonal(){
    // {IcosIdx, IcosIdx}
    IcosIdxNbrForPlusDiagonal.insert( { 7, 10 } );
    IcosIdxNbrForPlusDiagonal.insert( { 10, 1 } );
    IcosIdxNbrForPlusDiagonal.insert( { 1, 2 } );
    IcosIdxNbrForPlusDiagonal.insert( { 2, 4 } );
    IcosIdxNbrForPlusDiagonal.insert( { 4, 7 } );

    IcosIdxNbrForPlusDiagonal.insert( { 3, 5 } );
    IcosIdxNbrForPlusDiagonal.insert( { 5, 8 } );
    IcosIdxNbrForPlusDiagonal.insert( { 8, 9 } );
    IcosIdxNbrForPlusDiagonal.insert( { 9, 11 } );
    IcosIdxNbrForPlusDiagonal.insert( { 11, 3 } );
  }


  // fill sites = this->**Sites
  void ReadIcosCoords( std::vector<VE>& sites,
                       const int n_refine ) const {
    std::cout << "# reading simplicial points" << std::endl;
    std::ifstream file(dir+"pts_n"+std::to_string(n_refine)+".dat");

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      double v1, v2, v3;
      iss >> v1;
      iss >> v2;
      iss >> v3;
      sites.push_back( VE(v1, v2, v3) );
    }
  }

  void SetGridCoord2R3Coord(){
    for(int Patch=0; Patch<Grid::IcosahedralPatches; Patch++){
      const int ix0 = Patch2IcosIdx.at(Patch);
      const VE x0 = IcosSites[ix0];

      const int ix1 = IcosIdxNbrForPlusX.at(ix0);
      const VE x1 = IcosSites[ix1];
      const int ix2 = IcosIdxNbrForPlusDiagonal.at(ix0);
      const VE x2 = IcosSites[ix2];
      {
        //      x2
        //      /|
        //     / |
        //    /  |
        //   /   |
        //  /    |
        // x0----x1
        for(int x=0; x<L; x++){
          const double rx = 1.0*x/L;
          for(int y=0; y<=x; y++){
            const double ry = 1.0*y/L;
            const VE z = x0 + rx*(x1-x0) + ry*(x2-x1);
            const VE r = z/z.norm();
            GridCoord2R3Coord.insert( { GridCoord({x,y,Patch}), r } );
          }}
      }
      const int ix3 = IcosIdxNbrForPlusY.at(ix0);
      const VE x3 = IcosSites[ix3];
      {
        // x3----x2
        // |     /
        // |    /
        // |   /
        // |  /
        // | /
        // x0
        for(int y=0; y<L; y++){
          const double ry = 1.0*y/L;
          for(int x=0; x<y; x++){
            const double rx = 1.0*x/L;
            const VE z = x0 + rx*(x2-x3) + ry*(x3-x0);
            const VE r = z/z.norm();
            GridCoord2R3Coord.insert( { GridCoord({x,y,Patch}), r } );
          }}
      }
    }
    assert( GridCoord2R3Coord.size()==RefinedIcosSites.size()-2 ); // no N/S
  }

  void SetGridCoord2Idx(){
    for(auto it = GridCoord2R3Coord.begin(); it != GridCoord2R3Coord.end(); it++){
      const GridCoord x_Grid = it->first;
      const VE x_R3 = it->second;
      for(Idx idx=0; idx<RefinedIcosSites.size(); idx++){
        const VE x_prime = RefinedIcosSites[idx];

        if( (x_prime-x_R3).norm()<1.0e-8 ) {
          GridCoord2Idx.insert( { x_Grid, idx } );
          break;
        }

        if(idx==RefinedIcosSites.size()-1){
          std::cout << "!!!! not found !!!: "
                    << "x_Grid = " << x_Grid[0]
                    << " " << x_Grid[1]
                    << " " << x_Grid[2] << std::endl;
          // << "x_R3 = " << x_R3 << std::endl
        }
      }
    }

    // std::cout << "GridCoord2Idx.size() = " << GridCoord2Idx.size() << std::endl;
    // std::cout << "RefinedIcosSites.size() = " << RefinedIcosSites.size() << std::endl;
    assert( GridCoord2Idx.size()==RefinedIcosSites.size()-2 ); // no N/S
  }

  void SetIdx2GridCoord(){
    for(auto it = GridCoord2Idx.begin(); it != GridCoord2Idx.end(); it++){
      Idx2GridCoord.insert( {it->second, it->first} );
    }
  }

  void SetIdxPoles(){
    { // N
      VE x_R3; x_R3 << 0,0,1;
      for(Idx idx=0; idx<RefinedIcosSites.size(); idx++){
        const VE x_prime = RefinedIcosSites[idx];
        if( (x_prime-x_R3).norm()<1.0e-10 ) {
          IdxN = idx;
          break;
        }
      }
    }
    { // S
      VE x_R3; x_R3 << 0,0,-1;
      for(Idx idx=0; idx<RefinedIcosSites.size(); idx++){
        const VE x_prime = RefinedIcosSites[idx];
        if( (x_prime-x_R3).norm()<1.0e-10 ) {
          IdxS = idx;
          break;
        }
      }
    }
  }



};


using namespace std;
using namespace Grid;




const int MyNd=3;

class IcosahedralGimpl
{
public:
  typedef LatticeLorentzColourMatrix  GaugeField;
  typedef LatticeColourMatrix         GaugeLinkField;
  typedef LatticeDoubledGaugeField    DoubledGaugeField;
  typedef LatticeComplex              ComplexField;
};
  
template< class Gimpl>
class IcosahedralSupport {
public:
  // Move to inherit types macro as before
  typedef typename Gimpl::GaugeField GaugeField;
  typedef typename Gimpl::GaugeLinkField GaugeLinkField;
  typedef typename Gimpl::ComplexField ComplexField;
  typedef typename Gimpl::DoubledGaugeField DoubledGaugeField;
  //
  GridBase *VertexGrid;
  GridBase *EdgeGrid;

  IcosahedralStencil FaceStencil;
  IcosahedralStencil NNee; // edge neighbours with edge domain
  IcosahedralStencil NNev; // vertex neighbours but in edge domain
  IcosahedralStencil NNvv; // vertex neighbours with vertex domain
  IcosahedralStencil NNve; // edge neighbours with vertex domain
  
  IcosahedralSupport(GridBase *_VertexGrid,GridBase *_EdgeGrid)
    : FaceStencil (_EdgeGrid,_VertexGrid),
      NNee(_EdgeGrid,_VertexGrid),
      NNev(_EdgeGrid,_VertexGrid),
      NNve(_EdgeGrid,_VertexGrid),
      NNvv(_EdgeGrid,_VertexGrid),
      VertexGrid(_VertexGrid), EdgeGrid(_EdgeGrid)
  {
    FaceStencil.FaceStencil();
    std::cout << "NNee"<<std::endl;
    // true/false is "vertexInput, VertexOutput"
    NNee.NearestNeighbourStencil(false,false);// edge input + output      ; used by face stencil
    NNev.NearestNeighbourStencil(true,false); // vertex input, edge ouput ; used by gauge transform
    NNvv.NearestNeighbourStencil(true,true);  // vertex input + output    ; used by Laplacian
    NNve.NearestNeighbourStencil(false,true); // edge input, vertex output; used by double store
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // Gauge Link field GT is the gauge transform and lives on the VERTEX field
  ////////////////////////////////////////////////////////////////////////////////////
  void ForwardTriangles(GaugeField &Umu,LatticeComplex &plaq1,LatticeComplex &plaq2)
  {
    {
      autoView(Umu_v,Umu,AcceleratorRead);
      autoView(plaq1_v,plaq1,AcceleratorWrite);
      autoView(plaq2_v,plaq2,AcceleratorWrite);
      autoView(stencil_v,FaceStencil,AcceleratorRead);

      accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{
	
	auto Lx = Umu_v(ss)(IcosahedronPatchX);
	auto Ly = Umu_v(ss)(IcosahedronPatchY);
	auto Ld = Umu_v(ss)(IcosahedronPatchDiagonal);

	// for trace [ U_x(z) U_y(z+\hat x) adj(U_d(z)) ]
	{
	  auto SE1 = stencil_v.GetEntry(0,ss);
	  auto doAdj = SE1->_adjoint;
	  auto pol   = SE1->_polarisation;
	  auto s1    = SE1->_offset;
	  auto L1 = Umu_v(s1)(pol);
	  if(doAdj)
	    L1 = adj(L1);

	  coalescedWrite(plaq1_v[ss](),trace(Lx*L1*adj(Ld) ) );
	}

	// for trace [  U_y(z) U_x(z+\hat y) adj(U_d(z))   ]
	{
	  auto SE2 = stencil_v.GetEntry(1,ss);
	  auto doAdj = SE2->_adjoint;
	  auto pol   = SE2->_polarisation;
	  auto s2    = SE2->_offset;
	  auto L2 = Umu_v(s2)(pol);
	  if(doAdj)
	    L2 = adj(L2);
	  coalescedWrite(plaq2_v[ss](),trace(Ly*L2*adj(Ld) ) );
	}
      });
    }
  }

  // Staples for gauge force
  void IcosahedralStaples(GaugeField &Umu,
			  GaugeLinkField &stapleXY,
			  GaugeLinkField &stapleYX,
			  GaugeLinkField &stapleXD,
			  GaugeLinkField &stapleDX,
			  GaugeLinkField &stapleYD,
			  GaugeLinkField &stapleDY)
  {
      autoView(Umu_v,Umu,AcceleratorRead);
      autoView(stapleXY_v,stapleXY,AcceleratorWrite);
      autoView(stapleYX_v,stapleYX,AcceleratorWrite);
      autoView(stapleXD_v,stapleXD,AcceleratorWrite);
      autoView(stapleDX_v,stapleDX,AcceleratorWrite);
      autoView(stapleYD_v,stapleYD,AcceleratorWrite);
      autoView(stapleDY_v,stapleDY,AcceleratorWrite);
      autoView(stencil_v,NNee,AcceleratorRead);

      const int np = NNee._npoints;
      
      const int ent_Xp = 0;
      const int ent_Yp = 1;
      const int ent_Xm = 3;
      const int ent_Ym = 4;

      accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{

	// Three forward links from this site
	auto Lx = Umu_v(ss)(IcosahedronPatchX);
	auto Ly = Umu_v(ss)(IcosahedronPatchY);
	auto Ld = Umu_v(ss)(IcosahedronPatchDiagonal);

	///////////////////////////////////////////////////////////////////
	// Terms for the staple orthog to PlusDiagonal
	///////////////////////////////////////////////////////////////////
	// adj( U_y(z+\hat x)) adj(U_x(z)) 
	{
	  auto SE1 = stencil_v.GetEntry(ent_Xp,ss);
	  auto doAdj = SE1->_adjoint;
	  auto pol   = SE1->_polarisation;
	  auto s1    = SE1->_offset;
	  auto Ly_at_xp = Umu_v(s1)(pol);
	  if(doAdj)
	    Ly_at_xp = adj(Ly_at_xp);
	  coalescedWrite(stapleXY_v[ss](),adj(Ly_at_xp)*adj(Lx)  );
	}
	// adj( U_y(z) )  adj(U_x(z+\hat y))
	{
	  auto SE2 = stencil_v.GetEntry(ent_Yp,ss);
	  auto doAdj = SE2->_adjoint;
	  auto pol   = SE2->_polarisation;
	  auto s2    = SE2->_offset;
	  auto Lx_at_yp = Umu_v(s2)(pol);
	  if(doAdj)
	    Lx_at_yp = adj(Lx_at_yp);

	  coalescedWrite(stapleYX_v[ss](),adj(Lx_at_yp)*adj(Ly) );
	}
	///////////////////////////////////////////////////////////////////
	// Terms for the staple covering Xp : Dp Yp(x++) and Yp(y--) Dp(y--)
	///////////////////////////////////////////////////////////////////
	// U_y(z+\hat x)*adj(U_d(z)) 
	{
	  auto SE1 = stencil_v.GetEntry(ent_Xp,ss);
	  auto doAdj = SE1->_adjoint;
	  auto pol   = SE1->_polarisation;
	  auto s1    = SE1->_offset;
	  auto Ly_at_xp = Umu_v(s1)(pol);
	  if(doAdj)
	    Ly_at_xp = adj(Ly_at_xp);
	  coalescedWrite(stapleDY_v[ss](), Ly_at_xp *adj(Ld));
	}
	//  adj(U_d(z-\hat y)) U_y(z-\hat y)
	{
	  auto SE2 = stencil_v.GetEntry(ent_Ym,ss);
	  auto doAdj = SE2->_adjoint;
	  auto pol   = SE2->_polarisation;
	  auto s2    = SE2->_offset;
	  int pol1 = IcosahedronPatchY;
	  int pol2 = IcosahedronPatchDiagonal;
	  if ( pol != IcosahedronPatchDiagonal ) {
	    pol1 = IcosahedronPatchDiagonal;
	    pol2 = IcosahedronPatchX;
	  }
	  auto Ly_at_ym = Umu_v(s2)(pol1);
	  auto Ld_at_ym = Umu_v(s2)(pol2);
	  coalescedWrite(stapleYD_v[ss](),adj(Ld_at_ym)*Ly_at_ym );
	}
	///////////////////////////////////////////////////////////////////
	// Terms for the staple covering Yp : Dp Xp(y++) and Xp(x--) Dp(x--)
	///////////////////////////////////////////////////////////////////
	//  U_x(z+\hat y)adj(U_d(z))
	{
	  auto SE1 = stencil_v.GetEntry(ent_Yp,ss);
	  auto doAdj = SE1->_adjoint;
	  auto pol   = SE1->_polarisation;
	  auto s1    = SE1->_offset;
	  auto Lx_at_yp = Umu_v(s1)(pol);
	  if(doAdj)
	    Lx_at_yp = adj(Lx_at_yp);
	  coalescedWrite(stapleDX_v[ss](),Lx_at_yp*adj(Ld) );
	}

	//  adj(U_d(z-\hat x))U_x(z-\hat x)
	{
	  auto SE2 = stencil_v.GetEntry(ent_Xm,ss);
	  auto doAdj = SE2->_adjoint;
	  auto pol   = SE2->_polarisation;
	  auto s2    = SE2->_offset;
	  int pol1 = IcosahedronPatchX;
	  int pol2 = IcosahedronPatchDiagonal;
	  if ( pol != IcosahedronPatchDiagonal ) {
	    pol1 = IcosahedronPatchDiagonal;
	    pol2 = IcosahedronPatchY;
	  }
	  auto Ly = Umu_v(ss)(IcosahedronPatchY);
	  auto Lx_at_xm = Umu_v(s2)(pol1);
	  auto Ld_at_xm = Umu_v(s2)(pol2);
	  coalescedWrite(stapleXD_v[ss](),adj(Ld_at_xm)*Lx_at_xm );
	}
      });
  }
  
  template<class MatterField>
  void Laplacian(MatterField &in,MatterField &out)
  {
    autoView(out_v,out,AcceleratorWrite);
    autoView(in_v,in,AcceleratorRead);
    autoView(stencil_v,NNvv,AcceleratorRead);

    const int np = NNvv._npoints;

    const int ent_Xp = 0;
    const int ent_Yp = 1;
    const int ent_Dp = 2;
    const int ent_Xm = 3;
    const int ent_Ym = 4;
    const int ent_Dm = 5;
    
    accelerator_for(ss,VertexGrid->oSites(),vComplex::Nsimd(),{
	  
      auto SE = stencil_v.GetEntry(ent_Xp,ss);
      uint64_t xp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Yp,ss);
      uint64_t yp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Dp,ss);
      uint64_t dp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Xm,ss);
      uint64_t xm_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Ym,ss);
      uint64_t ym_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Dm,ss);
      uint64_t dm_idx = SE->_offset;
      int missingLink = SE->_missing_link;
      
      auto i    = in_v(ss)();
      auto inxp = in_v(xp_idx)();
      auto inyp = in_v(yp_idx)();
      auto indp = in_v(dp_idx)();
      auto inxm = in_v(xm_idx)();
      auto inym = in_v(ym_idx)();
      auto indm = in_v(dm_idx)();

      auto o    = out_v(ss)();

      if ( missingLink ) {
	o = (1.0/5.0)*(inxp+inyp+indp+inxm+inym)-o;
      } else {
	o = (1.0/6.0)*(inxp+inyp+indp+inxm+inym+indm)-o;
      }
      
      coalescedWrite(out_v[ss](),o);
      });
  }


  template<class MatterField>
  void CovariantLaplacian(MatterField &in,MatterField &out,DoubledGaugeField &Uds)
  {
    autoView(out_v,out,AcceleratorWrite);
    autoView(in_v,in,AcceleratorRead);
    autoView(U_v,Uds,AcceleratorRead);
    autoView(stencil_v,NNvv,AcceleratorRead);

    const int np = NNvv._npoints;

    const int ent_Xp = 0;
    const int ent_Yp = 1;
    const int ent_Dp = 2;
    const int ent_Xm = 3;
    const int ent_Ym = 4;
    const int ent_Dm = 5;
    
    accelerator_for(ss,VertexGrid->oSites(),vComplex::Nsimd(),{
	  
      auto SE = stencil_v.GetEntry(ent_Xp,ss);
      uint64_t xp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Yp,ss);
      uint64_t yp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Dp,ss);
      uint64_t dp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Xm,ss);
      uint64_t xm_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Ym,ss);
      uint64_t ym_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Dm,ss);
      uint64_t dm_idx = SE->_offset;
      int missingLink = SE->_missing_link;
      
      auto i    = in_v(ss)();
      auto inxp = in_v(xp_idx)();
      auto inyp = in_v(yp_idx)();
      auto indp = in_v(dp_idx)();
      auto inxm = in_v(xm_idx)();
      auto inym = in_v(ym_idx)();
      auto indm = in_v(dm_idx)();

      inxp = U_v(ss)(0)*inxp;
      inyp = U_v(ss)(1)*inyp;
      indp = U_v(ss)(2)*indp;
      inxm = U_v(ss)(3)*inxm;
      inym = U_v(ss)(4)*inym;
      
      auto o    = i;

      if ( missingLink ) {
	o = (1.0/5.0)*(inxp+inyp+indp+inxm+inym)-i;
      } else {
	indm = U_v(ss)(5)*indm;
	o = (1.0/6.0)*(inxp+inyp+indp+inxm+inym+indm)-i;
      }
      
      coalescedWrite(out_v[ss](),o);
      });
  }
  
  void DoubleStore(GaugeField &U,DoubledGaugeField &Uds)
  {
    assert(U.Grid()==EdgeGrid);
    assert(Uds.Grid()==VertexGrid);
    autoView(Uds_v,Uds,AcceleratorWrite);
    autoView(U_v,U,AcceleratorRead);
    autoView(stencil_v,NNvv,AcceleratorRead);

    // Vertex result
    // Edge valued input
    // Might need an extra case (?)
    const int np = NNvv._npoints;

    const int ent_Xm = 3;
    const int ent_Ym = 4;
    const int ent_Dm = 5;
    
    accelerator_for(ss,VertexGrid->CartesianOsites(),vComplex::Nsimd(),{
	
      // Three local links
      {
	auto Lx = U_v(ss)(IcosahedronPatchX);
	auto Ly = U_v(ss)(IcosahedronPatchY);
	auto Ld = U_v(ss)(IcosahedronPatchDiagonal);
	coalescedWrite(Uds_v[ss](IcosahedronPatchX),Lx);
	coalescedWrite(Uds_v[ss](IcosahedronPatchY),Ly);
	coalescedWrite(Uds_v[ss](IcosahedronPatchDiagonal),Ld);
      }	
      // Three backwards links
      {
	auto SE = stencil_v.GetEntry(ent_Xm,ss);
	auto pol = SE->_polarisation;
	auto s   = SE->_offset;
	int pol1 = IcosahedronPatchX;
	//
	// Xm link is given diagonal unless HemiPatch changes in Northern hemisphere
	// But here for double storing we need the reverse link so
	// return either the Xplus or Diag plus direction
	//
	if ( pol != IcosahedronPatchDiagonal ) {
	  pol1 = IcosahedronPatchDiagonal;
	}
	auto Lx_at_xm = U_v(s)(pol1);
	Lx_at_xm = adj(Lx_at_xm);
	coalescedWrite(Uds_v[ss](3),Lx_at_xm );
      }
      {
	auto SE   = stencil_v.GetEntry(ent_Ym,ss);
	auto pol  = SE->_polarisation;
	auto s    = SE->_offset;
	//
	// Ym link is given diagonal unless HemiPatch changes in the Southern hemisphere
	// But here for double storing we need the reverse link so
	// return either the Yplus or Diag plus direction
	// 
	int pol1 = IcosahedronPatchY;
	if ( pol != IcosahedronPatchDiagonal ) {
	  pol1 = IcosahedronPatchDiagonal;
	}
	auto Ly_at_ym = U_v(s)(pol1);
	Ly_at_ym = adj(Ly_at_ym);
	coalescedWrite(Uds_v[ss](4),Ly_at_ym );
      }
      int missingLink;
      {      
	auto SE = stencil_v.GetEntry(ent_Dm,ss);
	auto s  = SE->_offset;
	// Dm link given the correct value for this use case
	auto pol= SE->_polarisation;
	missingLink = SE->_missing_link;
	if ( ! missingLink ) {
	  auto Ld_at_dm  =  U_v(s)(pol);
	  Ld_at_dm = adj(Ld_at_dm);
	  coalescedWrite(Uds_v[ss](5),Ld_at_dm );
	}
      }
    });
    auto pole_sites = VertexGrid->oSites() - VertexGrid->CartesianOsites();
    auto pole_offset= VertexGrid->CartesianOsites();
    
    accelerator_for(ss,pole_sites,vComplex::Nsimd(),{
	for(int p=0;p<HemiPatches;p++){ // Neighbours of each pole
	auto SE = stencil_v.GetEntry(p,pole_offset+ss);
	auto s  = SE->_offset;
	auto pol= SE->_polarisation;
	auto Link =  adj(U_v(s)(pol));
	coalescedWrite(Uds_v[pole_offset+ss](p),Link);
      }
    });
  }
  void  GaugeTransform(GaugeLinkField &gt, GaugeField &Umu)
  {
    autoView(Umu_v,Umu,AcceleratorWrite);
    autoView(g_v,gt,AcceleratorRead);
    autoView(stencil_v,NNev,AcceleratorRead);

    const int np = NNev._npoints;
    
    const int ent_Xp = 0;
    const int ent_Yp = 1;
    const int ent_Dp = 2;

    accelerator_for(ss,EdgeGrid->oSites(),vComplex::Nsimd(),{

      // Three forward links from this site
      auto Lx = Umu_v(ss)(IcosahedronPatchX);
      auto Ly = Umu_v(ss)(IcosahedronPatchY);
      auto Ld = Umu_v(ss)(IcosahedronPatchDiagonal);

      auto SE = stencil_v.GetEntry(ent_Xp,ss);
      uint64_t xp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Yp,ss);
      uint64_t yp_idx = SE->_offset;

      SE = stencil_v.GetEntry(ent_Dp,ss);
      uint64_t dp_idx = SE->_offset;

      auto g  = g_v(ss)();
      auto gx = g_v(xp_idx)();
      auto gy = g_v(yp_idx)();
      auto gd = g_v(dp_idx)();

      auto lx = Umu_v(ss)(IcosahedronPatchX);
      auto ly = Umu_v(ss)(IcosahedronPatchY);
      auto ld = Umu_v(ss)(IcosahedronPatchDiagonal);

      lx = g*lx*adj(gx);
      ly = g*ly*adj(gy);
      ld = g*ld*adj(gd);

      coalescedWrite(Umu_v[ss](IcosahedronPatchX),lx);
      coalescedWrite(Umu_v[ss](IcosahedronPatchY),ly);
      coalescedWrite(Umu_v[ss](IcosahedronPatchDiagonal),ld);
      });
  }

};


int main (int argc, char ** argv)
{
  // ====================================================
  std::cout << "@@@@@@ SET ICOSAHEDRAL GRID @@@@@@" << std::endl;
  // ====================================================

  Grid_init(&argc,&argv);

  // const int L=8;
  // const int Npatch = IcosahedralPatches;

  // Put SIMD all in time direction
  Coordinate latt_size = GridDefaultLatt();
  Coordinate simd_layout({1,1,vComplex::Nsimd(),1});
  Coordinate mpi_layout = GridDefaultMpi();

  std::cout << GridLogMessage << " mpi "<<mpi_layout<<std::endl;
  std::cout << GridLogMessage << " simd "<<simd_layout<<std::endl;
  std::cout << GridLogMessage << " latt "<<latt_size<<std::endl;
  GridCartesianCrossIcosahedron EdgeGrid(latt_size,simd_layout,mpi_layout,IcosahedralEdges);
  std::cout << GridLogMessage << " Created edge grid "<<std::endl;
  GridCartesianCrossIcosahedron VertexGrid(latt_size,simd_layout,mpi_layout,IcosahedralVertices);
  std::cout << GridLogMessage << " Created vertex grid "<<std::endl;

  // ====================================================
  std::cout << "@@@@@@ SET FIELDS AND ZERO INIT @@@@@@" << std::endl;
  // ====================================================

  LatticeLorentzColourMatrix Umu(&EdgeGrid);
  LatticeLorentzColourMatrix Umuck(&EdgeGrid);
  // LatticeIcosahedralLorentzColourMatrix Umu(&EdgeGrid);
  // LatticeIcosahedralLorentzColourMatrix Umuck(&EdgeGrid);
  LatticeComplex Phi(&VertexGrid);
  std::cout << GridLogMessage << " Created two fields "<<std::endl;

  Phi = Zero();
  Umu = Zero();
  std::cout << GridLogMessage << " Zeroed two fields "<<std::endl;


  // ====================================================
  std::cout << "!!!!!!!!!!!! RANDOM GAUGE TRANSFORM !!!!!!!!!!!!!!" << std::endl;
  // ====================================================

  Complex one (1.0);
  Phi = one;
  Umu = one;

  std::cout << GridLogMessage << "Creating face stencil"<<std::endl;
  IcosahedralSupport<IcosahedralGimpl> Support(&VertexGrid,&EdgeGrid);

  std::cout << GridLogMessage << " Calling Test Geometry "<<std::endl;
  Support.FaceStencil.TestGeometry();


  // Random gauge xform
  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG  vRNG(&EdgeGrid);   vRNG.SeedFixedIntegers(seeds);
  
  LatticeColourMatrix g(&VertexGrid);
  LatticeReal gr(&VertexGrid);
  LatticeComplex gc(&VertexGrid);
  //  gr = Zero();
  gaussian(vRNG,gr);
  Complex ci(0.0,1.0);
  gc = toComplex(gr);
  g=one;
  g = g * exp(ci*gc);
 
  std::cout << GridLogMessage << "****************************************"<<std::endl;
  std::cout << GridLogMessage << " Check plaquette is gauge invariant "<<std::endl;
  std::cout << GridLogMessage << "****************************************"<<std::endl;
  std::cout << GridLogMessage << " applying gauge transform"<<std::endl;
  Support.GaugeTransform    (g,Umu);
  std::cout << GridLogMessage << " applied gauge transform "<<std::endl;

  LatticeComplex plaq1(&EdgeGrid);
  LatticeComplex plaq2(&EdgeGrid);

  std::cout << GridLogMessage << " recalculating plaquette "<<std::endl;
  Support.ForwardTriangles(Umu,plaq1,plaq2);
  std::cout << GridLogMessage << " plaq1 "<< norm2(plaq1)<<std::endl;
  std::cout << GridLogMessage << " plaq2 "<< norm2(plaq2)<<std::endl;




  // ====================================================
  std::cout << "!!!!!!!!!!!! TYPE CONVERSION IMPLEMENTATION !!!!!!!!!!!!!!" << std::endl;
  // ====================================================


  // ************************************************************** //
  // Conversion from Grid::LatticeLorentzColourMatrix to Gauge U_NM //
  // ************************************************************** //

  auto* grid = &EdgeGrid;
  int nd = grid->Nd();
  int L  = grid->LocalDimensions()[0];
  assert (grid->LocalDimensions()[0] == grid->LocalDimensions()[1]);
  int T = grid->GlobalDimensions()[2];

  assert( Comp::Nt==T );
  assert( Comp::N_REFINE==L );

  SphereUtils util(L);
  IcosahedralStencil NNee(&EdgeGrid, &VertexGrid); // edge neighbours with edge domain

  using Base=S2Simp;
  Base base(Comp::N_REFINE);
  std::cout << "# lattice set. " << std::endl;
  using Gauge=GaugeExt<Base, Comp::Nt, Comp::is_compact>;
  Gauge U_NM(base);

  using Action=U1WilsonExt<Base>;

  for(int x=0; x<L; x++){
    for(int y=0; y<L; y++){
       for(int t=0; t<T; t++){
         for(int Patch=0; Patch<IcosahedralPatches; Patch++){
           // const int HemiPatch = Patch%HemiPatches;
           // const int north = Patch/HemiPatches;
           // const int south = 1-north;

           Coordinate n({x,y,t,Patch});
           const Idx in = util.GridCoord2Idx.at( GridCoord{x,y,Patch} );
           const auto Un = peekSite(Umu, n);

           int isPoleY;
           int isPoleX;

           // ~~~~~~~~~~~~~
           // loop over m = n + mu (mu=X,Y,D)
           // ~~~~~~~~~~~~~
           Coordinate m;
           Idx im;
           ComplexD cplx;

           // X
           auto Unmu = peekIndex<LorentzIndex>(Un, IcosahedronPatchX);
           NNee.GetNbrForPlusX       (grid,n,m,isPoleX);
           if(!isPoleX) im = util.GridCoord2Idx.at( GridCoord{m[IcosahedronPatchX],m[IcosahedronPatchY],m[IcosahedronTime]} );
           else im = util.IdxS;
           {
             cplx = peekIndex<ColourIndex>( Unmu, 0, 0 );
             const Link ell{in,im};
             const Idx il = base.map2il.at(ell);
             const int sign = base.map2sign.at(ell);
             U_NM.sp( t, il ) = sign*std::arg(cplx);
           }

           // Y
           Unmu = peekIndex<LorentzIndex>(Un, IcosahedronPatchY);
           NNee.GetNbrForPlusY       (grid,n,m,isPoleY);
           if(!isPoleY) im = util.GridCoord2Idx.at( GridCoord{m[IcosahedronPatchX],m[IcosahedronPatchY],m[IcosahedronTime]} );
           else im = util.IdxN;
           {
             cplx = peekIndex<ColourIndex>( Unmu, 0, 0 );
             const Link ell{in,im};
             const Idx il = base.map2il.at(ell);
             const int sign = base.map2sign.at(ell);
             U_NM.sp( t, il ) = sign*std::arg(cplx);
           }

           // D
           Unmu = peekIndex<LorentzIndex>(Un, IcosahedronPatchDiagonal);
           NNee.GetNbrForPlusDiagonal(grid,n,m);
           im = util.GridCoord2Idx.at( GridCoord{m[IcosahedronPatchX],m[IcosahedronPatchY],m[IcosahedronTime]} );
           cplx = peekIndex<ColourIndex>( Unmu, 0, 0 );

           const Link ell{in,im};
           const Idx il = base.map2il.at(ell);
           const int sign = base.map2sign.at(ell);
           U_NM.sp( t, il ) = sign*std::arg(cplx);


           // ~~~~~~~~~~~~~
           // mu = t
           // ~~~~~~~~~~~~~
           Unmu = peekIndex<LorentzIndex>(Un, IcosahedronTime);
           cplx = peekIndex<ColourIndex>( Unmu, 0, 0 );
           U_NM.tp( t, in ) = std::arg(cplx);
         }}}}


  // ************************************************************** //
  // CHECK //
  // ************************************************************** //

  const double gsq = 1.0;
  const double at = 16/Comp::Nt;
  Action SW( gsq, at, base );
  std::cout << "SW = " << SW(U_NM) << std::endl;


  std::vector<double> tmp(U_NM.Nt, 0.0);

  // spatial
  for(int s=0; s<U_NM.Nt; s++){
    for(Idx i=0; i<base.faces.size(); i++) {
      const Face& face = base.faces[i];

      double angle = 0.0;
      if constexpr(U_NM.is_compact) angle += - SW.beta_s[i] * ( std::cos( U_NM.plaquette_angle(s, face) ) - 1.0);
      else angle += 0.5*SW.beta_s[i] * std::pow( U_NM.plaquette_angle(s, face), 2 );
      std::cout << "s = " << s << "; i = " << i << "; angle = " << angle << std::endl;
    }
  }

  // temporal
  for(int s=0; s<U_NM.Nt ; s++){
    for(const Link& link : base.links) {
      const Idx il = base.map2il.at(link);

      double angle = 0.0;
      if constexpr(U_NM.is_compact) angle += - SW.beta_t[il]*(std::cos( U_NM.plaquette_angle(s, link) ) - 1.0);
      else angle += 0.5*SW.beta_t[il] * std::pow( U_NM.plaquette_angle(s, link), 2 );
      std::cout << "s = " << s << "; il = " << il << "; angle = " << angle << std::endl;
    }
  }




  // ************************************************************** //
  // Spin connection
  // ************************************************************** //

  using WilsonDirac=DiracExt<Base, DiracS2Simp>;

  // LatticeSpinMatrixD MatPs(&EdgeGrid);
  // MatPs = Zero();
  // std::cout << "MatPs = " << std::endl
  //           << MatPs << std::endl;

  const double m = 0.0;
  const double r = 1.0;
  const double M5 = 0.0;
  WilsonDirac DW(base, m, r, M5, at);
  const auto& bd = DW.bd;

  {
    Eigen::MatrixXcd id = Eigen::MatrixXcd::Identity(4,4);
    std::array<Eigen::MatrixXcd,4> gammas;
    {
      for(auto& gamma : gammas) gamma = Eigen::MatrixXcd::Zero(4,4);
      const Complex             I(0., 1.), mI(0., -1.);
      gammas[0](0, 3) = I;
      gammas[0](1, 2) = I;
      gammas[0](2, 1) = mI;
      gammas[0](3, 0) = mI;
      //
      gammas[1](0, 3) = -1.;
      gammas[1](1, 2) = 1.;
      gammas[1](2, 1) = 1.;
      gammas[1](3, 0) = -1.;
      //
      gammas[2](0, 2) = I;
      gammas[2](1, 3) = mI;
      gammas[2](2, 0) = mI;
      gammas[2](3, 1) = I;
      //
      gammas[3](0, 2) = 1.;
      gammas[3](1, 3) = 1.;
      gammas[3](2, 0) = 1.;
      gammas[3](3, 1) = 1.;
    }


    Eigen::MatrixXcd eDotGamma = [&](const Idx ix, const Idx iy) { // located at x
      const double al = DW.alpha.at(Link{ix,iy});
      return std::cos(al)*gammas[1] + std::sin(al)*gammas[2];
    };


    for(int s=0; s<U_NM.Nt; s++){
      for(Idx ix=0; ix<base.n_sites; ix++){
        Idx counter = 4*base.counter_accum.back()*s + 4*base.counter_accum[ix];
        for(const Idx iy : base.nns[ix]){
          const Idx il = base.map2il.at(Link{ix,iy});
          Eigen::MatrixXcd MatP = 0.5*( -r * id + eDotGamma(ix, iy) );
          const MS MatOmega = bd.Omega(ix, iy);
        }
      }
    }
  }


  Grid_finalize();
}
































  // template<class GaugeNM>
  // void Grid2NMGauge( GaugeNM& U_NM, const Grid::LatticeLorentzColourMatrix& Umu,
  //                    const int L, const int T,
  //                    Grid::IcosahedralStencil NNee ) const {
  //   // using namespace Grid;

  //   const auto& base = U_NM.base;
  //   auto* grid = NNee._grid;

  //   for(int x=0; x<L; x++){
  //     for(int y=0; y<L; y++){
  //       for(int t=0; t<T; t++){
  //         for(int Patch=0; Patch<Grid::IcosahedralPatches; Patch++){
  //           // const int HemiPatch = Patch%HemiPatches;
  //           // const int north = Patch/HemiPatches;
  //           // const int south = 1-north;

  //           Grid::Coordinate n({x,y,t,Patch});
  //           const Idx in = this->GridCoord2Idx.at( GridCoord{x,y,Patch} );
  //           const auto Un = Grid::peekSite(Umu, n);

  //           int isPoleY;
  //           int isPoleX;

  //           // ~~~~~~~~~~~~~
  //           // loop over m = n + mu (mu=X,Y,D)
  //           // ~~~~~~~~~~~~~
  //           Grid::Coordinate m;
  //           Idx im;
  //           Grid::ComplexD cplx;

  //           // X
  //           auto Unmu = Grid::peekIndex<Grid::LorentzIndex,decltype(Un)>(Un, Grid::IcosahedronPatchX);
  //           NNee.GetNbrForPlusX       (grid,n,m,isPoleX);
  //           if(!isPoleX) im = this->GridCoord2Idx.at( GridCoord{m[Grid::IcosahedronPatchX],m[Grid::IcosahedronPatchY],m[Grid::IcosahedronTime]} );
  //           else im = this->IdxS;
  //           {
  //             cplx = Grid::peekIndex<Grid::ColourIndex>( Unmu, 0, 0 );
  //             const Link ell{in,im};
  //             const Idx il = base.map2il.at(ell);
  //             const int sign = base.map2sign.at(ell);
  //             U_NM.sp( t, il ) = sign*std::arg(cplx);
  //           }

  //           // // Y
  //           // Unmu = Grid::peekIndex<Grid::LorentzIndex>(Un, Grid::IcosahedronPatchY);
  //           // NNee.GetNbrForPlusY       (grid,n,m,isPoleY);
  //           // if(!isPoleY) im = this->GridCoord2Idx.at( GridCoord{m[Grid::IcosahedronPatchX],m[Grid::IcosahedronPatchY],m[Grid::IcosahedronTime]} );
  //           // else im = this->IdxN;
  //           // {
  //           //   cplx = Grid::peekIndex<Grid::ColourIndex>( Unmu, 0, 0 );
  //           //   const Link ell{in,im};
  //           //   const Idx il = base.map2il.at(ell);
  //           //   const int sign = base.map2sign.at(ell);
  //           //   U_NM.sp( t, il ) = sign*std::arg(cplx);
  //           // }

  //           // // D
  //           // Unmu = Grid::peekIndex<Grid::LorentzIndex>(Un, Grid::IcosahedronPatchDiagonal);
  //           // NNee.GetNbrForPlusDiagonal(grid,n,m);
  //           // im = this->GridCoord2Idx.at( GridCoord{m[Grid::IcosahedronPatchX],m[Grid::IcosahedronPatchY],m[Grid::IcosahedronTime]} );
  //           // cplx = Grid::peekIndex<Grid::ColourIndex>( Unmu, 0, 0 );

  //           // const Link ell{in,im};
  //           // const Idx il = base.map2il.at(ell);
  //           // const int sign = base.map2sign.at(ell);
  //           // U_NM.sp( t, il ) = sign*std::arg(cplx);


  //           // // ~~~~~~~~~~~~~
  //           // // mu = t
  //           // // ~~~~~~~~~~~~~
  //           // Unmu = Grid::peekIndex<Grid::LorentzIndex>(Un, Grid::IcosahedronTime);
  //           // cplx = Grid::peekIndex<Grid::ColourIndex>( Unmu, 0, 0 );
  //           // U_NM.tp( t, in ) = std::arg(cplx);
  //         }}}}

  // }
