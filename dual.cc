#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <limits>
#include <stdfloat>


#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>



// using Double = double;
using Double = std::float64_t;
using Idx = std::int32_t;
// using Idx=std::size_t;
// using V2=Eigen::Matrix<Double, 1, 2>;
// constexpr int EDIM = 3;
// using V3=Eigen::Matrix<Double, 1, 3>;
using V3=Eigen::Matrix<Double, 3, 1>;
// using V3 = Eigen::Vector3d;
using M3=Eigen::Matrix<Double, 3, 3>;
// using M3 = Eigen::Matrix3d;//
// using V3=Eigen::Matrix<Double, 1, 3>;
using Complex = std::complex<double>;




#include "geodesic.h"
#include "integral.h"


#include "icosahedron.h"
#include "lattices.h"
#include "groups.h"

// using namespace Geodesic

#include "connection.h"


// using Link = std::array<Idx,2>; // <Idx,Idx>;
// using Face = std::vector<Idx>;
// using VD=Eigen::Matrix<Double, 2, 1>;
// using VE=Eigen::Matrix<Double, 3, 1>;


Double TOL=1.0e-14;

constexpr Double TOLLOOSE=1.0e-6;

const int limit = 4000; // 1000;
// Double epsabs = 1.0e-15; // 0.;
Double epsabs = 1.0e-13; // 0.;
Double epsrel = 1.0e-13; // TOLLOOSE;
int key = 5;


std::string dir = "./data/";


int main(int argc, char* argv[]){
  int L=128;
  if(argc==2) L = atoi(argv[1]);


  RefinedIcosahedron lattice(L);
  FullIcosahedralGroup Ih( "multtablemathematica.dat",
                           3, 19, 60 );
  Rotation rot;
  RefinedIcosahedronDual dual(lattice);

  // Double tht = 1.0/L;
  // V3 tht;
  // tht << -0.6/L, 0.0, 0.0;
  // Eigen::Matrix<Double, 3, 3> rot_;
  // rot << std::cos(tht), 0.0, std::sin(tht),
  //   0.0, 1.0, 0.0,
  //   -std::sin(tht), 0.0, std::cos(tht);
  // V3 deltan << ;

  // std::vector<VE> simp_sites;
  // std::vector<VE> dual_sites;
  // std::vector<Link> links;
  // std::vector<std::vector<Idx>> nns;
  // std::vector<Face> faces;

  std::vector<double> omegas;
  std::vector<Geodesic::Sol> sols;

  for(Idx if1=0; if1<dual.NVertices(); if1++){
    const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
    // assert(f1[2]==0);
    if(f1[3]!=XZ) continue;
    const V3 rf1 = dual.vertices[if1];

    for(int df=0; df<3; df++){
      FaceCoords f2;
      dual.shift( f2, f1, df );
      const Idx if2 = dual.idx( f2 );
      const V3 rf2 = dual.vertices[if2];

      assert( std::abs( rf1.norm() -1.0)<1.0e-14 );
      assert( std::abs( rf2.norm() -1.0)<1.0e-14 );

      using namespace Geodesic;

      Pt x1( rf1 );
      Pt x2( rf2 );
      // std::cout << "debug. x1 = " << x1.x << std::endl;
      // std::cout << "debug. x1 = " << x1.xi << std::endl;
      // std::cout << "debug. x2 = " << x2.x << std::endl;
      // std::cout << "debug. x2 = " << x2.xi << std::endl;

      Sol sol = SolveGeodesics( x1, x2 );
          gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
          double result, error;

          if(!sol.is_split){
            F f1 = [&](Double s){ return sol.Dphi(s) * std::cos( sol.theta(s) ); };
            gsl_function F;
            F.function = &unwrap;
            F.params = &f1;
            gsl_integration_qag(&F, 0., sol.ell, epsabs, epsrel, limit,
                                key, w, &result, &error);
          }
          else{
            double result1, result2;
            double error1, error2;
            {
              F f1 = [&](Double s){ return sol.Dphi(s, 0) * std::cos( sol.theta(s, 0) ); };
              gsl_function F;
              F.function = &unwrap;
              F.params = &f1;
              gsl_integration_qag(&F, 0.0, sol.sE, epsabs, epsrel, limit,
                                  key, w, &result1, &error1);
            }
            {
              F f1 = [&](Double s){ return sol.Dphi(s, 1) * std::cos( sol.theta(s, 1) ); };
              gsl_function F;
              F.function = &unwrap;
              F.params = &f1;
              gsl_integration_qag(&F, sol.sE, sol.ell, epsabs, epsrel, limit,
                                  key, w, &result2, &error2);
            }
            result = result1+result2;
            error=error1+error2;
          }
          omegas.push_back( result );
          sols.push_back( sol );
    }}

  DualSpinStructure spin(dual);

  {
    const double phi0 = M_PI/48.0;
    std::vector<DualSpinStructure::DirectedLink> section;

    // for(Idx il=0; il<links.size(); il++){
    //   Link link = links[il];
    //   const double phi1 = Pt(sites[link[0]]).xi[1];
    //   const double phi2 = Pt(sites[link[1]]).xi[1];
    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      // assert(f1[2]==0);
      if(f1[3]!=XZ) continue;
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df );
        const Idx if2 = dual.idx( f2 );
        const V3 rf2 = dual.vertices[if2];

        using namespace Geodesic;

        Pt x1( rf1 );
        Pt x2( rf2 );
        const double phi1 = x1.xi[1];
        const double phi2 = x2.xi[1];

        double dphi1 = phi1-phi0;
        double dphi2 = phi2-phi0;
        double dphi = phi1-phi2;

        if(dphi1 >= M_PI) dphi1 -= 2.0*M_PI;
        else if(dphi1 < -M_PI) dphi1 += 2.0*M_PI;
        if(dphi2 >= M_PI) dphi2 -= 2.0*M_PI;
        else if(dphi2 < -M_PI) dphi2 += 2.0*M_PI;
        if(dphi >= M_PI) dphi -= 2.0*M_PI;
        else if(dphi < -M_PI) dphi += 2.0*M_PI;

        if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < TOLLOOSE ){
          section.push_back( std::make_pair(if1, df) );
        }
      }}
    for(auto& link : section) omegas[spin.idx(link)] -= 2.0*M_PI;
  }


  std::map<const DualSpinStructure::DirectedLink, double> omega;
  for(Idx il=0; il<spin.NDirectedLinks(); il++){
    // const DualSpinStructure::DirectedLink link = ;
    omega.insert( { spin.idx2DirectedLink(il), omegas[il] } );
    omega.insert( { spin.idx2DirectedLinkReverse(il), -omegas[il] } );
  }


  std::vector<double> alpha0s, alpha1s;
  for( auto& sol : sols){
    using namespace Geodesic;
    using VD = V2;

    double alpha0, alpha1;
    if(!sol.is_split){
      const double theta0 = sol.theta(0);
      const double dtheta_ds0 = sol.Dtheta(0);
      const double dphi_ds0 = sol.Dphi(0);
      assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );

      VD e0( dtheta_ds0, std::sin(theta0)*dphi_ds0 );
      alpha0 = _acos( e0[0] );
      if(e0[1]<0) alpha0 *= -1.0;
      assert( std::abs(std::cos(alpha0)-e0[0]) + std::abs(std::sin(alpha0)-e0[1]) < TOLLOOSE );

      const double theta1 = sol.theta(sol.ell);
      const double dtheta_ds1 = sol.Dtheta(sol.ell);
      const double dphi_ds1 = sol.Dphi(sol.ell);
      assert( std::abs(dtheta_ds1*dtheta_ds1 + std::sin(theta1)*std::sin(theta1)*dphi_ds1*dphi_ds1-1.0)<TOLLOOSE );

      VD e1( dtheta_ds1, std::sin(theta1)*dphi_ds1 );
      alpha1 = _acos( e1[0] );
      if(e1[1]<0) alpha1 *= -1.0;
      assert( std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1]) < TOLLOOSE );
    }
    else{
      const double theta0 = sol.theta(0, 0);
      const double dtheta_ds0 = sol.Dtheta(0, 0);
      const double dphi_ds0 = sol.Dphi(0, 0);
      assert( std::abs(dtheta_ds0*dtheta_ds0 + std::sin(theta0)*std::sin(theta0)*dphi_ds0*dphi_ds0-1.0)<TOLLOOSE );

      VD e0( dtheta_ds0, std::sin(theta0)*dphi_ds0 );
      alpha0 = _acos( e0[0] );
      if(e0[1]<0) alpha0 *= -1.0;
      assert( std::abs(std::cos(alpha0)-e0[0]) + std::abs(std::sin(alpha0)-e0[1]) < TOLLOOSE );

      const double theta1 = sol.theta(sol.ell, 1);
      const double dtheta_ds1 = sol.Dtheta(sol.ell, 1);
      const double dphi_ds1 = sol.Dphi(sol.ell, 1);
      assert( std::abs(dtheta_ds1*dtheta_ds1 + std::sin(theta1)*std::sin(theta1)*dphi_ds1*dphi_ds1-1.0)<TOLLOOSE );

      VD e1( dtheta_ds1, std::sin(theta1)*dphi_ds1 );
      alpha1 = _acos( e1[0] );
      if(e1[1]<0) alpha1 *= -1.0;
      assert( std::abs(std::cos(alpha1)-e1[0]) + std::abs(std::sin(alpha1)-e1[1]) < TOLLOOSE );
    }
    alpha0s.push_back( alpha0 );
    alpha1s.push_back( Mod(alpha1 + M_PI) );
  }

  // std::map<const Link, double> alpha;
  // for(Idx il=0; il<links.size(); il++){
  //     Link link = links[il];
  //   alpha.insert( { link, alpha0s[il] } );
  //   alpha.insert( { Link{link[1],link[0]}, alpha1s[il] } );
  // }

  std::map<const DualSpinStructure::DirectedLink, double> alpha;
  for(Idx il=0; il<spin.NDirectedLinks(); il++){
    alpha.insert( { spin.idx2DirectedLink(il), alpha0s[il] } );
    alpha.insert( { spin.idx2DirectedLinkReverse(il), alpha1s[il] } );
  }

  // for(Idx ix=0; ix<sites.size(); ix++){
  // for(Idx if1=0; if1<dual.NVertices(); if1++){

  //   const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
  //   const V3 rf1 = dual.vertices[if1];

  //   // using namespace Geodesic;
  //   const auto e0x = Geodesic::Pt( rf1 ).e0();
  //   const auto e1x = Geodesic::Pt( rf1 ).e1();

  //   for(int df=0; df<3; df++){
  //     //    for(Idx iy : nns[ix]){
  //     FaceCoords f2;
  //     dual.shift( f2, f1, df );
  //     const Idx if2 = dual.idx( f2 );
  //     const V3 rf2 = dual.vertices[if2];

  //     const V3 diff = rf2-rf1;
  //     const double c0 = e0x.dot(diff);
  //     const double c1 = e1x.dot(diff);
  //     double al = std::atan(c1/c0);
  //     const int br = Geodesic::decide_branch( al-alpha.at( std::make_pair(if1, df) ) );
  //     al -= M_PI*br;
  //     // alpha[Link({ix, iy})] = al;
  //     alpha[ std::make_pair(if1, df) ] = al;
  //   }
  // }


  {
    std::clog << "# checking spin structure" << std::endl;
    // for(Idx ix=0; ix<sites.size(); ix++){
    //   for(Idx iy : nns[ix]){

    for(Idx il=0; il<spin.NDirectedLinks(); il++){
      // for(Idx if1=0; if1<dual.NVertices(); if1++){
      //   for(int df=0; df<3; df++){
      // omega.insert( { , omegas[il] } );
      // omega.insert( { spin.idx2DirectedLinkReverse(il), -omegas[il] } );
      DualSpinStructure::DirectedLink ell = spin.idx2DirectedLink(il);
      DualSpinStructure::DirectedLink ellR = spin.idx2DirectedLinkReverse(il);
      //    for(Idx iy : nns[ix]){
      // FaceCoords f2;
      // dual.shift( f2, f1, df );
      // const Idx if2 = dual.idx( f2 );
      // const V3 rf2 = dual.vertices[if2];

      const double alpha1 = alpha.at( ell );
      double alpha2 = alpha.at( ellR );
      double omega12 = omega.at( ell );

      double diff = (alpha2 + M_PI + omega12) - alpha1;
      std::cout << "diff: " << diff << std::endl;
      assert( Geodesic::isModdable(diff, 2.0*M_PI, TOL) );

      double om = alpha1 - (alpha2 + M_PI);
      std::cout << "diff = " << om-omega12 << std::endl;
      // const int br = decide_branch( om-omega12 );
      // om -= M_PI*br;
      // omega[Link({ix, iy})] = om;
    }
  }

  // {
  //   std::clog << "# checking deficits" << std::endl;
  //   int counter=0;
  //   for(auto& face : faces){
  //     int sign = 1;
  //     {
  //       VE x0 = simp_sites[counter];
  //       VE x1 = sites[face[0]];
  //       VE x2 = sites[face[1]];
  //       if((x1-x0).cross(x2-x1).dot(x0) < 0) sign = -1;
  //     }

  //     double sum = 0.0;
  //     for(int i=0; i<face.size(); i++){
  //       const Idx ix = face[i];
  //       const Idx j = (i+1)%face.size();
  //       const Idx jx = face[j];
  //       sum += omega.at( Link{ix,jx} );
  //     }
  //     sum *= sign;
  //     double mod = Mod(sum, 4.0*M_PI);
  //     if(mod>2.0*M_PI) mod -= 4.0*M_PI;
  //     assert( -1.5 * 4.0*M_PI/face.size() < mod && mod < 0.0 );
  //     counter++;
  //   }
  // }



  // {
  //   std::ofstream ofs(dir+"omega_dual_n"+std::to_string(n_refine)+".dat");
  //   ofs << std::scientific << std::setprecision(15);
  //   for(Idx il=0; il<links.size(); il++){
  //     Link link = links[il];
  //     ofs << std::setw(25) << link[0] << " ";
  //     ofs << std::setw(25) << link[1] << " ";
  //     ofs << std::setw(25) << omegas[il] << std::endl;
  //   }
  // }

  // {
  //   std::ofstream ofs(dir+"alpha_dual_n"+std::to_string(n_refine)+".dat");
  //   ofs << std::scientific << std::setprecision(15);
  //   for(Idx il=0; il<links.size(); il++){
  //     Link link = links[il];
  //     ofs << std::setw(25) << link[0] << " ";
  //     ofs << std::setw(25) << link[1] << " ";
  //     ofs << std::setw(25) << alpha0s[il] << std::endl;

  //     ofs << std::setw(25) << link[1] << " ";
  //     ofs << std::setw(25) << link[0] << " ";
  //     ofs << std::setw(25) << alpha1s[il] << std::endl;
  //   }
  // }

  // // dualtriangleareas
  // std::vector<Face> simp_faces;
  // {
  //   std::ifstream file(dir+"face_n"+std::to_string(n_refine)+".dat");

  //   std::string str;
  //   while (std::getline(file, str)){
  //     std::istringstream iss(str);
  //     Idx v;
  //     Face face;
  //     while( iss >> v ) {
  //       face.push_back( v );
  //     }
  //     simp_faces.push_back( face );
  //   }
  // }
  // std::vector<double> triangleareas;
  // for(const Face& face : simp_faces){
  //   double area = sphericalarea( Pt(simp_sites[face[0]]),
  //                                Pt(simp_sites[face[1]]),
  //                                Pt(simp_sites[face[2]]));
  //   triangleareas.push_back( area );
  // }
  // assert( simp_faces.size()==triangleareas.size() );
  // std::cout << sites.size() << std::endl;
  // std::cout << simp_faces.size() << std::endl;
  // std::cout << triangleareas.size() << std::endl;
  // assert( sites.size()==triangleareas.size() );

  // {
  //   std::ofstream ofs(dir+"dualtriangleareas_n"+std::to_string(n_refine)+".dat");
  //   ofs << std::scientific << std::setprecision(15);
  //   for(double elem : triangleareas) ofs << std::setw(25) << elem << std::endl;
  // }



  return 0;
}

