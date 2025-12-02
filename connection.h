#pragma once

#include <gsl/gsl_integration.h>
#include<utility>



using Func=std::function<double(const double)>;
double unwrap(double x, void *p) {
  auto fp = static_cast<Func*>(p);
  return (*fp)(x);
}


struct DualSpinStructure {
  using DirectedLink=std::pair<Idx,int>;
  using Face=std::vector<DirectedLink>;

  RefinedIcosahedronDual dual; // copy

  V3 ptilde(const V3& r1, const V3& r2, const double stilde) const {
    assert(0.<=stilde && stilde<=1.);
    assert( std::abs(r1.norm()-1.0)<1.0e-14 );
    assert( std::abs(r2.norm()-1.0)<1.0e-14 );
    return r1 + stilde*(r2-r1);
  }

  V3 p(const V3& r1, const V3& r2, const double stilde) const {
    const V3 tmp = ptilde(r1,r2,stilde);
    return tmp/tmp.norm();
  }

  V3 etilde(const V3& r1, const V3& r2, const double stilde) const {
    const V3 dp = r2-r1;
    const V3 p_ = p(r1,r2,stilde);
    const V3 ptilde_ = ptilde(r1,r2,stilde);
    const double tmp = p_.dot(dp);
    return (dp - tmp*p_)/ptilde_.norm();
  }

  V3 dtheta_dp(const V3& r) const {
    const double sqrt_xsq_ysq = std::sqrt( r(0)*r(0)+r(1)*r(1) );
    const double normsq = r.squaredNorm();

    const double dtheta_dx = r(0)*r(2)/ normsq / sqrt_xsq_ysq;
    const double dtheta_dy = r(1)*r(2)/ normsq / sqrt_xsq_ysq;
    const double dtheta_dz = -sqrt_xsq_ysq / normsq;
    V3 res;
    res << dtheta_dx, dtheta_dy, dtheta_dz;
    return res;
  }

  V3 dphi_dp(const V3& r) const {
    const double dphi_dx = -r(1)/(r(0)*r(0)+r(1)*r(1));
    const double dphi_dy = r(0)/(r(0)*r(0)+r(1)*r(1));
    V3 res;
    res << dphi_dx, dphi_dy, 0.0;
    return res;
  }

  // V3 dp_dstilde(const V3& r1, const V3& r2, const double stilde) const {
  //   const V3 dp = r2-r1;
  //   const V3 ptilde = ptilde(r1, r2, stilde);
  //   const double norm = ptilde.norm();

  //   const double innner = ptilde.dot(dp);
  // }

  V3 e(const V3& r1, const V3& r2, const double stilde) const {
    const V3 tmp = etilde(r1,r2,stilde);
    return tmp/tmp.norm();
  }

  V2 e2(const V3& r1, const V3& r2, const double stilde) const {
    const double xx = dtheta_dstilde(r1,r2,stilde);
    const double yy = std::sin(theta(p(r1,r2,stilde))) * dphi_dstilde(r1,r2,stilde);
    V2 res;
    res << xx, yy;
    return res;
  }

  inline double theta(const V3& r) const {
    assert( std::abs(r.norm()-1.0)<1.0e-14 );
    const double res = std::acos(r[2]);
    assert( 0<=res && res<=M_PI);
    return res;
  }

  inline double phi(const V3& r) const {
    assert( std::abs(r.norm()-1.0)<1.0e-14 );
    const double res = std::atan2(r[1], r[0]);
    assert( -M_PI<=res && res<=M_PI);
    return res;
  }


  double dtheta_dstilde(const V3& r1, const V3& r2, const double stilde) const {
    const V3 p_ = p(r1,r2,stilde);
    const V3 dp_dstilde_ = etilde(r1,r2,stilde);
    const V3 dtheta_dp_ = dtheta_dp( p_ );
    return dp_dstilde_.dot(dtheta_dp_);
  }

  double dphi_dstilde(const V3& r1, const V3& r2, const double stilde) const {
    const V3 p_ = p(r1,r2,stilde);
    const V3 dp_dstilde_ = etilde(r1,r2,stilde);
    const V3 dphi_dp_ = dphi_dp( p_ );
    return dp_dstilde_.dot(dphi_dp_);
  }

  inline double alpha(const V3& r1, const V3& r2, const double stilde) const {
    // const V3 p_ = p(r1,r2,stilde);
    // const double xx = dtheta_dstilde(r1,r2,stilde);
    // const double yy = std::sin(theta(p_)) * dphi_dstilde(r1,r2,stilde);
    // assert( std::abs(xx*xx+yy*yy-1.0)<1.0e-14 );
    const V2 v = e2(r1,r2,stilde);
    return std::atan2( v[1], v[0] );
  }





  std::vector<double> alphas;
  std::vector<double> omegas;

  std::vector<Face> faces;

  const double phi0;

  DualSpinStructure(const RefinedIcosahedronDual& dual_,
                    const double phi0_=M_PI/48.0)
    : dual( dual_ )
    , phi0(phi0_)
  {
    // Eigen::Matrix<double, 3, 3> rot;
    // rot << std::cos(tht), 0.0, std::sin(tht),
    //   0.0, 1.0, 0.0,
    //   -std::sin(tht), 0.0, std::cos(tht);
    // for(auto& elem : dual.vertices) elem = rot*elem;
    // rot << std::cos(tht), std::sin(tht), 0.0,
    //   -std::sin(tht), std::cos(tht), 0.0,
    //   0.0, 0.0, 1.0;
    // for(auto& elem : dual.vertices) elem = rot*elem;

    FillAlphas();
    FillOmegas();
    FillFaces();
  }

  inline Idx idx( const Idx iff, const int df ) const {
    return 3*iff + df;
  }

  inline Idx idx( const DirectedLink& ell ) const {
    return idx(ell.first, ell.second);
  }

  DirectedLink idx2DirectedLink( const Idx il ) const {
    const int df = il%3;
    const Idx iff = il/3;
    return std::make_pair(iff, df);
  }

  DirectedLink idx2DirectedLinkReverse( const Idx il ) const {
    const int df = il%3;
    const Idx iff = il/3;
    const FaceCoords f1 = dual.idx2FaceCoords(iff);

    FaceCoords f2;
    dual.shift(f2, f1, df);

    return std::make_pair(dual.idx(f2), df);
  }

  inline Idx NDirectedLinks() const {
    return 3*dual.NVertices();
  }

  void FillAlphas(){
    alphas.clear();
    alphas.resize( NDirectedLinks() );

    // Idx counter=0;
    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      // if(f1[3]!=XZ) continue;
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );
        const V3 rf2 = dual.vertices[if2];

        double alpha1 = alpha( rf1, rf2, 0.0 );

        // double alpha2f = alpha( rf1, rf2, 1.0 );
        // double alpha2b = alpha( rf2, rf1, 0.0 );

        // std::cout << "forward = " << e2(rf1,rf2,1.0).transpose() << std::endl;
        // std::cout << "backward = " << e2(rf2,rf1,0.0).transpose() << std::endl;
        // double diff = alpha2f - alpha2b - M_PI;

        // if(diff>M_PI) diff -= 2.0*M_PI;
        // else if(diff<=-M_PI) diff += 2.0*M_PI;
        // std::cout << "debug. diff = " << diff << std::endl;
        // assert( std::abs(diff)<1.0e-14 );

        // std::cout << "idx = " << idx(if1, df) << std::endl;
        // std::cout << "idx = " << idx(if2, df) << std::endl;
        alphas[idx(if1, df)] = alpha1;
        // alphas[idx(if2, df)] = alpha2b;

        // double omegaf = alpha2f - alpha1;
        // while(omegaf>M_PI) omegaf -= 2.0*M_PI;
        // while(omegaf<=-M_PI) omegaf += 2.0*M_PI;

        // const double phi1 = phi(rf1);
        // const double phi2 = phi(rf2);

        // double dphi1 = phi1-phi0;
        // double dphi2 = phi2-phi0;
        // double dphi = phi1-phi2;

        // if(dphi1 >= M_PI) dphi1 -= 2.0*M_PI;
        // else if(dphi1 < -M_PI) dphi1 += 2.0*M_PI;
        // if(dphi2 >= M_PI) dphi2 -= 2.0*M_PI;
        // else if(dphi2 < -M_PI) dphi2 += 2.0*M_PI;
        // if(dphi >= M_PI) dphi -= 2.0*M_PI;
        // else if(dphi < -M_PI) dphi += 2.0*M_PI;
        // if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < 1.0e-14 ) omegaf -= 2.0*M_PI;

        // omegas[idx(if1, df)] = omegaf;
        // omegas[idx(if2, df)] = -omegaf;
      }
    }
  }



  const int limit = 4000; // 1000;
  // Double epsabs = 1.0e-15; // 0.;
  double epsabs = 1.0e-13; // 0.;
  double epsrel = 1.0e-13; // TOLLOOSE;
  int key = 5;



  void FillOmegas(){
    // omegas.clear();
    // Idx counter=0;
    omegas.clear();
    omegas.resize( NDirectedLinks() );

    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      // if(f1[3]!=XZ) continue;
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );

        // debug.
        // {
        //   Idx il = idx(if1, df);
        //   DirectedLink ell = idx2DirectedLink( il );
        //   assert( ell.first==if1 );
        //   assert( ell.second==df );
        //   assert(il == idx(ell));
        // }
        // {
        //   Idx il = idx(if2, df);
        //   DirectedLink ell = idx2DirectedLink( il );
        //   assert( ell.first==if2 );
        //   assert( ell.second==df );
        //   assert(il == idx(ell));
        // }

        const V3 rf2 = dual.vertices[if2];

        double omega12;
        {
          gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
          double result, error;

          Func f1 = [&](const double stilde){
            const V3 p_ = p(rf1, rf2, stilde);
            const double dphi_dstilde_ = dphi_dstilde(rf1, rf2, stilde);
            // std::cout << "debug. dphi_dstilde = " << dphi_dstilde << std::endl;
            const double theta_ = theta( p_ );

            return -dphi_dstilde_ * std::cos( theta_ );
          };
          gsl_function F;
          F.function = &unwrap;
          F.params = &f1;
          gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit,
                              key, w, &result, &error);

          omega12 = result;
          gsl_integration_workspace_free (w);
        }

        // check
        {
          const double alpha_1 = alpha( rf1,rf2,0.0);
          const double alpha_2 = alpha( rf1,rf2,1.0);

          double diff = alpha_2 - omega12 - alpha_1;
          // std::cout << "debug. diff = " << diff << std::endl;
          if(diff>M_PI) diff -= 2.0*M_PI;
          else if(diff<=-M_PI) diff += 2.0*M_PI;
          assert( std::abs(diff)<1.0e-14 );
        }

        // alpha checks
        // const double alpha1 = alphas[ idx(if1, df) ];
        // const double alpha2 = alphas[ idx(if2, df) ];
        // double alpha1c = alpha( rf1, rf2, 0.0 );

        // double alpha2f = alpha( rf1, rf2, 1.0 );
        // double alpha2b = alpha( rf2, rf1, 0.0 );
        // assert( std::abs(alpha1-alpha1c)<1.0e-14 );
        // assert( std::abs(alpha2-alpha2b)<1.0e-14 );
        // assert( std::abs(alpha_2-alpha2f)<1.0e-14 );

        // {
        //   double diff = alpha2f - alpha2b - M_PI;

        //   if(diff>M_PI) diff -= 2.0*M_PI;
        //   else if(diff<=-M_PI) diff += 2.0*M_PI;
        //   std::cout << "debug2. diff = " << diff << std::endl;
        //   assert( std::abs(diff)<1.0e-14 );
        // }

        // double diff3 = M_PI + alpha2 - omega12 - alpha1;
        // std::cout << "debug. diff3 = " << diff3 << std::endl;
        // if(diff3>M_PI) diff3 -= 2.0*M_PI;
        // else if(diff3<=-M_PI) diff3 += 2.0*M_PI;
        // assert( std::abs(diff3)<1.0e-14 );

        // assert( counter==idx(if1,df) ); // debug
        // counter++;
        // assert( counter==omegas.size() ); // debug

        const double phi1 = phi(rf1);
        const double phi2 = phi(rf2);

        double dphi1 = phi1-phi0;
        double dphi2 = phi2-phi0;
        double dphi = phi1-phi2;

        if(dphi1 >= M_PI) dphi1 -= 2.0*M_PI;
        else if(dphi1 < -M_PI) dphi1 += 2.0*M_PI;
        if(dphi2 >= M_PI) dphi2 -= 2.0*M_PI;
        else if(dphi2 < -M_PI) dphi2 += 2.0*M_PI;
        if(dphi >= M_PI) dphi -= 2.0*M_PI;
        else if(dphi < -M_PI) dphi += 2.0*M_PI;
        if( std::abs(dphi1) + std::abs(dphi2) - std::abs(dphi) < 1.0e-14 ) omega12 -= 2.0*M_PI;

        omegas[idx(if1, df)] = omega12;
        // omegas[idx(if2, df)] = -omega12;
      }
    }
  }

  void FillFaces(){
    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1);
      if(f1[3]!=XZ) continue;

      Face face;
      FaceCoords f2 = f1;

      bool is_south = (f1[2]<5);
      // std::cout << "debug. ### f1 = " << std::endl;
      // for(auto elem : f1) std::cout << elem << " ";
      // std::cout << std::endl;

      int df = dual.nC;
      face.push_back(std::make_pair(if1,df));
      // std::cout << "debug. df = " << df << std::endl;

      // int x_tmp = f2[0];
      // int y_tmp = f2[1];
      int type_tmp = f2[3];
      int is_shift = dual.shift( f2, f1, df );
      Idx if2 = dual.idx(f2);

      // std::cout << "debug. f2 = " << std::endl;
      // for(auto elem : f2) std::cout << elem << " ";
      // std::cout << std::endl;

      int counter;
      for(counter=0; counter<6; counter++){
        // df to next point
        if(is_south){
          if( f2[3] != type_tmp ) df = (df+1)%3; // C->A, A->B, B->C
          else if( is_shift==1 ) df = df;
          else df = (df+2)%3;
        }
        else{
          if( f2[3] != type_tmp ) {
            df = (df+1)%3; // C->A, A->B, B->C
            // std::cout << "debug pt1" << std::endl;
          }
          else if( is_shift==1 ) {
            df = (df+2)%3;
            // std::cout << "debug pt2" << std::endl;
          }
          else {
            df = df;
            // std::cout << "debug pt3" << std::endl;
          }
        }
        face.push_back(std::make_pair(if2,df));
        // if( f2[0]==x_tmp && f2[1]==y_tmp ) df = (df+2)%3; // C->A, A->B, B->C
        // else df = (df+1)%3;

        // std::cout << "debug. df = " << df << std::endl;
        // x_tmp = f2[0];
        // y_tmp = f2[1];
        type_tmp = f2[3];
        FaceCoords f3;
        is_shift = dual.shift( f3, f2, df );
        f2 = f3;
        // std::cout << "debug. f2 = " << std::endl;
        // for(auto elem : f2) std::cout << elem << " ";
        // std::cout << std::endl;
        // std::cout << "debug. is_shift = " << is_shift << std::endl;
        if2 = dual.idx(f2);

        if(if2==if1) break;
      } // end for
      if(counter==6){
        // std::cout << "debug. face = " << std::endl;
        // for(auto elem : face) std::cout << elem << " ";
        // std::cout << std::endl;
        assert(false);
      }
      faces.push_back(face);
    }

    // south pole
    {
      Face face;
      for(int s=0; s<NPatches/2; s++){
        const FaceCoords f1{dual.L-1, 0, s, XZ};
        const Idx if1 = dual.idx(f1);
        const int df = dual.nA;
        face.push_back(std::make_pair(if1,df));
      }
      faces.push_back(face);
    }
    // north pole
    {
      Face face;
      for(int s=NPatches/2; s<NPatches; s++){
        const FaceCoords f1{0, dual.L-1, s, ZY};
        const Idx if1 = dual.idx(f1);
        const int df = dual.nA;
        face.push_back(std::make_pair(if1,df));
      }
      faces.push_back(face);
    }

  }

};
