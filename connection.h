#pragma once

#include <gsl/gsl_integration.h>



using Func=std::function<double(const double)>;
double unwrap(double x, void *p) {
  auto fp = static_cast<Func*>(p);
  return (*fp)(x);
}


struct DualSpinStructure {

  const RefinedIcosahedronDual& dual;

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
    const V2 v = e2(r1,r2,stilde);
    return std::atan2( v[1], v[0] );
  }





  std::vector<double> alphas;
  std::vector<double> omegas;

  const double TOL;
  const double phi0;

  DualSpinStructure(const RefinedIcosahedronDual& dual_,
                    const double TOL_=1.0e-12, // okay up to L~160
                    const double phi0_=M_PI/48.0
                    )
    : dual( dual_ )
    , TOL(TOL_)
    , phi0(phi0_)
  {
    for(const V3& r : dual.vertices){
      const double phi_ = phi(r);
      assert( std::norm(phi0-phi_)>1.0e-14 );
    }

    FillAlphas();
    FillOmegas();
    check_omega();
  }

  void check_omega() const {
    for(Idx i=0; i<dual.simplicial.NVertices(); i++){
      const DualFace face = dual.faces[i];
      double omega_sum = 0.0;
      for(auto elem : face) {
        omega_sum += omegas[ dual.linkidx( elem.first, elem.second ) ];
      }
      while(omega_sum>2.0*M_PI) omega_sum -= 4.0*M_PI;
      while(omega_sum<=-2.0*M_PI) omega_sum += 4.0*M_PI;

      assert( std::abs(omega_sum+dual.vols[i])<TOL );
    }
  }

  // DirectedLink idx2DirectedLinkReverse( const Idx il ) const {
  //   const int df = il%3;
  //   const Idx iff = il/3;
  //   const FaceCoords f1 = dual.idx2FaceCoords(iff);

  //   FaceCoords f2;
  //   dual.shift(f2, f1, df);

  //   return std::make_pair(dual.idx(f2), df);
  // }

  void FillAlphas(){
    alphas.clear();
    alphas.resize( dual.NDirectedLinks() );

    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );
        const V3 rf2 = dual.vertices[if2];

        alphas[dual.linkidx(if1, df)] = alpha( rf1, rf2, 0.0 );
      }
    }
  }



  const int limit = 5000; // 1000;
  double epsabs = 1.0e-15; // 0.;
  double epsrel = 1.0e-13; // TOLLOOSE;
  int key = 6;



  void FillOmegas(){
    omegas.clear();
    omegas.resize( dual.NDirectedLinks() );

    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );
        const V3 rf2 = dual.vertices[if2];

        double omega12;
        {
          gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
          double result, error;

          Func f1 = [&](const double stilde){
            const V3 p_ = p(rf1, rf2, stilde);
            const double dphi_dstilde_ = dphi_dstilde(rf1, rf2, stilde);
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

        {
          const double alpha_1 = alpha( rf1,rf2,0.0);
          const double alpha_2 = alpha( rf1,rf2,1.0);

          double diff = alpha_2 - omega12 - alpha_1;
          // std::cout << "debug. diff = " << diff << std::endl;
          if(diff>M_PI) diff -= 2.0*M_PI;
          else if(diff<=-M_PI) diff += 2.0*M_PI;
          assert( std::abs(diff)<TOL );
        }

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

        omegas[dual.linkidx(if1, df)] = omega12;
      }
    }
  }


};
