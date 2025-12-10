#pragma once

struct Sphere{
  Sphere(){}

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
};
