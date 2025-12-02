#pragma once


using Face=std::vector<Idx>;


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
    const V3 dr = r2-r1;
    const V3 p_ = p(r1,r2,stilde);
    const double tmp = p_.dot(dr);
    return dr - tmp*p_;
  }

  V3 e(const V3& r1, const V3& r2, const double stilde) const {
    const V3 tmp = etilde(r1,r2,stilde);
    return tmp/tmp.norm();
  }

  inline double theta(const V3& r) const {
    assert( std::abs(r.norm()-1.0)<1.0e-14 );
    return std::acos(r[2])/r.norm();
  }

  inline double phi(const V3& r) const {
    assert( std::abs(r.norm()-1.0)<1.0e-14 );
    return std::atan2(r[1], r[0]);
  }


  double dtheta_dstilde(const V3& r1, const V3& r2, const double stilde) const {
    const V3 p_ = p(r1,r2,stilde);
    const V3 e_ = e(r1,r2,stilde);
    return -e_[2]/std::sqrt(1.0-p_[2]*p_[2]);
  }


  double dphi_dstilde(const V3& r1, const V3& r2, const double stilde) const {
    const V3 p_ = p(r1,r2,stilde);
    const V3 e_ = e(r1,r2,stilde);
    return (p_[0]*e_[1]-p_[1]*e_[0])/(1.0-p_[2]*p_[2]);
  }

  inline double alpha(const V3& r1, const V3& r2, const double stilde) const {
    const V3 p_ = p(r1,r2,stilde);
    const double xx = dtheta_dstilde(r1,r2,stilde);
    const double yy = std::sin(theta(p_)) * dphi_dstilde(r1,r2,stilde);
    return std::atan2( yy, xx );
  }

  std::vector<double> alphas;
  std::vector<double> omegas;

  std::vector<Face> faces;

  DualSpinStructure(const RefinedIcosahedronDual& dual_)
    : dual(dual_)
  {
    FillAlphas();
    FillOmegas();
    FillFaces();
  }

  inline Idx idx( const Idx iff, const int df ) const {
    return 3*iff + df;
  }

  void FillAlphas(){
    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );
        const V3 rf2 = dual.vertices[if2];

        alphas.push_back( alpha( rf1, rf2, 0.0 ) );
      }
    }
  }

  void FillOmegas(){
    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );

        double alpha1 = alphas[ idx(if1, df) ];
        double alpha2 = alphas[ idx(if2, df) ] + M_PI;
        double omega = alpha2 - alpha1;
        while(omega>M_PI) omega -= 2.0*M_PI;
        while(omega<=-M_PI) omega += 2.0*M_PI;
        omegas.push_back( omega );
      }
    }
  }

  void FillFaces(){
    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1);
      if(f1[3]!=XZ) continue;

      Face face;
      FaceCoords f2 = f1;

      std::cout << "debug. ### f1 = " << std::endl;
      for(auto elem : f1) std::cout << elem << " ";
      std::cout << std::endl;

      face.push_back(if1);
      int df = dual.nC;

      // int x_tmp = f2[0];
      // int y_tmp = f2[1];
      int type_tmp = f2[3];
      int is_shift = dual.shift( f2, f1, df );
      Idx if2 = dual.idx(f2);

      std::cout << "debug. f2 = " << std::endl;
      for(auto elem : f2) std::cout << elem << " ";
      std::cout << std::endl;

      int counter;
      for(counter=0; counter<6; counter++){
        face.push_back(if2);
        // df to next point
        if( f2[3] != type_tmp ) df = (df+1)%3; // C->A, A->B, B->C
        else if( is_shift==1 ) df = df;
        else df = (df+2)%3;
        // if( f2[0]==x_tmp && f2[1]==y_tmp ) df = (df+2)%3; // C->A, A->B, B->C
        // else df = (df+1)%3;

        // x_tmp = f2[0];
        // y_tmp = f2[1];
        type_tmp = f2[3];
        is_shift = dual.shift( f2, f2, df );
        std::cout << "debug. f2 = " << std::endl;
        for(auto elem : f2) std::cout << elem << " ";
        std::cout << std::endl;
        if2 = dual.idx(f2);

        if(if2==if1) break;
      } // end for
      if(counter==6){
        std::cout << "debug. face = " << std::endl;
        for(auto elem : face) std::cout << elem << " ";
        std::cout << std::endl;
      }
      faces.push_back(face);
    }
  }

};
