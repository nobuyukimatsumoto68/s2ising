#pragma once

<<<<<<< HEAD

using Face=std::vector<Idx>;


struct DualSpinStructure {
=======
#include<utility>



struct DualSpinStructure {
  using DirectedLink=std::pair<Idx,int>;
  using Face=std::vector<DirectedLink>;

>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
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
<<<<<<< HEAD
    return std::acos(r[2])/r.norm();
=======
    const double res = std::acos(r[2]);
    assert( 0<=res && res<=M_PI);
    return res;
>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
  }

  inline double phi(const V3& r) const {
    assert( std::abs(r.norm()-1.0)<1.0e-14 );
<<<<<<< HEAD
    return std::atan2(r[1], r[0]);
=======
    const double res = std::atan2(r[1], r[0]);
    assert( -M_PI<=res && res<=M_PI);
    return res;
>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
  }


  double dtheta_dstilde(const V3& r1, const V3& r2, const double stilde) const {
    const V3 p_ = p(r1,r2,stilde);
    const V3 e_ = e(r1,r2,stilde);
    return -e_[2]/std::sqrt(1.0-p_[2]*p_[2]);
  }

<<<<<<< HEAD

=======
>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
  double dphi_dstilde(const V3& r1, const V3& r2, const double stilde) const {
    const V3 p_ = p(r1,r2,stilde);
    const V3 e_ = e(r1,r2,stilde);
    return (p_[0]*e_[1]-p_[1]*e_[0])/(1.0-p_[2]*p_[2]);
  }

  inline double alpha(const V3& r1, const V3& r2, const double stilde) const {
    const V3 p_ = p(r1,r2,stilde);
    const double xx = dtheta_dstilde(r1,r2,stilde);
    const double yy = std::sin(theta(p_)) * dphi_dstilde(r1,r2,stilde);
<<<<<<< HEAD
    return std::atan2( yy, xx );
  }

=======

    assert( std::abs(xx*xx+yy*yy-1.0)<1.0e-14 );

    return std::atan2( yy, xx );
  }





>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
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

<<<<<<< HEAD
  void FillAlphas(){
=======
  inline Idx idx( const DirectedLink& ell ) const {
    return idx(ell.first, ell.second);
  }

  void FillAlphas(){
    alphas.clear();
    Idx counter=0;
>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
      const V3 rf1 = dual.vertices[if1];

      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );
        const V3 rf2 = dual.vertices[if2];

        alphas.push_back( alpha( rf1, rf2, 0.0 ) );
<<<<<<< HEAD
=======
        assert( counter==idx(if1,df) ); // debug
        counter++;
        assert( counter==alphas.size() ); // debug
>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
      }
    }
  }

  void FillOmegas(){
<<<<<<< HEAD
    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };
=======
    omegas.clear();
    Idx counter=0;
    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1); //{ n[0],n[1],n[2],type };

>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
      for(int df=0; df<3; df++){
        FaceCoords f2;
        dual.shift( f2, f1, df);
        const Idx if2 = dual.idx( f2 );

<<<<<<< HEAD
        double alpha1 = alphas[ idx(if1, df) ];
        double alpha2 = alphas[ idx(if2, df) ] + M_PI;
        double omega = alpha2 - alpha1;
        while(omega>M_PI) omega -= 2.0*M_PI;
        while(omega<=-M_PI) omega += 2.0*M_PI;
        omegas.push_back( omega );
=======
        const double alpha1 = alphas[ idx(if1, df) ];
        const double alpha2 = alphas[ idx(if2, df) ];
        double omega12 = alpha1 - alpha2 - M_PI;
        // while(omega12>M_PI) omega12 -= 2.0*M_PI;
        // while(omega12<=-M_PI) omega12 += 2.0*M_PI;
        omegas.push_back( omega12 );

        assert( counter==idx(if1,df) ); // debug
        counter++;
        assert( counter==omegas.size() ); // debug
>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
      }
    }
  }

  void FillFaces(){
    for(Idx if1=0; if1<dual.NVertices(); if1++){
      const FaceCoords f1 = dual.idx2FaceCoords(if1);
      if(f1[3]!=XZ) continue;

      Face face;
      FaceCoords f2 = f1;

<<<<<<< HEAD
      std::cout << "debug. ### f1 = " << std::endl;
      for(auto elem : f1) std::cout << elem << " ";
      std::cout << std::endl;

      face.push_back(if1);
      int df = dual.nC;
=======
      bool is_south = (f1[2]<5);
      // std::cout << "debug. ### f1 = " << std::endl;
      // for(auto elem : f1) std::cout << elem << " ";
      // std::cout << std::endl;

      int df = dual.nC;
      face.push_back(std::make_pair(if1,df));
      // std::cout << "debug. df = " << df << std::endl;
>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41

      // int x_tmp = f2[0];
      // int y_tmp = f2[1];
      int type_tmp = f2[3];
      int is_shift = dual.shift( f2, f1, df );
      Idx if2 = dual.idx(f2);

<<<<<<< HEAD
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
=======
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
>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
        if2 = dual.idx(f2);

        if(if2==if1) break;
      } // end for
      if(counter==6){
<<<<<<< HEAD
        std::cout << "debug. face = " << std::endl;
        for(auto elem : face) std::cout << elem << " ";
        std::cout << std::endl;
      }
      faces.push_back(face);
    }
=======
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

>>>>>>> ac6cc0274ea3bd53ad375a8b93cedf4cffc02c41
  }

};
