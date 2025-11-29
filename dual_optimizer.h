#pragma once


using XiCoords = std::array<double,3>;
using ZetaCoords = std::array<double,3>;


struct DualLatticeOptimizerBase{
  RefinedIcosahedronDual& dual;
  const Icosahedron& icos;
  Orbits orbits;

  const BasePoints& basePoints;
  const BaseTypes& baseTypes;

  // both project onto the axes intersect with the center
  std::array<V3,3> xiHats;
  std::array<V3,2> xHats;
  std::array<V3,3> zetaHats;

  // pointers collecting all the coordinates
  std::vector<double> p;
  std::vector<Idx> if2pIdx;

  DualLatticeOptimizerBase
  (
   RefinedIcosahedronDual& dual_,
   const Icosahedron& icos_,
   const FullIcosahedralGroup& Ih_,
   const Rotation& rot_
   )
    : dual(dual_)
    , icos(icos_)
    , orbits( dual.vertices,
              dual.basePoints,
              dual.baseTypes,
              icos_, Ih_, rot_ )
    , basePoints( dual.basePoints )
    , baseTypes( dual.baseTypes )
  {
    Idx counter=0;
    if2pIdx.clear();
    for(Idx iff=0; iff<baseTypes.size(); iff++){
      const int type = baseTypes[iff];
      if(type==1) {
        if2pIdx.push_back(counter);
        counter++;
      }
      else if(type==2) {
        if2pIdx.push_back(counter);
        counter++;
      }
      else if(type==4) {
        if2pIdx.push_back(counter);
        counter+=2;
      }
      else if(type!=0) assert(false);
    }

    {
      const V3 pA = icos[icos.patches[0][0]];
      const V3 pB = icos[icos.patches[0][1]];
      const V3 pD = icos[icos.patches[0][3]];

      xiHats[0] = 0.5*(pB+pD)-pA;
      xiHats[0] /= xiHats[0].norm();
      xiHats[1] = pD-pB;
      xiHats[1] /= xiHats[1].norm();
      xiHats[2] = (pA+pB+pD)/3.0;
      xiHats[2] /= xiHats[2].norm();
      assert( xiHats[0].dot(xiHats[1])<1.0e-14 );
      assert( xiHats[1].dot(xiHats[2])<1.0e-14 );
      assert( xiHats[2].dot(xiHats[0])<1.0e-14 );

      xHats[0] = pB-pA;
      xHats[0] /= xHats[0].norm();
      xHats[1] = xiHats[1];

      zetaHats[0] = xHats[0];
      zetaHats[1] = pD - 0.5*(pA+pB);
      zetaHats[1] /= zetaHats[1].norm();
      zetaHats[2] = xiHats[2];
      assert( zetaHats[0].dot(zetaHats[1])<1.0e-14 );
      assert( zetaHats[1].dot(zetaHats[2])<1.0e-14 );
      assert( zetaHats[2].dot(zetaHats[0])<1.0e-14 );
    }

    p.clear();
    p.resize( counter );

    FillBaseXiCoords();
    FillBasePointsFromCoords();
    orbits.RefreshOrbits();
  }



  XiCoords XiProjector( const V3& r ) const {
    XiCoords res;
    for(int mu=0; mu<3; mu++) res[mu] = r.dot(xiHats[mu]);
    return res;
  }

  V3 embedXi( const XiCoords& xi ) const {
    V3 tmp = V3::Zero();
    for(int mu=0; mu<3; mu++) tmp += xi[mu] * xiHats[mu];
    tmp /= tmp.norm();
    return tmp;
  }

  ZetaCoords ZetaProjector( const V3& r ) const {
    ZetaCoords res;
    for(int mu=0; mu<3; mu++) res[mu] = r.dot(zetaHats[mu]);
    return res;
  }

  V3 embedZeta( const ZetaCoords& zeta ) const {
    V3 tmp = V3::Zero();
    for(int mu=0; mu<3; mu++) tmp += zeta[mu] * zetaHats[mu];
    tmp /= tmp.norm();
    return tmp;
  }


  bool is_type1( const V3& r ) const {
    const double inner = r.dot( zetaHats[0] );
    return (std::abs(inner)<1.0e-14);
  }

  bool is_type2( const V3& r ) const {
    const double inner = r.dot( xiHats[1] );
    return (std::abs(inner)<1.0e-14);
  }



  void FillBaseXiCoords(){
    for(Idx iff=0; iff<basePoints.size(); iff++){
      const int type = baseTypes[iff];
      const V3 r = dual.vertices[ basePoints[iff] ];

      if(type==1) {
        assert( is_type1( r ) );

        const double first = 0.0;
        const double second = r.dot( zetaHats[1] );
        p[ if2pIdx.at(iff) ] = second;

        const V3 check = embedZeta( ZetaCoords{ first, second, std::sqrt(1.0-first*first-second*second) } );
        assert( (r-check).norm()<1.0e-14 );
      }
      else if(type==2) {
        assert( is_type2( r ) );

        const double first = r.dot( xiHats[0] );
        const double second = 0.0;
        p[ if2pIdx.at(iff) ] = first;

        const V3 check = embedXi( XiCoords{ first, second, std::sqrt(1.0-first*first-second*second) } );
        assert( (r-check).norm()<1.0e-14 );
      }
      else if(type==4) {
        const double first = r.dot( xiHats[0] );
        const double second = r.dot( xiHats[1] );
        p[ if2pIdx.at(iff) ] = first;
        p[ if2pIdx.at(iff)+1 ] = second;

        const V3 check = embedXi( XiCoords{ first, second, std::sqrt(1.0-first*first-second*second) } );
        assert( (r-check).norm()<1.0e-14 );
      }
    }
  }


  void FillBasePointsFromCoords(){
    for(Idx iff=0; iff<basePoints.size(); iff++){
      const int type = baseTypes[iff];

      if(type==1) {
        const double first = 0.0;
        const double second = p[ if2pIdx.at(iff) ];
        dual.vertices[ basePoints[iff] ] = embedZeta( ZetaCoords{ first, second, std::sqrt(1.0-first*first-second*second) } );
      }
      else if(type==2) {
        const double first = p[ if2pIdx.at(iff) ];
        const double second = 0.0;
        dual.vertices[ basePoints[iff] ] = embedXi( XiCoords{ first, second, std::sqrt(1.0-first*first-second*second) } );
      }
      else if(type==4) {
        const double first = p[ if2pIdx.at(iff) ];
        const double second = p[ if2pIdx.at(iff)+1 ];
        dual.vertices[ basePoints[iff] ] =  embedXi( XiCoords{ first, second, std::sqrt(1.0-first*first-second*second) } );
      }
    }
  }

};




class DualLatticeAngleCostEvaluator : DualLatticeOptimizerBase{
public:

  // <CAB
  // https://www.johndcook.com/blog/spherical_trigonometry/
  double sphericalAngle( const V3& rC, const V3& rA, const V3& rB ) const {
    const double b = arcLength(rC, rA);
    const double c = arcLength(rA, rB);
    const double a = arcLength(rB, rC);

    const double numer = std::cos(a)-std::cos(b)*std::cos(c);
    const double denom = std::sin(b)*std::sin(c);
    return std::acos(numer/denom);
  }


  DualLatticeAngleCostEvaluator
  (
   RefinedIcosahedronDual& dual_,
   const Icosahedron& icos_,
   const FullIcosahedralGroup& Ih_,
   const Rotation& rot_
   )
    :  DualLatticeOptimizerBase(dual_, icos_, Ih_, rot_)
  {}


  double cost( const double tol = 1.0e-12) const {
    double sqrd = 0.0;

    for(const Idx iff : basePoints){

      const V3 rf = dual.vertices[iff];
      const FaceCoords nf = dual.idx2FaceCoords( iff );

      FaceCoords nfA, nfB, nfC;
      if(nf[3]==XZ){
        dual.shiftPA( nfA, nf );
        dual.shiftPB( nfB, nf );
        dual.shiftPC( nfC, nf );
      }
      else if(nf[3]==ZY){
        dual.shiftMA( nfA, nf );
        dual.shiftMB( nfB, nf );
        dual.shiftMC( nfC, nf );
      }
      else assert(false);

      const V3 rfA = dual.vertices[dual.idx(nfA)];
      const V3 rfB = dual.vertices[dual.idx(nfB)];
      const V3 rfC = dual.vertices[dual.idx(nfC)];

      const double A0B = sphericalAngle( rfA, rf, rfB );
      const double B0C = sphericalAngle( rfB, rf, rfC );
      const double C0A = sphericalAngle( rfC, rf, rfA );

      // std::cout << "debug. angles = " << A0B << " " << B0C << " " << C0A << std::endl;
      // std::cout << "debug. sum-2Pi = " << A0B+B0C+C0A-2.0*M_PI << std::endl;
      assert( std::abs( A0B+B0C+C0A-2.0*M_PI ) < tol );

      const double diff = A0B - B0C;
      sqrd += diff*diff;
    }
    return std::sqrt( sqrd/basePoints.size() );
  }

};
