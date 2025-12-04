#pragma once


constexpr int XZ = 0;
constexpr int ZY = 1;

using Vertices = std::vector<V3>;

using Coords = std::array<int, 3>; // (x,y,s)
using FaceCoords = std::array<int, 4>; // (x,y,s,XZ/ZY)


using BasePoints = std::vector<Idx>;
using BaseTypes = std::vector<int>;


#include<utility>
using Link=std::pair<Idx,Idx>;
using DirectedLink=std::pair<Idx,int>;
using DualFace=std::vector<DirectedLink>;



inline V3 project( const V3& x ) { return x/x.norm(); }

V3 circumcenter(const V3& r0, const V3& r1, const V3& r2){
  const V3 r10 = r1 - r0;
  const V3 r20 = r2 - r0;

  const V3 tmp1 = r10.squaredNorm() * r20 - r20.squaredNorm() * r10;
  const V3 cross = r10.cross(r20);
  const V3 numer = tmp1.cross(cross);
  const double denom = 2.0*cross.squaredNorm();

  return numer/denom + r0;
}

double arcLength( const V3& rA, const V3& rB ){
  const double AB = rA.dot(rB) / rA.norm() / rB.norm();
  return std::acos(AB);
}

double sphericalarea( const V3& x1, const V3& x2, const V3& x3 ){
  double a = arcLength( x2, x3 );
  double b = arcLength( x3, x1 );
  double c = arcLength( x1, x2 );
  double s = 0.5*(a+b+c);
  double tantan = std::tan(0.5*s) * std::tan(0.5*(s-a)) * std::tan(0.5*(s-b)) * std::tan(0.5*(s-c));
  return 4.0*std::atan(std::sqrt(tantan));
}




Vertices operator+ ( Vertices r1, const Vertices& r2 ){ // copy
  assert(r1.size()==r2.size());
  for(Idx i=0; i<r1.size(); i++) r1[i] += r2[i];
  return r1;
}

Vertices operator- ( Vertices r1, const Vertices& r2 ){ // copy
  assert(r1.size()==r2.size());
  for(Idx i=0; i<r1.size(); i++) r1[i] -= r2[i];
  return r1;
}

double squaredNorm ( const Vertices& r ){
  double res = 0.0;
  for(Idx i=0; i<r.size(); i++) res += r[i].squaredNorm();
  return res;
}


#include <set>

struct RefinedIcosahedron : Icosahedron{
  const int L;
  Vertices vertices;

  BasePoints basePoints;
  BaseTypes baseTypes;

  // -----------------------------

  RefinedIcosahedron( const int L_ )
    : L(L_)
  {
    FillVerticesSimpleRefinement();
    GetBasePoints();
  }

  // -----------------------------

  Idx NVertices() const { return 10*L*L+2; }
  Idx NLinks() const { return 30*L*L; }
  Idx NFaces() const { return 20*L*L; }

  std::vector<std::set<Idx>> nns;

  void set_nns(){
    for(Idx in=0; in<NVertices(); in++){
      std::set<Idx> nn;
      const Coords n = idx2Coords(in);
      Coords np;
      for(int mu=0; mu<6; mu++) {
        shift( np, n, mu );
        if(np[2]!=PatchIdxVoid) nn.insert(idx(np));
      }
    }
  }

  std::string print() const {
    std::stringstream ss;
    Idx counter=0;
    ss << std::scientific << std::setprecision(25);
    for(const V3& x : vertices){
      ss << counter << ": " << x.transpose() << std::endl;
      counter++;
    }
    return ss.str();
  }

  // Idx
  // first 10L^2 is for the regular points
  // for irregular points
  Idx idxS() const { return 10*L*L; }
  Idx idxN() const { return 10*L*L + 1; }

  Idx idx( const int xl0, const int xl1, const int s ) const {
    assert(s!=PatchIdxVoid);
    if(s==PatchIdxS) return idxS();
    else if(s==PatchIdxN) return idxN();
    else{
      assert( xl0 < L );
      assert( xl1 < L );
      return s*L*L + xl1*L + xl0;
    }
  }
  Idx idx( const Coords& n ) const { return idx(n[0], n[1], n[2]); }

  Coords idx2Coords( Idx i ) const { // copy i
    if(i==idxS()) return Coords{0,0,PatchIdxS};
    else if(i==idxN()) return Coords{0,0,PatchIdxN};
    else{
      const int xl0 = i%L;
      i /= L;
      const int xl1 = i%L;
      i /= L;
      const int s = i;
      assert( s<NPatches );

      return Coords{xl0, xl1, s};
    }
  }

  ////////


  // -----------------------------

  void FillVerticesSimpleRefinement() { // (pA, pB, pD), (pA, pD, pC): right-handed
    vertices.clear();
    vertices.resize( NVertices() );

    for(int s=0; s<NPatches; s++){
      const V3 pA = icos[patches[s][0]];
      const V3 pB = icos[patches[s][1]];
      const V3 pC = icos[patches[s][2]];
      const V3 pD = icos[patches[s][3]];

      // xz face
      //      pD
      //      /|
      //     / |
      //    /  |
      //   /   |
      //  /    |
      // pA----pB
      const V3 vAB = pB-pA;
      const V3 vBD = pD-pB;
      for(int i=0; i<L; i++){
        const double ri = 1.0*i/L;
        for(int j=0; j<=i; j++){
          const double rj = 1.0*j/L;
          const V3 r = pA + ri*vAB + rj*vBD;
          vertices[idx(i,j,s)] = project( r );
        }}

      // yz face
      // pC----pD
      // |     /
      // |    /
      // |   /
      // |  /
      // | /
      // pA
      const V3 vCD = pD-pC;
      const V3 vAC = pC-pA;
      for(int j=0; j<L; j++){
        const double rj = 1.0*j/L;
        for(int i=0; i<j; i++){
          const double ri = 1.0*i/L;
          const V3 r = pA + ri*vCD + rj*vAC;
          vertices[idx(i,j,s)] = project( r );
        }}
    }

    vertices[idxS()] = icos[0];
    vertices[idxN()] = icos[NIcosVertices-1];
  }


  int shiftPX( Coords& np,
               const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    // std::cout << "debug. s = " << s << std::endl;
    assert(0<=s && s<PatchIdxVoid);
    if(s==PatchIdxS) {
      np = Coords{ L-1,0,0 };
      return 0;
    }
    if(s==PatchIdxN) {
      np = Coords{ 0,L-1,5 };
      return 1;
    }

    int res = 0;
    if(x==L-1){ // steps over the patch bdy
      // res = +1;
      if(s<NPatches/2){ // southern hemisphere
        if(y==0) np = Coords{-1, -1, PatchIdxS};
        else np = Coords{L-y, 0, (s+1)%(NPatches/2)};
      }
      else np = Coords{0, y, (s-4)%(NPatches/2)}; // northern hemisphere
    }
    else np = Coords{x+1, y, s}; // within patch
    return res;
  }

  // only defined for regular points
  int shiftMX( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    // assert(0<=s && s<NPatches);
    assert(0<=s && s<PatchIdxVoid);
    if(s==PatchIdxS) {
      np = Coords{ L-1,0,3 };
      return 0;
    }
    if(s==PatchIdxN) {
      np = Coords{ 0,L-1,8 };
      return 1;
    }


    int res = 0;
    if(x==0){ // steps over the patch bdy
      // res = 0;
      if(s<NPatches/2) {
        // std::cout << "debug mx s." << std::endl; // used
        np = Coords{L-1, y, (s+4)%(NPatches/2)+NPatches/2}; // southern hemisphere
      }
      else {
        // std::cout << "debug mx n." << std::endl; // used
        res = 2;
        np = Coords{L-1-y, L-1, (s-1)%(NPatches/2)+NPatches/2}; // northern hemisphere
      }
    }
    else np = Coords{x-1, y, s}; // within patch
    return res;
  }

  // only defined for regular points
  int shiftPY( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    // assert(0<=s && s<NPatches);
    assert(0<=s && s<PatchIdxVoid);
    if(s==PatchIdxS) {
      np = Coords{ L-1,0,1 };
      return 2;
    }
    if(s==PatchIdxN) {
      np = Coords{ 0,L-1,6 };
      return 0;
    }


    int res = 0;
    if(y==L-1) { // steps over the patch bdy
      // res = +1;
      if(s<NPatches/2) np = Coords{x, 0, s+5}; // southern hemisphere
      else { // northern hemisphere
        if(x==0) np = Coords{-1, -1, PatchIdxN};
        else np = Coords{0, L-x, (s+1)%(NPatches/2)+NPatches/2};
      }
    }
    else np = Coords{x, y+1, s}; // within patch
    return res;
  }

  // only defined for regular points
  int shiftMY( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    // assert(0<=s && s<NPatches);
    assert(0<=s && s<PatchIdxVoid);
    if(s==PatchIdxS) {
      np = Coords{ L-1,0,4 };
      return 2;
    }
    if(s==PatchIdxN) {
      np = Coords{ 0,L-1,9 };
      return 0;
    }


    int res = 0;
    if(y==0){ // steps over the patch bdy
      res = +1;
      if(s<NPatches/2) {
        // std::cout << "debug my s." << std::endl; // used
        np = Coords{L-1, L-1-x, (s-1+NPatches/2)%(NPatches/2)}; // southern hemisphere
      }
      else {
        // std::cout << "debug my n." << std::endl; // used
        np = Coords{x, L-1, (s-5)%(NPatches/2)}; // northern hemisphere
      }
    }
    else np = Coords{x, y-1, s}; // within patch
    return res;
  }


  // only defined for regular points
  int shiftPZ( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    // assert(0<=s && s<NPatches);
    assert(0<=s && s<PatchIdxVoid);
    if(s==PatchIdxS) {
      np = Coords{ L-1,0,2 };
      return 1;
    }
    if(s==PatchIdxN) {
      np = Coords{ 0,L-1,7 };
      return 2;
    }


    int res = 0;
    if(s<NPatches/2){ // southern hemisphere
      // res = 1;
      if(x==L-1) np = Coords{L-1-y, 0, (s+1)%(NPatches/2)};
      else if(y==L-1) np = Coords{x+1, 0, s+5};
      else np = Coords{x+1, y+1, s};
    }
    else { // northern hemisphere
      if(y==L-1) np = Coords{0, L-1-x, (s+1)%(NPatches/2)+NPatches/2};
      else if(x==L-1) np = Coords{0, y+1, (s-4)%(NPatches/2)};
      else np = Coords{x+1, y+1, s};
    }
    return res;
  }

  // only defined for regular points
  int shiftMZ( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    // assert(0<=s && s<NPatches);
    assert(0<=s && s<PatchIdxVoid);
    if(s==PatchIdxS) {
      np = Coords{ 0,0,PatchIdxVoid };
      return 0;
    }
    if(s==PatchIdxN) {
      np = Coords{ 0,0,PatchIdxVoid };
      return 0;
    }


    int res = 0;
    if(s<NPatches/2){ // southern hemisphere
      if(x==0) {
        // res = 2;
        if(y==0) {
          // std::cout << "debug mz s1." << std::endl; // used
          res += 1;
          np = Coords{0, 0, PatchIdxVoid};
        }
        else {
          // std::cout << "debug mz s2." << std::endl; // used
          np = Coords{L-1, y-1, (s+4)%(NPatches/2)+NPatches/2};
        }
      }
      else if(y==0) {
        // std::cout << "debug mz s3." << std::endl; // used
        res += 1;
        np = Coords{L-1, L-x, (s-1+NPatches/2)%(NPatches/2)};
      }
      else np = Coords{x-1, y-1, s};
    }
    else { // northern hemisphere
      if(x==0) {
        // res = 2;
        if(y==0) {
          // std::cout << "debug mz n1." << std::endl; // used
          res += 1;
          np = Coords{0, 0, PatchIdxVoid};
        }
        else {
          // std::cout << "debug mz n2." << std::endl;
          np = Coords{L-y, L-1, (s-1)%(NPatches/2)+NPatches/2};
        }
      }
      else if(y==0) {
        // std::cout << "debug mz n3." << std::endl;
        res += 1;
        np = Coords{x-1, L-1, (s-5+NPatches/2)%(NPatches/2)};
      }
      else np = Coords{x-1, y-1, s};
    }
    return res;
  }

  const int nX=0;
  const int nY=1;
  const int nZ=2;

  int shift( Coords& np,
             const Coords& n,
             const int dn ) const {
    assert(0<=dn && dn<6);

    int res;
    if(dn==nX) res = shiftPX(np, n);
    else if(dn==nY) res = shiftPY(np, n);
    else if(dn==nZ) res = shiftPZ(np, n);
    else if(dn==nX+3) res = shiftMX(np, n);
    else if(dn==nY+3) res = shiftMY(np, n);
    else if(dn==nZ+3) res = shiftMZ(np, n);
    else assert(false);
    return res;
  }


  // basepoints
  void GetBasePoints(){
    assert( L==1 || L%2==0 );
    basePoints.clear();
    baseTypes.clear();
    constexpr int s=0;

    // type 0
    basePoints.push_back( idx(0,0,s) );
    baseTypes.push_back( 0 );

    if(L!=1){
      basePoints.push_back( idx(L/2,0,s) );
      baseTypes.push_back( 0 );

      // type 1
      for(int x=L/2+1; 2*(2*x-L)<=x; x+=2){
        basePoints.push_back( idx(x,2*x-L,s) );
        baseTypes.push_back( 1 );
      }

      // type 2
      for(int x=2; 3*x/2<=L; x+=2){
        basePoints.push_back( idx(x,x/2,s) );
        baseTypes.push_back( 2 );
      }

      // type 3
      for(int x=1; x<L/2; x++){
        basePoints.push_back( idx(x,0,s) );
        baseTypes.push_back( 3 );
      }

      // type 4
      for(int x=0; x<=L/2; x++){
        for(int y=0; y<=L/2; y++){
          if( y==0 || 0>=x-2*y || 2*x-y>=L ) continue;
          basePoints.push_back( idx(x,y,s) );
          baseTypes.push_back( 4 );
        }
      }
    }
  }


};




struct RefinedIcosahedronDual {
  const RefinedIcosahedron& simplicial;
  const int L;
  const Idx Nx;

  Vertices vertices;

  BasePoints basePoints;
  BaseTypes baseTypes;

  std::vector<DualFace> faces;

  const int nA=0;
  const int nB=1;
  const int nC=2;

  std::vector<double> ells;
  std::vector<double> link_volumes;

  double mean_ell;
  double mean_link_volume;
  double mean_dual_area;

  std::vector<double> vols; // hex vol
  double mean_vol;

  RefinedIcosahedronDual(const RefinedIcosahedron& simplicial_)
    : simplicial(simplicial_)
    , L(simplicial.L)
    , Nx( simplicial.idx(Coords{L-1,L-1,NPatches-1})+1 )
  {
    FillPointsCircumCenterDual();
    assert(vertices.size()==2*Nx);
    GetBasePoints();
    FillFaces();
    FillDualLinks();
    set_ell_link_volume();
    set_vol();
  }

  Idx NVertices() const { return simplicial.NFaces(); }
  Idx NLinks() const { return simplicial.NLinks(); }
  Idx NFaces() const { return simplicial.NVertices(); }

  Idx idx( const FaceCoords& n ) const {
    const Idx in = simplicial.idx( Coords{n[0], n[1], n[2]} );
    return 2*in+n[3];
  }

  FaceCoords idx2FaceCoords( Idx i ) const { // copy i
    const int type = i%2;
    i /= 2;
    const Coords n = simplicial.idx2Coords(i);
    return FaceCoords{n[0],n[1],n[2],type};
  }


  int shiftPA( FaceCoords& fp,
                const FaceCoords& f ) const {
    const int x = f[0];
    const int y = f[1];
    const int s = f[2];
    const int type = f[3];

    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);
    assert(type==XZ);

    Coords np;
    simplicial.shiftPX(np, Coords{x,y,s});

    int res=0;
    if(np[2]==PatchIdxS) { // rotate around SP
      assert(x==L-1);
      assert(y==0);
      fp = FaceCoords{x, y, (s+1)%(NPatches/2), XZ};
    }
    else if(s<NPatches/2 && x==L-1) {
      fp = FaceCoords{np[0]-1, np[1], np[2], XZ};
      res = 1;
    }
    else fp = FaceCoords{np[0], np[1], np[2], ZY};
    return res;
  }

  int shiftMA( FaceCoords& fp,
                const FaceCoords& f ) const {
    const int x = f[0];
    const int y = f[1];
    const int s = f[2];
    const int type = f[3];

    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);
    assert(type==ZY);

    Coords np;
    simplicial.shiftMX(np, Coords{x,y,s});

    if(NPatches/2<=s && x==0) fp = FaceCoords{np[0], np[1], np[2], ZY};
    else fp = FaceCoords{np[0], np[1], np[2], XZ};
    return 0;
  }

  int shiftPB( FaceCoords& fp,
                const FaceCoords& f ) const {
    const int x = f[0];
    const int y = f[1];
    const int s = f[2];
    const int type = f[3];

    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);
    assert(type==XZ);

    fp = FaceCoords{x, y, s, ZY};
    return 0;
  }

  int shiftMB( FaceCoords& fp,
                const FaceCoords& f ) const {
    const int x = f[0];
    const int y = f[1];
    const int s = f[2];
    const int type = f[3];

    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);
    assert(type==ZY);

    fp = FaceCoords{x, y, s, XZ};
    return 0;
  }

  int shiftPC( FaceCoords& fp,
                const FaceCoords& f ) const {
    const int x = f[0];
    const int y = f[1];
    const int s = f[2];
    const int type = f[3];

    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);
    assert(type==XZ);

    Coords np;
    simplicial.shiftMY(np, Coords{x,y,s});

    if(s<NPatches/2 && y==0) fp = FaceCoords{np[0], np[1], np[2], XZ};
    else fp = FaceCoords{np[0], np[1], np[2], ZY};
    return 0;
  }

  int shiftMC( FaceCoords& fp,
                const FaceCoords& f ) const {
    const int x = f[0];
    const int y = f[1];
    const int s = f[2];
    const int type = f[3];

    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);
    assert(type==ZY);

    Coords np;
    simplicial.shiftPY(np, Coords{x,y,s});

    int res = 0;
    if(np[2]==PatchIdxN) { // rotate around NP
      assert(y==L-1);
      assert(x==0);
      fp = FaceCoords{x, y, (s+1)%(NPatches/2)+(NPatches/2), ZY};
      res = 1;
    }
    else if(NPatches/2<=s && y==L-1) {
      fp = FaceCoords{np[0], np[1]-1, np[2], ZY};
      res = 1;
    }
    else {
      fp = FaceCoords{np[0], np[1], np[2], XZ};
    }
    return res;
  }


  int shift( FaceCoords& fp,
              const FaceCoords& f,
              const int df ) const {
    const int type = f[3];

    int res = 0;
    if(type==XZ){
      if(df==nA) res = shiftPA(fp, f);
      else if(df==nB) res = shiftPB(fp, f);
      else if(df==nC) res = shiftPC(fp, f);
      else assert(false);
    }
    else if(type==ZY){
      if(df==nA) res = shiftMA(fp, f);
      else if(df==nB) res = shiftMB(fp, f);
      else if(df==nC) res = shiftMC(fp, f);
      else assert(false);
    }
    else assert(false);
    return res;
  }



  // with circumcenterdual
  void FillPointsCircumCenterDual(){
    vertices.clear();

    for(int s=0; s<NPatches; s++){
      for(int y=0; y<L; y++){
        for(int x=0; x<L; x++){
          const Coords n{x,y,s};
          const Idx in = simplicial.idx(n);

          Coords nPX, nPY, nPZ;
          simplicial.shiftPX( nPX, n );
          simplicial.shiftPY( nPY, n );
          simplicial.shiftPZ( nPZ, n );
          const Idx inPX = simplicial.idx(nPX);
          const Idx inPY = simplicial.idx(nPY);
          const Idx inPZ = simplicial.idx(nPZ);

          const V3 vn = simplicial.vertices[in];
          const V3 vnPX = simplicial.vertices[inPX];
          const V3 vnPY = simplicial.vertices[inPY];
          const V3 vnPZ = simplicial.vertices[inPZ];

          // xz face
          V3 r = project( circumcenter(vn, vnPX, vnPZ) );
          assert(vertices.size()==idx(FaceCoords{x,y,s,XZ}));
          vertices.push_back( r );

          // zy face
          r = project( circumcenter(vn, vnPZ, vnPY ) );
          assert(vertices.size()==idx(FaceCoords{x,y,s,ZY}));
          vertices.push_back( r );
        }}}

    assert( vertices.size()==NVertices() );
  }

  // @@
  // basepoints
  void GetBasePoints(){
    basePoints.clear();
    baseTypes.clear();
    constexpr int s=0;

    for(int x=0; x<L; x++){
      for(int y=0; y<L; y++){
        for( const int type : std::vector<int>{XZ, ZY} ){

          // orthogonal coordinates (xp, yp)
          // xp == x; yp==y-x/2;
          // xp3 == 3*xp, yp2 == 2y
          int xp3, yp2;

          if(type==XZ){
            xp3 = 3*x + 2;
            yp2 = 2*y - x;
          }
          else { // type=ZY
            xp3 = 3*x + 1;
            yp2 = 2*y - x + 1;
          }

          if( 3*yp2+xp3>=0 && xp3-yp2-2*L<=0 && yp2<=0 ){ // if in fund domain
            const bool isType1 = (xp3-yp2-2*L==0);
            const bool isType2 = (yp2==0);
            const bool isType3 = (3*yp2+xp3==0);

            if(isType3) assert( false ); // no type 3

            basePoints.push_back( idx(FaceCoords{x, y, s, type}) );
            if(isType1&&isType2) baseTypes.push_back( 0 ); // type 0
            else if(isType1) baseTypes.push_back( 1 ); // type 1
            else if(isType2) baseTypes.push_back( 2 ); // type 2
            else baseTypes.push_back( 4 ); // type 4
          }
        }}}
  }


  inline Idx linkidx( const Idx iff, const int df ) const {
    return 3*iff + df;
  }

  inline Idx linkidx( const DirectedLink& ell ) const {
    return linkidx(ell.first, ell.second);
  }

  DirectedLink linkidx2DirectedLink( const Idx il ) const {
    const int df = il%3;
    const Idx iff = il/3;
    return std::make_pair(iff, df);
  }

  inline Idx NDirectedLinks() const {
    return 3*NVertices();
  }

  void FillFaces(){
    for(Idx if1=0; if1<NVertices(); if1++){
      const FaceCoords f1 = idx2FaceCoords(if1);
      if(f1[3]!=XZ) continue;

      DualFace face;
      FaceCoords f2 = f1;

      bool is_south = (f1[2]<5);

      int df = nC;
      face.push_back(std::make_pair(if1,df));

      int type_tmp = f2[3];
      int is_shift = shift( f2, f1, df );
      Idx if2 = idx(f2);

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
          }
          else if( is_shift==1 ) {
            df = (df+2)%3;
          }
          else {
            df = df;
          }
        }
        face.push_back(std::make_pair(if2,df));

        type_tmp = f2[3];
        FaceCoords f3;
        is_shift = shift( f3, f2, df );
        f2 = f3;
        if2 = idx(f2);

        if(if2==if1) break;
      } // end for
      if(counter==6){
        assert(false);
      }
      faces.push_back(face);
    }

    // south pole
    {
      DualFace face;
      for(int s=0; s<NPatches/2; s++){
        const FaceCoords f1{L-1, 0, s, XZ};
        const Idx if1 = idx(f1);
        const int df = nA;
        face.push_back(std::make_pair(if1,df));
      }
      faces.push_back(face);
    }
    // north pole
    {
      DualFace face;
      for(int s=NPatches/2; s<NPatches; s++){
        const FaceCoords f1{0, L-1, s, ZY};
        const Idx if1 = idx(f1);
        const int df = nA;
        face.push_back(std::make_pair(if1,df));
      }
      faces.push_back(face);
    }
  }


  std::map<const DirectedLink, const DirectedLink> linkSimp2Dual;

  void FillDualLinks(){
    linkSimp2Dual.clear();
    for(Idx in=0; in<simplicial.NVertices(); in++){
      const Coords n = simplicial.idx2Coords(in);
      const DualFace& face = faces[in];
      auto pf = face.begin();
      if(face.size()==6){
        // x=0
        linkSimp2Dual.insert( { DirectedLink{in,0}, *pf } ); pf++;
        // -y=4
        linkSimp2Dual.insert( { DirectedLink{in,4}, *pf } ); pf++;
        // -z=5
        linkSimp2Dual.insert( { DirectedLink{in,5}, *pf } ); pf++;
        // -x=3
        linkSimp2Dual.insert( { DirectedLink{in,3}, *pf } ); pf++;
        // y=1
        linkSimp2Dual.insert( { DirectedLink{in,1}, *pf } ); pf++;
        // z=2
        linkSimp2Dual.insert( { DirectedLink{in,2}, *pf } ); pf++;
      }
      else{
        // x=0
        linkSimp2Dual.insert( { DirectedLink{in,0}, *pf } ); pf++;
        // -y=4
        linkSimp2Dual.insert( { DirectedLink{in,4}, *pf } ); pf++;
        // // -z=5
        // linkSimp2Dual.insert( { DirectedLink{in,5}, *pf } ); pf++;
        // -x=3
        linkSimp2Dual.insert( { DirectedLink{in,3}, *pf } ); pf++;
        // y=1
        linkSimp2Dual.insert( { DirectedLink{in,1}, *pf } ); pf++;
        // z=2
        linkSimp2Dual.insert( { DirectedLink{in,2}, *pf } ); pf++;
      }
      // else assert(false);
    }
  }

  // void get_dual_link( FaceCoords& fA,
  //                     FaceCoords& fB,
  //                     const Coords n_,
  //                     const int mu_ ) const {

  //   Coords n;
  //   int mu;
  //   if(n_[2]==PatchIdxS){
  //     n = Coords{L-1,0,mu_};
  //     mu = 0;
  //   }
  //   else if(n_[2]==PatchIdxN){
  //     n = Coords{0,L-1,mu_+5};
  //     mu = 1;
  //   }
  //   else if(0<=mu_ && mu_<3){
  //     n=n_;
  //     mu = mu_;
  //   }
  //   else if(3<=mu_ && mu_<6){
  //     int dn = simplicial.shift(n, n_, mu_);
  //     // std::cout << "debug. " << dn << std::endl;
  //     if(n_[2]>=5) dn=4-dn; // @@@
  //     mu = (mu_+dn)%3;
  //   }
  //   else assert(false);

  //   if(mu==simplicial.nX) {
  //     fA = FaceCoords{n[0], n[1], n[2], XZ};
  //     shiftPC(fB, fA);
  //   }
  //   else if(mu==simplicial.nY) {
  //     fB = FaceCoords{n[0], n[1], n[2], ZY};
  //     shiftMA(fA, fB);
  //   }
  //   else if(mu==simplicial.nZ) {
  //     fA = FaceCoords{n[0], n[1], n[2], ZY};
  //     fB = FaceCoords{n[0], n[1], n[2], XZ};
  //   }
  //   // else if(mu==simplicial.nX+3) {
  //   //   fA = FaceCoords{n[0], n[1], n[2], XZ};
  //   //   shiftPC(fB, fA);
  //   // }
  //   // else if(mu==simplicial.nY+3) {
  //   //   fB = FaceCoords{n[0], n[1], n[2], ZY};
  //   //   shiftMA(fA, fB);
  //   // }
  //   // else if(mu==simplicial.nZ+3) {
  //   //   fA = FaceCoords{n[0], n[1], n[2], ZY};
  //   //   fB = FaceCoords{n[0], n[1], n[2], XZ};
  //   // }
  //   else {
  //     // std::cout << "debug. mu_ = " << mu_ << std::endl;
  //     std::cout << "debug. mu = " << mu << std::endl;
  //     assert( false );
  //   }
  // }

  void dual_simp_link( Coords& pA,
                       Coords& pB,
                       const FaceCoords& f1,
                       const int df
                       ) const {
    const Coords base{ f1[0], f1[1], f1[2] };
    if(f1[3]==XZ){
      if(df==nA){
        simplicial.shiftPX( pA, base );
        simplicial.shiftPZ( pB, base );
      }
      else if(df==nB){
        simplicial.shiftPZ( pA, base );
        pB = base;
      }
      else if(df==nC){
        pA = base;
        simplicial.shiftPX( pB, base );
      }
      else assert(false);
    }
    else if(f1[3]==ZY){
      if(df==nA){
        simplicial.shiftPY( pA, base );
        pB = base;
      }
      else if(df==nB){
        pA = base;
        simplicial.shiftPZ( pB, base );
      }
      else if(df==nC){
        simplicial.shiftPZ( pA, base );
        simplicial.shiftPY( pB, base );
      }
      else assert(false);
    }
    else assert(false);
  }


  void set_ell_link_volume(){
    ells.resize( NDirectedLinks() );
    link_volumes.resize( NDirectedLinks() );

    for(Idx il=0; il<NDirectedLinks(); il++){
      const DirectedLink ell = linkidx2DirectedLink(il);
      const Idx if1 = ell.first;
      const FaceCoords f1 = idx2FaceCoords( if1 );
      const int df = ell.second;

      FaceCoords f2;
      shift( f2, f1, df);
      const Idx if2 = idx(f2);

      Coords pA, pB;
      dual_simp_link(pA, pB, f1, df);
      const Idx iA = simplicial.idx(pA);
      const Idx iB = simplicial.idx(pB);

      const V3 x = vertices[if1];
      const V3 y = vertices[if2];
      ells[il] = arcLength(x,y);

      const double areaA = sphericalarea(simplicial.vertices[iA],x,y);
      const double areaB = sphericalarea(simplicial.vertices[iB],x,y);
      link_volumes[il] = (areaA + areaB);
    }

    mean_link_volume = 0.0;
    for(const double elem : link_volumes) {
      // std::cout << "debug. elem = " << elem << std::endl;
      mean_link_volume+=elem;
    }
    // std::cout << "debug. diff = " << mean_link_volume-8.0*M_PI << std::endl;
    assert( std::abs(mean_link_volume-8.0*M_PI)/std::sqrt(link_volumes.size())<1.0e-12 ); // double count
    mean_link_volume /= link_volumes.size();

    mean_ell = 0.0;
    for(const double elem : ells) mean_ell+=elem;
    mean_ell /= ells.size();
  }


  void set_vol(){
    vols.clear();
    assert( faces.size()==simplicial.NVertices() );
    for(Idx i=0; i<simplicial.NVertices(); i++){
      const V3 p = simplicial.vertices[i];

      const DualFace face = faces[i];
      double area = 0.0;
      for(const DirectedLink& link : face){
        const Idx if1 = link.first;
        const FaceCoords f1 = idx2FaceCoords( if1 );
        const V3 rf1 = vertices[if1];

        const int df = link.second;
        FaceCoords f2;
        shift( f2, f1, df);
        const Idx if2 = idx(f2);
        const V3 rf2 = vertices[if2];

        area += sphericalarea( p, rf1, rf2 );
      }
      vols.push_back( area );
    }
    mean_vol = 0.0;
    for(const double elem : vols) {
      // std::cout << "debug. elem = " << elem << std::endl;
      mean_vol+=elem;
    }
    assert( std::abs(mean_vol-4.0*M_PI)/std::sqrt(link_volumes.size())<1.0e-12 );
    mean_vol /= vols.size();
  }




};

