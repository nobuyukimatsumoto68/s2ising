#pragma once


constexpr int XZ = 0;
constexpr int ZY = 1;

using Vertices = std::vector<V3>;

using Coords = std::array<int, 3>; // (x,y,s)
using FaceCoords = std::array<int, 4>; // (x,y,s,XZ/ZY)




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







struct RefinedIcosahedron : Icosahedron{
  const int L;
  Vertices vertices;

  // -----------------------------

  RefinedIcosahedron( const int L_ )
    : L(L_)
  {
    FillVerticesSimpleRefinement();
  }

  // -----------------------------

  Idx NVertices() const { return 10*L*L+2; }
  Idx NLinks() const { return 30*L*L; }
  Idx NFaces() const { return 20*L*L; }

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
    if(s==PatchIdxS) return idxS();
    else if(s==PatchIdxN) return idxN();
    else{
      assert( xl0 < L );
      assert( xl1 < L );
      return s*L*L + xl1*L + xl0;
    }
  }
  Idx idx( const Coords& n ) const { return idx(n[0], n[1], n[2]); }


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

  // only defined for regular points
  void shiftPX( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);

    if(x==L-1){ // steps over the patch bdy
      if(s<NPatches/2){ // southern hemisphere
        if(y==0) np = Coords{-1, -1, PatchIdxS};
        else np = Coords{y, 0, (s+1)%(NPatches/2)};
      }
      else np = Coords{0, y, (s-4)%(NPatches/2)}; // northern hemisphere
    }
    else np = Coords{x+1, y, s}; // within patch
  }

  // only defined for regular points
  void shiftMX( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);

    if(x==0){ // steps over the patch bdy
      if(s<NPatches/2) np = Coords{L-1, y, (s+4)%(NPatches/2)+NPatches/2}; // southern hemisphere
      else np = Coords{L-1-y, L-1, (s-1)%(NPatches/2)+NPatches/2}; // northern hemisphere
    }
    else np = Coords{x-1, y, s}; // within patch
  }

  // only defined for regular points
  void shiftPY( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);

    if(y==L-1) { // steps over the patch bdy
      if(s<NPatches/2) np = Coords{x, 0, s+5}; // southern hemisphere
      else { // northern hemisphere
        if(x==0) np = Coords{-1, -1, PatchIdxN};
        else np = Coords{0, L-x, (s+1)%(NPatches/2)+NPatches/2};
      }
    }
    else np = Coords{x, y+1, s}; // within patch
  }

  // only defined for regular points
  void shiftMY( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);

    if(y==0){ // steps over the patch bdy
      if(s<NPatches/2) np = Coords{L-1, L-1-x, (s-1+NPatches/2)%(NPatches/2)}; // southern hemisphere
      else np = Coords{x, L-1, (s-5)%(NPatches/2)}; // northern hemisphere
    }
    else np = Coords{x, y-1, s}; // within patch
  }


  // only defined for regular points
  void shiftPZ( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);

    if(s<NPatches/2){ // southern hemisphere
      if(x==L-1) np = Coords{L-1-y, 0, (s+1)%(NPatches/2)};
      else if(y==L-1) np = Coords{x+1, 0, s+5};
      else np = Coords{x+1, y+1, s};
    }
    else { // northern hemisphere
      if(y==L-1) np = Coords{0, L-1-x, (s+1)%(NPatches/2)+NPatches/2};
      else if(x==L-1) np = Coords{0, y+1, (s-4)%(NPatches/2)};
      else np = Coords{x+1, y+1, s};
    }
  }

  // only defined for regular points
  void shiftMZ( Coords& np,
                const Coords& n ) const {
    const int x = n[0];
    const int y = n[1];
    const int s = n[2];
    assert(0<=x && x<L);
    assert(0<=y && y<L);
    assert(0<=s && s<NPatches);

    if(s<NPatches/2){ // southern hemisphere
      if(x==0) {
        if(y==0) np = Coords{0, 0, PatchIdxVoid};
        else np = Coords{L-1, y-1, (s+4)%(NPatches/2)+NPatches/2};
      }
      else if(y==0) np = Coords{L-1, L-x, (s-1+NPatches/2)%(NPatches/2)};
      else np = Coords{x-1, y-1, s};
    }
    else { // northern hemisphere
      if(x==0) {
        if(y==0) np = Coords{0, 0, PatchIdxVoid};
        else np = Coords{L-y, L-1, (s-1)%(NPatches/2)+NPatches/2};
      }
      else if(y==0) np = Coords{x-1, L-1, (s-5+NPatches/2)%(NPatches/2)};
      else np = Coords{x-1, y-1, s};
    }
  }

  // @@
  // basepoints


};


struct RefinedIcosahedronDual {
  const RefinedIcosahedron& simplicial;
  const int L;
  const Idx Nx;

  Vertices vertices;

  RefinedIcosahedronDual(const RefinedIcosahedron& simplicial_)
    : simplicial(simplicial_)
    , L(simplicial.L)
    , Nx( simplicial.idx(Coords{L-1,L-1,NPatches-1})+1 )
  {
    FillDualPoints();
    assert(vertices.size()==2*Nx);
  }

  Idx NVertices() const { return simplicial.NFaces(); }
  Idx NLinks() const { return simplicial.NLinks(); }
  Idx NFaces() const { return simplicial.NVertices(); }

  Idx idx( const FaceCoords& n ) const {
    const Idx in = simplicial.idx( Coords{n[0], n[1], n[2]} );
    return 2*in+n[3];
  }


  void shiftPA( FaceCoords& fp,
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

    if(np[2]==PatchIdxS) { // rotate around SP
      assert(x==L-1);
      assert(y==0);
      fp = FaceCoords{x, y, (s+1)%(NPatches/2), XZ};
    }
    else if(s<NPatches/2 && x==L-1) fp = FaceCoords{np[0]-1, np[1], np[2], XZ};
    else fp = FaceCoords{np[0], np[1], np[2], ZY};
  }

  void shiftMA( FaceCoords& fp,
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
  }

  void shiftPB( FaceCoords& fp,
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
  }

  void shiftMB( FaceCoords& fp,
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
  }

  void shiftPC( FaceCoords& fp,
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
  }

  void shiftMC( FaceCoords& fp,
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

    // std::cout << "n = " << x << " " << y << " " << s << std::endl;
    // std::cout << "np = " << np[0] << " " << np[1] << " " << np[2] << std::endl;

    if(np[2]==PatchIdxN) { // rotate around SP
      assert(y==L-1);
      assert(x==0);
      fp = FaceCoords{x, y, (s+1)%(NPatches/2)+(NPatches/2), ZY};
    }
    else if(NPatches/2<=s && y==L-1) fp = FaceCoords{np[0], np[1]-1, np[2], ZY};
    else fp = FaceCoords{np[0], np[1], np[2], XZ};
  }


  void FillDualPoints(){
    vertices.clear();

    Idx counter=0;
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

};

