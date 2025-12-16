#pragma once


#include <random>

template<Idx N>
class Ising{
public:
  DualLoop<N>& loops;
  const Fermion& D;
  const DualSpinStructure& spin;
  const RefinedIcosahedronDual& dual;
  const RefinedIcosahedron& lattice;

  std::map<const DirectedLink, const double> ts_aug;
  std::map<const DirectedLink, const double> betas;
  // std::map<const Link, const double> As;

  Ising( DualLoop<N>& loops_ )
    : loops(loops_)
    , D(loops.D)
    , spin(D.spin)
    , dual(D.dual)
    , lattice(dual.simplicial)
  {
    FillAugmentedTs();
    FillBetas();
  }

  void FillAugmentedTs(){
    ts_aug.clear();
    for(Idx if0=0; if0<dual.NVertices(); if0++){
      // if(!D.is_triv){
      const FaceCoords f0 = dual.idx2FaceCoords( if0 );
      FaceCoords fA, fB, fC;
      dual.shift( fA, f0, 0 );
      dual.shift( fB, f0, 1 );
      dual.shift( fC, f0, 2 );
      const Idx ifA = dual.idx(fA);
      const Idx ifB = dual.idx(fB);
      const Idx ifC = dual.idx(fC);

      double uA0B=loops.us.at(Corner{ifA,if0,ifB});
      double uB0C=loops.us.at(Corner{ifB,if0,ifC});
      double uC0A=loops.us.at(Corner{ifC,if0,ifA});

      if( std::abs(uB0C)<1.0e-15 ) uB0C = 1.0;
      const double tA = std::sqrt( uC0A*uA0B/uB0C );
      const double tB = std::sqrt( uA0B*uB0C/uC0A );
      const double tC = std::sqrt( uB0C*uC0A/uA0B );

      ts_aug.insert( { DirectedLink{if0,0}, tA } );
      ts_aug.insert( { DirectedLink{if0,1}, tB } );
      ts_aug.insert( { DirectedLink{if0,2}, tC } );
      // }
      // else{
      //   ts_aug.insert( { DirectedLink{if0,0}, 0.0 } );
      //   ts_aug.insert( { DirectedLink{if0,1}, 0.0 } );
      //   ts_aug.insert( { DirectedLink{if0,2}, 0.0 } );
      // }
    }
  }

  void FillBetas(){
    for(Idx il=0; il<dual.NDirectedLinks(); il++) {
      const DirectedLink link = dual.directedlinkidx2DirectedLink(il);
      const Idx if1 = link.first;
      const FaceCoords f1 = dual.idx2FaceCoords(if1);
      const int df12 = link.second;
      FaceCoords f2; dual.shift(f2, f1, df12);
      const Idx if2 = dual.idx(f2);
      const int df21 = dual.getDirection.at(Link{if2,if1});

      const double t12 = ts_aug[DirectedLink{if1,df12}];
      const double t21 = ts_aug[DirectedLink{if2,df21}];
      const double t = t12*t21;

      betas.insert( {DirectedLink{if1,df12}, std::atanh(t)} );
      betas.insert( {DirectedLink{if2,df21}, std::atanh(t)} );

      // const double b12 = std::atanh(t12);
      // const double b21 = std::atanh(t21);
      // const double asq = std::cosh(b12+b21)*std::cosh(b12-b21);
      // As.insert( {DirectedLink{if1,df12}, std::sqrt(asq)} );
      // As.insert( {DirectedLink{if2,df21}, std::sqrt(asq)} );
    }
  }

  double eval_loop( const Idx q ) {
    loops.set(q);

    double w_prod = 1.;
    for(const auto& loop : loops){
      double w = 1.;
      for(int i=0; i<(int)loop.size(); i++){
        const Idx if1 = loop[i];
        const Idx if2 = loop[(i+1)%loop.size()];
        const int df12 = dual.getDirection.at(Link{if1,if2});
        w *= std::tanh( betas.at( DirectedLink{if1,df12} ) );
      }
      w_prod *= w;
    }
    return w_prod;
  }

};


// template<Idx N>
// class IsingField{
// public:
//   const Ising<N>& ising;
//   const RefinedIcosahedronDual& dual;

//   std::vector<int> s;

//   IsingField( const Ising<N>& ising_ )
//     : ising(ising_)
//     , dual(ising.loop.D.dual)
//     , s(dual.NVertices())
//   {}

// };




template<Idx N, Idx N2>
class Spin {
public:
  using Config = std::bitset<N2>;

  Ising<N>& ising;
  const RefinedIcosahedronDual& dual;

private:
  std::vector<bool> s;

public:
  std::mt19937 gen;
  std::uniform_real_distribution<double> d01D; // (0.0, 1.0); // (1, 6);
  std::uniform_int_distribution<int> d01I; // (0, 1); // (1, 6);
  std::uniform_int_distribution<Idx> d0N; // (0, Lx*Ly-1); // (1, 6);
  double dist01(){ return d01D(gen); }
  Idx dist0N(){ return d0N(gen); }
  int dist01I(){ return d01I(gen); }

  int operator[](const Idx i) const {
    const int si = s[i];
    return 2*si-1;
  }


  Spin() = delete;

  Spin( Ising<N>& ising_, const int seed=0 )
    : ising(ising_)
    , dual(ising.loops.D.dual)
    , s( dual.NVertices() )
    , d01D(0.0, 1.0)
    , d01I(0, 1)
    , d0N(0, dual.NVertices()-1)
  {
    assert( N2==dual.NVertices() );
    std::mt19937 gen0( seed );
    gen.seed( gen0() );
    //
    set1();
  }

  void set1() {
    for(Idx i=0; i<s.size(); i++){
      s[i] = 1;
    }
  }

  void set0() {
    for(Idx i=0; i<s.size(); i++){
      s[i] = 0;
    }
  }

  void random() {
    for(Idx i=0; i<s.size(); i++){
      s[i] = dist01I();
    }
  }

  void set( const Idx q ){
    assert(0<=q && q<std::pow(2,N2));
    const Config c = Config(q);
    for(Idx i=0; i<s.size(); i++){
      s[i] = c[i];
    }
  }

  void flip(){
    for(Idx i=0; i<s.size(); i++){
      s[i] = !s[i];
    }
  }

  void set_from_loop( const Idx q ){
    assert(0<=q && q<std::pow(2,N));
    ising.loops.set(q);

    set0();
    for(const auto& loop : ising.loops){
      for(const Idx if1 : loop){
        s[if1] = true;
      }
    }
  }


  std::string print() const {
    std::stringstream ss;
    for(Idx i=0; i<s.size(); i++) ss << s[i];
    return ss.str();
  }

  double S() const {
    double res = 0.0;

    for(Idx il=0; il<2*dual.NLinks(); il++){
      const DirectedLink ell = dual.directedlinkidx2DirectedLink( il );
      const Idx if1 = ell.first;
      FaceCoords f = dual.idx2FaceCoords( ell.first );
      // if(f[3]!=XZ) continue;

      FaceCoords f2;
      dual.shift( f2, dual.idx2FaceCoords(ell.first), ell.second );
      const Idx if2 = dual.idx(f2);
      if(if1>=if2) continue;

      const double beta = ising.betas.at( ell );
      const int prod = (*this)[ell.first]*(*this)[if2];
      assert( prod*prod == 1 );
      res +=  beta * prod;
    }
    return res;
  }


  double weight() const {
    double res = 1.0;

    for(Idx il=0; il<2*dual.NLinks(); il++){
      const DirectedLink ell = dual.directedlinkidx2DirectedLink( il );
      const Idx if1 = ell.first;
      FaceCoords f = dual.idx2FaceCoords( if1 );
      // if(f[3]!=XZ) continue;

      FaceCoords f2;
      dual.shift( f2, f, ell.second );
      const Idx if2 = dual.idx(f2);
      if(if1>=if2) continue;

      const double beta = ising.betas.at( ell );
      const int prod = (*this)[ell.first]*(*this)[if2];
      assert( prod*prod == 1 );
      // std::cout << "debug. " << ell.first << " " << if2 << " tanh = " << std::tanh(beta) << std::endl;
      res *=  1.0 + std::tanh(beta) * prod;
    }
    return res;
  }



  void ckpoint( const std::string& str ) const {
    std::ofstream of( str, std::ios::out | std::ios::binary | std::ios::trunc);
    if(!of) assert(false);

    // bool tmp = 0;
    // for(Idx i=0; i<Lx*Ly; i++){
    for(Idx i=0; i<s.size(); i++){
      of.write( (char*) &s+i, sizeof(bool) );
    }
    of.close();
  }

  void read( const std::string& str ) {
    std::ifstream ifs( str, std::ios::in | std::ios::binary );
    if(!ifs) assert(false);

    // bool tmp;
    // for(Idx i=0; i<Lx*Ly; ++i){
    for(Idx i=0; i<s.size(); i++){
      ifs.read((char*) &s+i, sizeof(bool) );
      // (*this)[i] = tmp;
    }
    ifs.close();
  }



};


// void heatbath( Spin& s ){
//   // omp not allowed
//   for(Idx i=0; i<Lx*Ly; i++){
//     if( !is_site(i) ) continue;

//     double henv = 0;
//     for(int mu=0; mu<SIX; mu++){
//       if( !is_link(i,mu) ) continue;
//       Idx j;
//       cshift( j, i, mu );
//       henv += Beta[mu%THREE]*s[j];
//     }

//     const double p = std::exp(2.0*henv);
//     const double r = dist01();
//     if( r<p/(1.0+p) ) s[i] = 1;
//     else s[i] = -1;
//   }
// }


// void wolff( Spin& s ){
//   std::vector<int> is_cluster(Lx*Ly, 0);
//   std::stack<Idx> stack_idx;

//   Idx init = dist0N();
//   while( !is_site(init) ) init = dist0N();

//   is_cluster[init] = 1;
//   stack_idx.push(init);

//   while( stack_idx.size() != 0 ){

//     const Idx p = stack_idx.top();
//     stack_idx.pop();
//     s[p] = -s[p]; // flip when visited

//     for(int mu = 0; mu < SIX; mu++){
//       if( !is_link(p,mu) ) continue;
//       Idx q;
//       cshift(q, p, mu);
//       if( s[q] == s[p] || is_cluster[q]==1 ) continue; // s[x]*sR[y]<0 or y in c

//       const double r = dist01();
//       if( r < std::exp(-2.0 * Beta[mu%THREE]) ) continue; // reject

//       is_cluster[q] = 1;
//       stack_idx.push(q);
//     }
//   }
// }


