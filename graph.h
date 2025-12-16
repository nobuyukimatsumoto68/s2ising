#pragma once


#include <bitset>
#include <set>


using Corner = std::array<Idx, 3>;


template<Idx N0,Idx N1,Idx N2>
class Graph{
public:
  using SimpConfig = std::bitset<N0>;
  using LinkConfig = std::bitset<N1>;
  using DualConfig = std::bitset<N2>;

  // using LinkGraph = std::vector<Idx>;
  const Fermion& D;
  const DualSpinStructure& spin;
  const RefinedIcosahedronDual& dual;
  const RefinedIcosahedron& lattice;

  SimpConfig s_sites;
  LinkConfig s_links;
  DualConfig s_faces;

  Graph
  (
   const Fermion& D_
   )
    : D(D_)
    , spin(D.spin)
    , dual(D.dual)
    , lattice(dual.simplicial)
  {
  }

  // Idx operator()() const { return q; }

  bool site(const Idx i) const {
    assert(i<N0);
    return s_sites[i];
  }
  bool link(const Idx i) const {
    assert(i<N1);
    return s_links[i];
  }
  bool face(const Idx i) const {
    assert(i<N2);
    return s_faces[i];
  }

  // Config config() const { return c; }
  // std::string printLoops() const {
  //   std::stringstream ss;
  //   bool is_first = true;
  //   for(auto loop : loops){
  //     if(is_first) is_first = false;
  //     else ss << std::endl;
  //     for(auto elem : loop) ss << elem << " ";
  //   }
  //   return ss.str();
  // }

  void set_from_siteconfig( const Idx q ){
    s_sites = SimpConfig(q);
    s_links = LinkConfig(0);
    s_faces = FaceConfig(0);

    for(Idx il=0; il<dual.NLinks(); il++){
      const DirectedLink ell = dual.link2DirectedLink(il);
    }
  }
  void set_from_linkconfig( const Idx q ){};
  void set_from_faceconfig( const Idx q ){};



  double eval_prod(){
    Complex w_prod = 1.;
    for(const auto& loop : loops){
      double factor = 1.;
      SpinMatrix mat = SpinMatrix::Identity();
      for(int i=0; i<(int)loop.size(); i++){
        const Idx if1 = loop[i];
        const Idx if2 = loop[(i+1)%loop.size()];

        const FaceCoords f1 = dual.idx2FaceCoords(if1);
        const int df2 = dual.getDirection.at(Link{if1,if2});

        factor *= D.kappas[dual.directedlinkidx(if1,df2)];
        factor /= D.mus[if1];
        mat = mat * spin.P( f1, df2 );
      }
      const Complex w = factor*mat.trace();
      w_prod *= w;
    }
    assert(std::abs(w_prod.imag())<1.0e-13);
    return std::pow(-1,loops.size()) * w_prod.real();
  }

  // define u012 and rewrite
  double eval(){
    double w_prod = 1.;
    for(const auto& loop : loops){
      double w = 1.;
      for(int i=0; i<(int)loop.size(); i++){
        const Idx if0 = loop[(i+loop.size()-1)%loop.size()];
        const Idx if1 = loop[i];
        const Idx if2 = loop[(i+1)%loop.size()];
        w *= us.at(Corner{if0,if1,if2});
      }
      w_prod *= w;
    }
    return w_prod;
  }

  // define u012 and rewrite
  double eval_(){
    double w_prod = 1.;
    for(const auto& loop : loops){
      double w = 1.;
      for(int i=0; i<(int)loop.size(); i++){
        const Idx if0 = loop[(i+loop.size()-1)%loop.size()];
        const Idx if1 = loop[i];
        const Idx if2 = loop[(i+1)%loop.size()];

        const int df2 = dual.getDirection.at(Link{if1,if2});
        const int df0 = dual.getDirection.at(Link{if1,if0});

        double al10 = spin.alphas[dual.directedlinkidx(if1, df0)];
        double al12 = spin.alphas[dual.directedlinkidx(if1, df2)];
        double dalpha = (al10+M_PI)-al12;
        w *= D.kappas[dual.directedlinkidx(if1,df2)];
        w /= D.mus[if1];
        while(dalpha>M_PI) dalpha -= 2.0*M_PI;
        while(dalpha<-M_PI) dalpha += 2.0*M_PI;
        w *= std::cos( dalpha/2 );
      }
      // w *= -1.0;
      w_prod *= w;
    }
    return w_prod;
  }

};
