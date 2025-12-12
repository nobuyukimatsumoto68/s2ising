#pragma once


#include <bitset>
#include <set>


template<Idx N>
class DualLoops {
public:
  using Config = std::bitset<N>;
  using Loop = std::vector<Idx>;

  const Fermion& D;
  const DualSpinStructure& spin;
  const RefinedIcosahedronDual& dual;
  const RefinedIcosahedron& lattice;

  Idx q;
  Config c;

  std::vector<Loop> loops;
  auto begin(){ return loops.begin(); }
  auto end(){ return loops.end(); }
  auto begin() const { return loops.begin(); }
  auto end() const { return loops.end(); }

  DualLoops
  (
   const Fermion& D_
   )
    : D(D_)
    , spin(D.spin)
    , dual(D.dual)
    , lattice(dual.simplicial)
  {}

  Idx operator()() const { return q; }

  bool operator[](const Idx i) const {
    assert(i<N);
    return c[i];
  }

  Config config() const { return c; }

  std::string printLoops() const {
    std::stringstream ss;
    bool is_first = true;
    for(auto loop : loops){
      if(is_first) is_first = false;
      else ss << std::endl;
      for(auto elem : loop) ss << elem << " ";
    }
    return ss.str();
  }


  void set( const Idx q_ ){
    q = q_;
    //    c.
    c = Config( q );

    loops.clear();

    // 1:up, 0:down
    std::set<Idx> ups;
    for(Idx in=0; in<N; in++) if(c[in]==1) ups.insert(in);

    // for(auto elem : ups) std::cout << "debug. elem = " << elem << std::endl;

    std::set<Link> links_in_loops;
    for(const Idx in : ups){
      const Coords n = lattice.idx2Coords(in);

      for(int mu=0; mu<6; mu++){
        Coords np;
        lattice.shift( np, n, mu );
        if(np[2]==PatchIdxVoid) continue;
        const Idx inp = lattice.idx(np);

        if( !ups.contains(inp) ){ // (n, np) : broken
          // std::cout << "debug. broken. in: " << in << " inp: " << inp << std::endl;
          const DirectedLink duallink = dual.linkSimp2Dual.at( DirectedLink{in,mu} );
          const FaceCoords fA = dual.idx2FaceCoords( duallink.first );
          const int df = duallink.second;
          // std::cout << "debug. broken. dual: "
          //           << fA[0] << " " << fA[1] << " " << fA[2] << " " << fA[3] << " / "
          //           << duallink.second << std::endl;
          FaceCoords fB;
          dual.shift( fB, fA, df );

          const Idx ifA = dual.idx(fA);
          const Idx ifB = dual.idx(fB);
          links_in_loops.insert( Link{ifA, ifB} );
        }
      }
    }

    // for(auto em : links_in_loops) std::cout << "debug. em = " << em.first << " " << em.second << std::endl;
    // return;

    // group links
    while(links_in_loops.size()!=0){
      // initialize a loop
      Loop loop;
      const Link ell = *links_in_loops.begin();
      const Idx if_end = ell.first, if_init = ell.second;

      // retrieve a loop starting from if_init; coming back to if_end
      Idx if_current = if_init;
      loop.push_back(if_current);

      // std::cout << "debug. pt1" << std::endl;

      while(*std::prev(loop.end())!=if_end){
        for(auto em=std::next(links_in_loops.begin()); em!=links_in_loops.end(); em++){
          // std::cout << "debug. em = " << em->first << " " << em->second << std::endl;
          if(em->first==if_current){
            if_current = em->second;
            loop.push_back(if_current);
            links_in_loops.erase(em);
            break;
          }
          else if(em->second==if_current){
            if_current = em->first;
            loop.push_back(if_current);
            links_in_loops.erase(em);
            break;
          }
        } // end for
      }
      links_in_loops.erase(ell);
      loops.push_back(loop);
      // std::cout << "debug. loops.size() == " << loops.size() << std::endl;
    }
  }

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

        factor *= D.kappas[dual.linkidx(if1,df2)];
        factor /= D.mus[if1];
        mat = mat * spin.P( f1, df2 );
      }
      const Complex w = factor*mat.trace();
      w_prod *= w;
    }
    assert(std::abs(w_prod.imag())<1.0e-13);
    return w_prod.real();
  }

  double eval(){
    double w_prod = 1.;
    for(const auto& loop : loops){
      double w = 1.;
      for(int i=0; i<(int)loop.size(); i++){
        const Idx if0 = loop[(i+loop.size()-1)%loop.size()];
        const Idx if1 = loop[i];
        const Idx if2 = loop[(i+1)%loop.size()];

        const int df2 = dual.getDirection.at(Link{if1,if2});
        const int df0 = dual.getDirection.at(Link{if1,if0});

        double al10 = spin.alphas[dual.linkidx(if1, df0)];
        double al12 = spin.alphas[dual.linkidx(if1, df2)];
        double dalpha = (al10+M_PI)-al12;
        w *= D.kappas[dual.linkidx(if1,df2)];
        w /= D.mus[if1];
        while(dalpha>M_PI) dalpha -= 2.0*M_PI;
        while(dalpha<-M_PI) dalpha += 2.0*M_PI;
        w *= std::cos( dalpha/2 );
      }
      w *= -1.0;
      w_prod *= w;
    }
    return w_prod;
  }

};
