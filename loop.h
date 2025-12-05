#pragma once


#include <bitset>
#include <set>


template<Idx N>
class DualLoops {
public:
  using Config = std::bitset<N>;
  using Loop = std::vector<Idx>;

  const RefinedIcosahedronDual& dual;
  const RefinedIcosahedron& lattice;

  Idx q;
  Config c;

  std::vector<Loop> loops;
  // auto begin(){ return loops.begin(); }
  // auto end(){ return loops.end(); }
  // auto begin() const { return loops.begin(); }
  // auto end() const { return loops.end(); }

  DualLoops
  (
   const RefinedIcosahedronDual& dual_
   )
    : dual(dual_)
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
    for(auto loop : loops){
      for(auto elem : loop) ss << elem << " ";
      ss << std::endl;
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

    std::set<Link> links_in_loops;
    for(const Idx in : ups){
      const Coords n = lattice.idx2Coords(in);

      for(int mu=0; mu<6; mu++){
        Coords np;
        lattice.shift( np, n, mu );
        if(np[2]==PatchIdxVoid) continue;
        const Idx inp = lattice.idx(np);

        if( !ups.contains(inp) ){ // (n, np) : broken
          const DirectedLink duallink = dual.linkSimp2Dual.at( DirectedLink{in,mu} );
          const FaceCoords fA = dual.idx2FaceCoords( duallink.first );
          const int df = duallink.second;
          FaceCoords fB;
          dual.shift( fB, fA, df );

          const Idx ifA = dual.idx(fA);
          const Idx ifB = dual.idx(fB);
          links_in_loops.insert( Link{ifA, ifB} );
        }
      }
    }

    // group links
    while(links_in_loops.size()!=0){
      // initialize a loop
      Loop loop;
      const Link ell = *links_in_loops.begin();
      const Idx if_end = ell.first, if_init = ell.second;

      // retrieve a loop starting from if_init; coming back to if_end
      Idx if_current = if_init;
      loop.push_back(if_current);

      while(*std::prev(loop.end())!=if_end){
        for(auto em=std::next(links_in_loops.begin()); em!=links_in_loops.end(); em++){
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
    }
  }
};
