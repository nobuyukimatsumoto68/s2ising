#pragma once

template<Idx N>
class Ising{
public:
  DualLoop<N>& loops;
  const Fermion& D;
  const DualSpinStructure& spin;
  const RefinedIcosahedronDual& dual;
  const RefinedIcosahedron& lattice;

  std::map<const DirectedLink, const double> ts_aug;
  std::map<const Link, const double> betas;
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
      const FaceCoords f0 = dual.idx2FaceCoords( if0 );
      FaceCoords fA, fB, fC;
      dual.shift( fA, f0, 0 );
      dual.shift( fB, f0, 1 );
      dual.shift( fC, f0, 2 );
      const Idx ifA = dual.idx(fA);
      const Idx ifB = dual.idx(fB);
      const Idx ifC = dual.idx(fC);

      const double uA0B=loops.us.at(Corner{ifA,if0,ifB});
      const double uB0C=loops.us.at(Corner{ifB,if0,ifC});
      const double uC0A=loops.us.at(Corner{ifC,if0,ifA});

      const double tA = std::sqrt( uC0A*uA0B/uB0C );
      const double tB = std::sqrt( uA0B*uB0C/uC0A );
      const double tC = std::sqrt( uB0C*uC0A/uA0B );

      ts_aug.insert( { DirectedLink{if0,0}, tA } );
      ts_aug.insert( { DirectedLink{if0,1}, tB } );
      ts_aug.insert( { DirectedLink{if0,2}, tC } );
    }
  }

  void FillBetas(){
    for(Idx il=0; il<dual.NDirectedLinks(); il++) {
      const DirectedLink link = dual.linkidx2DirectedLink(il);
      const Idx if1 = link.first;
      const FaceCoords f1 = dual.idx2FaceCoords(if1);
      const int df12 = link.second;
      FaceCoords f2; dual.shift(f2, f1, df12);
      const Idx if2 = dual.idx(f2);
      const int df21 = dual.getDirection.at(Link{if2,if1});

      const double t12 = ts_aug[DirectedLink{if1,df12}];
      const double t21 = ts_aug[DirectedLink{if2,df21}];
      const double t = t12*t21; // @@@ WHY SQRT?
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
