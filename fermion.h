#pragma once

template<class Lattice>
struct FermionVector {
  Lattice& lattice;
  std::vector<Complex> v;

  Fermion(const Lattice& lattice_)
    : lattice(lattice_)
    , v(2*lattice.NVertices(), 0.0)
  {}

  Idx idx( const Idx i, const int spin ) const {
    assert( 0<=spin && spin<=1);
    return 2*i+spin;
  }

  template<class CoordType>
  Idx idx( const CoordType x, const int spin ) const {
    return idx( lattice.idx(x), spin );
  }

  Complex operator[](const Idx i) const { return v[i]; }
  Complex& operator[](const Idx i) { return v[i]; }

  void insertMultBlock( const Idx i,
                        const SpinMatrix& mat,
                        const FermionVector& x,
                        const Idx j){
    Complex& y0 = this->idx(i,0);
    Complex& y1 = this->idx(i,1);

    const Complex x0 = x.idx(j,0);
    const Complex x1 = x.idx(j,1);

    y0 = mat(0,0)*x0 + mat(0,1)*x1;
    y1 = mat(1,0)*x0 + mat(1,1)*x1;
  }

};




template<class V>
struct Fermion {
  const DualSpinStructure& spin;
  const RefinedIcosahedronDual& dual;

  Fermion( const DualSpinStructure& spin_ )
    : spin(spin_)
    , dual(spin.dual)
  {}

  V operator()( const V& v ) const {
    V res(dual);

    for(Idx i=0; i<dual.NVertices(); i++){
      const FaceCoord f = dual.idx2FaceCoords(i);
      for(int dmu=; dmu<3; dmu++){
        const int mu = 3*(f[3]==ZY) + dmu;
        FaceCoord fp;
        dual.shift(fp,f,mu);
        const Idx j=dual.idx(fp);

        const SpinMatrix mat = - kappa * spin.P(f,df);
        res.insertMultBlock(i, mat, v, j);
      }
      mu@@;
    }

    return res;
  }

  set_kappa;

  set_mu;


};
