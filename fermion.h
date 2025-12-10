#pragma once

template<class Lattice>
struct FermionVector {
  const Lattice& lattice;
  std::vector<Complex> v;

  FermionVector(const Lattice& lattice_)
    : lattice(lattice_)
    , v(2*lattice.NVertices(), 0.0)
  {}

  FermionVector operator=(const FermionVector& other) {
    if (this == &other) return *this;
    assert(&(this->lattice) == &(other.lattice));
    this->v = other.v;
    return *this;
  }

  std::string print() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(25);
    for(const Complex c : v) ss << c.real() << " " << c.imag() << std::endl;
    return ss.str();
  }

  Complex* data() { return v.data(); }

  Idx idx( const Idx i, const int spin ) const {
    assert( 0<=spin && spin<=1);
    return 2*i+spin;
  }

  template<class CoordType>
  Idx idx2( const CoordType x, const int spin ) const {
    return idx( lattice.idx(x), spin );
  }

  Complex operator[](const Idx i) const { return v[i]; }
  Complex& operator[](const Idx i) { return v[i]; }

  Complex operator()( const Idx i, const int spin ) const {
    return v[idx(i,spin)];
  }

  Complex& operator()( const Idx i, const int spin ) {
    return v[idx(i,spin)];
  }

  void set_zero() { for(Complex& elem : v) elem = 0.0; }


  void addMultBlock( const Idx i,
                        const SpinMatrix& mat,
                        const FermionVector& x,
                        const Idx j){
    (*this)(i,0) += mat(0,0)*x(j,0) + mat(0,1)*x(j,1);
    (*this)(i,1) += mat(1,0)*x(j,0) + mat(1,1)*x(j,1);
  }

};




struct Fermion {
  const DualSpinStructure& spin;
  const RefinedIcosahedronDual& dual;


  std::vector<double> kappas; // dual.linkidx
  std::vector<double> mus; // dual.idx


  Fermion( const DualSpinStructure& spin_ )
    : spin(spin_)
    , dual(spin.dual)
  {
    set_kappas();
    set_mus();
  }

  using V = FermionVector<RefinedIcosahedronDual>;
  void operator()( V& res, const V& v ) const {
    res.set_zero();
    for(Idx i=0; i<dual.NVertices(); i++){
      const FaceCoords f = dual.idx2FaceCoords(i);
      for(int df=0; df<3; df++){
        FaceCoords fp;
        dual.shift(fp,f,df);
        const Idx j=dual.idx(fp);

        const SpinMatrix mat = - kappas[dual.linkidx(i,df)] * spin.P(f,df);
        res.addMultBlock(i, mat, v, j);
      }
      const SpinMatrix mat = 0.5 * mus[i] * spin.sigma[0];
      res.addMultBlock(i, mat, v, i);
    }
  }

  Eigen::MatrixXcd matrix() const {
    const Idx size = 2*dual.NVertices();
    Eigen::MatrixXcd res(size, size);

    FermionVector phi(dual), Dphi(dual);
    for(Idx i=0; i<dual.NVertices(); i++){
      for(int s=0; s<2; s++){
        phi.set_zero(); Dphi.set_zero();
        phi( i, s ) = 1.0;
        (*this)(Dphi, phi);

        res.block( 0, phi.idx(i,s), size, 1 ) = Eigen::Map<Eigen::MatrixXcd>( Dphi.data(), size, 1 );
      }
    }
    return res;
  }


  void set_kappas(){
    kappas.clear();
    kappas.resize( dual.NDirectedLinks() );
    for(Idx il=0; il<dual.NDirectedLinks(); il++) {
      kappas[il] = dual.link_volumes[il] / dual.ells[il];
    }
  }

  void set_mus(){
    mus.clear();
    mus.resize( dual.NVertices() );

    for(Idx i=0; i<dual.NVertices(); i++){
      for(int df=0; df<3; df++) mus[i] += kappas[dual.linkidx(i,df)];
    }
  }


};
