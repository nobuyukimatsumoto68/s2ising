#pragma once



constexpr int NIh = 120;



struct Rotation {
  std::array<M3, 3> J;

  Rotation(){
    // eps_{1ij}
    J[0] << 0., 0., 0.,
      0., 0., 1.,
      0., -1., 0.;
    // eps_{2ij}
    J[1] << 0., 0., -1.,
      0., 0., 0.,
      1., 0., 0.;
    // eps_{3ij}
    J[2] << 0., 1., 0.,
      -1., 0., 0.,
      0., 0., 0.;
  }

  M3 operator()( const V3& theta ) const {
    M3 expn = M3::Zero();
    for(int i=0; i<3; i++){
      expn -= theta[i]*J[i];
    }
    return expn.exp();
  }

};





struct FullIcosahedralGroup {

  using MultTable = std::vector<std::vector<int>>;
  using Rep = std::vector<Eigen::MatrixXcd>;

  MultTable tab;
  const int is;
  const int it;
  const int iu;

  std::vector<std::string> words;

  Rep reg;

  FullIcosahedralGroup( const std::string& tablefile,
                        const int is_,
                        const int it_,
                        const int iu_
                        )
    : is(is_)
    , it(it_)
    , iu(iu_)
  {
    ReadMultTable( tablefile );
    GetRegRepFromMultTable();

    const int ist = tab[is][it];

    {
      auto check = reg[is]*reg[is];
      double norm = ( check-Eigen::MatrixXcd::Identity(NIh,NIh) ).norm();
      assert( norm<1.0e-12 );
    }
    {
      auto check = reg[it]*reg[it]*reg[it];
      double norm = ( check-Eigen::MatrixXcd::Identity(NIh,NIh) ).norm();
      assert( norm<1.0e-12 );
    }
    {
      auto check = reg[ist]*reg[ist]*reg[ist]*reg[ist]*reg[ist];
      double norm = ( check-Eigen::MatrixXcd::Identity(NIh,NIh) ).norm();
      assert( norm<1.0e-12 );
    }
    {
      auto check = reg[iu]*reg[is]*reg[iu];
      double norm = ( check-reg[is] ).norm();
      assert( norm<1.0e-12 );
    }
    {
      auto check = reg[iu]*reg[it]*reg[iu];
      double norm = ( check-reg[it] ).norm();
      assert( norm<1.0e-12 );
    }

    GetIndependentWords();
  }


  void ReadMultTable( const std::string& filename ){
    tab.clear();

    // std::cout << "# reading simplicial points" << std::endl;
    std::ifstream file(filename);

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      std::vector<int> row;
      for(int i=0; i<NIh; i++){
        int j; iss >> j;
        row.push_back(j);
      }
      tab.push_back( row );
    }
    assert( tab.size()==NIh );
  }


  void GetRegRepFromMultTable(){
    reg.clear();
    reg.resize(NIh);
    for(int g=0; g<NIh; g++){
      reg[g] = Eigen::MatrixXcd::Zero( NIh, NIh );
      for(int i=0; i<NIh; i++){
        reg[g]( tab[g][i], i ) = 1.0;
      }
    }
  }


  void GenerateRep(Rep& rep,
                   const Eigen::MatrixXcd& ms,
                   const Eigen::MatrixXcd& mt
                   ){
    const int nd=ms.cols();

    rep.clear();
    rep.resize(NIh);

    const auto mu = Eigen::MatrixXcd::Identity( nd, nd );

    for(int ig=0; ig<NIh; ig++){
      rep[ig] = Eigen::MatrixXcd::Identity( nd, nd );
      for(const char y : words[ig]) rep[ig] = rep[ig] * replace( y, ms, mt, mu );
    }

    // check
    for(int i=0; i<NIh; i++){
      for(int g=0; g<NIh; g++){
        assert( (rep[ tab[g][i] ] - rep[g] * rep[i]).norm()<1.0e-12 );
      }}
  }


  int replace( const char x ) const {
    if(x=='s') return is;
    else if(x=='t') return it;
    else if(x=='u') return iu;
    else {
      assert(false);
      return 0;
    }
  }

  Eigen::MatrixXcd replace( const char x,
              const Eigen::MatrixXcd& ms,
              const Eigen::MatrixXcd& mt,
              const Eigen::MatrixXcd& mu ) const {
    if(x=='s') return ms;
    else if(x=='t') return mt;
    else if(x=='u') return mu;
    else {
      assert(false);
      return Eigen::MatrixXcd::Zero(ms.rows(),ms.cols());
    }
  }

  M3 replace( const char x,
              const M3& ms,
              const M3& mt,
              const M3& mu ) const {
    if(x=='s') return ms;
    else if(x=='t') return mt;
    else if(x=='u') return mu;
    else {
      assert(false);
      return M3::Zero();
    }
  }

  void GetIndependentWords(){
    std::vector<std::string> words_tmp;
    std::vector<int> elems;

    words_tmp.push_back("");
    elems.push_back(0);

    words_tmp.push_back("u");
    elems.push_back(iu);

    int counter=0;
    bool flag=true;

    while(flag){
      bool flag2=false;

      for( const char x : {'s', 't'} ){
        std::string word = words_tmp[counter] + x;
        std::vector<int> trsl; // (word.size());
        for(const char y : word) trsl.push_back( replace(y) );

        int elem = 0;
        for(int i=0; i<trsl.size(); i++){
          elem = tab[elem][trsl[i]];
        }

        if(std::find(elems.begin(), elems.end(), elem) == elems.end()) {
          words_tmp.push_back(word);
          elems.push_back(elem);
          flag2 = flag2 || true;
        }
      } // end for x

      flag = flag2 || (counter!=words_tmp.size()-1);
      counter++;
    }

    words.clear();
    words.resize(words_tmp.size());

    for(int i=0; i<elems.size(); i++){
      words[elems[i]] = words_tmp[i];
    }
  }

  std::string print() const {
    std::stringstream ss;
    Idx counter=0;
    for(auto word : words){
      ss << counter << ": " << word << std::endl;
      counter++;
    }
    return ss.str();
  }

  M3 rotation( const int i,
               const M3& ms,
               const M3& mt ) const {
    assert( i<NIh );
    const M3 mu = -M3::Identity();

    M3 mat = M3::Identity();
    for(const char y : words[i]) mat = mat*replace(y, ms, mt, mu);

    return mat;
  }

};




struct Orbits {
  Vertices& vertices; // reference
  const Icosahedron& icos;
  const FullIcosahedralGroup& Ih;
  const Rotation& rot;
  const Idx nVertices;

  M3 ms, mt;
  std::vector<std::vector<Idx>> orbits;

  Orbits(Vertices& vertices_,
         const Icosahedron& icos_,
         const FullIcosahedralGroup& Ih_,
         const Rotation& rot_,
         const std::array<int, 2>& link=std::array<int, 2>{0,1},
         const int ipatch=0 )
    : vertices(vertices_)
    , icos(icos_)
    , Ih(Ih_)
    , rot(rot_)
    , nVertices(vertices.size())
  {
    {
      // s
      const V3 p0 = icos[link[0]];
      const V3 p1 = icos[link[1]];
      V3 theta = project( p0+p1 );
      theta *= M_PI;
      ms = rot(theta);
      assert( (ms*ms-M3::Identity()).norm()<1.0e-13 );
    }
    {
      // t
      const V3 pa = icos[icos.patches[ipatch][0]];
      const V3 pb = icos[icos.patches[ipatch][1]];
      const V3 pd = icos[icos.patches[ipatch][3]];
      V3 theta = project( pa+pb+pd );
      theta *= 2.0*M_PI/3.0;

      mt = rot(theta);
      assert( (mt*mt*mt-M3::Identity()).norm()<1.0e-13 );
    }
    {
      const M3 mst = ms*mt;
      assert( (mst*mst*mst*mst*mst-M3::Identity()).norm()<1.0e-13 );
    }

    GetOrbits();
  }

  bool is_identical( const V3& x, const V3& y, const double tol=1.0e-12 ) const {
    const double norm = (x-y).norm();
    return norm < tol;
  }

  void GetOrbits(){
    orbits.clear();
    for(const V3& x : vertices){
      std::vector<Idx> orbit;

      for(int ig=0; ig<NIh; ig++){
        const M3 mg = Ih.rotation(ig, ms, mt);
        const V3 gx = mg*x;

        for(int iy=0; iy<nVertices; iy++){
          const V3 y = vertices[iy];
          if( is_identical(gx, y) ) {
            orbit.push_back(iy);
            break;
          }
        }
      }

      // std::cout << "debug. x = " << x.transpose() << std::endl;
      // std::cout << "debug. orbit.size() = " << orbit.size() << std::endl;
      assert( orbit.size()==NIh );
      orbits.push_back(orbit);
    }
  }

  void RefreshOrbits( const std::vector<Idx>& basepts ){
    for(const Idx in : basepts){
      const V3 r = vertices[in];
      const auto orbit = orbits[in];

      for(int ig=0; ig<NIh; ig++){
        const M3 mg = Ih.rotation(ig, ms, mt);
        const V3 gr = mg*r;
        vertices[orbit[ig]] = gr;
      }
    }
  }

  std::string print() const {
    std::stringstream ss;
    Idx counter=0;
    for(const auto& orbit : orbits){
      ss << counter << ": ";
      for(const auto& ix : orbit){
        ss << ix << " ";
      }
      ss << std::endl;
      counter++;
    }
    return ss.str();
  }
};
