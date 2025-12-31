#pragma once


#include <functional>


template<typename T1, typename T2>
class Jackknife {
public:
  // std::vector<T2> config;
  std::unique_ptr<T2[]> config;

  // std::vector<T1> bin_avg;
  std::vector<T1> jack_avg;
  T1 mean;
  T1 var;

  Idx N;
  Idx nbins;
  Idx binsize;

  // Jackknife()
  // {
  //   config.clear();
  // }

  Jackknife( const Idx N_ )
    : N(N_)
      // , config(N)
  {
    config = std::make_unique<T2[]>(N);
  }


  // Idx size() const { return config.size(); }
  Idx size() const { return N; }

  // void meas( const T2& w ) { config.push_back(w); }
  void meas( const Idx i, const T2& w ) { config[i] = w; }

  T1 jk_avg( const int i,
             const std::function<T1(const std::vector<T2>&)> f ) const {
    // std::vector<T2> vr = config;
    // vr.erase( vr.begin()+N, vr.end() );
    std::vector<T2> vr (&config[0], &config[0]+N );
    vr.erase( vr.begin()+i*binsize, vr.begin()+(i+1)*binsize );
    assert( vr.size()==(nbins-1)*binsize );
    return f( vr );
  }

  void init( const int binsize_ ){
    this->binsize = binsize_;
    this->nbins = N/binsize;
    this->N = binsize*nbins;
    jack_avg.clear();
    jack_avg.resize(nbins);
  }

  void init( const int binsize_, const int nbins_ ){
    this->binsize = binsize_;
    this->nbins = nbins_;
    this->N = binsize*nbins;
    jack_avg.clear();
    jack_avg.resize(nbins);
  }


  void finalize( const std::function<T1(const T1&)> square,
                 const T1& zero ){
    assert(jack_avg.size()==nbins);

    this->mean = zero;
    this->var = zero;

    for(int i=0; i<nbins; i++) mean += jack_avg[i];
    mean /= nbins;
    for(int i=0; i<nbins; i++) var += square(jack_avg[i] - mean);
    var *= 1.0*(nbins-1)/nbins;
  }

  void do_all( const std::function<T1(const std::vector<T2>&)> f,
               const std::function<T1(const T1&)> square,
               const T1& zero,
               const int binsize_ ) {
    init( binsize_ );
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(nparallel)
// #endif
    // jack_avg.clear();
    // jack_avg.resize(nbins);
    for(int i=0; i<nbins; i++) jack_avg[i]= jk_avg( i, f );

    finalize(square, zero);
    // this->mean = zero;
    // this->var = zero;
    // for(int i=0; i<nbins; i++) mean += jack_avg[i];
    // mean /= nbins;
    // for(int i=0; i<nbins; i++) var += square(jack_avg[i] - mean);
    // var *= 1.0*(nbins-1)/nbins;
  }


  void write( const T1& dat, const std::string& filepath, const Idx size ) const {
    std::ofstream of( filepath, std::ios::out | std::ios::binary | std::ios::trunc);
    if(!of) assert(false);

    double tmp = 0;
    for(Idx i=0; i<size; i++){
      tmp = dat[i];
      of.write( (char*) &tmp, sizeof(double) );
    }
    of.close();
  }

  void read( T1& dat, const std::string& filepath, const Idx size ) {
    std::ifstream ifs( filepath, std::ios::in | std::ios::binary );
    if(!ifs) assert(false);

    double tmp;
    for(Idx i=0; i<size; i++){
      ifs.read((char*) &tmp, sizeof(double) );
      dat[i] = tmp;
    }
    ifs.close();
  }


};



template<typename T1, typename T2>
class JackknifeSimp {
public:
  // std::vector<T2> config;
  std::unique_ptr<T2[]> config;

  std::vector<T1> bin_avg;
  std::vector<T1> jack_avg;
  T1 mean;
  T1 var;

  Idx N;
  Idx nbins;
  Idx binsize;

  // JackknifeSimp()
  // {
  //   //config.clear();
  // }

  JackknifeSimp( const Idx N_ )
    : N(N_)
      // , config(N)
  {
    config = std::make_unique<T2[]>(N);
  }


  // Idx size() const { return config.size(); }
  Idx size() const { return N; }

  // void meas( const T2& w ) { config.push_back(w); }
  void meas( const Idx i, const T2& w ) { config[i] = w; }

  T1 get_bin_avg( const int i,
                  const std::function<T1(const std::vector<T2>&)>& mean ) const {
    // std::vector<T2> vr (config.begin()+i*binsize, config.begin()+(i+1)*binsize );
    std::vector<T2> vr (&config[0]+i*binsize, &config[0]+(i+1)*binsize );
    assert( vr.size()==binsize );
    return mean( vr );
  }


  T1 jk_avg( const int i,
             const std::function<T1(const std::vector<T1>&)>& mean ) const {
    std::vector<T1> vr = bin_avg;
    vr.erase( vr.begin()+i );
    assert( vr.size()==(nbins-1) );
    return mean( vr );
  }

  void init( const int binsize_ ){
    this->binsize = binsize_;
    // this->nbins = config.size()/binsize;
    this->nbins = N/binsize;
    this->N = binsize*nbins;
    jack_avg.clear();
    bin_avg.clear();
    bin_avg.resize(nbins);
    jack_avg.resize(nbins);
  }

  void init( const int binsize_, const int nbins_ ){
    this->binsize = binsize_;
    this->nbins = nbins_;
    this->N = binsize*nbins;
    jack_avg.clear();
    bin_avg.clear();
    bin_avg.resize(nbins);
    jack_avg.resize(nbins);
  }


  void finalize( const std::function<T1(const T1&)> square,
                 const T1& zero ){
    assert(jack_avg.size()==nbins);

    this->mean = zero;
    this->var = zero;

    for(int i=0; i<nbins; i++) mean += jack_avg[i];
    mean /= nbins;
    for(int i=0; i<nbins; i++) var += square(jack_avg[i] - mean);
    var *= 1.0*(nbins-1)/nbins;
  }

//   void do_all( const std::function<T1(const std::vector<T2>&)> f,
//                const std::function<T1(const T1&)> square,
//                const T1& zero,
//                const int binsize_ ) {
//     init( binsize_ );
// // #ifdef _OPENMP
// // #pragma omp parallel for num_threads(nparallel)
// // #endif
//     // jack_avg.clear();
//     // jack_avg.resize(nbins);
//     // for(int i=0; i<nbins; i++) jack_avg[i]= jk_avg( i, f );

//     finalize(square, zero);
//     // this->mean = zero;
//     // this->var = zero;
//     // for(int i=0; i<nbins; i++) mean += jack_avg[i];
//     // mean /= nbins;
//     // for(int i=0; i<nbins; i++) var += square(jack_avg[i] - mean);
//     // var *= 1.0*(nbins-1)/nbins;
//   }


  void write( const T1& dat, const std::string& filepath, const Idx size ) const {
    std::ofstream of( filepath, std::ios::out | std::ios::binary | std::ios::trunc);
    if(!of) assert(false);

    double tmp = 0;
    for(Idx i=0; i<size; i++){
      tmp = dat[i];
      of.write( (char*) &tmp, sizeof(double) );
    }
    of.close();
  }

  void read( T1& dat, const std::string& filepath, const Idx size ) {
    std::ifstream ifs( filepath, std::ios::in | std::ios::binary );
    if(!ifs) assert(false);

    double tmp;
    for(Idx i=0; i<size; i++){
      ifs.read((char*) &tmp, sizeof(double) );
      dat[i] = tmp;
    }
    ifs.close();
  }


};
