#pragma once


#include <functional>


template<typename T1, typename T2>
class Jackknife {
public:
  std::vector<T2> config;

  // std::vector<T1> bin_avg;
  std::vector<T1> jack_avg;
  T1 mean;
  T1 var;

  Idx N;
  Idx nbins;
  Idx binsize;

  Jackknife()
  {
    config.clear();
  }

  Jackknife( const Idx N_ )
    : N(N_)
    , config(N)
  {
  }

  Idx size() const { return config.size(); }

  void meas( const T2& w ) { config.push_back(w); }
  void meas( const Idx i, const T2& w ) { config[i] = w; }

  void do_it( const std::function<T1(const std::vector<T2>&)> f,
              const std::function<T1(const T1&)> square,
              const T1& zero,
              const int binsize_ ) {
    this->binsize = binsize_;
    this->nbins = config.size()/binsize;
    this->N = binsize*nbins;


// #ifdef _OPENMP
// #pragma omp parallel for num_threads(nparallel)
// #endif
    jack_avg.clear();
    jack_avg.resize(nbins);
    for(int i=0; i<nbins; i++){
      std::clog << "# debug. ibin = " << i << std::endl;
      std::vector<T2> vr = config;
      vr.erase( vr.begin()+N, vr.end() );
      vr.erase( vr.begin()+i*binsize, vr.begin()+(i+1)*binsize );
      assert( vr.size()==(nbins-1)*binsize );
      jack_avg[i]= f( vr );
    }

    this->mean = zero;
    this->var = zero;

    for(int i=0; i<nbins; i++) mean += jack_avg[i];
    mean /= nbins;
    for(int i=0; i<nbins; i++) var += square(jack_avg[i] - mean);
    var *= 1.0*(nbins-1)/nbins;
  }
};
