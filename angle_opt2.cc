#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>
#include <complex>

// #include <hdf5.h>
// #include "highfive/H5File.hpp"
// #include <highfive/H5File.hpp>
// #include <highfive/H5DataSet.hpp>
// #include <highfive/H5DataSpace.hpp>

#include <algorithm>
#include <filesystem>


#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>


using Idx=std::size_t;
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
using Complex = std::complex<double>;

#include "icosahedron.h"
#include "lattices.h"
#include "groups.h"
#include "dual_optimizer.h"


#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>





double my_f (const gsl_vector *v, void *params){
  // double x, y;
  DualLatticeAngleCostEvaluator *opt = (DualLatticeAngleCostEvaluator *)params;

  for(Idx i=0; i<opt->size(); i++) opt->p[i] = gsl_vector_get(v, i);
  opt->update();

  return opt->cost();
  // x = gsl_vector_get(v, 0);
  // y = gsl_vector_get(v, 1);
  // return p[2] * (x - p[0]) * (x - p[0]) + p[3] * (y - p[1]) * (y - p[1]) + p[4];
}



// int main(int argc, char* argv[]){
int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------
  // const int L = 16;
  // const int L = 32;
  const int L = 48;
  RefinedIcosahedron lattice(L);

  FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  Rotation rot;

  RefinedIcosahedronDual dual(lattice);
  Orbits orbits(dual.vertices,
                dual.basePoints,
                dual.baseTypes,
                lattice, Ih, rot);

  DualLatticeAngleCostEvaluator opt(dual, lattice, Ih, rot);

  // std::cout << "p:" << std::endl;
  // for(Idx i=0; i<opt.size(); i++) std::cout << std::setw(15) << opt.p[i];
  // std::cout << std::endl;
  // std::cout << "cost = " << opt.cost() << std::endl;

  // const double orig = opt.p[0];
  // const int nticks = 100;
  // const double eps = 1.8*std::abs(orig)/nticks;
  // for(int i=-nticks/2; i<=nticks/2; i++){
  //   opt.p[0] = orig + i*eps;
  //   opt.update();
  //   std::cout << "p = " << opt.p[0] << " cost = " << opt.cost() << std::endl;
  // }

  const int iter_max = 100*L*L;

  {
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    // const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *s = NULL;

    // variables
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point */
    x = gsl_vector_alloc ( opt.size() );

    for(Idx i=0; i<opt.size(); i++) gsl_vector_set (x, i, opt.p[i]);

    /* Set initial step sizes */
    ss = gsl_vector_alloc ( opt.size() );
    gsl_vector_set_all (ss, 0.01*lattice.len/L);

    /* Initialize method and iterate */
    minex_func.n = opt.size();
    minex_func.f = my_f;
    minex_func.params = &opt; // par;

    s = gsl_multimin_fminimizer_alloc (T, opt.size());
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do{
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status) break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-8);

      std::cout << "p:" << std::endl;
      for(Idx i=0; i<opt.size(); i++) std::cout << std::setw(15) << opt.p[i];
      std::cout << std::endl;

      if (status == GSL_SUCCESS) printf ("converged to minimum at\n");

      printf ("%15d f() = %7.15f size = %.15f\n",
              iter,
              // gsl_vector_get (s->x, 0),
              // gsl_vector_get (s->x, 1),
              s->fval, size);
    } while (status == GSL_CONTINUE && iter < iter_max);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
  }





  return 0;
}
