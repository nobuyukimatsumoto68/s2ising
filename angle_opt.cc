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

#include <omp.h>




double my_f (const gsl_vector *v, void *params){
  DualLatticeAngleCostEvaluator *opt = (DualLatticeAngleCostEvaluator *)params;

  for(Idx i=0; i<opt->size(); i++) opt->p[i] = gsl_vector_get(v, i);
  opt->update();

  return opt->cost();
}


/* The gradient of f, df = (df/dx, df/dy). */
void my_df (const gsl_vector *v, void *params,
            gsl_vector *df)
{
  const double eps = 1.0e-6;
  DualLatticeAngleCostEvaluator *opt = (DualLatticeAngleCostEvaluator *)params;

  for(Idx i=0; i<opt->size(); i++) opt->p[i] = gsl_vector_get(v, i);
  opt->update();

  for(Idx i=0; i<opt->size(); i++){
    opt->p[i] = gsl_vector_get(v, i) + eps;
    opt->update();
    double plus = opt->cost();
    opt->p[i] = gsl_vector_get(v, i) - eps;
    opt->update();
    double minus = opt->cost();
    opt->p[i] = gsl_vector_get(v, i);

    const double deriv = (plus-minus)/(2.0*eps);
    gsl_vector_set(df, i, deriv);
  }
}

/* Compute both f and df together. */
void my_fdf (const gsl_vector *x, void *params,
             double *f, gsl_vector *df)
{
  *f = my_f(x, params);
  my_df(x, params, df);
}


// int main(int argc, char* argv[]){
int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------
  // const int L = 16;
  // const int L = 32;
  const int L = 48;

  std::cout << "debug. pt1" << std::endl;
  RefinedIcosahedron lattice(L);

  std::cout << "debug. pt2" << std::endl;
  FullIcosahedralGroup Ih( "multtablemathematica.dat", 3, 19, 60 );
  Rotation rot;

  std::cout << "debug. pt3" << std::endl;
  RefinedIcosahedronDual dual(lattice);

  std::cout << "debug. pt4" << std::endl;
  Orbits orbits(dual.vertices,
                dual.basePoints,
                dual.baseTypes,
                lattice, Ih, rot);

  std::cout << "debug. pt5" << std::endl;
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
    size_t iter = 0;
    int status;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

    // variables
    gsl_vector *ss, *x;
    // gsl_multimin_function minex_func;
    gsl_multimin_function_fdf my_func;
    // gsl_multimin_fdfminimizer_vector_bfgs2 my_func;

    my_func.n = opt.size();
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = &opt;

    /* Starting point */
    x = gsl_vector_alloc ( opt.size() );
    for(Idx i=0; i<opt.size(); i++) gsl_vector_set (x, i, opt.p[i]);

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, opt.size());

    gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

    do
      {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status) break;

        status = gsl_multimin_test_gradient (s->gradient, 1e-3);

        if (status == GSL_SUCCESS) printf ("Minimum found at:\n");

        printf ("%5d %10.5f\n", iter,
                // gsl_vector_get (s->x, 0),
                // gsl_vector_get (s->x, 1),
                s->f);

      }
    while (status == GSL_CONTINUE && iter < iter_max);


    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);
  }





  return 0;
}
