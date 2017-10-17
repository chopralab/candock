#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <vector>
#include <tuple>
#include <exception>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "helper/debug.hpp"
#include "helper/error.hpp"
#include "powerfit.hpp"

namespace Score {

        struct data {
                size_t  n;
                size_t  p;
                double *x;
                double *y;
        };

        int ex_rep_f (const gsl_vector *params, void *data, gsl_vector *f) {

                size_t  n = ( (struct data *) data)->n;
                size_t  p = ( (struct data *) data)->p;
                double *x = ( (struct data *) data)->x;
                double *y = ( (struct data *) data)->y;

                double sigma = gsl_vector_get (params, 0);
                double b = gsl_vector_get (params, 1);

                for (size_t i = 0; i < n; i++) {
                        /* Model Yi =sigma/(i)^p + b */
                        double t = x[i];
                        double Yi = sigma / std::pow (t, p) + b;
                        gsl_vector_set (f, i, Yi - y[i]);
                }

                return GSL_SUCCESS;
        }

        int ex_rep_df (const gsl_vector *params, void *data, gsl_matrix *J) {
                const size_t n = ( (struct data *) data)->n;
                const size_t p = ( (struct data *) data)->p;
                double *x = ( (struct data *) data)->x;

                for (size_t i = 0; i < n; i++) {
                        /* Jacobian matrix J(i,j) = dfi / dxj, */
                        /* where fi = (Yi - yi)/sigma[i],      */
                        /*       Yi = sigman/(i)^p + b  */
                        /* and the xj are the parameters (sigma,b) */
                        double t = x[i];
                        gsl_matrix_set (J, i, 0, 1.0 / std::pow (t,p));
                        gsl_matrix_set (J, i, 1, 1.0);
                }

                return GSL_SUCCESS;
        }

        void ex_rep_callback (const size_t iter, void *unused, const gsl_multifit_nlinear_workspace *w) {
#ifndef NDEBUG
                gsl_vector *f      = gsl_multifit_nlinear_residual (w);
                gsl_vector *params = gsl_multifit_nlinear_position (w);
#endif
                double rcond;

                /* compute reciprocal condition number of J(x) */
                gsl_multifit_nlinear_rcond (&rcond, w);

                dbgmsg ("iter: " << iter
                        << " sigma = " << gsl_vector_get (params, 0)
                        << ", b = " << gsl_vector_get (params, 1)
                        << ", conj(J) = " << 1.0 / rcond
                        << ", |f(x)| = " << gsl_blas_dnrm2 (f));
        }

        std::tuple<double, double, double> fit_power_function (std::vector<double> x, std::vector<double> y,
                        size_t power, std::pair<double,double> guess,
                        size_t max_iter) {
                if (x.size() != y.size()) {
                        throw std::exception();
                }

                // number of parameters
                const size_t param_count = 2;

                const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
                gsl_multifit_nlinear_workspace *w;
                gsl_multifit_nlinear_fdf fdf;
                gsl_multifit_nlinear_parameters fdf_params =
                        gsl_multifit_nlinear_default_parameters();

                data d = { x.size(), power, x.data(), y.data() };
                /* define the function to be minimized */
                fdf.f = ex_rep_f;
                fdf.df = ex_rep_df;   /* set to NULL for finite-difference Jacobian */
                fdf.fvv = NULL;     /* not using geodesic acceleration */
                fdf.n = x.size();
                fdf.p = param_count;
                fdf.params = &d;

                gsl_vector *f;
                gsl_matrix *J;
                gsl_matrix *covar = gsl_matrix_alloc (param_count, param_count);

                double x_init[2] = { guess.first, guess.second }; /* starting values */

                gsl_vector_view params = gsl_vector_view_array (x_init, param_count);
                double chisq, chisq0;
                int status, info;

                const double xtol = 1e-12;
                const double gtol = 1e-12;
                const double ftol = 0.0;

                /* allocate workspace with default parameters */
                w = gsl_multifit_nlinear_alloc (T, &fdf_params, x.size(), param_count);

                /* initialize solver with starting point and weights */
                gsl_multifit_nlinear_init (&params.vector, &fdf, w);

                /* compute initial cost function */
                f = gsl_multifit_nlinear_residual (w);

                gsl_blas_ddot (f, f, &chisq0);

                /* solve the system with a maximum of 20 iterations */
#ifndef _NDEBUG
                status = gsl_multifit_nlinear_driver (max_iter, xtol, gtol, ftol,
                                                      ex_rep_callback, NULL, &info, w);
#else
                status = gsl_multifit_nlinear_driver (max_iter, xtol, gtol, ftol,
                                                      NULL, NULL, &info, w);

#endif
                /* compute covariance of best fit parameters */
                J = gsl_multifit_nlinear_jac (w);

                gsl_multifit_nlinear_covar (J, 0.0, covar);

                /* compute final cost */
                gsl_blas_ddot (f, f, &chisq);


                dbgmsg ("summary from method '" << gsl_multifit_nlinear_name (w) << "/" << gsl_multifit_nlinear_trs_name (w));
                dbgmsg ("number of iterations: " << gsl_multifit_nlinear_niter (w));
                dbgmsg ("function evaluations: " << fdf.nevalf);
                dbgmsg ("Jacobian evaluations: " << fdf.nevaldf);
                dbgmsg ("reason for stopping: ");
                
                if (info == 1) {
                        dbgmsg("small step size");
                } else {
                        dbgmsg("small gradient\n");
                }

                dbgmsg ("initial |f(x)| = " << std::sqrt (chisq0));
                dbgmsg ("final   |f(x)| = " << std::sqrt (chisq));
                dbgmsg ("dof            = " << x.size() - param_count);

                double dof = x.size() - param_count;
                double c = std::max (1.0, std::sqrt (chisq / dof));

                double sigma = gsl_vector_get (w->x, 0);
                double b     = gsl_vector_get (w->x, 1);

                dbgmsg ("chisq/dof = " << chisq / dof);
                dbgmsg ("sigma  = " << sigma << " +/- " << c *std::sqrt (gsl_matrix_get (covar,0,0)));
                dbgmsg ("b      = " << b     << " +/- " << c *std::sqrt (gsl_matrix_get (covar,1,1)));

                dbgmsg ("status = " << gsl_strerror (status) << "\n");

                gsl_multifit_nlinear_free (w);
                gsl_matrix_free (covar);

                if (status != GSL_SUCCESS) {
                        throw Error ("Fitting to powerfunction failed");
                }

                return std::make_tuple (sigma,b,c);
        }

        std::tuple<double, double, double> fit_range_power_function (std::vector<double> x, std::vector<double> y) {

                double best_sigma, best_b, best_p;
                double best_error = HUGE_VAL;
#ifndef NDEBUG
                double best_guess_sigma, best_guess_b;
#endif
                // try a range of coefficients to get the best fit
                for (double b_guess = -2.00; b_guess <= 2; b_guess += 0.50)
                        for (size_t p_guess = 9; p_guess <= 13; p_guess += 1)
                                for (double s_guess = 1e+5; s_guess < 1e+10; s_guess *= 2) {
                                        double sigma, b, chi;

                                        try {
                                                std::tie (sigma, b, chi) = fit_power_function (x,y,p_guess,
                                                                           std::make_pair (s_guess,b_guess),
                                                                           100);
                                        } catch (Error &e) {
                                                continue;
                                        }

                                        if (chi < best_error) {
                                                best_sigma = sigma;
                                                best_b = b;
                                                best_error = chi;
#ifndef NDEBUG
                                                best_guess_sigma = s_guess;
                                                best_guess_b = b_guess;
#endif
                                                best_p = p_guess;
                                        }
                                }

                dbgmsg ("Best error: " << best_error);
                dbgmsg ("Best  guess: (" << best_guess_sigma << ", " << best_guess_b << ")");

                if (best_error == HUGE_VAL) {
                        throw Error ("warning : could not fit repulsion term, zero everything");
                }

                return std::make_tuple (best_sigma, best_b, best_p);
        }

        std::tuple<double, double, double> fit_range_power_function_fast (std::vector<double> x, std::vector<double> y) {

                double best_sigma, best_b, best_p;
                double best_error = HUGE_VAL;
#ifndef NDEBUG
                double best_guess_sigma, best_guess_b;
#endif
                // try a range of coefficients to get the best fit
                double b_guess = -1.00;
                size_t p_guess = 12;

                for (double s_guess = 1e+5; s_guess < 1e+10; s_guess *= 2) {
                        double sigma, b, chi;

                        try {
                                std::tie (sigma, b, chi) = fit_power_function (x,y,p_guess,
                                                           std::make_pair (s_guess,b_guess),
                                                           100);
                        } catch (Error &e) {
                                continue;
                        }

                        if (chi < best_error) {
                                best_sigma = sigma;
                                best_b = b;
                                best_error = chi;
#ifndef NDEBUG
                                best_guess_sigma = s_guess;
                                best_guess_b = b_guess;
#endif
                                best_p = p_guess;
                        }
                }

                dbgmsg ("Best error: " << best_error);
                dbgmsg ("Best  guess: (" << best_guess_sigma << ", " << best_guess_b << ")");

                if (best_error == HUGE_VAL) {
                        throw InterpolationError ("warning : could not fit repulsion term, zero everything");
                }

                return std::make_tuple (best_sigma, best_b, best_p);
        }

}
