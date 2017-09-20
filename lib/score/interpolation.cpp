#include "helper/inout.hpp"
#include "helper/help.hpp"
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <assert.h>
#include "interpolation.hpp"

namespace Score {
	vector<double> Interpolation::derivative(const vector<double> &y, const double step) {

		assert(y.size() > 0);
		vector<double> deriva;

		deriva.push_back((y[1] - y[0]) / step); // first

		for (size_t i = 1; i < y.size() - 1; ++i) {
			deriva.push_back((y[i + 1] - y[i - 1]) / (2 * step)); // symmetric difference quotient
		}

		deriva.push_back(0); // last
		
		return deriva;
	}
	
	/**
	 * Interpolate through every point
	 * 
	 */
	
	vector<double> Interpolation::interpolate(const vector<double> &dataX, const vector<double> &dataY, const double step) {

		assert(dataX.size() > 0 && dataX.size() == dataY.size());
		vector<double> pot;
		
		gsl_interp_accel *acc = gsl_interp_accel_alloc();
		gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, dataX.size());
		gsl_spline_init(spline, &dataX[0], &dataY[0], dataX.size());
		
		for (double xi = dataX.front(); xi <= dataX.back(); xi += step) {
			dbgmsg("x = " << xi << " yeval = " << gsl_spline_eval(spline, xi, acc));
			pot.push_back(gsl_spline_eval(spline, xi, acc));
		}

		gsl_spline_free(spline);
		gsl_interp_accel_free(acc);

		return pot;
	}

	/**
	 * B-spline interpolation to get smooth curve
	 * ------------------------------------------
	 * k : number of data points to fit
	 * ncoeff : number of fit coefficients
	 */

	vector<double> Interpolation::interpolate_bspline(const vector<double> &dataX, const vector<double> &dataY, const double step, const size_t k, const size_t ncoeffs) {

		assert(dataX.size() > 0 && dataX.size() == dataY.size());
                
		vector<double> pot;

		const size_t n = dataX.size();
                const size_t ncoeffs2 = (ncoeffs < n)? ncoeffs : n;
		const size_t nbreak = ncoeffs2 + 2 - k; // nbreak = ncoeffs2 + 2 - k

		gsl_rng_env_setup();

		gsl_bspline_workspace *bw = gsl_bspline_alloc(k, nbreak); // allocate a cubic bspline workspace (k = 4)
		gsl_vector *B = gsl_vector_alloc(ncoeffs2);
		gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
		gsl_vector *c = gsl_vector_alloc(ncoeffs2), *w = gsl_vector_alloc(n);
		gsl_vector *x = gsl_vector_alloc(n), *y = gsl_vector_alloc(n);
		gsl_matrix *X = gsl_matrix_alloc(n, ncoeffs2), *cov = gsl_matrix_alloc(ncoeffs2, ncoeffs2);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(n, ncoeffs2);
		double chisq;
	
		for (size_t i = 0; i < dataX.size(); ++i) {
			const double xi = dataX[i];
			const double yi = dataY[i];
			const double sigma = 0.1 * yi;
		
			gsl_vector_set(x, i, xi);
			gsl_vector_set(y, i, yi);
			gsl_vector_set(w, i, 1.0 / (sigma * sigma));
		}

		/* use uniform breakpoints on [0, 15] */
		gsl_bspline_knots_uniform(dataX.front(), dataX.back(), bw);
	
		/* construct the fit matrix X */
		for (size_t i = 0; i < n; ++i) {
			double xi = gsl_vector_get(x, i);
			
			/* compute B_j(xi) for all j */
			gsl_bspline_eval(xi, B, bw);
			
			/* fill in row i of X */
			for (size_t j = 0; j < ncoeffs2; ++j) {
				double Bj = gsl_vector_get(B, j);
				gsl_matrix_set(X, i, j, Bj);
			}
	    }

		/* do the fit */
		gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);
#ifdef NDEBUG
		gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
#else
                double tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
		dbgmsg("chisq/dof = " << chisq / (n - ncoeffs) << ", Rsq = " << 1.0 - chisq / tss << "\n");
#endif
		/* output the smoothed curve */
		for (double xi = dataX.front(); xi <= dataX.back(); xi += step) {
			double yi, yerr;		
			gsl_bspline_eval(xi, B, bw);
			gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
			pot.push_back(yi);
		}
		
		gsl_rng_free(r);
		gsl_bspline_free(bw);
		gsl_vector_free(B);
		gsl_vector_free(x);
		gsl_vector_free(y);
		gsl_matrix_free(X);
		gsl_vector_free(c);
		gsl_vector_free(w);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
	
		return pot;
	}

}
