#include "candock/score/interpolation.hpp"
#include <assert.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>
#include <functional>
#include "candock/helper/help.hpp"
#include "candock/helper/inout.hpp"
using namespace std;

namespace candock {

namespace score {
namespace Interpolation {
vector<double> derivative(const vector<double>& y, const double step) {
    assert(y.size() > 0);
    vector<double> deriva;

    deriva.push_back((y[1] - y[0]) / step);  // first

    for (size_t i = 1; i < y.size() - 1; ++i) {
        deriva.push_back((y[i + 1] - y[i - 1]) /
                         (2 * step));  // symmetric difference quotient
    }

    deriva.push_back(0);  // last

    return deriva;
}

/**
 * Interpolate through every point
 *
 */

vector<double> interpolate(const vector<double>& dataX,
                           const vector<double>& dataY, const double step) {
    assert(dataX.size() > 0 && dataX.size() == dataY.size());
    vector<double> pot;

    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, dataX.size());
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

BSplineFit::BSplineFit(const vector<double>& dataX, const vector<double>& dataY,
                       const size_t k, const size_t ncoeffs)
    : n(dataX.size()),
      ncoeffs2((ncoeffs < n) ? ncoeffs : n),
      nbreak(ncoeffs2 + 2 - k),  // nbreak = ncoeffs2 + 2 - k
      bw(gsl_bspline_alloc(
          k, nbreak)),  // allocate a cubic bspline workspace (k = 4)
      B(gsl_vector_alloc(ncoeffs2)),
      c(gsl_vector_alloc(ncoeffs2)),
      w(gsl_vector_alloc(n)),
      x(gsl_vector_alloc(n)),
      y(gsl_vector_alloc(n)),
      X(gsl_matrix_alloc(n, ncoeffs2)),
      cov(gsl_matrix_alloc(ncoeffs2, ncoeffs2)),
      mw(gsl_multifit_linear_alloc(n, ncoeffs2))

{
    assert(dataX.size() > 0 && dataX.size() == dataY.size());

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
    dbgmsg("chisq/dof = " << chisq / (n - ncoeffs)
                          << ", Rsq = " << 1.0 - chisq / tss << "\n");
#endif
}

vector<double> BSplineFit::interpolate_bspline(const double start,
                                               const double end,
                                               const double step) {
    vector<double> pot;

    /* output the smoothed curve */
    for (double xi = start; xi <= end; xi += step) {
        double yi, yerr;
        gsl_bspline_eval(xi, B, bw);
        gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
        pot.push_back(yi);
    }

    return pot;
}

BSplineFit::~BSplineFit() {
    gsl_bspline_free(bw);
    gsl_vector_free(B);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_vector_free(w);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(mw);
}
}
}
}
