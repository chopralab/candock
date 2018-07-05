#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include <vector>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

namespace candock {

namespace score {
        namespace Interpolation {
                std::vector<double> derivative(const std::vector<double> &y, const double step);
                std::vector<double> interpolate(const std::vector<double> &dataX, const std::vector<double> &dataY, const double step);

                class BSplineFit {
                        const size_t n;
                        const size_t ncoeffs2;
                        const size_t nbreak;

                        gsl_bspline_workspace *bw;
                        gsl_vector *B;
                        gsl_vector *c;
                        gsl_vector *w;
                        gsl_vector *x;
                        gsl_vector *y;
                        gsl_matrix *X;
                        gsl_matrix *cov;
                        gsl_multifit_linear_workspace *mw;

                public:
                        BSplineFit (const std::vector<double> &dataX, const std::vector<double> &dataY, const size_t k=4, const size_t ncoeff=24 );
                        virtual ~BSplineFit();
                        std::vector<double> interpolate_bspline( const double start, const double end, const double step );
                };
        }
}

}

#endif
