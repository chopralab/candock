#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include <vector>

namespace Score {
        namespace Interpolation {
                std::vector<double> derivative(const std::vector<double> &y, const double step);
                std::vector<double> interpolate(const std::vector<double> &dataX, const std::vector<double> &dataY, const double step);
                std::vector<double> interpolate_bspline(const std::vector<double> &dataX, const std::vector<double> &dataY, const double step, const size_t k=4, const size_t ncoeff=24);
        }
}

#endif
