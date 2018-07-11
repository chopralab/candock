#ifndef PCA_H
#define PCA_H
#include <tuple>

namespace candock {

namespace geometry {
typedef std::tuple<gsl_matrix*, gsl_matrix*, gsl_vector*> PrincipalComponents;
PrincipalComponents pca(const gsl_matrix* data, unsigned int L);
}
}

#endif
