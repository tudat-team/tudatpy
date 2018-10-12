#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace interpolators
{

template class LagrangeInterpolator< double, Eigen::VectorXd >;
template class LagrangeInterpolator< double, Eigen::Vector6d >;
template class LagrangeInterpolator< double, Eigen::MatrixXd >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class LagrangeInterpolator< Time, Eigen::VectorXd, long double >;
template class LagrangeInterpolator< Time, Eigen::Vector6d, long double >;
template class LagrangeInterpolator< Time, Eigen::MatrixXd, long double >;

template class LagrangeInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 > >;
template class LagrangeInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic, 6 > >;
template class LagrangeInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic,  Eigen::Dynamic > >;

template class LagrangeInterpolator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 >, long double >;
template class LagrangeInterpolator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 6 >, long double >;
template class LagrangeInterpolator< Time, Eigen::Matrix< long double, Eigen::Dynamic,  Eigen::Dynamic >, long double >;
#endif

} // namespace interpolators

} // namespace tudat
