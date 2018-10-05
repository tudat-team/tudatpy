#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"

namespace tudat
{
namespace interpolators
{

template class CubicSplineInterpolator< double, Eigen::VectorXd >;
template class CubicSplineInterpolator< double, Eigen::Vector6d >;
template class CubicSplineInterpolator< double, Eigen::MatrixXd >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class CubicSplineInterpolator< Time, Eigen::VectorXd, long double >;
template class CubicSplineInterpolator< Time, Eigen::Vector6d, long double >;
template class CubicSplineInterpolator< Time, Eigen::MatrixXd, long double >;

template class CubicSplineInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 > >;
template class CubicSplineInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic, 6 > >;
template class CubicSplineInterpolator< double, Eigen::Matrix< long double, Eigen::Dynamic,  Eigen::Dynamic > >;

template class CubicSplineInterpolator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 >, long double >;
template class CubicSplineInterpolator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 6 >, long double >;
template class CubicSplineInterpolator< Time, Eigen::Matrix< long double, Eigen::Dynamic,  Eigen::Dynamic >, long double >;
#endif

} // namespace interpolators
} // namespace tudat

