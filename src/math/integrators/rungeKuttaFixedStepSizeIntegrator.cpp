#include "tudat/math/integrators/rungeKuttaFixedStepSizeIntegrator.h"

namespace tudat
{
namespace numerical_integrators
{

template class RungeKuttaFixedStepSizeIntegrator < double, Eigen::VectorXd, Eigen::VectorXd >;
template class RungeKuttaFixedStepSizeIntegrator < double, Eigen::Vector6d, Eigen::Vector6d >;
template class RungeKuttaFixedStepSizeIntegrator < double, Eigen::MatrixXd, Eigen::MatrixXd >;

} // namespace numerical_integrators
} // namespace tudat

