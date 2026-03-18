#ifndef TUDAT_CALYLEY_KLEIN_PARAMETERS_H
#define TUDAT_CALYLEY_KLEIN_PARAMETERS_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <complex>

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace basic_mathematics
{

void convertQuaterionToCayleyKleinParameters(
        const Eigen::Quaterniond quaternion, std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB );

void convert323EulerAnglesToCayleyKleinParameters(
        const double firstZRotation, const double yRotation, const double secondZRotation,
        std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB );

void convert313EulerAnglesToCayleyKleinParameters(
        const double firstZRotation, const double xRotation, const double secondZRotation,
        std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB );

} // namespace basic_mathematics

} // namespace tudat

#endif // TUDAT_CALYLEY_KLEIN_PARAMETERS_H
