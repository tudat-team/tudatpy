/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120209    K. Kumar          File created.
 *
 *    References
 */

#include <cmath>

#include "Tudat/Astrodynamics/Gravitation/gravitationalAccelerationModel.h"

namespace tudat
{
namespace astrodynamics
{
namespace acceleration_models
{

//! Compute gravitational acceleration.
Eigen::Vector3d computeGravitationalAcceleration(
        const double universalGravitationalParameter,
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double massOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    return computeGravitationalAcceleration(
                positionOfBodySubjectToAcceleration,
                universalGravitationalParameter * massOfBodyExertingAcceleration,
                positionOfBodyExertingAcceleration );
}

//! Compute gravitational acceleration.
Eigen::Vector3d computeGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    return -gravitationalParameterOfBodyExertingAcceleration
            * ( positionOfBodySubjectToAcceleration - positionOfBodyExertingAcceleration )
            / std::pow( ( positionOfBodySubjectToAcceleration
                          - positionOfBodyExertingAcceleration ).norm( ), 3.0 );
}

} // namespace acceleration_models
} // namespace astrodynamics
} // namespace tudat
