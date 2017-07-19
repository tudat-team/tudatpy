#ifndef TUDAT_RELATIVISTIVTIMECONVERSION_H
#define TUDAT_RELATIVISTIVTIMECONVERSION_H

#include <Eigen/Core>

#include <Tudat/Basics/basicTypedefs.h>

namespace tudat
{

namespace relativity
{

double calculateFirstCentralBodyProperTimeRateDifference(
        const double relativeSpeed, const double gravitationalScalarPotential,
        const double equivalencePrincipleLpiViolationParameter = 0.0 );

double calculateFirstCentralBodyProperTimeRateDifference(
        const Eigen::Vector6d relativeStateVector, const double centralBodyGravitationalParameter,
        const double equivalencePrincipleLpiViolationParameter = 0.0 );

}

}

#endif // TUDAT_RELATIVISTIVTIMECONVERSION_H
