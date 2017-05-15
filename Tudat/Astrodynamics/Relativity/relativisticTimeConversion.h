#ifndef TUDAT_RELATIVISTIVTIMECONVERSION_H
#define TUDAT_RELATIVISTIVTIMECONVERSION_H

#include <Eigen/Core>

#include <Tudat/Basics/basicTypedefs.h>

namespace tudat
{

namespace relativity
{

double calculateFirstCentralBodyProperTimeRateDifference(
        const double relativeSpeed, const double gravitationalScalarPotential );

double calculateFirstCentralBodyProperTimeRateDifference(
        const Eigen::Vector6d relativeStateVector, const double centralBodyGravitationalParameter );

}

}

#endif // TUDAT_RELATIVISTIVTIMECONVERSION_H
