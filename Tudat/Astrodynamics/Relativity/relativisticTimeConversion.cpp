#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Relativity/relativisticTimeConversion.h"

namespace tudat
{

namespace relativity
{

double calculateFirstCentralBodyProperTimeRateDifference(
        const double relativeSpeed, const double gravitationalScalarPotential )
{
    return ( -( 0.5 * relativeSpeed * relativeSpeed + gravitationalScalarPotential )  ) *
            physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
}

double calculateFirstCentralBodyProperTimeRateDifference(
        const Eigen::Vector6d relativeStateVector, const double centralBodyGravitationalParameter )
{
    return calculateFirstCentralBodyProperTimeRateDifference(
                relativeStateVector.segment( 3, 3 ).norm( ), centralBodyGravitationalParameter /
                relativeStateVector.segment( 0, 3 ).norm( ) );
}

}

}

