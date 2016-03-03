#include <stdexcept>

#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"

namespace tudat
{
namespace aerodynamics
{

ExponentialAtmosphere::ExponentialAtmosphere(
        const BodiesWithPredefinedExponentialAtmospheres bodyWithPredefinedExponentialAtmosphere )
{
    switch( bodyWithPredefinedExponentialAtmosphere )
    {
    case earth:
        // Set local variables for Earth exponential atmosphere. Based on  lecture notes
        // Rocket Motion by Prof. Ir. B.A.C. Ambrosius, November 2009.

        // Set scale height.
        scaleHeight_ = 7.200e3;

        //Set density at zero altitude.
        densityAtZeroAltitude_ = 1.225;

        //Set atmosphere temperature.
        constantTemperature_ = 246.0;

        // Set specific gas constant.
        specificGasConstant_ = physical_constants::SPECIFIC_GAS_CONSTANT_AIR;

        ratioOfSpecificHeats_ = 1.4;

        break;

    default:
        throw std::runtime_error(
                    "Error when making exponential atmosphere, predefined atmosphere" +
                    boost::lexical_cast< std::string >(
                        bodyWithPredefinedExponentialAtmosphere ) + "not recognized." );
    }
}

}

}
