/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/firstOrderRelativisticLightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute partial derivative of 1st order relativistic correction w.r.t. gravitational parameter.
double getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtSingleGravitationalParameter(
        const double singleBodyLightTimeCorrection, const double bodyGravitationalParameter )
{
    return singleBodyLightTimeCorrection / bodyGravitationalParameter;
}

//! Function to compute partial derivative of 1st order relativistic correction w.r.t. PPN parameter gamma.
double getPartialOfFirstOrderRelativisticLightTimeCorrectionWrtPpnParameterGamma(
        const double totalLightTimeCorrection, const double ppnParameterGamma )
{
    return totalLightTimeCorrection / ( ppnParameterGamma + 1.0 );
}



}

}
