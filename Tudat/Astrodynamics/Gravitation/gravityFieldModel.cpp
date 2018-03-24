/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <boost/make_shared.hpp>
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"

namespace tudat
{
namespace gravitation
{

//! Set predefined central gravity field settings.
boost::shared_ptr< GravityFieldModel > getPredefinedCentralGravityField(
    BodiesWithPredefinedCentralGravityFields bodyWithPredefinedCentralGravityField )
{
    double gravitationalParameter = 0.0;

    // Select body with prefined central gravity field.
    switch( bodyWithPredefinedCentralGravityField )
    {
    case sun:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter = 1.32712440018e20;

        break;

    case mercury:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter = 2.203289218e13;

        break;

    case venus:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter = 3.2485504415e14;

        break;

    case earth:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter = 3.9859383624e14;

        break;

    case moon:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter = 4.903686391e12;

        break;

    case mars:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter = 4.2828018915e13;

        break;

    case jupiter:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.3, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter = 1.2668579374e17;

        break;

    case saturn:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.3, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter = 3.793100511400001e16;

        break;

    case uranus:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.3, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter = 5.793943348799999e15;

        break;

    case neptune:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.3, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter = 6.834733937e15;

        break;

    default:

        std::string errorMessage = "Desired predefined central gravity field " +
                std::to_string( bodyWithPredefinedCentralGravityField ) +
                " does not exist";
        throw std::runtime_error( errorMessage );
    }
    return boost::make_shared< GravityFieldModel >( gravitationalParameter );
}

} // namespace gravitation
} // namespace tudat
