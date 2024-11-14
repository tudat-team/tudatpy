/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/propagators/nBodyGaussModifiedEquinoctialStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the Gauss planetary equations for modified equinictial elements
Eigen::Vector6d computeGaussPlanetaryEquationsForModifiedEquinoctialElements(
        const Eigen::Vector6d& osculatingModifiedEquinoctialElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter )
{
    using namespace orbital_element_conversions;

    double semiLatusRectrum = osculatingModifiedEquinoctialElements( semiParameterIndex );
    double angularMomentumPerUnitGravitationalParameter =
            std::sqrt( semiLatusRectrum / centralBodyGravitationalParameter );

    double sineTrueLongitude = std::sin( osculatingModifiedEquinoctialElements( trueLongitudeIndex ) );
    double cosineTrueLongitude = std::cos( osculatingModifiedEquinoctialElements( trueLongitudeIndex ) );

    double parameterF = osculatingModifiedEquinoctialElements( fElementIndex );
    double parameterG = osculatingModifiedEquinoctialElements( gElementIndex );
    double parameterH = osculatingModifiedEquinoctialElements( hElementIndex );
    double parameterK = osculatingModifiedEquinoctialElements( kElementIndex );

    double parameterW = 1.0 + parameterF * cosineTrueLongitude + parameterG * sineTrueLongitude;
    double parameterSSquared = 1.0 + parameterH * parameterH + parameterK * parameterK;

    Eigen::Vector6d stateDerivative;
    stateDerivative( 0 ) = 2.0 * semiLatusRectrum / parameterW * angularMomentumPerUnitGravitationalParameter *
          accelerationsInRswFrame( 1 );

    double recurringTermInFGTerms = ( parameterH * sineTrueLongitude - parameterK * cosineTrueLongitude ) / parameterW;
    stateDerivative( 1 ) = angularMomentumPerUnitGravitationalParameter * (

                sineTrueLongitude * accelerationsInRswFrame( 0 ) +

                ( ( parameterW + 1.0 ) * cosineTrueLongitude + parameterF ) / parameterW * accelerationsInRswFrame( 1 ) -

                parameterG * recurringTermInFGTerms * accelerationsInRswFrame( 2 ) );

    stateDerivative( 2 ) = angularMomentumPerUnitGravitationalParameter * (

                -cosineTrueLongitude * accelerationsInRswFrame( 0 ) +

                ( ( parameterW + 1.0 ) * sineTrueLongitude + parameterG ) / parameterW * accelerationsInRswFrame( 1 ) +

                parameterF * recurringTermInFGTerms * accelerationsInRswFrame( 2 ) );

    double recurringTermInHKJTerms = angularMomentumPerUnitGravitationalParameter * parameterSSquared /
            ( 2.0 * parameterW );

    stateDerivative( 3 ) = recurringTermInHKJTerms * cosineTrueLongitude * accelerationsInRswFrame( 2 );
    stateDerivative( 4 ) = recurringTermInHKJTerms * sineTrueLongitude * accelerationsInRswFrame( 2 );

    stateDerivative( 5 ) =
            std::sqrt( semiLatusRectrum * centralBodyGravitationalParameter ) * parameterW * parameterW / (
                semiLatusRectrum  * semiLatusRectrum ) +
            angularMomentumPerUnitGravitationalParameter * recurringTermInFGTerms * accelerationsInRswFrame( 2 );

    return stateDerivative;
}

//template class NBodyGaussModifiedEquinictialStateDerivative< double, double >;

} // namespace propagators

} // namespace tudat
