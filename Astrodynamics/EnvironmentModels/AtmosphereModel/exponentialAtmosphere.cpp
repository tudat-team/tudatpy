/*! \file exponentialAtmosphere.cpp
 *    Source file that defines the exponential atmosphere model included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/AtmosphereModel/
 *    Version           : 4
 *    Check status      : Unchecked
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 24 February, 2011
 *    Last modified     : 29 June, 2011
 *
 *    References
 *
 *    Notes
 *      The accuracy of this model could be increased by implementing different
 *      values for the scale height and temperature for different altitudes
 *      (e.g., lower, middle and upper atmosphere).
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      110224    F.M. Engelen      File created.
 *      110324    J. Melman         Added overloaded get functions.
 *      110427    F.M. Engelen      Changed input parameter to altitude, longitude and latitude.
 *      110629    F.M. Engelen      Added predefined feature.
 *      110705    F.M. Engelen      Changed to passing by reference. Changed reference values.
 */

// Include statements.
#include "exponentialAtmosphere.h"

//! Default constructor.
ExponentialAtmosphere::ExponentialAtmosphere( )
{
    scaleHeight_ = -0.0;
    densityAtZeroAltitude_ = -0.0;
    constantTemperature_ = -0.0;
}

//! Default destructor.
ExponentialAtmosphere::~ExponentialAtmosphere( )
{
}

//! Set predefined exponential atmosphere settings.
void ExponentialAtmosphere::setPredefinedExponentialAtmosphere(
        BodiesWithPredefinedExponentialAtmospheres bodyWithPredefindExponentialAtmosphere )
{
    switch( bodyWithPredefindExponentialAtmosphere )
    {
    case earth:
        // Set local variables for Earth exponential atmosphere. Based on  lecture notes
        // Rocket Motion by Prof. Ir. B.A.C. Ambrosius, November 2009.

        // Set scale height.
        this->setScaleHeight( 7.200e3 );

        //Set density at zero altitude.
        this->setDensityAtZeroAltitude( 1.225 );

        //Set atmosphere temperature.
        this->setConstantTemperature( 246.0 );

        break;

    default:
        std::cerr << "This is not a body with a predefined exponential atmophere." << std::endl;
    }
}

//! Set scale altitude.
void ExponentialAtmosphere::setScaleHeight( const double& scaleHeight )
{
    scaleHeight_ = scaleHeight;
}

//! Get scale altitude.
double ExponentialAtmosphere::getScaleHeight( )
{
    return scaleHeight_;
}

//! Set density at zero altitude.
void ExponentialAtmosphere::setDensityAtZeroAltitude( const double& densityAtZeroAltitude )
{
    densityAtZeroAltitude_ = densityAtZeroAltitude;
}

//! Get density at zero altitude.
double ExponentialAtmosphere::getDensityAtZeroAltitude( )
{
    return densityAtZeroAltitude_;
}

//! Set constant temperature.
void ExponentialAtmosphere::setConstantTemperature( const double& constantTemperature )
{
    constantTemperature_ = constantTemperature;
}

//! Get constant temperature.
double ExponentialAtmosphere::getConstantTemperature( )
{
    return constantTemperature_;
}

//! Get local density.
double ExponentialAtmosphere::getDensity( const double& altitude )
{
    return densityAtZeroAltitude_ * exp( - altitude / scaleHeight_ );
}


//! Get local density in the general way.
double ExponentialAtmosphere::getDensity( const double& altitude,
                                          const double& longitude,
                                          const double& latitude,
                                          const double& time )
{
    return getDensity( altitude );
}

//! Get local pressure.
double ExponentialAtmosphere::getPressure( const double& altitude )
{
    return getDensity( altitude ) * PhysicalConstants::SPECIFIC_GAS_CONSTANT_AIR *
            constantTemperature_;
}

//! Get local pressure in the general way.
double ExponentialAtmosphere::getPressure( const double& altitude,
                                           const double& longitude,
                                           const double& latitude,
                                           const double& time )
{
    return getPressure( altitude );
}

//! Get local temperature.
double ExponentialAtmosphere::getTemperature( const double& altitude )
{
    return constantTemperature_;
}

//! Get local temperature in the general way.
double ExponentialAtmosphere::getTemperature( const double& altitude,
                                              const double& longitude,
                                              const double& latitude,
                                              const double& time )
{
    return constantTemperature_;
}

// End of file.
