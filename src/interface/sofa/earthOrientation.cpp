/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/interface/sofa/earthOrientation.h"
#include "tudat/basics/timeType.h"



namespace tudat
{

namespace sofa_interface
{

//! Function to calculate CIP and CIO locator according to requested IAU conventions
std::pair< Eigen::Vector2d, double > getPositionOfCipInGcrs(
        const double terrestrialTime, const double julianDaysEpochShift,
        const basic_astrodynamics::IAUConventions precessionNutationTheory )
{
    // Declare Sofa function return arguments (by reference)
    double xAngle, yAngle;
    double originLocator;

    // Check for IAU convention and retrieve requested values.
    switch( precessionNutationTheory )
    {
    case basic_astrodynamics::iau_2000_a:
        iauXys00a( julianDaysEpochShift, terrestrialTime / physical_constants::JULIAN_DAY,
                   &xAngle, &yAngle, &originLocator );
        break;

    case basic_astrodynamics::iau_2000_b:
        iauXys00b( julianDaysEpochShift, terrestrialTime / physical_constants::JULIAN_DAY,
                   &xAngle, &yAngle, &originLocator );
        break;

    case basic_astrodynamics::iau_2006:
        iauXys06a( julianDaysEpochShift, terrestrialTime / physical_constants::JULIAN_DAY,
                   &xAngle, &yAngle, &originLocator );
        break;
    default:
        throw std::runtime_error( "Warning, precession nutation theory selection not recongnized" );

    }

    // Set and return requested values.
    Eigen::Vector2d cioPosition;
    cioPosition << xAngle, yAngle;
    return std::pair< Eigen::Vector2d, double >( cioPosition, originLocator );
}

//! Function to calculate GMST according to requested IAU conventions
double calculateGreenwichMeanSiderealTime(
        const double terrestrialTime, const double universalTime1,
        const double referenceJulianDay, const basic_astrodynamics::IAUConventions iauConvention )
{
    // Declare GMST variable
    double gmst = TUDAT_NAN;

    // Check for IAU convention and retrieve requested GMST
    switch( iauConvention )
    {
    case basic_astrodynamics::iau_2000_a:
        gmst = iauGmst00( referenceJulianDay, universalTime1 / physical_constants::JULIAN_DAY,
                          referenceJulianDay, terrestrialTime / physical_constants::JULIAN_DAY );
        break;

    case basic_astrodynamics::iau_2000_b:
        gmst = iauGmst00( referenceJulianDay, universalTime1 / physical_constants::JULIAN_DAY,
                          referenceJulianDay, terrestrialTime / physical_constants::JULIAN_DAY );
        break;

    case basic_astrodynamics::iau_2006:
        gmst = iauGmst06( referenceJulianDay, universalTime1 / physical_constants::JULIAN_DAY,
                          referenceJulianDay, terrestrialTime / physical_constants::JULIAN_DAY );
        break;
    default:
        throw std::runtime_error( "Warning, IAU convention for GMST calculation not recognized" );

    }

    return gmst;

}


double calculateEquationOfEquinoxes(
		const double terrestrialTime, const double referenceJulianDay,
		const basic_astrodynamics::IAUConventions iauConvention )
{
	double equationOfEquinoxes = TUDAT_NAN;

	switch( iauConvention )
	{
		case basic_astrodynamics::iau_2000_a:
			equationOfEquinoxes = iauEe00a( referenceJulianDay, terrestrialTime / physical_constants::JULIAN_DAY );
			break;
		case basic_astrodynamics::iau_2000_b:
			equationOfEquinoxes = iauEe00b( referenceJulianDay, terrestrialTime / physical_constants::JULIAN_DAY );
			break;
		case basic_astrodynamics::iau_2006:
			equationOfEquinoxes = iauEe06a( referenceJulianDay, terrestrialTime / physical_constants::JULIAN_DAY );
			break;
		default:
			throw std::runtime_error( "Warning, IAU convention for EoE calculation not recognized" );
	}

	return equationOfEquinoxes;
}

Eigen::Matrix3d getPrecessionNutationMatrix( const double terrestrialTime, const double referenceJulianDay )
{
	double pnm[3][3];

//	std::cout << "Ref JD: " << referenceJulianDay << "\nTT in Julian days: " << terrestrialTime / physical_constants::JULIAN_DAY << std::endl;

	iauPnm80( referenceJulianDay, terrestrialTime / physical_constants::JULIAN_DAY, pnm );

	return ( Eigen::Matrix3d( ) << pnm[ 0 ][ 0 ], pnm[ 0 ][ 1 ], pnm[ 0 ][ 2 ],
			pnm[ 1 ][ 0 ], pnm[ 1 ][ 1 ], pnm[ 1 ][ 2 ],
			pnm[ 2 ][ 0 ], pnm[ 2 ][ 1 ], pnm[ 2 ][ 2 ] ).finished( );
}

void getPrecessionAngles( double &zeta, double &z, double &theta, const double terrestrialTime, const double referenceJulianDay )
{
	iauPrec76( referenceJulianDay, terrestrialTime, referenceJulianDay, 0, &zeta, &z, &theta );
}

//! Function to calculate ERA (earth rotation angle)
double calculateEarthRotationAngle( const double ut1, const double julianDaysEpochShift )
{
    return iauEra00( julianDaysEpochShift, ut1 / physical_constants::JULIAN_DAY );
}

//! Function to calculate ERA (earth rotation angle) in double precision
template< >
double calculateEarthRotationAngleTemplated< double >(
        const double currentUt1 )
{
    return calculateEarthRotationAngle( currentUt1, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}

//! Function to calculate ERA (earth rotation angle) in Time precision
template< >
double calculateEarthRotationAngleTemplated< Time >(
        const Time currentUt1 )
{

    int hoursSinceEpoch = currentUt1.getFullPeriods( );
    int fullDaysSinceEpoch = currentUt1.getFullPeriods( ) / 24;
    int hoursIntoCurrentDay = hoursSinceEpoch - 24 * fullDaysSinceEpoch;


    return calculateEarthRotationAngle( currentUt1.getSecondsIntoFullPeriod( ) + hoursIntoCurrentDay * 3600.0,
                                        basic_astrodynamics::JULIAN_DAY_ON_J2000 + fullDaysSinceEpoch );
}

//! Function to compute the rotation matrix from GCRS to J2000 at epoch
Eigen::Matrix3d getFrameBias(
        const double julianDaysSinceReference,
        const basic_astrodynamics::IAUConventions precessionNutationTheory,
        const double referenceJulianDay )
{
    double rb[3][3], rp[3][3], rbp[3][3];
    if( precessionNutationTheory == basic_astrodynamics::iau_2006 )
    {
        iauBp06( referenceJulianDay, julianDaysSinceReference, rb, rp, rbp );
    }
    else
    {
        iauBp00( referenceJulianDay, julianDaysSinceReference, rb, rp, rbp );
    }
    return ( Eigen::Matrix3d( )<< rb[ 0 ][ 0 ], rb[ 0 ][ 1 ], rb[ 0 ][ 2 ],
            rb[ 1 ][ 0 ], rb[ 1 ][ 1 ], rb[ 1 ][ 2 ],
            rb[ 2 ][ 0 ], rb[ 2 ][ 1 ], rb[ 2 ][ 2 ] ).finished( );

}

}

}
