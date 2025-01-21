/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/interface/sofa/earthOrientation.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "boost/algorithm/string.hpp"

namespace tudat
{


namespace ephemerides
{

Eigen::Matrix3d getRotationMatrixFromTemeToJ2000( const double epochSinceJ2000  )
{
    // First, rotate to the True Of Date (TOD) frame.
    double equationOfEquinoxes = sofa_interface::calculateEquationOfEquinoxes( epochSinceJ2000 );

    // Rotate around pole (z-axis)
    Eigen::AngleAxisd rotationObject = Eigen::AngleAxisd( equationOfEquinoxes, Eigen::Vector3d::UnitZ( ) );
    Eigen::Matrix3d rotationMatrix1 = rotationObject.toRotationMatrix( );

    // Now that we have our state vector in the TOD frame, we need to obtain the combined precession + nutation matrix from Sofa
    // (according to the 1976/1980 model)
    Eigen::Matrix3d precessionNutationMatrix = sofa_interface::getPrecessionNutationMatrix( epochSinceJ2000 );
    Eigen::Matrix3d rotationMatrix2 = precessionNutationMatrix.transpose( );

    return rotationMatrix2 * rotationMatrix1;
}

Eigen::Matrix3d getRotationMatrixFromJ2000ToTeme( const double epochSinceJ2000  )
{
    return getRotationMatrixFromTemeToJ2000( epochSinceJ2000 ).transpose( );
}


	TleEphemeris::TleEphemeris( const std::string &referenceFrameOrigin, const std::string &referenceFrameOrientation,
							   const std::shared_ptr< Tle > tle_ptr, const bool useSDP ) :
							   Ephemeris( referenceFrameOrigin, referenceFrameOrientation )
	{
		if( referenceFrameOrigin != "Earth" )
		{
			throw std::runtime_error( "Error: TleEphemeris only supports an Earth-centered reference frame." );
		}
		if( useSDP )
		{
			throw std::runtime_error( "TLE SDP propagator requested, which is not yet implemented." );
		}
		useSDP_ = useSDP;
		tle_ = tle_ptr;
		// Create table with ephemeris data
		// Create interpolator based on table

	}

    Eigen::Vector6d TleEphemeris::getCartesianStateInTemeFrame( double secondsSinceEpoch )
    {
        double ttSinceEpoch = secondsSinceEpoch - sofa_interface::getTDBminusTT( secondsSinceEpoch, Eigen::Vector3d::Zero( ));
        double utcSecondSinceEpoch = sofa_interface::convertTTtoUTC( ttSinceEpoch );

        // Call Spice interface to retrieve the spacecraft's state from the TLE in the True Equator, Mean Equinox frame (see
        // Vallado: Fundamentals of Astrodynamics and Applications 4th ed. (2013)). This frame is idiosyncratic in nature and
        // therefore needs to be converted to an intermediate standard reference frame.
        return spice_interface::getCartesianStateFromTleAtEpoch( utcSecondSinceEpoch, tle_ );
    }

	Eigen::Vector6d TleEphemeris::getCartesianState( double secondsSinceEpoch )
	{
        // Compute TEME state
        Eigen::Vector6d cartesianStateAtEpochTEME = getCartesianStateInTemeFrame( secondsSinceEpoch );

        // Compute rotation to J2000
        Eigen::Matrix3d rotationToJ2000 = getRotationMatrixFromTemeToJ2000( secondsSinceEpoch );

        Eigen::Vector6d cartesianStateOutput = Eigen::Vector6d::Zero( );
        if( referenceFrameOrientation_ == "J2000" )
		{
            cartesianStateOutput.segment( 0, 3 ) = rotationToJ2000 * cartesianStateAtEpochTEME.segment( 0, 3 );
            cartesianStateOutput.segment( 3, 3 ) = rotationToJ2000 * cartesianStateAtEpochTEME.segment( 3, 3 );
		}
		else if( referenceFrameOrientation_ == "ECLIPJ2000" )
		{
			Eigen::Matrix3d itrsToEclipticQuaternion = spice_interface::getRotationFromJ2000ToEclipJ2000( );

            cartesianStateOutput.segment( 0, 3 ) = itrsToEclipticQuaternion * rotationToJ2000 * cartesianStateAtEpochTEME.segment( 0, 3 );
            cartesianStateOutput.segment( 3, 3 ) = itrsToEclipticQuaternion * rotationToJ2000 * cartesianStateAtEpochTEME.segment( 3, 3 );
		}
		else
		{
			throw std::runtime_error( "TLE state conversion to target frame " + referenceFrameOrientation_ + " is currently unsupported." );
		}

		return cartesianStateOutput;
	}

	Tle::Tle( const std::string& lines )
	{
		// First of all, the TLE lines to be checked for validity: they shall not be empty, contain more than 69 characters
		// or have an invalid checksum.
		// Although unlikely, TLEs can become corrupted due to (un)intentional changes or data corruption during transmission.

		// Check if string is not empty
		if( lines.empty() )
		{
			throw std::runtime_error( "Error: TLE class was instantiated with string object, but string is empty." );
		}
		std::vector< std::string > tleLines;
		boost::algorithm::split( tleLines, lines, boost::is_any_of( "\n\r" ) );
		if( tleLines.size() != 2 )
		{
			throw std::runtime_error( "Error: TLE class was instantiated with string object, but string does not contain two 2 lines." );
		}
		// Check line length
		for( std::string line : tleLines )
		{
			if( line.length( ) != 69 )
			{
				std::cout << line << std::endl;
				std::string s( 1, line.at( 0 ) );
				throw std::runtime_error("Error: TLE class was instantiated with string object, but line " + s +
				" contains an invalid number of characters (counted " + std::to_string( line.length( ) ) + " characters)" );
			}
			// Checksum
			int checksum = 0;
			for( unsigned int i = 0; i < line.length( ) - 1; i++ )
			{
				char c = line[ i ];
				if( c == '-' )
				{
					checksum++;
				}
				else if( std::isdigit( c ) )
				{
					checksum += static_cast< int >( c ) - 48;
				}
				// Else do nothing
			}
			checksum %= 10;
			if(checksum == static_cast< int >( line[ 68 ] ) - 48 )
			{
				// std::cout << "TLE checksum verified." << std::endl;
			}
			else
			{
				throw std::runtime_error( "TLE checksum invalid!" );
			}
		}

		// Parse TLE
		std::string line1 = tleLines.at( 0 );
		std::string line2 = tleLines.at( 1 );

		int epochYear = std::stoi( line1.substr( 18, 2 ) );
		double epochDayFraction = std::stod( line1.substr( 20, 12 ) );

		// Convert to seconds since J2000
		// TLE day number starts with a 1, so a day fraction of 1.0 would mean Jan 1st, 0:00. Hence, we have to subtract 1.5 days
		// in seconds to obtain number of seconds since Jan 1st, noon (as dictated by J2000).
		if( epochYear < 57 )
		{
			epochYear += 2000;
		}
		else
		{
			epochYear += 1900;
		}
		// TLE day numbering starts with 1, whereas Tudat assumes January 1st to be number 0
        epoch_ = basic_astrodynamics::DateTime( epochYear, 1, 1, 0, 0, 0.0 ).epoch< Time >( ).getSeconds< double >( ) + ( epochDayFraction - 1.0 ) * physical_constants::JULIAN_DAY;

		double bStar = std::stod( line1.substr( 53, 6 ) );
		double bStarExp = std::stod( line1.substr( 59, 2 ) );
		bStar_ = bStar * std::pow( 10, bStarExp - 5 );

		// Convert angles to radians
		inclination_ = unit_conversions::convertDegreesToRadians( std::stod( line2.substr( 8, 7 ) ) );
		rightAscension_ = unit_conversions::convertDegreesToRadians( std::stod( line2.substr( 17, 8 ) ) );

		std::string eccentricityString = "0." + line2.substr( 26, 7 );
		eccentricity_ = std::stod( eccentricityString );

		argOfPerigee_ = unit_conversions::convertDegreesToRadians( std::stod( line2.substr( 34, 8 ) ) );
		meanAnomaly_ = unit_conversions::convertDegreesToRadians( std::stod( line2.substr( 43, 8 ) ) );

		// Convert mean motion to radians per min (as opposed to rev/day as given in TLE)
		double meanMotionRevDay = std::stod( line2.substr( 52, 11 ) );
		meanMotion_ = meanMotionRevDay * 2.0 * mathematical_constants::PI / ( 60 * 24 );
	}

	Tle::Tle( const std::string& tleLine1, const std::string& tleLine2 ): Tle( std::string( tleLine1 + "\n" + tleLine2 ) )
	{ }

    Tle::Tle( const double *spiceElements )
	{
		for( double &spiceElement : spiceElements_ )
		{
			spiceElement = *spiceElements;
			spiceElements++;
		}
	}

}
}
