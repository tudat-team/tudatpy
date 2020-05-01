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

#include "Tudat/Astrodynamics/Ephemerides/tleEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/External/SofaInterface/earthOrientation.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"
#include "boost/algorithm/string.hpp"

namespace tudat
{

namespace ephemerides
{

	TleEphemeris::TleEphemeris( const std::string &referenceFrameOrigin, const std::string &referenceFrameOrientation,
							   const std::shared_ptr< Tle > tle_ptr, const bool useSDP ) :
							   Ephemeris( referenceFrameOrigin, referenceFrameOrientation )
	{
		useSDP_ = useSDP;
		tle_ = tle_ptr;
		// Create table with ephemeris data
		// Create interpolator based on table

	}

	Eigen::Vector6d TleEphemeris::getCartesianState( double secondsSinceEpoch )
	{
		// Call Spice interface to retrieve the spacecraft's state from the TLE in the True Equator, Mean Equinox frame (see
		// Vallado: Fundamentals of Astrodynamics and Applications 4th ed. (2013)). This frame is idiosyncratic in nature and
		// therefore needs to be converted to an intermediate standard reference frame.
		const Eigen::Vector6d cartesianStateAtEpochTEME =
				spice_interface::getCartesianStateFromTleAtEpoch( secondsSinceEpoch, tle_ );

		std::cout << "State at TLE epoch: " << spice_interface::getCartesianStateFromTleAtEpoch( tle_->getEpoch( ), tle_ ) << std::endl;

		Eigen::Vector3d positionTEME = cartesianStateAtEpochTEME.head( 3 );
		std::cout << "Position in TEME:\n" << positionTEME << std::endl;
		Eigen::Vector3d velocityTEME = cartesianStateAtEpochTEME.tail( 3 );

		// First, rotate to the True Of Date (TOD) frame.
		double equationOfEquinoxes = sofa_interface::calculateEquationOfEquinoxes( secondsSinceEpoch );
		std::cout << "Eq of the equinoxes: " << equationOfEquinoxes << std::endl;

		// Rotate around pole (z-axis)
		Eigen::AngleAxisd rotationObject = Eigen::AngleAxisd( equationOfEquinoxes, Eigen::Vector3d::UnitZ( ) );
		Eigen::Vector3d positionTOD = rotationObject.toRotationMatrix( ) * positionTEME;
		Eigen::Vector3d velocityTOD = rotationObject.toRotationMatrix( ) * velocityTEME;

		std::cout << "Position in TOD:\n" << positionTOD << std::endl;
		std::cout << "Velocity in TOD:\n" << velocityTOD << std::endl;

		double zeta, z, theta;
		sofa_interface::getPrecessionAngles(zeta, z, theta, secondsSinceEpoch );
		std::cout << "Zeta: " << zeta << " , z: " << z << " , theta: " << theta << std::endl;

		// Now that we have our state vector in the TOD frame, we need to obtain the combined precession + nutation matrix from Sofa
		// (according to the 1976/1980 model)
		Eigen::Matrix3d precessionNutationMatrix = sofa_interface::getPrecessionNutationMatrix( secondsSinceEpoch );
		// Multiply by inverted matrix to get to J2000
		Eigen::Vector3d  positionJ2000 = precessionNutationMatrix.transpose( ) * positionTOD;
		Eigen::Vector3d  velocityJ2000 = precessionNutationMatrix.transpose( ) * velocityTOD;

		// Get hour angle (theta GMST) to rotate to PEF
		// double thetaGmst = sofa_interface::calculateGreenwichMeanSiderealTime( timeUTC, timeUT1,
		// 		tle_->getEpoch(), basic_astrodynamics::iau_2000_b );

		// TODO: convert state from TEME frame to frame set by the ephemeris settings
		if( referenceFrameOrientation_ == "J2000" )
		{
			Eigen::Vector6d stateJ2000;
			stateJ2000 << positionJ2000, velocityJ2000;
			return stateJ2000;
		}
		else if( referenceFrameOrientation_ == "ECLIPJ2000" )
		{
			Eigen::Quaterniond itrsToEclipticQuaternion = spice_interface::computeRotationQuaternionBetweenFrames(
					"J2000", "ECLIPJ2000", secondsSinceEpoch );
			Eigen::Vector3d positionEclipJ2000 = itrsToEclipticQuaternion * positionJ2000;
			Eigen::Vector3d velocityEclipJ2000 = itrsToEclipticQuaternion * velocityJ2000;
			Eigen::Vector6d stateEclipJ2000;
			stateEclipJ2000 << positionEclipJ2000, velocityEclipJ2000;
			return stateEclipJ2000;
		}
		else
		{
			throw std::runtime_error( "TLE state conversion to target frame " + referenceFrameOrientation_ + " is currently unsupported." );
		}

		return cartesianStateAtEpochTEME;
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
		boost::algorithm::split( tleLines, lines, boost::is_any_of( "\n" ) );
		if( tleLines.size() != 2 )
		{
			throw std::runtime_error( "Error: TLE class was instantiated with string object, but string contains more than 2 lines." );
		}
		// Check line length
		for( std::string& line : tleLines )
		{
			if( line.length() != 69 )
			{
				throw std::runtime_error("Error: TLE class was instantiated with string object, but one or more lines contain an invalid "
							 "number of characters." );
			}
			// TODO: checksum verification
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
				std::cout << "TLE checksum verified." << std::endl;
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
		std::cout << "Epoch year: " << epochYear << "\nDay fraction: " << epochDayFraction << std::endl;

		// Convert to seconds since J2000
		// TLE day number starts with a 1, so a day fraction of 1.0 would mean Jan 1st, 0:00. Hence, we have to subtract 1.5 days
		// in seconds to obtain number of seconds since Jan 1st, noon (as dictated by J2000).
		// TODO: fix calculation (does not account for leap days!)
		if( epochYear < 57 )
		{
			epochYear += 2000;
		}
		else
		{
			epochYear += 1900;
		}
		// TLE day numbering starts with 1, whereas Tudat assumes January 1st to be number 0
		boost::gregorian::date date = basic_astrodynamics::convertYearAndDaysInYearToDate( epochYear, std::floor( epochDayFraction ) - 1 );
		auto jdSinceJ2000 = basic_astrodynamics::calculateJulianDaySinceEpoch< double >( date, epochDayFraction );
		epoch_ = (epochYear - 2000) * physical_constants::JULIAN_YEAR + ( epochDayFraction - 1.5 ) * physical_constants::JULIAN_DAY;

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

	Tle::Tle( const std::string& tleLine1, const std::string& tleLine2 )
	{

	}

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
