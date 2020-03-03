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

	Eigen::Vector6d TleEphemeris::getCartesianState( const double secondsSinceEpoch )
	{
		const Eigen::Vector6d cartesianStateAtEpoch =
				spice_interface::getCartesianStateFromTleAtEpoch( secondsSinceEpoch, tle_ );

		return cartesianStateAtEpoch;
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

		double epochYear = std::stod( line1.substr( 18, 2 ) );
		double epochDayFraction = std::stod( line1.substr( 20, 12 ) );
		std::cout << "Epoch year: " << epochYear << "\nDay fraction: " << epochDayFraction << std::endl;

		// Convert to seconds since J2000
		// TLE day number starts with a 1, so a day fraction of 1.0 would mean Jan 1st, 0:00. Hence, we have to subtract 1.5 days
		// in seconds to obtain number of seconds since Jan 1st, noon (as dictated by J2000).
		epoch_ = epochYear * physical_constants::JULIAN_YEAR + ( epochDayFraction - 1.5 ) * physical_constants::JULIAN_DAY;

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
