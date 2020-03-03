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

#ifndef TUDATBUNDLE_TLEEPHEMERIS_H
#define TUDATBUNDLE_TLEEPHEMERIS_H

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

namespace tudat
{

namespace ephemerides
{

class Tle
{
public:

	//! Default (empty) constructor
	Tle( ) = default;

	//! Cosntructor to create TLE object from a string containing the two lines for the relevant satellite.
	explicit Tle( const std::string& tleLines );

	Tle( const std::string& tleLine1, const std::string& tleLine2 );

	//explicit Tle( const std::string& satelliteName, const std::string& url = "https://celestrak.com/NORAD/elements/active.txt" ){ };

	explicit Tle( const double *spiceElements );

	//! Constructor
	Tle( const double epoch, const double bStar, const double inclination, const double rightAscension,
		 const double eccentricity, const double argOfPerigee, const double meanAnomaly,
		 const double meanMotion ):
			epoch_( epoch ), bStar_( bStar ), inclination_( inclination ), rightAscension_( rightAscension ),
			eccentricity_( eccentricity ), argOfPerigee_( argOfPerigee ), meanAnomaly_( meanAnomaly ),
			meanMotion_( meanMotion ) { };

	double getEpoch() const
	{
		return epoch_;
	}

	double getBStar() const
	{
		return bStar_;
	}

	double getInclination() const
	{
		return inclination_;
	}

	double getRightAscension() const
	{
		return rightAscension_;
	}

	double getEccentricity() const
	{
		return eccentricity_;
	}

	double getArgOfPerigee() const
	{
		return argOfPerigee_;
	}

	double getMeanAnomaly() const
	{
		return meanAnomaly_;
	}

	double getMeanMotion() const
	{
		return meanMotion_;
	}

private:

	double epoch_;
	double bStar_;
	double inclination_;
	double rightAscension_;
	double eccentricity_;
	double argOfPerigee_;
	double meanAnomaly_;
	double meanMotion_;

	// Declare two-line elements array that is compatible with CSpice
	// I.e. in units and order
	double spiceElements_[ 10 ];

};

class TleEphemeris : public Ephemeris
{
public:

	//! Constructor
	/*!
	 *
	 * @param referenceFrameOrigin
	 * @param referenceFrameOrientation
	 * @param useSDP
	 */
	TleEphemeris( const std::string& referenceFrameOrigin,
				  const std::string& referenceFrameOrientation,
				  const std::shared_ptr< Tle > tle, const bool useSDP = false);

	Eigen::Vector6d getCartesianState(const double secondsSinceEpoch) override;

private:

	std::shared_ptr< Tle > tle_;
	bool useSDP_;

};

}
}
#endif //TUDATBUNDLE_TLEEPHEMERIS_H
