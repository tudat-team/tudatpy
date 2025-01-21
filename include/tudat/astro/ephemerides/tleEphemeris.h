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

#ifndef TUDAT_TLEEPHEMERIS_H
#define TUDAT_TLEEPHEMERIS_H

#include "tudat/astro/ephemerides/ephemeris.h"

namespace tudat
{

namespace ephemerides
{

Eigen::Matrix3d getRotationMatrixFromTemeToJ2000( const double epochSinceJ2000  );

Eigen::Matrix3d getRotationMatrixFromJ2000ToTeme( const double epochSinceJ2000  );

//! Class holding data for one set of two-line elements (TLE) for a specific Earth-orbiting satellite.
/*!
 * Class holding data for one set of two-line elements (TLE) for a specific Earth-orbiting satellite. It is valid at its epoch plus/minus
 * approximately 1 week, depending on the required accuracy of the predictions. The class offers different options for construction:
 * from one string holding both lines; from two strings, each holding one line; from an array of doubles containing Spice-compatible
 * elements; and an explicit constructor from the elements themselves.
 * This class is used in conjunction with the TleEphemeris class, which provides position and velocity predictions based on the SGP and SDP
 * propagators which are used for TLE propagation.
 */
class Tle
{
public:

	//! Default (empty) constructor
	Tle( ) = default;

	//! Constructor to create TLE object from a string containing the two lines for the relevant satellite.
	/*!
	 * Class constructor to generate a TLE object based on a string containing the two element lines for the relevant satellite.
	 * @param tleLines String containing the TLE lines
	 */
	explicit Tle( const std::string& tleLines );

	/*!
	 * Class constructor to generate a TLE object from two separate string objects, each containg a TLE line.
	 * @param tleLine1 String containing the first TLE line.
	 * @param tleLine2 String containing the second TLE line.
	 */
	Tle( const std::string& tleLine1, const std::string& tleLine2 );

	//explicit Tle( const std::string& satelliteName, const std::string& url = "https://celestrak.com/NORAD/elements/active.txt" ){ };

	/*!
	 * Class constructor to generate a TLE object directly from Spice-compatible elements, given as a pointer/array. Mostly intended for debugging
	 * purposes.
	 * @param spiceElements Array of doubles containing the Spice-compatible elements (see Spice documentation for more information).
	 */
	explicit Tle( const double *spiceElements );

	/*!
	 * Class constructor to generate a TLE object directly from standard two-line elements. Intended for debugging purposes.
	 * @param epoch TLE set epoch in seconds from J2000.
	 * @param bStar B-Star coefficient.
	 * @param inclination Inclination of the orbit in radians.
	 * @param rightAscension Right ascension of the orbit in radians.
	 * @param eccentricity Eccentricty of the orbit.
	 * @param argOfPerigee Argument of perigee of the orbit in radians.
	 * @param meanAnomaly Mean anomaly in radians.
	 * @param meanMotion Mean motion in radians/minute (!).
	 */
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

//! Ephemeris derived class that calculates the Cartesian position as a function of time assuming
//! a two-line elements based orbit.
/*!
 *  Ephemeris derived class that calculates the Cartesian position as a function of time assuming
 *  a two-line elements based orbit. The class supports state propagation for Earth-orbiting
 *  satellites, given the set of two-line elements (TLE) that form the parameters for the SGP
 *  and SDP models. Currently, the class supports target reference frames with their origin at
 *  the center of the Earth, and either J2000 or ECLIPJ2000 orientation.
 */
class TleEphemeris : public Ephemeris
{
public:

	//! Class constructor.
	/*!
	 * Class constructor for a TLE ephemeris based on a set of two-line elements.
	 * @param referenceFrameOrigin Origin of the target frame of reference (default: Earth). Currently (v1.0) any value other than Earth is
	 * unsupported.
	 * @param referenceFrameOrientation Orientation of the target frame of reference (default: J2000).
	 * @param tle Smart pointer to a TLE object from which the ephemeris is to be generated.
	 * @param useSDP Boolean specifying whether the near-Earth SGP (false) or deep-space (SDP) model (true) should be used.
	 */
	TleEphemeris( const std::string& referenceFrameOrigin = "Earth",
				  const std::string& referenceFrameOrientation = "J2000",
				  const std::shared_ptr< Tle > tle = nullptr, const bool useSDP = false);

    Eigen::Vector6d getCartesianStateInTemeFrame( double secondsSinceEpoch );

    //! Function to get state from ephemeris.
	/*!
	 *  Returns state from ephemeris at given time, given a two-line element set using either SGP4 or SDP4.
	 *  \param secondsSinceEpoch Seconds since J2000 epoch at which ephemeris is to be evaluated.
	 *  \return TLE-derived orbit Cartesian state at given time.
	 */
	Eigen::Vector6d getCartesianState(const double secondsSinceEpoch) override;

private:

	std::shared_ptr< Tle > tle_;
	bool useSDP_;

};

}
}
#endif //TUDAT_TLEEPHEMERIS_H
