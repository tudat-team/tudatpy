/*! \file approximatePlanetPositionsDataContainer.cpp
 *    This source file contains the definition of a data container class for
 *    data extracted from the JPL "Approximate Positions of Major Planets"
 *    ephemeris for a specific planet ( http://ssd.jpl.nasa.gov/?planet_pos ).
 *    The ephemeris file used is for the period 3000 BC to 3000 AD.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 24 February, 2011
 *    Last modified     : 24 February, 2011
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the
 *          Major Planets, http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf,
 *          last accessed: 24 February, 2011.
 *
 *    Notes
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
 *      110224    K. Kumar          First creation of code.
 */

// Include statements.
#include "approximatePlanetPositionsDataContainer.h"

// Using declarations.
using std::endl;

//! Default constructor.
ApproximatePlanetPositionsDataContainer::
        ApproximatePlanetPositionsDataContainer( )
            : semiMajorAxis_( -0.0 ), eccentricity_( -0.0 ),
              inclination_( -0.0 ), meanLongitude_( -0.0 ),
              longitudeOfPerihelion_( -0.0 ),
              longitudeOfAscendingNode_( -0.0 ),
              rateOfChangeOfSemiMajorAxis_( -0.0 ),
              rateOfChangeOfEccentricity_( -0.0 ),
              rateOfChangeOfInclination_( -0.0 ),
              rateOfChangeOfMeanLongitude_( -0.0 ),
              rateOfChangeOfLongitudeOfPerihelion_( -0.0 ),
              rateOfChangeOfLongitudeOfAscendingNode_( -0.0 ),
              additionalTermB_( -0.0 ), additionalTermC_( -0.0 ),
              additionalTermS_( -0.0 ), additionalTermF_( -0.0 )
{
}

//! Default destructor.
ApproximatePlanetPositionsDataContainer::
        ~ApproximatePlanetPositionsDataContainer( )
{
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          ApproximatePlanetPositionsDataContainer&
                          approximatePlanetPositionsDataContainer )
{
    stream << "This is an ApproximatePlanetPositionsDataContainer object. "
           << endl;
    stream << "The data corresponds to the table entry for "
           << approximatePlanetPositionsDataContainer.planetName_
           << endl;
    stream << "The semi-major axis in AU is set to: "
           << approximatePlanetPositionsDataContainer.semiMajorAxis_
           << endl;
    stream << "The eccentricity in radians is set to: "
           << approximatePlanetPositionsDataContainer.eccentricity_
           << endl;
    stream << "The inclination in degrees is set to: "
           << approximatePlanetPositionsDataContainer.inclination_
           << endl;
    stream << "The mean longitude in degrees is set to: "
           << approximatePlanetPositionsDataContainer.meanLongitude_
           << endl;
    stream << "The longitude of perihelion in degrees is set to: "
           << approximatePlanetPositionsDataContainer
                .longitudeOfPerihelion_ << endl;
    stream << "The longitude of the ascending node in degrees is set to: "
           << approximatePlanetPositionsDataContainer
                .longitudeOfAscendingNode_<< endl;
    stream << "The rate of change of semi-major axis AU per century is set "
           << "to: " << approximatePlanetPositionsDataContainer
                          .rateOfChangeOfSemiMajorAxis_ << endl;
    stream << "The rate of change of eccentricity in radians per century is "
           << "set to: " << approximatePlanetPositionsDataContainer
                              .rateOfChangeOfEccentricity_ << endl;
    stream << "The rate of change of inclination in degrees per century is "
           << "set to: " << approximatePlanetPositionsDataContainer
                              .rateOfChangeOfInclination_ << endl;
    stream << "The rate of change of mean longitude in degrees per century is "
           << "set to: " << approximatePlanetPositionsDataContainer
                              .rateOfChangeOfMeanLongitude_ << endl;
    stream << "The rate of change of longitude of perihelion in degrees per "
           << "century is set to: "
           << approximatePlanetPositionsDataContainer
                .rateOfChangeOfLongitudeOfPerihelion_ << endl;
    stream << "The rate of change of longitude of ascending node in degrees "
           << "per century is set to: "
           << approximatePlanetPositionsDataContainer
                .rateOfChangeOfLongitudeOfAscendingNode_ << endl;

    // Check if additional terms are defined for outer planets.
    if ( approximatePlanetPositionsDataContainer.additionalTermB_ != -0.0 )
    {
        stream << "The additional term B is set to: "
               << approximatePlanetPositionsDataContainer
                    .additionalTermB_ << endl;
        stream << "The additional term C is set to: "
               << approximatePlanetPositionsDataContainer
                    .additionalTermC_ << endl;
        stream << "The additional term S is set to: "
               << approximatePlanetPositionsDataContainer
                    .additionalTermS_ << endl;
        stream << "The additional term F is set to: "
               << approximatePlanetPositionsDataContainer
                    .additionalTermF_ << endl;
    }

    // Return stream.
    return stream;
}

// End of file.
