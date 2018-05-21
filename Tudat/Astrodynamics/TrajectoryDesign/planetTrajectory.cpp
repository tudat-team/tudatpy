/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120606    P. Musegaas       First creation of code.
 *      120827    P. Musegaas       Adaptation to own ephemeris type.
 *
 *    References
 *
 */

// The ephemeris extractor now takes planet ephemeris as only means of extracting the data. This
//      could be improved.

#include <cmath>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "exportTrajectory.h"
#include "planetTrajectory.h"

//#include "Temporary/temporaryFunctions.h"

namespace tudat
{
namespace spaceTrajectories
{

//! Return a vector of positions and times from ephemeris data for a certain epoch and duration.
void returnPlanetTrajectory(const ephemerides::EphemerisPointer &ephemerisPtr,
                             const double centralBodyGravitationalParameter,
                             const double startingEpochMJD2000,
                             const double duration,
                             const double maximumTimeStep,
                             std::vector < Eigen::Vector3d >& positionVector,
                             std::vector < double >& timeVector,
                             const double startingTime )
{
    // Calculate the number of intermediate steps required and the corresponding time step.
    int numberOfSteps = std::ceil( duration / maximumTimeStep );
    double timeStep = duration / numberOfSteps;

    // Resize the position and time vectors. Note that the "+ 1" comes from the end points.
    positionVector.resize( numberOfSteps + 1 );
    timeVector.resize( numberOfSteps + 1 );

    // Initiate variables to temporarily store the ephemeris data.
    Eigen::Vector6d temporaryKeplerianElements, temporaryCartesianElements;

    // Initiate time variable in JD for extracting data, and a time variable in seconds to store
    // the time in the vector.
    double MJD2000 = 51544.5;
    double timeMJD2000 = startingEpochMJD2000 + MJD2000;
    double time = 0.0;

    boost::shared_ptr<ephemerides::ApproximatePlanetPositions> approxEphemerisPtr =
            boost::dynamic_pointer_cast<ephemerides::ApproximatePlanetPositions>(ephemerisPtr);
    double timeJD2000;
    double timeSecondsSinceEpoch;
    // Start the for loop to obtain positions along the trajectory.
    for ( int counter = 0; counter < numberOfSteps + 1; counter++ )
    {
        timeJD2000 = basic_astrodynamics::convertModifiedJulianDayToJulianDay(timeMJD2000);
        timeSecondsSinceEpoch = basic_astrodynamics::convertJulianDayToSecondsSinceEpoch(timeJD2000);
        // Get Keplerian state of the planet at the corresponding time.
        temporaryKeplerianElements = ( *approxEphemerisPtr ).getKeplerianStateFromEphemeris( timeSecondsSinceEpoch );

        // Convert into cartesian elements.
        temporaryCartesianElements = orbital_element_conversions::
                convertKeplerianToCartesianElements( temporaryKeplerianElements,
                                                     centralBodyGravitationalParameter );

        // Store the position and time
        positionVector[ counter ] = temporaryCartesianElements.segment( 0, 3 );
        timeVector[ counter ] = time + startingTime;

        // Update the times.
        time += timeStep;
        timeMJD2000 += timeStep / physical_constants::JULIAN_DAY;
    }
}

void returnSingleRevolutionPlanetTrajectory(const ephemerides::EphemerisPointer &ephemerisPtr,
        const double centralBodyGravitationalParameter,
        const double startingEpochMJD2000,
        const double maximumTimeStep,
        std::vector < Eigen::Vector3d >& positionVector,
        std::vector < double >& timeVector,
        const double startingTime )
{
    // Initiate variables to temporarily store the ephemeris data.
    Eigen::Vector6d temporaryKeplerianElements, initialCartesianElements;

    // Convert MJD2000 to seconds since epoch
    double startingEpochMJD = startingEpochMJD2000 + 51544.5;
    double timeJD2000 = basic_astrodynamics::convertModifiedJulianDayToJulianDay(startingEpochMJD2000);
    double timeSecondsSinceEpoch = basic_astrodynamics::convertJulianDayToSecondsSinceEpoch(timeJD2000);

    // Get Cartesian state of the planet at the corresponding time.
    initialCartesianElements = ephemerisPtr->getCartesianState( timeSecondsSinceEpoch );

    // Compute duration (one revolution).
    double period = 2 * mathematical_constants::PI / std::sqrt( centralBodyGravitationalParameter /
            std::pow( temporaryKeplerianElements( 0 ), 3 ) );

    // Compute position and time vector, by pasing computed data to other function.
    returnTrajectory( initialCartesianElements, centralBodyGravitationalParameter, period,
                      maximumTimeStep, positionVector, timeVector, startingTime );
}

} // namespace spaceTrajectories
} // namespace tudat

