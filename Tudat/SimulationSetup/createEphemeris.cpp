/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      150501    D. Dirkx          Ported from personal code
 *
 *    References
 *
 *    Notes
 *
 */

#include <boost/lambda/lambda.hpp>
#include <boost/lexical_cast.hpp>

#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositionsCircularCoplanar.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/SimulationSetup/createEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;

//! Function to create a tabulated ephemeris using data from Spice.
boost::shared_ptr< Ephemeris > createTabulatedEphemerisFromSpice(
        const std::string& body,
        const double initialTime,
        const double endTime,
        const double timeStep,
        const std::string& observerName,
        const std::string& referenceFrameName )
{
    using namespace interpolators;

    std::map< double, basic_mathematics::Vector6d > timeHistoryOfState;

    // Calculate state from spice at given time intervals and store in timeHistoryOfState.
    double currentTime = initialTime;
    while( currentTime < endTime )
    {
        timeHistoryOfState[ currentTime ] = spice_interface::getBodyCartesianStateAtEpoch(
                    body, observerName, referenceFrameName, "none", currentTime );
        currentTime += timeStep;
    }

    // Create interpolator.
    boost::shared_ptr< LagrangeInterpolator< double, basic_mathematics::Vector6d > > interpolator =
            boost::make_shared< LagrangeInterpolator< double, basic_mathematics::Vector6d > >(
                timeHistoryOfState, 6, huntingAlgorithm,
                lagrange_cubic_spline_boundary_interpolation );

    // Create ephemeris and return.
    return boost::make_shared< TabulatedCartesianEphemeris< > >(
                interpolator, observerName, referenceFrameName );
}
//! Function to create a ephemeris model.
boost::shared_ptr< ephemerides::Ephemeris > createBodyEphemeris(
        const boost::shared_ptr< EphemerisSettings > ephemerisSettings,
        const std::string& bodyName )
{
    // Declare return object.
    boost::shared_ptr< ephemerides::Ephemeris > ephemeris;

    // Check which type of ephemeris model is to be created.
    switch( ephemerisSettings->getEphemerisType( ) )
    {
    case direct_spice_ephemeris:
    {
        // Check consistency of type and class.
        boost::shared_ptr< DirectSpiceEphemerisSettings > directEphemerisSettings =
                boost::dynamic_pointer_cast< DirectSpiceEphemerisSettings >( ephemerisSettings );
        if( directEphemerisSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected direct spice ephemeris settings for body " + bodyName );
        }
        else
        {
            // Create corresponding ephemeris object.
            ephemeris = boost::make_shared< SpiceEphemeris >(
                        bodyName,
                        directEphemerisSettings->getFrameOrigin( ),
                        directEphemerisSettings->getCorrectForStellarAbberation( ),
                        directEphemerisSettings->getCorrectForLightTimeAbberation( ),
                        directEphemerisSettings->getConvergeLighTimeAbberation( ),
                        directEphemerisSettings->getFrameOrientation( ) );
        }
        break;
    }
    case interpolated_spice:
    {
        // Check consistency of type and class.
        boost::shared_ptr< InterpolatedSpiceEphemerisSettings > interpolatedEphemerisSettings =
                boost::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >(
                    ephemerisSettings );
        if( interpolatedEphemerisSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected tabulated spice ephemeris settings for body " + bodyName );
        }
        else
        {
            // Since only the barycenters of planetary systems are included in the standard DE
            // ephemerides, append 'Barycenter' to body name.
            std::string inputName;
            inputName = bodyName;
            if( bodyName == "Mercury" || bodyName == "Venus" || bodyName == "Mars" ||
                    bodyName == "Jupiter"  || bodyName == "Saturn" ||
                    bodyName == "Uranus" || bodyName == "Neptune" )
            {
                inputName += " Barycenter";
                std::cerr<<"Warning, position of "<<bodyName<<" taken as baycenter of that body's "
                        <<"planetary system."<<std::endl;
            }

            // Create corresponding ephemeris object.
            ephemeris = createTabulatedEphemerisFromSpice(
                        inputName,
                        interpolatedEphemerisSettings->getInitialTime( ),
                        interpolatedEphemerisSettings->getFinalTime( ),
                        interpolatedEphemerisSettings->getTimeStep( ),
                        interpolatedEphemerisSettings->getFrameOrigin( ),
                        interpolatedEphemerisSettings->getFrameOrientation( ) );
        }
        break;
    }
    case tabulated_ephemeris:
    {
        // Check consistency of type and class.
        boost::shared_ptr< TabulatedEphemerisSettings > tabulatedEphemerisSettings =
                boost::dynamic_pointer_cast< TabulatedEphemerisSettings >( ephemerisSettings );
        if( tabulatedEphemerisSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected tabulated ephemeris settings for body " + bodyName );
        }
        else
        {
            // Create corresponding ephemeris object.
            ephemeris = boost::make_shared< TabulatedCartesianEphemeris< > >(
                        boost::make_shared<
                        interpolators::LagrangeInterpolator< double, basic_mathematics::Vector6d > >
                        ( tabulatedEphemerisSettings->getBodyStateHistory( ), 6,
                          interpolators::huntingAlgorithm,
                          interpolators::lagrange_cubic_spline_boundary_interpolation ),
                        tabulatedEphemerisSettings->getFrameOrigin( ),
                        tabulatedEphemerisSettings->getFrameOrientation( ) );
        }
        break;
    }
    case approximate_planet_positions:
    {
        // Check consistency of type and class.
        boost::shared_ptr< ApproximatePlanetPositionSettings > approximateEphemerisSettings =
                boost::dynamic_pointer_cast< ApproximatePlanetPositionSettings >(
                    ephemerisSettings );
        if( approximateEphemerisSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected approximate ephemeris settings for body " + bodyName );
        }
        else
        {
            // Create corresponding ephemeris object.
            if( approximateEphemerisSettings->getUseCircularCoplanarApproximation( ) )
            {
                ephemeris = boost::make_shared< ApproximatePlanetPositionsCircularCoplanar >(
                            approximateEphemerisSettings->getBodyIdentifier( ) );
            }
            else
            {
                ephemeris = boost::make_shared< ApproximatePlanetPositions >(
                            approximateEphemerisSettings->getBodyIdentifier( ) );
            }
        }
        break;
    }
    default:
    {
        throw std::runtime_error(
                    "Error, did not recognize ephemeris model settings type " +
                    boost::lexical_cast< std::string >( ephemerisSettings->getEphemerisType( ) ) );
    }
    }
    return ephemeris;

}

}

}
