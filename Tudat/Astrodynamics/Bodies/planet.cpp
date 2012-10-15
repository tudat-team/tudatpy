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
 *      110623    K. Kumar          Creation of code.
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *
 *    References
 *
 */

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Bodies/planet.h"

namespace tudat
{
namespace bodies
{

//! Set predefined planet settings.
void Planet::setPredefinedPlanetSettings( PredefinedPlanets predefinedPlanet )
{
    using ephemerides::Ephemeris;
    using ephemerides::ApproximatePlanetPositions;
    using astrodynamics::gravitation::CentralGravityField;

    // Select predefined planet.
    switch( predefinedPlanet )
    {
    case sun:

        // Set predefined Sun central gravity field.
        gravityFieldModel_ = boost::make_shared< CentralGravityField >( CentralGravityField::sun );

        break;

    case mercury:

        // Set predefined Mercury central gravity field.
        gravityFieldModel_ = boost::make_shared< CentralGravityField >(
                    CentralGravityField::mercury );

        // Set Mercury as planet for ephemeris.
        ephemeris_ = boost::make_shared< ApproximatePlanetPositions >(
                    ApproximatePlanetPositions::mercury );
        break;

    case venus:

        // Set predefined Venus central gravity field.
        gravityFieldModel_ = boost::make_shared< CentralGravityField >(
                    CentralGravityField::venus );

        // Set Venus as planet for ephemeris.
        ephemeris_ = boost::make_shared< ApproximatePlanetPositions >(
                    ApproximatePlanetPositions::venus );
        break;

    case earth:

        // Set predefined Earth central gravity field.
        gravityFieldModel_ = boost::make_shared< CentralGravityField >(
                    CentralGravityField::earth );

        // Set Earth as planet for ephemeris.
        ephemeris_ = boost::make_shared< ApproximatePlanetPositions >(
                    ApproximatePlanetPositions::earthMoonBarycenter );

        break;

    case moon:

        // Set predefined Moon central gravity field.
        gravityFieldModel_ = boost::make_shared< CentralGravityField >(
                    CentralGravityField::moon );

        // Set Moon as body for ephemeris.
        ephemeris_ = boost::make_shared< ApproximatePlanetPositions >(
                    ApproximatePlanetPositions::earthMoonBarycenter );
        break;

    case mars:

        // Set predefined Mars central gravity field.
        gravityFieldModel_ = boost::make_shared< CentralGravityField >(
                    CentralGravityField::mars );

        // Set Mars as planet for ephemeris.
        ephemeris_ = boost::make_shared< ApproximatePlanetPositions >(
                    ApproximatePlanetPositions::mars );
        break;

    case jupiter:

        // Set predefined Jupiter central gravity field.
        gravityFieldModel_ = boost::make_shared< CentralGravityField >(
                    CentralGravityField::jupiter );

        // Set Jupiter as planet for ephemeris.
        ephemeris_ = boost::make_shared< ApproximatePlanetPositions >(
                    ApproximatePlanetPositions::jupiter );
        break;

    case saturn:

        // Set predefined Saturn central gravity field.
        gravityFieldModel_ = boost::make_shared< CentralGravityField >(
                    CentralGravityField::saturn );

        // Set Saturn as planet for ephemeris.
        ephemeris_ = boost::make_shared< ApproximatePlanetPositions >(
                    ApproximatePlanetPositions::saturn );
        break;

    case uranus:

        // Set predefined Uranus central gravity field.
        gravityFieldModel_ = boost::make_shared< CentralGravityField >(
                    CentralGravityField::uranus );

        // Set Uranus as planet for ephemeris.
        ephemeris_ = boost::make_shared< ApproximatePlanetPositions >(
                    ApproximatePlanetPositions::uranus );
        break;

    case neptune:

        // Set predefined Neptune central gravity field.
        gravityFieldModel_ = boost::make_shared< CentralGravityField >(
                    CentralGravityField::neptune );

        // Set Neptune as planet for ephemeris.
        ephemeris_ = boost::make_shared< ApproximatePlanetPositions >(
                    ApproximatePlanetPositions::neptune );
        break;

    default:

        // Print cerr statement.
        std::cerr << "Desired predefined planet does not exist." << std::endl;
    };
}

} // namespace bodies
} // namespace tudat
