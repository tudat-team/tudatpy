/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      110623    K. Kumar          First creation of code.
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
