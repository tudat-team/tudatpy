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
 *
 *    References
 *
 */

#ifndef TUDAT_PLANET_H
#define TUDAT_PLANET_H

#include "Tudat/Astrodynamics/Bodies/celestialBody.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityField.h"
#include "Tudat/Astrodynamics/Bodies/Ephemeris/approximatePlanetPositions.h"

namespace tudat
{

//! Planet class.
/*!
 * Class that contains planet properties.
 */
class Planet : public CelestialBody
{
public:

    //! Predefined planets
    /*!
     * Predefined planets.
     */
    enum PredefinedPlanets
    { sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune };

    //! Set predefined planet settings.
    /*!
     * Sets predefined planet settings.
     * \param predefinedPlanet Predefined planet.
     */
    void setPredefinedPlanetSettings( PredefinedPlanets predefinedPlanet );

protected:

private:

    //! Predefined central gravity field.
    /*!
     * Predefined central gravity field.
     */
    CentralGravityField predefinedCentralGravityField_;

    //! Approximate planet positions ephemeris.
    /*!
     * JPL apprximate planet positions ephemeris.
     */
    ApproximatePlanetPositions approximatePlanetPositions_;
};

} // namespace tudat

#endif // TUDAT_PLANET_H
