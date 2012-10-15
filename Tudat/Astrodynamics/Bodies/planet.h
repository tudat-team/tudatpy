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

#ifndef TUDAT_PLANET_H
#define TUDAT_PLANET_H

#include "Tudat/Astrodynamics/Bodies/celestialBody.h"
#include "Tudat/Astrodynamics/Bodies/Ephemeris/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityField.h"

namespace tudat
{
namespace bodies
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
    {
        sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune
    };

    //! Set predefined planet settings.
    /*!
     * Sets predefined planet settings.
     * \param predefinedPlanet Predefined planet.
     */
    void setPredefinedPlanetSettings( PredefinedPlanets predefinedPlanet );

protected:

private:
};

} // namespace bodies
} // namespace tudat

#endif // TUDAT_PLANET_H
