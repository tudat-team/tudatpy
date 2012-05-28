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
 *      100906    J. Melman         Creation of code.
 *      100910    J. Melman         No more switch statement and enum.
 *      100929    K. Kumar          Minor comment changes.
 *      110112    K. Kumar          Modified to use GravityFieldModel; corrected path.
 *      110113    K. Kumar          Added getGravityFieldModel( ) function.
 *      110310    K. Kumar          Added ephemeris; added missing destructor.
 *      120322    D. Dirkx          Moved Ephemeris member from here to Body
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.*
 *    References
 *
 */

#ifndef TUDAT_CELESTIAL_BODY_H
#define TUDAT_CELESTIAL_BODY_H

#include <iostream>

#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"
#include "Tudat/Astrodynamics/Bodies/body.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"

namespace tudat
{
namespace bodies
{

//! Celestial body class.
/*!
 * Celestial body class.
 */
class CelestialBody : public Body
{
public:

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~CelestialBody( ) { }

    //! Set gravity field model.
    /*!
     * Sets the gravity field model.
     * \param gravityFieldModel Pointer to a gravity field model.
     */
    void setGravityFieldModel(
            boost::shared_ptr< astrodynamics::gravitation::GravityFieldModel > gravityFieldModel )
    {
        gravityFieldModel_ = gravityFieldModel;
    }

    //! Set atmosphere model.
    /*!
     * Sets the atmosphere model.
     * \param atmosphereModel Pointer to an atmosphere model.
     */
    void setAtmosphereModel( boost::shared_ptr< AtmosphereModel > atmosphereModel )
    {
        atmosphereModel_ = atmosphereModel;
    }

    //! Get gravity field model.
    /*!
     * Returns the gravity field model.
     * \return Pointer to gravity field model.
     */
    boost::shared_ptr< astrodynamics::gravitation::GravityFieldModel > getGravityFieldModel( )
    {
        return gravityFieldModel_;
    }

    //! Get atmosphere model.
    /*!
     * Returns the atmosphere model.
     * \return Pointer to the atmosphere model.
     */
    boost::shared_ptr< AtmosphereModel > getAtmospheremodel( ) { return atmosphereModel_; }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param celestialBody Celestial body.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, CelestialBody& celestialBody )
    {
        stream << "This is a CelestialBody object. The gravitational parameter is set to: "
               << celestialBody.getGravityFieldModel( )->getGravitationalParameter( ) << std::endl;
        return stream;
    }

protected:

    //! Pointer to gravity field model.
    /*!
     * Pointer to gravity field model.
     */
    boost::shared_ptr< astrodynamics::gravitation::GravityFieldModel > gravityFieldModel_;

    //! Pointer to atmosphere model.
    /*!
     * Pointer to atmosphere model.
     */
    boost::shared_ptr< AtmosphereModel > atmosphereModel_;

private:
};

} // namespace bodies
} // namespace tudat

#endif // TUDAT_CELESTIAL_BODY_H
