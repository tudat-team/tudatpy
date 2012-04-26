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
 *      100906    J. Melman         First creation of code.
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
