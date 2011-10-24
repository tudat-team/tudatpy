/*! \file celestialBody.h
 *    This file contains a class that describes celestial bodies with all their characteristics.
 *
 *    Path              : /Astrodynamics/Bodies/CelestialBodies/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 6 September, 2010
 *    Last modified     : 10 March, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110113    K. Kumar          Added getGravityFieldModel() function.
 *      110310    K. Kumar          Added ephemeris; added missing destructor.
 */

#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H

// Include statements.
#include <iostream>
#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/EnvironmentModels/GravityFieldModel/gravityFieldModel.h"
#include "Astrodynamics/EnvironmentModels/AtmosphereModel/atmosphereModel.h"
#include "Astrodynamics/States/cartesianElements.h"
#include "Astrodynamics/States/ephemeris.h"

//! Celestial body class.
/*!
 * Celestial body class.
 */
class CelestialBody : public Body
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    CelestialBody( ) : pointerToEphemeris_( NULL ), pointerToGravityFieldModel_( NULL ),
        pointerToAtmosphereModel_( NULL ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~CelestialBody( ) { }

    //! Set ephemeris.
    /*!
     * Sets the ephemeris.
     * \param pointerToEphemeris Pointer to ephemeris.
     */
    void setEphemeris( Ephemeris* pointerToEphemeris )
    { pointerToEphemeris_ = pointerToEphemeris; }

    //! Set gravity field model.
    /*!
     * Sets the gravity field model.
     * \param pointerToGravityFieldModel Pointer to gravity field model.
     */
    void setGravityFieldModel( GravityFieldModel* pointerToGravityFieldModel )
    { pointerToGravityFieldModel_ = pointerToGravityFieldModel; }

    //! Set atmosphere model.
    /*!
     * Sets the atmosphere model.
     * \param pointerToAtmosphereModel Pointer to an atmosphere model.
     */
    void setAtmosphereModel( AtmosphereModel* pointerToAtmosphereModel )
    { pointerToAtmosphereModel_ = pointerToAtmosphereModel; }

    //! Get gravitational parameter.
    /*!
     * Returns the gravitational parameter in meter^3 per second^2.
     * \return Gravitational parameter.
     */
    const double getGravitationalParameter( ) const
    { return pointerToGravityFieldModel_->getGravitationalParameter( ); }

    //! Get state from ephemeris at given Julian date.
    /*!
     * Returns the state of the celestial body from the defined ephemeris at
     * the given Julian date in Cartesian elements.
     * \return Pointer to Cartesian elements.
     */
    CartesianElements* getStateFromEphemeris( const double& julianDate )
    { return pointerToEphemeris_->getStateFromEphemeris( julianDate ); }

    //! Get ephemeris.
    /*!
     * Gets the ephemeris.
     * \return Pointer to ephemeris.
     */
    Ephemeris* getEphemeris( ) { return pointerToEphemeris_; }

    //! Get gravity field model.
    /*!
     * Returns the gravity field model.
     * \return Gravity field model.
     */
    GravityFieldModel* getGravityFieldModel( ) { return pointerToGravityFieldModel_; }

    //! Get atmosphere model.
    /*!
     * Returns the atmosphere model.
     * \return Pointer to the atmosphere model.
     */
    AtmosphereModel* getAtmospheremodel( ) { return pointerToAtmosphereModel_; }

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
               << celestialBody.getGravitationalParameter( ) << std::endl;
        return stream;
    }

protected:

    //! Pointer to ephemeris.
    /*!
     * Pointer to ephemeris.
     */
    Ephemeris* pointerToEphemeris_;

    //! Pointer to gravity field model.
    /*!
     * Pointer to gravity field model.
     */
    GravityFieldModel* pointerToGravityFieldModel_;

    //! Pointer to atmosphere model.
    /*!
     * Pointer to atmosphere model.
     */
    AtmosphereModel* pointerToAtmosphereModel_;
};

#endif // CELESTIAL_BODY_H

// End of file.
