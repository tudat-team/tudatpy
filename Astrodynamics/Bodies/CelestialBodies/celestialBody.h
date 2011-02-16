/*! \file celestialBody.h
 *    This file contains a class that describes celestial bodies
 *    with all their characteristics.
 *
 *    Path              : /Astrodynamics/Bodies/CelestialBodies/
 *    Version           : 5
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
 *    Last modified     : 13 January, 2011
 *
 *    References
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
 *      YYMMDD    author        comment
 *      100906    J. Melman     First creation of code.
 *      100910    J. Melman     No more switch statement and enum.
 *      100929    K. Kumar      Minor comment changes.
 *      110112    K. Kumar      Modified to use GravityFieldModel; corrected
 *                              path.
 *      110113    K. Kumar      Added getGravityFieldModel() function.
 */

#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H

// Include statements.
#include "body.h"
#include "gravityFieldModel.h"

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
    CelestialBody( );

    //! Set gravity field model.
    /*!
     * Sets the gravity field model.
     * \param pointerToGravityFieldModel Pointer to gravity field model.
     */
    void setGravityFieldModel( GravityFieldModel* pointerToGravityFieldModel);

    //! Get gravitational parameter.
    /*!
     * Returns the gravitational parameter in meter^3 per second^2.
     * \return Gravitational parameter.
     */
    const double getGravitationalParameter( ) const;

    //! Get gravity field model.
    /*!
     * Returns the gravity field model.
     * \return Gravity field model.
     */
     GravityFieldModel* getGravityFieldModel( );

protected:

    //! Gravity field parameters.
    GravityFieldModel* pointerToGravityFieldModel_;
};

#endif // CELESTIAL_BODY_H

// End of file.
