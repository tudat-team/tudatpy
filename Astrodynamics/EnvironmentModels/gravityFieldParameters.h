/*! \file gravityFieldParameters.h
 *    This file contains a class that define a gravity field
 *    with all their characteristics.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/
 *    Version           : 2
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
 *    Date created      : 20 September, 2010
 *    Last modified     : 29 September, 2010
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modifiation is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    author        comment
 *      100920    J. Melman     First creation of code.
 *      100929    K. Kumar      Minor comment changes.
 */

#ifndef GRAVITY_FIELD_PARAMETERS_H
#define GRAVITY_FIELD_PARAMETERS_H

//! Gravity field parameters class.
/*!
 * Gravity field parameters class.
 */
class GravityFieldParameters
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    GravityFieldParameters( );

    //! Customized constructor.
    /*!
     * Customized constructor.
     */
    GravityFieldParameters( const double& gravitationalParameter,
                            const double& referenceRadius );

    //! Sets the gravitational parameter.
    /*!
     * Defines the gravitational parameter in meter^3 per second^2.
     * \param gravitationalParameter.
     */
    void setGravitationalParameter( const double& gravitationalParameter );

    //! Gets the gravitational parameter.
    /*!
     * Returns the gravitational parameter in meter^3 per second^2.
     * \return gravitationalParameter
     */
    const double getGravitationalParameter( ) const;

    //! Sets the reference radius.
    /*!
     * Defines the reference radius used for the spherical harmonics expansion
     * in meters.
     *  \param referenceRadius
     */
    void setReferenceRadius( const double& referenceRadius );

    //! Gets the reference radius.
    /*!
     * Returns the reference radius used for the spherical harmonics expansion
     * in meters.
     * \return referenceRadius
     */
    const double getReferenceRadius( ) const;

protected:

    //! Gravitational parameter.
    /*!
     * The gravitational parameter in meter^3 per second^2.
     */
    double gravitationalParameter_;

    //! Reference radius.
    /*!
     * The reference radius used for the spherical harmonics expansion in
     * meters.
     */
    double referenceRadius_;

};

#endif // GRAVITY_FIELD_PARAMETERS_H

// End of file.
