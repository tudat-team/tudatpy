/*! \file body.h
 *    This file contains a class that describes bodies, both
 *    celestial bodies and vehicles, with all their characteristics.
 *
 *    Path              : /Astrodynamics/Bodies/
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
 *      100920    J. Melman     First creation of code.
 *      100929    K. Kumar      Minor comment changes.
 */

#ifndef BODY_H
#define BODY_H

//! Body base class.
/*!
 * Body base class.
 */
class Body
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     */
    Body( );

    //! Sets the shape model.
    /*!
     *  Defines the shape model.
     */
    void setShapeModel( );
};

#endif // BODY_H

// End of file.
