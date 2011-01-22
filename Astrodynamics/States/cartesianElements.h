/*! \file cartesianElements.h
 *    This header file contains the Cartesian elements class included in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 20 Checked, 2010
 *    Last modified     : 02 December, 2010
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
 *      101020    K. Kumar      First creation of code.
 *      101026    K. Kumar      Added position and velocity vectors, modified
 *                              existing comments.
 *      101201    E. Iorfida    Modified punctuation.
 *      101202    J. Melman     Corrected some Doxygen comments, changed
 *                              Cartesian into Cartesian. Private variables
 *                              are obsolete now, since everything is done
 *                              with state_ from base class.
 */

#ifndef CARTESIANELEMENTS_H
#define CARTESIANELEMENTS_H

// Include statements.
#include "orbitalElements.h"

//! Cartesian elements class.
/*!
 * Cartesian elements class.
 */
class CartesianElements : public OrbitalElements
{
public:\

    //! Default constructor.
    /*!
     * Default constructor.
     */
    CartesianElements( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CartesianElements( );

    //! Set Cartesian element: x.
    /*!
     * This function sets the Cartesian element: x.
     * \param cartesianElementX Cartesian element: x.
     */
    void setCartesianElementX( const double& cartesianElementX );

    //! Set Cartesian element: y.
    /*!
     * This function sets the Cartesian element: y.
     * \param cartesianElementY Cartesian element: y.
     */
    void setCartesianElementY( const double& cartesianElementY );

    //! Set Cartesian element: z.
    /*!
     * This function sets the Cartesian element: z.
     * \param cartesianElementZ Cartesian element: z.
     */
    void setCartesianElementZ( const double& cartesianElementZ );

    //! Set Cartesian element: xDot.
    /*!
     * This function sets the Cartesian element: xDot.
     * \param cartesianElementXDot Cartesian element: xDot.
     */
    void setCartesianElementXDot( const double& cartesianElementXDot );

    //! Set Cartesian element: yDot.
    /*!
     * This function sets the Cartesian element: yDot.
     * \param cartesianElementYDot Cartesian element: yDot.
     */
    void setCartesianElementYDot( const double& cartesianElementYDot );

    //! Set Cartesian element: zDot.
    /*!
     * This function sets the Cartesian element: zDot.
     * \param cartesianElementZDot Cartesian element: zDot.
     */
    void setCartesianElementZDot( const double& cartesianElementZDot );

    //! Get Cartesian element: x.
    /*!
     * This function returns the Cartesian element: x.
     * \return Cartesian element: x.
     */
    double& getCartesianElementX( );

    //! Get Cartesian element: y.
    /*!
     * This function returns the Cartesian element: y.
     * \return Cartesian element: y.
     */
    double& getCartesianElementY( );

    //! Get Cartesian element: z.
    /*!
     * This function returns the Cartesian element: z.
     * \return Cartesian element: z.
     */
    double& getCartesianElementZ( );

    //! Get Cartesian element: xDot.
    /*!
     * This function returns the Cartesian element: xDot.
     * \return Cartesian element: xDot.
     */
    double& getCartesianElementXDot( );

    //! Get Cartesian element: yDot.
    /*!
     * This function returns the Cartesian element: yDot.
     * \return Cartesian element: yDot.
     */
    double& getCartesianElementYDot( );

    //! Get Cartesian element: zDot.
    /*!
     * This function returns the Cartesian element: zDot.
     * \return Cartesian element: zDot.
     */
    double& getCartesianElementZDot( );

    //! Set position vector
    /*!
     * This function sets the position vector.
     * \param positionVector Position vector.
     */
    void setPositionVector ( Vector3d& positionVector );

    //! Set velocity vector
    /*!
     * This function sets the velocity vector.
     * \param velocityVector Velocity vector.
     */
    void setVelocityVector ( Vector3d& velocityVector );

    //! Get position vector.
    /*!
     * This function returns the position vector.
     * \return Position vector.
     */
    Vector3d& getPositionVector( );

    //! Get velocity vector.
    /*!
     * This function returns the velocity vector.
     * \return Velocity vector.
     */
    Vector3d& getVelocityVector( );

protected:

private:

    //! Position vector.
    /*!
     * Position vector.
     */
    Vector3d positionVector_;

    //! Velocity vector.
    /*!
     * Velocity vector.
     */
    Vector3d velocityVector_;
};

#endif // CARTESIANELEMENTS_H

// End of file.
