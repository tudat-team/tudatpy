/*! \file cartesianElements.h
 *    This header file contains the Cartesian elements class included in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 20 October, 2010
 *    Last modified     : 4 February, 2011
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
 *      101020    K. Kumar          First creation of code.
 *      101026    K. Kumar          Added position and velocity vectors, modified existing
 *                                  comments.
 *      101201    E. Iorfida        Modified punctuation.
 *      101202    J. Melman         Corrected some Doxygen comments, changed Cartesian into
 *                                  Cartesian. Private variables are obsolete now, since everything
 *                                  is done with state_ from base class.
 *      110110    K. Kumar          Minor modifications.
 *      110131    K. Kumar          Moved code to CartesianPositionElements  class; added
 *                                  inheritance.
 *      110204    K. Kumar          Removed "vector" from naming.
 */

#ifndef CARTESIANELEMENTS_H
#define CARTESIANELEMENTS_H

// Include statements.
#include <iostream>
#include "Astrodynamics/States/cartesianPositionElements.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Cartesian elements class.
/*!
 * Cartesian elements class.
 */
class CartesianElements : public CartesianPositionElements
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    CartesianElements( ) : position_( Vector3d::Zero( ) ) { state.setZero( 6 ); }

    //! Set Cartesian element: xDot.
    /*!
     * Sets the Cartesian element: xDot.
     * \param cartesianElementXDot Cartesian element: xDot.
     */
    void setCartesianElementXDot( const double& cartesianElementXDot )
    { state( 3 ) = cartesianElementXDot; }

    //! Set Cartesian element: yDot.
    /*!
     * Sets the Cartesian element: yDot.
     * \param cartesianElementYDot Cartesian element: yDot.
     */
    void setCartesianElementYDot( const double& cartesianElementYDot )
    { state( 4 ) = cartesianElementYDot; }

    //! Set Cartesian element: zDot.
    /*!
     * Sets the Cartesian element: zDot.
     * \param cartesianElementZDot Cartesian element: zDot.
     */
    void setCartesianElementZDot( const double& cartesianElementZDot )
    { state( 5 ) = cartesianElementZDot; }

    //! Set position.
    /*!
     * Sets the position.
     * \return position Position.
     */
    void setPosition( Vector3d& position ) { state.segment( 0, 3 ) = position; }

    //! Set velocity.
    /*!
     * Sets the velocity.
     * \return velocity Velocity.
     */
    void setVelocity( Vector3d& velocity ) { state.segment( 3, 3 ) = velocity; }

    //! Get Cartesian element: xDot.
    /*!
     * Returns the Cartesian element: xDot.
     * \return Cartesian element: xDot.
     */
    double& getCartesianElementXDot( ) { return state( 3 ); }

    //! Get Cartesian element: yDot.
    /*!
     * Returns the Cartesian element: yDot.
     * \return Cartesian element: yDot.
     */
    double& getCartesianElementYDot( ) { return state( 4 ); }

    //! Get Cartesian element: zDot.
    /*!
     * Returns the Cartesian element: zDot.
     * \return Cartesian element: zDot.
     */
    double& getCartesianElementZDot( ) { return state( 5 ); }

    //! Get position.
    /*!
     * Returns the position.
     * \return Position.
     */
    Vector3d& getPosition( ) { position_ = state.segment( 0, 3 ); return position_; }

    //! Get velocity.
    /*!
     * Returns the velocity.
     * \return Velocity.
     */
    Vector3d& getVelocity( ) { velocity_ = state.segment( 3, 3 ); return velocity_; }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param cartesianElements CartesianElements object.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, CartesianElements& cartesianElements )
    { stream << "The state is set to: " << cartesianElements.state << std::endl; return stream; }

protected:

private:

    //! Position.
    /*!
     * Position.
     */
    Vector3d position_;

    //! Velocity.
    /*!
     * Velocity.
     */
    Vector3d velocity_;
};

}

#endif // CARTESIANELEMENTS_H

// End of file.
