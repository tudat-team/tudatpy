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
 *      110207    K. Kumar          Creation of code.
 *
 *    References
 *
 */

#ifndef TUDAT_CARTESIAN_VELOCITY_ELEMENTS_H
#define TUDAT_CARTESIAN_VELOCITY_ELEMENTS_H

#include <iostream>
#include "Tudat/Astrodynamics/States/state.h"

namespace tudat
{
namespace astrodynamics
{
namespace states
{

//! Cartesian position elements class.
/*!
 * This class contains the Cartesian position elements.
 */
class CartesianVelocityElements : public State
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    CartesianVelocityElements( ) { state.setZero( 3 ); }

    //! Set Cartesian element: xDot.
    /*!
     * Sets the Cartesian element: xDot.
     * \param cartesianElementXDot Cartesian element: xDot.
     */
    void setCartesianElementXDot( const double cartesianElementXDot )
    {
        state( 0 ) = cartesianElementXDot;
    }

    //! Set Cartesian element: yDot.
    /*!
     * Sets the Cartesian element: yDot.
     * \param cartesianElementYDot Cartesian element: yDot.
     */
    void setCartesianElementYDot( const double cartesianElementYDot )
    {
        state( 1 ) = cartesianElementYDot;
    }

    //! Set Cartesian element: zDot.
    /*!
     * Sets the Cartesian element: zDot.
     * \param cartesianElementZDot Cartesian element: zDot.
     */
    void setCartesianElementZDot( const double cartesianElementZDot )
    {
        state( 2 ) = cartesianElementZDot;
    }

    //! Get Cartesian element: xDot.
    /*!
     * Returns the Cartesian element: xDot.
     * \return Cartesian element: xDot.
     */
    double getCartesianElementXDot( ) { return state( 0 ); }

    //! Get Cartesian element: yDot.
    /*!
     * Returns the Cartesian element: yDot.
     * \return Cartesian element: yDot.
     */
    double getCartesianElementYDot( ) { return state( 1 ); }

    //! Get Cartesian element: zDot.
    /*!
     * Returns the Cartesian element: zDot.
     * \return Cartesian element: zDot.
     */
    double getCartesianElementZDot( ) { return state( 2 ); }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param cartesianVelocityElements CartesianVelocityElements object.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     CartesianVelocityElements& cartesianVelocityElements )
    {
        stream << "The state is set to: " << cartesianVelocityElements.state << std::endl;
        return stream;
    }

protected:

private:
};

} // namespace states
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_CARTESIAN_VELOCITY_ELEMENTS_H
