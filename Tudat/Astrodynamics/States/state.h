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
 *      101026    K. Kumar          Creation of code.
 *      101110    K. Kumar          Added setState( ) function for matrices.
 *      101202    J. Melman         Changed path and blank line removed.
 *      110110    K. Kumar          Changed matrices to vectors.
 *      110202    K. Kumar          Removed set/get functions for state and made state a public
 *                                  variable for functionality.
 *      110204    K. Kumar          Note added about public state variable.
 *      110207    K. Kumar          Added ostream overload.
 *      111123    B. Tong Minh      Added custom constructors.
 *      120509    K. Kumar          Modified length-initializer to initialize element values to
 *                                  zero.
 *
 *    References
 *
 *    The variable state is public as opposed to being protected and accessed by set and get
 *    functions to simplify the use of Eigen functions for vectors in the rest of the code,
 *    e.g. setZero( ), segment( ).
 *
 */

#ifndef TUDAT_STATE_H
#define TUDAT_STATE_H

#include <Eigen/Core>

#include <iostream>

namespace tudat
{
namespace astrodynamics
{
namespace states
{

//! State class.
/*!
 * Definition of State class. This is a base class for all state classes
 * defined within Tudat.
 */
class State
{
public:

    //! Default constructor.
    /*!
     * Default constructor, leaves the internal state vector uninitialized.
     */
    State( ) { }

    //! Constructor with state vector.
    /*!
     * Custom constructor ,which sets the internal state vector to the specified vector.
     * \param stateToSet State vector to set.
     */
    State( const Eigen::VectorXd& stateToSet ) : state( stateToSet ) { }

    //! Constructor with initial length.
    /*!
     * Custom constructor, which sets the internal state vector to the specified length, with all
     * elements set to zero.
     * \param stateLength Length of the state vector.
     */
    State( int stateLength ) : state( Eigen::VectorXd::Zero( stateLength ) ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~State( ) { }

    //! State.
    /*!
     * State.
     */
    Eigen::VectorXd state;

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param stateObject State object.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, State& stateObject )
    {
        stream << "The state is set to: " << stateObject.state << std::endl;
        return stream;
    }

protected:

private:
};

} // namespace states
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_STATE_H
