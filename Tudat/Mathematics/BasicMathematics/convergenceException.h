/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      120317    S. Billemont      File created.
 *      120317    S. Billemont      Added optional error message to the constructor.
 *      130305    S. Billemont      Renaming of message variables for GCC.
 *      130307    K. Kumar          Added missing Doxygen comments to exception struct.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_CONVERGENCE_EXEPTION_H
#define TUDAT_CONVERGENCE_EXEPTION_H

#include <iostream>
#include <stdexcept>

#include <boost/exception/all.hpp>

namespace tudat
{
namespace basic_mathematics
{

//! Exception indicating that an algorithm could not complete convergence to a solution.
/*!
 * You can add extra information using;
 * void f()
 * {
 *   throw ConvergenceExeption() << my_info(42);
 * }
 *
 * \sa Boost.Exception.
 */
struct ConvergenceException : public virtual boost::exception, public virtual std::exception
{
public:

    //! Constructor that sets exception message to default string.
    ConvergenceException( )
        : message( "Failed to converge to a solution.")
    { }

    //! Constructor that sets exception message based on input character array.
    ConvergenceException( const char* errMessage )
        : message( errMessage )
    { }

    //! Constructor that sets exception message based on input string.
    ConvergenceException( const std::string& errMessage )
        : message( errMessage.c_str( ) )
    { }

    //! Return what the exception message stored is.
    const char* what( ) const throw( ) { return message; }

protected:

private:

    //! Stored exception message.
    const char* message;
};

} // namespace basic_mathematics
} // namespace tudat

#endif  // TUDAT_CONVERGENCE_EXEPTION_H
