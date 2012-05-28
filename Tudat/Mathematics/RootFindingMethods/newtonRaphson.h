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
 *      101111    E. Iorfida        Creation of code.
 *      101116    E. Iorfida        Added setFunction and execute.
 *      101121    E. Iorfida        Added Doxygen comments.
 *      110111    E. Iorfida        Deleted useless lines, and modified punctuation.
 *      110111    K. Kumar          Changed variable and function names to be more descriptive;
 *                                  added "End of file."
 *      110114    K. Kumar          Removed circular code dependency.
 *      110119    K. Kumar          Updated code to work with adaptor and abstract base
 *                                  implementation so that pointer-to-member functions are not
 *                                  required; changed filename.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *
 *    References
 *
 */

#ifndef TUDAT_NEWTON_RAPHSON_H
#define TUDAT_NEWTON_RAPHSON_H

#include <iostream>

#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonBase.h"
#include "Tudat/Mathematics/RootFindingMethods/rootFinder.h"

namespace tudat
{

//! Newton-Raphson class.
/*!
 * Implementation of Newton-Raphson class in Tudat.
 */
class NewtonRaphson : public RootFinder
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    NewtonRaphson( ) : pointerToNewtonRaphsonBase_( NULL ) { }

    //! Set adaptor class for Newton-Raphson.
    /*!
     * Sets the adaptor class for the Newton-Raphson method, which serves as a
     * method to communicate with the class containing the mathematical
     * functions used for the Newton-Rapshon method. Objects of the
     * NewtonRaphsonAdaptor class need to be passed as the argument to this
     * function.
     * \param pointerToNewtonRaphsonBase Polymorphic pointer to adaptor class.
     */
    void setNewtonRaphsonAdaptor( NewtonRaphsonBase* pointerToNewtonRaphsonBase );

    //! Execute Newton-Raphson method.
    /*!
     * This function executes the Newton-Raphson root-finding method.
     */
    void execute( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param newtonRaphson Newton-Raphson.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, NewtonRaphson& newtonRaphson );

protected:

private:

    //! Polymorphic pointer to Newton-Raphson abstract base class.
    /*!
     * Polymorphic pointer to Newton-Raphson abstract base class.
     */
    NewtonRaphsonBase* pointerToNewtonRaphsonBase_;
};

} // namespace tudat

#endif // TUDAT_NEWTON_RAPHSON_H
