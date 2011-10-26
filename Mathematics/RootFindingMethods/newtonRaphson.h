/*! \file newtonRaphson.h
 *    Header file that defines the Newton-Raphson method implemented in Tudat.
 *
 *    Path              : /Mathematics/RootFindingMethods/
 *    Version           : 8
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 11 November, 2010
 *    Last modified     : 5 September, 2011
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
 *      101111    E. Iorfida        First creation of code.
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
 */

#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

// Include statements.
#include <iostream>
#include "Mathematics/RootFindingMethods/newtonRaphsonBase.h"
#include "Mathematics/RootFindingMethods/rootFinder.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
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

}

#endif // NEWTONRAPHSON_H

// End of file.
