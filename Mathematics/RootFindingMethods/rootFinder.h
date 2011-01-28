/*! \file rootFinder.h
 *    This header file contains a base class for all root-finder algorithms
 *    classes in Tudat.
 *
 *    Path              : /Mathematics/RootFindingMethods/
 *    Version           : 7
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
 *    Last modified     : 24 January, 2011
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
 *      101111    E. Iorfida    First creation of code.
 *      101116    E. Iorfida    Added set-get functions.
 *      101121    E. Iorfida    Added Doxygen comments.
 *      110111    E. Iorfida    Deleted useless lines, and modified
 *                              punctuation.
 *      110111    K. Kumar      Renamed functions and variables to be
 *                              more accurate and descriptive; added
 *                              "End of file." comment; added example of
 *                              pointer-to-member-function for LambertTargeter.
 *      110119    K. Kumar      Updated code to work with adaptor and abstract
 *                              base implementation so that pointer-to-member
 *                              functions are not required.
 *      110124    E. Iorfida    Added set/get functions for maximum number of
 *                              iterations.
 */

#ifndef ROOTFINDER_H
#define ROOTFINDER_H

// Include statements.
#include <ctime>
#include <iostream>
#include "basicMathematicsFunctions.h"
#include "linearAlgebra.h"
#include "rootFinderBase.h"

// Using declarations.
using std::cerr;
using std::endl;

//! Root-finder class.
/*!
 * This class serves as a base class for all root-finder algorithms included in
 * Tudat.
 */
class RootFinder
{
public:

    // Definition of typedef.
    // Functions to which root-finding methods are applied can be passed as
    // pointers to global functions. ( Polymorphic ) Pointers to RootFinder
    // objects must be used for global functions to be usable.
    typedef double ( *pointerToDoubleTakingFunction )( double& );

    //! Default constructor.
    /*!
     * Default constructor.
     */
    RootFinder( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~RootFinder( );

    //! Set initial guess of the root of mathematical function.
    /*!
     * Sets initial guess of the root of the mathematical function.
     * \param initialGuessOfRoot Initial guess of root of mathematical
     *          function.
     */
    void setInitialGuessOfRoot( const double& initialGuessOfRoot );

    //! Set maximum number of iterations.
    /*!
     * Sets maximum number of iterations for root-finding method.
     * \param maximumNumberOfIterations Maximum number of iterations.
     */
    void setMaximumNumberOfIterations( const unsigned int&
                                       maximumNumberOfIterations );

    //! Set tolerance.
    /*!
     * Sets tolerance for root-finding method.
     * \param tolerance Tolerance.
     */
    void setTolerance( const double& tolerance );

    //! Set pointer to mathematical function.
    /*!
     * Sets a pointer to the mathematical function ( global ) to which the
     * root-finding method is applied.
     * \param globalFunction Pointer to global mathematical function.
     */
    void setMathematicalFunction( pointerToDoubleTakingFunction
                                  globalFunction );

    //! Set pointer to first-derivative mathematical function.
    /*!
     * Sets a pointer to the first-derivative mathematical function ( global )
     * to which the root-finding method is applied.
     * \param globalFirstDerivativeFunction Pointer to global first-derivative
     *           mathematical function.
     */
    void setFirstDerivativeMathematicalFunction(
            pointerToDoubleTakingFunction globalFirstDerivativeFunction );

    //! Get maximum number of iterations.
    /*!
     * Returns the maximum number of iterations.
     * \return Number of iterations.
     */
    unsigned int& getMaximumNumberOfIterations( );

    //! Get tolerance.
    /*!
     * Returns the tolerance.
     * \return Tolerance.
     */
    double& getTolerance( );

    //! Get root of mathematical function.
    /*!
     * Returns the computed root of the mathmatical function.
     * \return Computed root of the mathematical function.
     */
    double& getComputedRootOfFunction( );

    //! Execute.
    /*!
     * This is called to execute a root-finder method.
     */
    virtual void execute( ) = 0;

protected:

    //! Maximum number of iterations.
    /*!
     * Maximum number of iterations.
     */
    unsigned int maximumNumberOfIterations_;

    //! Initial guess of the root of mathematical function.
    /*!
     * Initial guess of the root of mathematical function.
     */
    double initialGuessOfRoot_;

    //! Current value of the root of mathematical function.
    /*!
     * Current value of the root of mathematical function.
     */
    double currentValueOfRoot_;

    //! Next value of the root of mathematical function.
    /*!
     * Next value of the root of mathematical function.
     */
    double nextValueOfRoot_;

    //! Tolerance.
    /*!
     * Maximum allowed difference between the next value and the current value
     * of the root of the mathematical function.
     */
    double tolerance_;

    //! Pointer to global function.
    /*!
     * Pointer to global mathematical function to which the root-finding method
     * is applied.
     */
    pointerToDoubleTakingFunction pointerToGlobalFunction_;

    //! Pointer to global first-derivative function.
    /*!
     * Pointer to global first-derivative mathematical function to which the
     * root-finding method is applied.
     */
    pointerToDoubleTakingFunction pointerToGlobalFirstDerivativeFunction_;

    //! Pointer to root-finder abstract base class.
    /*!
     * Pointer to root-finder abstract base class.
     */
    RootFinderBase* pointerToRootFinderBase_;

private:
};

#endif // ROOTFINDER_H

// End of file.
