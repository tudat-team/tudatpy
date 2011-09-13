/*! \file rootFinder.cpp
 *    This source file contains a base class for all root-finder algorithms
 *    classes in Tudat.
 *
 *    Path              : /Mathematics/RootFindingMethods/
 *    Version           : 5
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
 *    Date created      : 21 November, 2010
 *    Last modified     : 24 January, 2011
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
 *      YYMMDD    Author            comment
 *      101121    E. Iorfida        First creation of code.
 *      110111    E. Iorfida        Deleted useless lines, and modified
 *                                  punctuation.
 *      110111    K. Kumar          Minor comment changes; "End of file."
 *                                  comment added; changes in .h file updated.
 *      110119    K. Kumar          Updated code to work with adaptor and
 *                                  abstract base implementation so that
 *                                  pointer-to-member functions are not
 *                                  required.
 *      110124    E. Iorfida        Added set/get functions for maximum number
 *                                  of iterations.
 */

// Include statements.
#include "rootFinder.h"

//! Default constructor.
RootFinder::RootFinder( ) : maximumNumberOfIterations_( 100 ),
                            initialGuessOfRoot_( -0.0 ),
                            currentValueOfRoot_( -0.0 ),
                            nextValueOfRoot_( -0.0 ),
                            tolerance_( 1.0e-12 ),
                            pointerToGlobalFunction_( NULL ),
                            pointerToGlobalFirstDerivativeFunction_( NULL )
{
}

//! Default destructor.
RootFinder::~RootFinder( )
{
}

//! Set initial guess of the root of mathematical function.
void RootFinder::setInitialGuessOfRoot( const double& initialGuessOfRoot )
{
    initialGuessOfRoot_ = initialGuessOfRoot;
    currentValueOfRoot_ = initialGuessOfRoot;
}

//! Set maximum number of iterations.
void RootFinder::setMaximumNumberOfIterations( const unsigned int&
                                               maximumNumberOfIterations )
{
    maximumNumberOfIterations_ = maximumNumberOfIterations;
}

//! Set tolerance.
void RootFinder::setTolerance( const double &tolerance )
{
    tolerance_ = tolerance;
}

//! Set pointer to mathematical function.
void RootFinder::setMathematicalFunction( pointerToDoubleTakingFunction
                                          globalFunction )
{
    pointerToGlobalFunction_ = globalFunction;
}

//! Set pointer to first-derivative mathematical function.
void RootFinder::setFirstDerivativeMathematicalFunction(
        pointerToDoubleTakingFunction globalFirstDerivativeFunction )
{
    pointerToGlobalFirstDerivativeFunction_ = globalFirstDerivativeFunction;
}

//! Get maximum number of iterations.
unsigned int& RootFinder::getMaximumNumberOfIterations( )
{
    return maximumNumberOfIterations_;
}

//! Get tolerance.
double& RootFinder::getTolerance( )
{
    return tolerance_;
}

//! Get root of mathematical function.
double& RootFinder::getComputedRootOfFunction( )
{
    return nextValueOfRoot_;
}

// End of file.
