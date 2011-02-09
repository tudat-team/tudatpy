/*! \file newtonRaphson.cpp
 *    Source file of the Newton-Raphson method implemented in Tudat.
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
 *      101118    E. Iorfida    Added algorithm implementation.
 *      101121    E. Iorfida    Modified algorithm and added Doxygen comments.
 *      101216    E. Iorfida    Modified punctuation.
 *      110111    E. Iorfida    Deleted useless lines, and modified
 *                              punctuation.
 *      110111    K. Kumar      Updated based on changes to .h file; changed
 *                              do-while loop to for loop to prevent infinite
 *                              loop; added "End of file." comment.
 *      110119    K. Kumar      Updated code to work with adaptor and abstract
 *                              base implementation so that pointer-to-member
 *                              functions are not required; changed filename.
 *      110124    K. Kumar      Added cerr statement for when loop does not
 *                              converge.
 */

// Include statements.
#include "newtonRaphson.h"

//! Default constructor.
NewtonRaphson::NewtonRaphson( ) : pointerToNewtonRaphsonBase_( NULL )
{
}

//! Default destructor.
NewtonRaphson::~NewtonRaphson( )
{
}

//! Set adaptor class for Newton-Raphson.
void NewtonRaphson::setNewtonRaphsonAdaptor( NewtonRaphsonBase*
                                             pointerToNewtonRaphsonBase )
{
    pointerToNewtonRaphsonBase_ = pointerToNewtonRaphsonBase;
    pointerToGlobalFunction_ = NULL;
    pointerToGlobalFirstDerivativeFunction_ = NULL;
}

//! Execute Newton-Raphson method.
void NewtonRaphson::execute( )
{
    // Algorithm implementation.
    // First step.
    nextValueOfRoot_ = initialGuessOfRoot_;

    for ( unsigned int i = 0; i < maximumNumberOfIterations_; i++ )
    {
        // Update current value of root.
        currentValueOfRoot_ = nextValueOfRoot_;

        // Compute next value of root using the following algorithm:
        // x_n+1 = x_n - F(x_n)/F'(x_n).

        // Check if necessary global mathematical functions are set.
        if ( pointerToGlobalFunction_
             && pointerToGlobalFirstDerivativeFunction_ )
        {
            // Use pointers to global functions.
            nextValueOfRoot_ = currentValueOfRoot_
                               - ( *pointerToGlobalFunction_ )(
                                       currentValueOfRoot_ )
                               / ( *pointerToGlobalFirstDerivativeFunction_ )(
                                       currentValueOfRoot_ );
        }

        // Else check if class with member mathematical functions is set.
        else if ( pointerToNewtonRaphsonBase_ )
        {
            // Use pointers to member functions.
            nextValueOfRoot_ = currentValueOfRoot_
                               - pointerToNewtonRaphsonBase_
                               ->computeFunction( currentValueOfRoot_ )
                               /  pointerToNewtonRaphsonBase_
                               ->computeFirstDerivativeFunction(
                                       currentValueOfRoot_ );
        }

        // Check if difference between successive iterations of the
        // root-finding method satisfies the set tolerance and break from loop
        // if it is satisfied.
        if ( mathematics::computeAbsoluteValue( nextValueOfRoot_
                                                - currentValueOfRoot_ )
            <= tolerance_ )
        {
            break;
        }

        // If the end of the loop is reached before the tolerance is satisfied,
        // print cerr statement to indicate to the user that there method did
        // not converge.
        if ( i == maximumNumberOfIterations_ - 1 )
        {
            // Reset values of root to initialized values.
            currentValueOfRoot_ = -0.0;
            nextValueOfRoot_ = -0.0;

            cerr << "The Newton-Raphson algorithm did not converge after "
                 << maximumNumberOfIterations_ << " iterations." << endl;
        }
    }
}

// End of file.
