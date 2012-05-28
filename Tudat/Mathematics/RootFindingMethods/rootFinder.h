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
 *      101116    E. Iorfida        Added set-get functions.
 *      101121    E. Iorfida        Added Doxygen comments.
 *      110111    E. Iorfida        Deleted useless lines, and modified punctuation.
 *      110111    K. Kumar          Renamed functions and variables to be more accurate and
 *                                  descriptive; added "End of file." comment; added example of
 *                                  pointer-to-member-function for LambertTargeter.
 *      110119    K. Kumar          Updated code to work with adaptor and  abstract base
 *                                  implementation so that pointer-to-member functions are not
 *                                  required.
 *      110124    E. Iorfida        Added set/get functions for maximum number of iterations.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *
 *    References
 *
 */

#ifndef TUDAT_ROOT_FINDER_H
#define TUDAT_ROOT_FINDER_H

#include <iostream>

#include "Tudat/Mathematics/RootFindingMethods/rootFinderBase.h"

namespace tudat
{

//! Root-finder class.
/*!
 * This class serves as a base class for all root-finder algorithms included in
 * Tudat.
 */
class RootFinder
{
public:

    //! Definition of typedef.
    /*!
     * Functions to which root-finding methods are applied can be passed as
     * pointers to global functions. ( Polymorphic ) Pointers to RootFinder
     * objects must be used for global functions to be usable.
     */
    typedef double ( *pointerToDoubleTakingFunction )( double& );

    //! Default constructor.
    /*!
     * Default constructor.
     */
    RootFinder( ) : maximumNumberOfIterations_( 100 ), initialGuessOfRoot_( -0.0 ),
        currentValueOfRoot_( -0.0 ), nextValueOfRoot_( -0.0 ), tolerance_( 1.0e-12 ),
        pointerToGlobalFunction_( NULL ), pointerToGlobalFirstDerivativeFunction_( NULL ),
        pointerToRootFinderBase_( NULL ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~RootFinder( ) { }

    //! Set initial guess of the root of mathematical function.
    /*!
     * Sets initial guess of the root of the mathematical function.
     * \param initialGuessOfRoot Initial guess of root of mathematical
     *          function.
     */
    void setInitialGuessOfRoot( double initialGuessOfRoot )
    { initialGuessOfRoot_ = initialGuessOfRoot; currentValueOfRoot_ = initialGuessOfRoot; }

    //! Set maximum number of iterations.
    /*!
     * Sets maximum number of iterations for root-finding method.
     * \param maximumNumberOfIterations Maximum number of iterations.
     */
    void setMaximumNumberOfIterations( const unsigned int& maximumNumberOfIterations ) 
    { maximumNumberOfIterations_ = maximumNumberOfIterations; }

    //! Set tolerance.
    /*!
     * Sets tolerance for root-finding method.
     * \param tolerance Tolerance.
     */
    void setTolerance( double tolerance ) { tolerance_ = tolerance; }

    //! Set pointer to mathematical function.
    /*!
     * Sets a pointer to the mathematical function ( global ) to which the
     * root-finding method is applied.
     * \param globalFunction Pointer to global mathematical function.
     */
    void setMathematicalFunction( pointerToDoubleTakingFunction globalFunction )
    { pointerToGlobalFunction_ = globalFunction; }

    //! Set pointer to first-derivative mathematical function.
    /*!
     * Sets a pointer to the first-derivative mathematical function ( global )
     * to which the root-finding method is applied.
     * \param globalFirstDerivativeFunction Pointer to global first-derivative
     *           mathematical function.
     */
    void setFirstDerivativeMathematicalFunction( pointerToDoubleTakingFunction
                                                 globalFirstDerivativeFunction )
    { pointerToGlobalFirstDerivativeFunction_ = globalFirstDerivativeFunction; }

    //! Get maximum number of iterations.
    /*!
     * Returns the maximum number of iterations.
     * \return Number of iterations.
     */
    unsigned int& getMaximumNumberOfIterations( ) { return maximumNumberOfIterations_; }

    //! Get tolerance.
    /*!
     * Returns the tolerance.
     * \return Tolerance.
     */
    double& getTolerance( ) { return tolerance_; }

    //! Get root of mathematical function.
    /*!
     * Returns the computed root of the mathmatical function.
     * \return Computed root of the mathematical function.
     */
    double& getComputedRootOfFunction( ) { return nextValueOfRoot_; }

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

} // namespace tudat

#endif // TUDAT_ROOT_FINDER_H
