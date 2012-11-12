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
 *      110527    L. van der Ham    Creation of code.
 *      110602    L. van der Ham    Made LibrationPoints class and added comments.
 *      110707    L. van der Ham    Make code compatible with Tudat revision 114.
 *      110629    L. van der Ham    Modifications according to comments first code check.
 *      110710    K. Kumar          Removed duplicated code; modified libration point
 *                                  compute-functions; changed filename and class; L1 and L2
 *                                  function coefficient discrepancies spotted.
 *      110712    L. van der Ham    Changed L1, L2 and L3 function coefficients.
 *      110927    L. van der Ham    Reverted to full equations of motion for determination location
 *                                  of colinear libration points.
 *      111027    K. Kumar          Moved 1-line functions to header file.
 *      120307    K. Kumar          Moved file.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *      van der Ham, L. Interplanetary trajectory design using dynamical systems theory,
 *          MSc thesis, Delft University of Technology, Delft, The Netherlands, 2012.
 *      Mireles James, J.D. Celestial Mechanics Notes Set 4: The Circular Restricted Three Body
 *          Problem, 2006, http://www.math.utexas.edu/users/jjames/hw4Notes.pdf,
 *          last accessed: 18th May, 2012.
 *
 *    WARNING: There seems to be a bug in the computation of the L3 location!
 *
 */

#include <cmath>
#include <iostream>

#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Compute location of Lagrange libration point.
void LibrationPoint::computeLocationOfLibrationPoint(
    LagrangeLibrationPoints lagrangeLibrationPoint )
{
    using std::pow;
    using std::sqrt;

    // Newton-Raphson method implementation.
    // Set the class that contains the functions needed for Newton-Raphson.
    newtonRaphsonAdaptorForLibrationPoint_.setClass( this );

    // Set functions for Newton-Raphson based on collinear libration point passed as input
    // parameter, or computed locations directly of equilateral libration points.
    switch( lagrangeLibrationPoint )
    {
    case L1:

        // Set mathematical functions for determination of libration point L1.
        newtonRaphsonAdaptorForLibrationPoint_.setPointerToFunction(
                    &LibrationPoint::computeL1LocationFunction_ );
        newtonRaphsonAdaptorForLibrationPoint_.setPointerToFirstDerivativeFunction(
                    &LibrationPoint::computeL1FirstDerivativeLocationFunction_ );

        // Set the adaptor for Newton-Raphson method.
        newtonRaphson_->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptorForLibrationPoint_ );

        // Set initial guess for root of equation of motion for L1.
        newtonRaphson_->setInitialGuessOfRoot( 1.0 );

        // Compute root of equation for L1 via Newton-Raphson method.
        newtonRaphson_->execute( );

        // Set position vector of L1 in Cartesian elements.
        positionOfLibrationPoint_.x( ) = newtonRaphson_->getComputedRootOfFunction( );
        positionOfLibrationPoint_.y( ) = 0.0;
        positionOfLibrationPoint_.z( ) = 0.0;

        break;

    case L2:

        // Set mathematical functions for determination of libration point L2.
        newtonRaphsonAdaptorForLibrationPoint_.setPointerToFunction(
                    &LibrationPoint::computeL2LocationFunction_ );
        newtonRaphsonAdaptorForLibrationPoint_.setPointerToFirstDerivativeFunction(
                    &LibrationPoint::computeL2FirstDerivativeLocationFunction_ );

        // Set the adaptor for Newton-Raphson method.
        newtonRaphson_->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptorForLibrationPoint_ );

        // Set initial guess for root of equation of motion in L2.
        newtonRaphson_->setInitialGuessOfRoot( 1.0 );

        // Compute root of equation for L2 via Newton-Raphson method.
        newtonRaphson_->execute( );

        // Set position vector of L2 in Cartesian elements.
        positionOfLibrationPoint_.x( ) = newtonRaphson_->getComputedRootOfFunction( );
        positionOfLibrationPoint_.y( ) = 0.0;
        positionOfLibrationPoint_.z( ) = 0.0;

        break;

    case L3:

        // Set mathematical functions for determination of libration point L3.
        newtonRaphsonAdaptorForLibrationPoint_.setPointerToFunction(
                    &LibrationPoint::computeL3LocationFunction_ );
        newtonRaphsonAdaptorForLibrationPoint_.setPointerToFirstDerivativeFunction(
                    &LibrationPoint::computeL3FirstDerivativeLocationFunction_ );

        // Set the adaptor for Newton-Raphson method.
        newtonRaphson_->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptorForLibrationPoint_ );

        // Set initial guess for root of equation of motion in L3.
        newtonRaphson_->setInitialGuessOfRoot( -1.0 );

        // Compute root of equation for L3 via Newton-Raphson method.
        newtonRaphson_->execute( );

        // Set position vector of L3 in Cartesian elements.
        positionOfLibrationPoint_.x( ) = newtonRaphson_->getComputedRootOfFunction( );
        positionOfLibrationPoint_.y( ) = 0.0;
        positionOfLibrationPoint_.z( ) = 0.0;

        break;

    case L4:

        // Set position vector of L4 in Cartesian elements.
        positionOfLibrationPoint_.x( ) = 0.5 - massParameter_;
        positionOfLibrationPoint_.y( ) = 0.5 * sqrt( 3.0 );
        positionOfLibrationPoint_.z( ) = 0.0;

        break;

    case L5:

        // Set position vector of L5 in Cartesian elements.
        positionOfLibrationPoint_.x( ) = 0.5 - massParameter_;
        positionOfLibrationPoint_.y( ) = -0.5 * sqrt( 3.0 );
        positionOfLibrationPoint_.z( ) = 0.0;

        break;

    default:

        std::cerr << "The Lagrange libration point requested does not exist." << std::endl;
    };
}

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat
