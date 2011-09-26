/*! \file librationPoint.cpp
 *    Source file that defines the computation of the location of a Lagrange
 *    libration point in a Circular Restricted Three-Body Problem (CRTBP).
 *
 *    Path              : /Astrodynamics/MissionSegments/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : L. van der Ham
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.vanderHam@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 7 June, 2011
 *    Last modified     : 12 July, 2011
 *
 *    References
 *      van der Ham, L. TBD.
 *      Mireles James, J.D., "Celestial Mechanics Notes Set 4: The Circular
 *          Restricted Three Body Problem",
 *          http://www.math.utexas.edu/users/jjames/celestMech, 2006.
 *
 *    Notes
 *      Determination of x-location of colinear libration points
 *      via 5th-order polynomials is based on method given in (James, 2006).
 *      However, coefficients of these polynomials are erroneous in
 *      (James, 2006). Derivation of the correct coefficients is given in
 *      (van der Ham, TBD).
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
 *      110527    L. van der Ham    First creation of code.
 *      110602    L. van der Ham    Made LibrationPoints class and added comments.
 *      110707    L. van der Ham    Make code compatible with Tudat revision 114.
 *      110629    L. van der Ham    Modifications according to comments first code check.
 *      110710    K. Kumar          Removed duplicated code; modified libration
 *                                  point compute-functions; changed filename
 *                                  and class; L1 and L2 function coefficient
 *                                  discrepancies spotted.
 *      110712    L. van der Ham    Changed L1, L2 and L3 function coefficients.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Astrodynamics/ForceModels/gravitationalForceModel.h"
#include "Astrodynamics/MissionSegments/LibrationPoints/librationPoint.h"
#include "Mathematics/basicMathematicsFunctions.h"

// Using declarations.
using std::cerr;
using std::endl;
using std::pow;
using std::sqrt;

//! Compute mass parameter.
void LibrationPoint::computeMassParameter( )
{
    // Calculate dimensionless mass parameter for bodies 1 and 2.
    massParameter_ = pointerToSecondaryCelestialBody_->getGravitationalParameter( )
            / ( pointerToPrimaryCelestialBody_->getGravitationalParameter( )
                + pointerToSecondaryCelestialBody_->getGravitationalParameter( ) );
}

//! Compute location of Lagrange libration point.
void LibrationPoint::computeLocationOfLibrationPoint(
    LagrangeLibrationPoints lagrangeLibrationPoint )
{
    // Compute mass parameter squared.
    massParameterSquared_ = pow( massParameter_, 2.0 );

    // Compute ( 1 - mass parameter )^2.
    oneMinusMassParameterSquared_ = pow( 1.0 - massParameter_, 2.0 );

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
        pointerToNewtonRaphson_->setNewtonRaphsonAdaptor(
                    &newtonRaphsonAdaptorForLibrationPoint_ );

        // Set initial guess for root of 5th-order polynomial.
        pointerToNewtonRaphson_->setInitialGuessOfRoot( 1.0 );

        // Compute root of equation for L1 via Newton-Raphson method.
        pointerToNewtonRaphson_->execute( );

        // Set position vector of L1 in Cartesian elements.
        positionOfLibrationPoint_.setCartesianElementX(
                    pointerToNewtonRaphson_->getComputedRootOfFunction( ) );
        positionOfLibrationPoint_.setCartesianElementY( 0.0 );
        positionOfLibrationPoint_.setCartesianElementZ( 0.0 );

        break;

    case L2:

        // Set mathematical functions for determination of libration point L2.
        newtonRaphsonAdaptorForLibrationPoint_.setPointerToFunction(
                    &LibrationPoint::computeL2LocationFunction_ );
        newtonRaphsonAdaptorForLibrationPoint_.setPointerToFirstDerivativeFunction(
                    &LibrationPoint::computeL2FirstDerivativeLocationFunction_ );

        // Set the adaptor for Newton-Raphson method.
        pointerToNewtonRaphson_->setNewtonRaphsonAdaptor(
                    &newtonRaphsonAdaptorForLibrationPoint_ );

        // Set initial guess for root of 5th-order polynomial.
        pointerToNewtonRaphson_->setInitialGuessOfRoot( 1.0 );

        // Compute root of equation for L2 via Newton-Raphson method.
        pointerToNewtonRaphson_->execute( );

        // Set position vector of L1 in Cartesian elements.
        positionOfLibrationPoint_.setCartesianElementX(
                    pointerToNewtonRaphson_->getComputedRootOfFunction( ) );
        positionOfLibrationPoint_.setCartesianElementY( 0.0 );
        positionOfLibrationPoint_.setCartesianElementZ( 0.0 );

        break;

    case L3:

        // Set mathematical functions for determination of libration point L3.
        newtonRaphsonAdaptorForLibrationPoint_.setPointerToFunction(
                    &LibrationPoint::computeL3LocationFunction_ );
        newtonRaphsonAdaptorForLibrationPoint_.setPointerToFirstDerivativeFunction(
                    &LibrationPoint::computeL3FirstDerivativeLocationFunction_ );

        // Set the adaptor for Newton-Raphson method.
        pointerToNewtonRaphson_->setNewtonRaphsonAdaptor(
                    &newtonRaphsonAdaptorForLibrationPoint_ );

        // Set initial guess for root of 5th-order polynomial.
        pointerToNewtonRaphson_->setInitialGuessOfRoot( -1.0 );

        // Compute root of equation for L1 via Newton-Raphson method.
        pointerToNewtonRaphson_->execute( );

        // Set position vector of L1 in Cartesian elements.
        positionOfLibrationPoint_.setCartesianElementX(
                    pointerToNewtonRaphson_->getComputedRootOfFunction( ) );
        positionOfLibrationPoint_.setCartesianElementY( 0.0 );
        positionOfLibrationPoint_.setCartesianElementZ( 0.0 );

        break;

    case L4:

        // Set position vector of L4 in Cartesian elements.
        positionOfLibrationPoint_.setCartesianElementX( 0.5 - massParameter_ );
        positionOfLibrationPoint_.setCartesianElementY( 0.5 * sqrt( 3.0 ) );
        positionOfLibrationPoint_.setCartesianElementZ( 0.0 );

        break;

    case L5:

        // Set position vector of L5 in Cartesian elements.
        positionOfLibrationPoint_.setCartesianElementX( 0.5 - massParameter_ );
        positionOfLibrationPoint_.setCartesianElementY( -0.5 * sqrt( 3.0 ) );
        positionOfLibrationPoint_.setCartesianElementZ( 0.0 );

        break;

    default:

        cerr << "The Lagrange libration point requested does not exist." << endl;
    };
}

//! Compute 5th-order polynomial for location of L1 libration point.
double LibrationPoint::computeL1LocationFunction_( double& xLocationEstimate )
{
    // Coefficients of the 5th-order polynomial for L1 (James, 2006).
    double a_ = 2.0 * ( ( 2.0 * massParameter_ ) - 1.0 );
    double b_ = oneMinusMassParameterSquared_ - 4.0 * massParameter_
                * ( 1.0 - massParameter_ ) + massParameterSquared_;
    double c_ = 2.0 * massParameter_ * ( 1.0 - massParameter_ )
                * ( 1.0 - ( 2.0 * massParameter_ ) ) - 1.0;
    double d_ = massParameterSquared_ * oneMinusMassParameterSquared_
                + 2.0 * ( massParameterSquared_ + oneMinusMassParameterSquared_ )
                -4.0 * massParameterSquared_;
    double e_ = -pow( massParameter_, 3.0 ) - pow( 1.0 - massParameter_, 3.0 );

    // Return value of the location function.
    return pow( xLocationEstimate, 5.0 ) + a_ * pow( xLocationEstimate, 4.0 )
            + b_ * pow( xLocationEstimate, 3.0 ) +  c_ * pow( xLocationEstimate, 2.0 )
            + d_ * xLocationEstimate + e_;
}

//! Compute derivative of 5th-order polynomial for location of L1 libration point.
double LibrationPoint::computeL1FirstDerivativeLocationFunction_( double& xLocationEstimate )
{
    // Coefficients of first-derivative of the 5th-order polynomial for L1 (James, 2006).
    double a_ = 2.0 * ( ( 2.0 * massParameter_ ) - 1.0 );
    double b_ = oneMinusMassParameterSquared_ - 4.0 * massParameter_
                * ( 1.0 - massParameter_ ) + massParameterSquared_;
    double c_ = 2.0 * massParameter_ * ( 1.0 - massParameter_ )
                * ( 1.0 - ( 2.0 * massParameter_ ) ) - 1.0 ;
    double d_ = massParameterSquared_ * oneMinusMassParameterSquared_
                + 2.0 * ( massParameterSquared_ + oneMinusMassParameterSquared_ )
                -4.0 * massParameterSquared_;

    // Return value of first-deriavtive of the location function.
    return 5.0 * pow( xLocationEstimate, 4.0 ) + 4.0 * a_ * pow( xLocationEstimate, 3.0 )
            + 3.0 * b_ * pow( xLocationEstimate, 2.0 ) + 2.0 * c_ * xLocationEstimate + d_;
}

//! Compute 5th-order polynomial for location of L2 libration point.
double LibrationPoint::computeL2LocationFunction_( double& xLocationEstimate )
{
    // Coefficients of the 5th-order polynomial for L2 (James, 2006).
    double a_ = 2.0 * ( 2.0 * massParameter_ - 1.0 );
    double b_ = oneMinusMassParameterSquared_ - 4.0 * massParameter_
                * ( 1.0 - massParameter_ ) + massParameterSquared_;
    double c_ = 2.0 * massParameter_ * ( 1.0 - massParameter_ )
                * ( 1.0 -  2.0 * massParameter_ ) - 1.0 + 2.0 * massParameter_;
    double d_ = massParameterSquared_ * oneMinusMassParameterSquared_
                + 2.0 * ( massParameterSquared_  - oneMinusMassParameterSquared_ )
                + 4.0 * massParameterSquared_ - 8.0 * massParameter_ + 4.0;
    double e_ = - pow( 1.0 - massParameter_, 3.0 ) + pow( massParameter_, 3.0 ) ;

    // Return value of the location function.
    return pow( xLocationEstimate, 5.0 ) + a_ * pow( xLocationEstimate, 4.0 )
            + b_ * pow( xLocationEstimate, 3.0 ) +  c_ * pow( xLocationEstimate, 2.0 )
            + d_ * xLocationEstimate + e_;
}

//! Compute derivative of 5th-order polynomial for location of L2 libration point.
double LibrationPoint::computeL2FirstDerivativeLocationFunction_( double& xLocationEstimate )
{
    // Coefficients of the first-derivative of the 5th-order polynomial for L2 (James, 2006).
    double a_ = 2.0 * ( 2.0 * massParameter_ - 1.0 );
    double b_ = oneMinusMassParameterSquared_ - 4.0 * massParameter_
                * ( 1.0 - massParameter_ ) + massParameterSquared_;
    double c_ = 2.0 * massParameter_ * ( 1.0 - massParameter_ )
                * ( 1.0 -  2.0 * massParameter_ ) - 1.0 + 2.0 * massParameter_;
    double d_ = massParameterSquared_ * oneMinusMassParameterSquared_
                + 2.0 * ( massParameterSquared_  - oneMinusMassParameterSquared_ )
                + 4.0 * massParameterSquared_ - 8.0 * massParameter_ + 4.0;

    // Return value of first-deriavtive of the location function.
    return 5.0 * pow( xLocationEstimate, 4.0 ) + 4.0 * a_ * pow( xLocationEstimate, 3.0 )
            + 3.0 * b_ * pow( xLocationEstimate, 2.0 ) + 2.0 * c_ * xLocationEstimate + d_;
}

//! Compute 5th-order polynomial for location of L3 libration point.
double LibrationPoint::computeL3LocationFunction_( double& xLocationEstimate )
{
    // Coefficients of the 5th-order polynomial for L3 (James, 2006).
    double a_ = 2.0 * ( 2.0 * massParameter_ - 1.0 );
    double b_ = oneMinusMassParameterSquared_ - 4.0 * massParameter_
                * ( 1.0 - massParameter_ ) + massParameterSquared_;
    double c_ = 2.0 * massParameter_ * ( 1.0 - massParameter_ )
                * ( 1.0 -  2.0 * massParameter_ ) + 1.0 - 2.0 * massParameter_;
    double d_ = massParameterSquared_ * oneMinusMassParameterSquared_
                + 2.0 * ( massParameterSquared_ + oneMinusMassParameterSquared_ )
                -8.0 * massParameterSquared_ + 8.0 * massParameter_ - 4.0;
    double e_ = pow ( 1.0 - massParameter_, 3.0 ) - pow( massParameter_, 3.0 );

    // Return value of the location function.
    return pow( xLocationEstimate, 5.0 ) + a_ * pow( xLocationEstimate, 4.0 )
            + b_ * pow( xLocationEstimate, 3.0 ) +  c_ * pow( xLocationEstimate, 2.0 )
            + d_ * xLocationEstimate + e_;
}

//! Compute derivative of 5th-order polynomial for location of L3 libration point.
double LibrationPoint::computeL3FirstDerivativeLocationFunction_( double& xLocationEstimate )
{
    // Coefficients of the first-derivative of the 5th-order polynomial for L3 (James, 2006).
    double a_ = 2.0 * ( 2.0 * massParameter_ - 1.0 );
    double b_ = oneMinusMassParameterSquared_ - 4.0 * massParameter_
                * ( 1.0 - massParameter_ ) + massParameterSquared_;
    double c_ = 2.0 * massParameter_ * ( 1.0 - massParameter_ )
                * ( 1.0 -  2.0 * massParameter_ ) + 1.0 - 2.0 * massParameter_;
    double d_ = massParameterSquared_ * oneMinusMassParameterSquared_
                + 2.0 * ( massParameterSquared_ + oneMinusMassParameterSquared_ )
                -8.0 * massParameterSquared_ + 8.0 * massParameter_ - 4.0;

    // Return value of first-deriavtive of the location function.
    return 5.0 * pow( xLocationEstimate, 4.0 ) + 4.0 * a_ * pow( xLocationEstimate, 3.0 )
            + 3.0 * b_ * pow( xLocationEstimate, 2.0 ) + 2.0 * c_ * xLocationEstimate + d_;
}

// End of file.
