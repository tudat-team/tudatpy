/*    Copyright (c) 2010 Delft University of Technology.
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
 *      110905    L. van der Ham    First creation of code.
 *      110922    L. van der Ham    Added computation of Jacobian energy constant.
 *      111107    L. van der Ham    Removed option planar.
 *      111110    K. Kumar          Minor comment and naming modifications; error removed from
 *                                  z-direction equation.
 *      120221    L. van der Ham    Removed computation Jacobian energy.
 *      120309    K. Kumar          Updated code to latest Tudat standards; updated
 *                                  computeStateDerivative() function.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 */

#include <cmath>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>

#include "Tudat/Astrodynamics/Gravitation/stateDerivativeCircularRestrictedThreeBodyProblem.h"

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Compute state derivative.
Eigen::VectorXd StateDerivativeCircularRestrictedThreeBodyProblem::
computeStateDerivative( const double time, const Eigen::VectorXd& cartesianState )
{
    using std::pow;
    using orbital_element_conversions::xPositionIndex;
    using orbital_element_conversions::yPositionIndex;
    using orbital_element_conversions::zPositionIndex;
    using orbital_element_conversions::xVelocityIndex;
    using orbital_element_conversions::yVelocityIndex;
    using orbital_element_conversions::zVelocityIndex;

    // Compute relative position vectors.
    double normDistanceToPrimaryBodyCubed_ = pow( pow( cartesianState( xPositionIndex )
                                                       + massParameter_, 2.0 )
                                           + pow( cartesianState( yPositionIndex ), 2.0 )
                                           + pow( cartesianState( zPositionIndex ), 2.0 ), 1.5 );

    double normDistanceToSecondaryBodyCubed_ = pow( pow( cartesianState( xPositionIndex )
                                                     - ( 1.0 - massParameter_ ), 2.0 )
                                             + pow( cartesianState( yPositionIndex ), 2.0 )
                                             + pow( cartesianState ( zPositionIndex ), 2.0 ),
                                                    1.5 );

    // Compute derivative of state.
    Eigen::VectorXd stateDerivative( 6 );
    unsigned int xAccelerationIndex = 3;
    unsigned int yAccelerationIndex = 4;
    unsigned int zAccelerationIndex = 5;

    stateDerivative.segment( xPositionIndex, 3 ) = cartesianState.segment( xVelocityIndex, 3 );

    stateDerivative( xAccelerationIndex ) = cartesianState( xPositionIndex )
            - ( ( 1.0 - massParameter_ ) / normDistanceToPrimaryBodyCubed_ )
            * ( cartesianState( xPositionIndex ) + massParameter_ )
            - ( massParameter_ / normDistanceToSecondaryBodyCubed_ )
            * ( cartesianState( xPositionIndex ) - ( 1.0 - massParameter_ ) )
            + 2.0 * cartesianState( yVelocityIndex );
    stateDerivative( yAccelerationIndex ) = cartesianState( yPositionIndex )
            * ( 1.0 - ( ( 1.0 - massParameter_ ) / normDistanceToPrimaryBodyCubed_ )
                - ( massParameter_ / normDistanceToSecondaryBodyCubed_ ) )
            - 2.0 * cartesianState( xVelocityIndex );
    stateDerivative( zAccelerationIndex ) = -cartesianState( zPositionIndex )
            * ( ( ( 1.0 - massParameter_ ) / normDistanceToPrimaryBodyCubed_ )
                + ( massParameter_ / normDistanceToSecondaryBodyCubed_ ) );

    // Return computed state derivative.
    return stateDerivative;
}

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat
