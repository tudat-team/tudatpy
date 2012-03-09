/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      110716    L. van der Ham    First creation of code.
 *      110804    L. van der Ham    Changed use of approximatePlanetPositionsCircularCoplanar.
 *      111108    K. Kumar          Simplified conversion functions; changed "normalized" to
 *                                  "dimensionless".
 *      120307    K. Kumar          Updated Doxygen comments; made functions const-correct;
 *                                  updated to satisfy Tudat protocols; moved file.
 *
 *    References
 *        Wakker, K.F.,"Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 *    Position dimension-scale is distance between the two massive bodies in the CRTBP.
 *    Time dimension-scale is based on orbital period of 2*pi.
 *
 */

#include <cmath>

#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"

namespace tudat
{
namespace circular_restricted_three_body_problem
{

//! Convert dimensionless Cartesian state to dimensional units.
Eigen::VectorXd convertDimensionlessCartesianStateToDimensionalUnits(
        const Eigen::VectorXd &dimensionlessCartesianState,
        const double gravitationalParameterOfPrimaryBody, const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries )
{
    // Declare dimensional Cartesian state.
    Eigen::VectorXd dimensionalCartesianState( 6 );

    // Convert position to dimensional units.
    dimensionalCartesianState.segment( 0, 3 )
            = dimensionlessCartesianState.segment( 0, 3 ) * distanceBetweenPrimaries;

    // Convert velocity to dimensional units.
    dimensionalCartesianState.segment( 3, 3 )
            = dimensionlessCartesianState.segment( 3, 3 )
            * std::sqrt( ( gravitationalParameterOfPrimaryBody
                           + gravitationalParameterOfSecondaryBody ) / distanceBetweenPrimaries );

    // Return state in dimensional units.
    return dimensionalCartesianState;
}

//! Convert dimensionless time to dimensional units.
double convertDimensionlessTimeToDimensionalTime(
        const double timeInDimensionlessUnits, const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries )
{
    return timeInDimensionlessUnits * std::sqrt( std::pow( distanceBetweenPrimaries, 3.0 )
                                                 / ( gravitationalParameterOfPrimaryBody
                                                     + gravitationalParameterOfSecondaryBody ) );
}

} // namespace circular_restricted_three_body_problem
} // namespace tudat
