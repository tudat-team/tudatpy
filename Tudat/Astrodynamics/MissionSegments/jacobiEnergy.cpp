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
 *      110829    L. van der Ham    First creation of code.
 *      120221    L. van der Ham    Update of code and namespace.
 *      120227    K. Kumar          Updated code to new Tudat N Commandments; added enum to
 *                                  replace "magic numbers"; renamed file.
 *      120306    K. Kumar          Removed erroneous Doxygen comments.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 */

#include <cmath>

#include "Tudat/Astrodynamics/MissionSegments/jacobiEnergy.h"

namespace tudat
{
namespace circular_restricted_three_body_problem
{

//! Compute Jacobi energy.
double computeJacobiEnergy( const double massParameter,
                            const Eigen::VectorXd& stateInNormalizedUnits )
{
    using std::sqrt;
    using std::pow;

    double distanceToPrimaryBody = sqrt( pow( massParameter
                                              + stateInNormalizedUnits(
                                                  normalizedXPositionIndex ),  2.0 )
                                         + pow( stateInNormalizedUnits(
                                                    normalizedYPositionIndex ), 2.0 )
                                         + pow( stateInNormalizedUnits(
                                                    normalizedZPositionIndex ), 2.0 ) );

    double distanceToSecondaryBody = sqrt( pow( - ( 1.0 - massParameter
                                                    - stateInNormalizedUnits(
                                                        normalizedXPositionIndex ) ), 2.0 )
                                           + pow( stateInNormalizedUnits(
                                                      normalizedYPositionIndex ), 2.0 )
                                           + pow( stateInNormalizedUnits(
                                                      normalizedZPositionIndex ), 2.0 ) );

    double velocitySquared = pow( stateInNormalizedUnits( normalizedXVelocityIndex ), 2.0 ) +
            pow( stateInNormalizedUnits( normalizedYVelocityIndex ), 2.0 ) +
            pow( stateInNormalizedUnits( normalizedZVelocityIndex ), 2.0 );

    return pow( stateInNormalizedUnits( normalizedXPositionIndex ), 2.0 )
            + pow( stateInNormalizedUnits( normalizedYPositionIndex ), 2.0 )
            + 2.0 * ( 1.0 - massParameter ) / distanceToPrimaryBody
            + 2.0 * massParameter / distanceToSecondaryBody - velocitySquared;
}

} // namespace circular_restricted_three_body_problem
} // tudat
