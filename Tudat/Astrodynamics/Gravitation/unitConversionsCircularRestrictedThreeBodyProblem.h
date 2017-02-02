/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      110717    L. van der Ham    File created.
 *      111108    K. Kumar          Simplified conversion functions; changed "normalized" to
 *                                  "dimensionless".
 *      120307    K. Kumar          Updated Doxygen comments; made functions const-correct;
 *                                  updated to satisfy Tudat protocols; moved file.
 *      130121    K. Kumar          Updated VectorXd to Vector6d.
 *
 *    References
 *        Wakker, K.F.,"Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 *    Notes
 *      Position dimension-scale is distance between the two massive bodies in the CRTBP.
 *      Time dimension-scale is based on orbital period of 2*pi.
 *
 */

#ifndef TUDAT_UNIT_CONVERSIONS_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
#define TUDAT_UNIT_CONVERSIONS_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Convert dimensionless Cartesian state to dimensional units.
/*!
 * Convert dimensionless Cartesian state to dimensional units.
 * \param dimensionlessCartesianState Dimensionless Cartesian state vector.
 *          Order is important!
 *          dimensionlessCartesianState( 0 ) = x-position coordinate,                           [m]
 *          dimensionlessCartesianState( 1 ) = y-position coordinate,                           [m]
 *          dimensionlessCartesianState( 2 ) = z-position coordinate,                           [m]
 *          dimensionlessCartesianState( 3 ) = x-velocity coordinate,                         [m/s]
 *          dimensionlessCartesianState( 4 ) = y-velocity coordinate,                         [m/s]
 *          dimensionlessCartesianState( 5 ) = z-velocity coordinate.                         [m/s]
 * \param gravitationalParameterOfPrimaryBody Gravitational parameter of primary body.   [m^3 s^-2]
 * \param gravitationalParameterOfSecondaryBody Gravitational parameter of secondary body.
 *                                                                                       [m^3 s^-2]
 * \param distanceBetweenPrimaries Distance between primaries.                                  [m]
 * \return Dimensional Cartesian state.
 */
Eigen::VectorXd convertDimensionlessCartesianStateToDimensionalUnits(
        const Eigen::Vector6d& dimensionlessCartesianState,
        const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries );

//! Convert dimensionless time to dimensional units.
/*!
 * Convert dimensionless time to dimensional units.
 * \param timeInDimensionlessUnits Time in normalized units.
 * \param gravitationalParameterOfPrimaryBody Gravitational parameter of primary body.   [m^3 s^-2]
 * \param gravitationalParameterOfSecondaryBody Gravitational parameter of secondary body.
 *                                                                                       [m^3 s^-2]
 * \param distanceBetweenPrimaries Distance between primaries.                                  [m]
 * \return Dimensional time.                                                                    [s]
 */
double convertDimensionlessTimeToDimensionalTime(
        const double timeInDimensionlessUnits, const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries );

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace tudat

#endif // TUDAT_UNIT_CONVERSIONS_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
