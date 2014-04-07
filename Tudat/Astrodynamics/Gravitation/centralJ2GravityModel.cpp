/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      121105    K. Kumar          File created from content in other files.
 *      121210    D. Dirkx          Added function implementations for class.
 *
 *    References
 *
 *    Notes
 *
 */

#include <cmath>

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h"

namespace tudat
{
namespace gravitation
{

//! Compute gravitational acceleration due to J2.
Eigen::Vector3d computeGravitationalAccelerationDueToJ2(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double equatorialRadiusOfBodyExertingAcceleration,
        const double j2CoefficientOfGravityField,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    // Set constant values reused for optimal computation of acceleration components.
    const double distanceBetweenBodies = ( positionOfBodySubjectToAcceleration
                                           - positionOfBodyExertingAcceleration ).norm( );

    const double preMultiplier = -gravitationalParameterOfBodyExertingAcceleration
            / std::pow( distanceBetweenBodies, 4.0 ) * 1.5 * j2CoefficientOfGravityField
            * equatorialRadiusOfBodyExertingAcceleration
            * equatorialRadiusOfBodyExertingAcceleration;

    const double scaledZCoordinate = ( positionOfBodySubjectToAcceleration.z( )
                                       - positionOfBodyExertingAcceleration.z( ) )
            / distanceBetweenBodies;

    const double scaledZCoordinateSquared = scaledZCoordinate * scaledZCoordinate;

    const double factorForXAndYDirections = ( 1.0 - 5.0 * scaledZCoordinateSquared )
            / distanceBetweenBodies;

    // Compute components of acceleration due to J2-effect.
    Eigen::Vector3d gravitationalAccelerationDueToJ2 = Eigen::Vector3d::Constant( preMultiplier );

    gravitationalAccelerationDueToJ2( basic_astrodynamics::xCartesianPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.x( )
                 - positionOfBodyExertingAcceleration.x( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ2( basic_astrodynamics::yCartesianPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.y( )
                 - positionOfBodyExertingAcceleration.y( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ2( basic_astrodynamics::zCartesianPositionIndex )
            *= ( 3.0 - 5.0 * scaledZCoordinateSquared ) * scaledZCoordinate;

    return gravitationalAccelerationDueToJ2;
}

//! Get gravitational acceleration.
Eigen::Vector3d CentralJ2GravitationalAccelerationModel::getAcceleration( )
{
    // Sum and return constituent acceleration terms.
    return computeGravitationalAcceleration(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->positionOfBodyExertingAcceleration )
            + computeGravitationalAccelerationDueToJ2(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->equatorialRadius,
                this->j2GravityCoefficient,
                this->positionOfBodyExertingAcceleration );
}

} // namespace gravitation
} // namespace tudat
