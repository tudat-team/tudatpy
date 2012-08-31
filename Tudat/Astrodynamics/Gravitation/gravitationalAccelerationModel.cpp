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
 *      120209    K. Kumar          File created.
 *
 *    References
 *      Wakker.
 *      Melman, J. PhD Thesis, 2012.
 *
 */

#define DEGREE 1
#define COEFFICIENTVALUE

#include <cmath>
#include <stdexcept>

#include <boost/exception/all.hpp>

#include "Tudat/Astrodynamics/Gravitation/gravitationalAccelerationModel.h"

namespace tudat
{
namespace astrodynamics
{
namespace acceleration_models
{

//! Compute gravitational acceleration.
Eigen::Vector3d computeGravitationalAcceleration(
        const double universalGravitationalConstant,
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double massOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    return computeGravitationalAcceleration(
                positionOfBodySubjectToAcceleration,
                universalGravitationalConstant * massOfBodyExertingAcceleration,
                positionOfBodyExertingAcceleration );
}

//! Compute gravitational acceleration.
Eigen::Vector3d computeGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    return -gravitationalParameterOfBodyExertingAcceleration
            * ( positionOfBodySubjectToAcceleration - positionOfBodyExertingAcceleration )
            / std::pow( ( positionOfBodySubjectToAcceleration
                          - positionOfBodyExertingAcceleration ).norm( ), 3.0 );
}

//! Compute gravitational acceleration due to J2.
Eigen::Vector3d computeGravitationalAccelerationDueToJ2(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double j2CoefficientOfGravityField,
        const double effectiveRadiusOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    // Set constant values reused for optimal computation of acceleration components.
    const double distanceBetweenBodies = ( positionOfBodySubjectToAcceleration
                                           - positionOfBodyExertingAcceleration ).norm( );

    const double preMultiplier = -gravitationalParameterOfBodyExertingAcceleration
            / std::pow( distanceBetweenBodies, 4.0 ) * 1.5 * j2CoefficientOfGravityField
            * effectiveRadiusOfBodyExertingAcceleration
            * effectiveRadiusOfBodyExertingAcceleration;

    const double scaledZCoordinate = ( positionOfBodySubjectToAcceleration.z( )
                                       - positionOfBodyExertingAcceleration.z( ) )
            / distanceBetweenBodies;

    const double scaledZCoordinateSquared = scaledZCoordinate * scaledZCoordinate;

    const double factorForXAndYDirections = ( 1.0 - 5.0 * scaledZCoordinateSquared )
            / distanceBetweenBodies;

    // Compute components of acceleration due to J2-effect.
    Eigen::Vector3d gravitationalAccelerationDueToJ2 = Eigen::Vector3d::Constant( preMultiplier );

    gravitationalAccelerationDueToJ2( cartesianXPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.x( )
                 - positionOfBodyExertingAcceleration.x( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ2( cartesianYPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.y( )
                 - positionOfBodyExertingAcceleration.y( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ2( cartesianZPositionIndex )
            *= ( 3.0 - 5.0 * scaledZCoordinateSquared ) * scaledZCoordinate;

    return gravitationalAccelerationDueToJ2;
}

Eigen::Vector3d computeGravitationalAccelerationDueToJ3(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double j3CoefficientOfGravityField,
        const double effectiveRadiusOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    // Set constant values reused for optimal computation of acceleration components.
    const double distanceBetweenBodies = ( positionOfBodySubjectToAcceleration
                                           - positionOfBodyExertingAcceleration ).norm( );

    const double preMultiplier = -gravitationalParameterOfBodyExertingAcceleration
            / std::pow( distanceBetweenBodies, 5.0 ) * 2.5 * j3CoefficientOfGravityField
            * std::pow( effectiveRadiusOfBodyExertingAcceleration, 3.0 );

    const double scaledZCoordinate = ( positionOfBodySubjectToAcceleration.z( )
                                       - positionOfBodyExertingAcceleration.z( ) )
            / distanceBetweenBodies;

    const double scaledZCoordinateSquared = scaledZCoordinate * scaledZCoordinate;

    const double factorForXAndYDirections = ( 3.0 - 7.0 * scaledZCoordinateSquared )
            * scaledZCoordinate / distanceBetweenBodies;

    // Compute components of acceleration due to J3-effect.
    Eigen::Vector3d gravitationalAccelerationDueToJ3 = Eigen::Vector3d::Constant( preMultiplier );

    gravitationalAccelerationDueToJ3( cartesianXPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.x( )
                 - positionOfBodyExertingAcceleration.x( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ3( cartesianYPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.y( )
                 - positionOfBodyExertingAcceleration.y( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ3( cartesianZPositionIndex )
            *= ( -0.6 + 6.0 * scaledZCoordinateSquared
                 - 7.0 * scaledZCoordinateSquared * scaledZCoordinateSquared );

    return gravitationalAccelerationDueToJ3;
}

//! Compute gravitational acceleration due to J4.
Eigen::Vector3d computeGravitationalAccelerationDueToJ4(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const double j4CoefficientOfGravityField,
        const double effectiveRadiusOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    // Set constant values reused for optimal computation of acceleration components.
    const double distanceBetweenBodies = ( positionOfBodySubjectToAcceleration
                                           - positionOfBodyExertingAcceleration ).norm( );

    const double preMultiplier = gravitationalParameterOfBodyExertingAcceleration
            / std::pow( distanceBetweenBodies, 6.0 ) * 4.375 * j4CoefficientOfGravityField
            * std::pow( effectiveRadiusOfBodyExertingAcceleration, 4.0 );

    const double scaledZCoordinate = ( positionOfBodySubjectToAcceleration.z( )
                                       - positionOfBodyExertingAcceleration.z( ) )
            / distanceBetweenBodies;

    const double scaledZCoordinateSquared = scaledZCoordinate * scaledZCoordinate;

    const double scaledZCoordinateToPower4 = scaledZCoordinateSquared * scaledZCoordinateSquared;

    const double factorForXAndYDirections = ( 3.0 / 7.0 - 6.0 * scaledZCoordinateSquared
                                              + 9.0 * scaledZCoordinateToPower4 )
            / distanceBetweenBodies;

    // Compute components of acceleration due to J4-effect.
    Eigen::Vector3d gravitationalAccelerationDueToJ4 = Eigen::Vector3d::Constant( preMultiplier );

    gravitationalAccelerationDueToJ4( cartesianXPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.x( )
                 - positionOfBodyExertingAcceleration.x( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ4( cartesianYPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.y( )
                 - positionOfBodyExertingAcceleration.y( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ4( cartesianZPositionIndex )
            *= ( 15.0 / 7.0 - 10.0 * scaledZCoordinateSquared + 9.0 * scaledZCoordinateToPower4 )
            * scaledZCoordinate;

    return gravitationalAccelerationDueToJ4;
}

//! Compute gravitational acceleration zonal sum.
Eigen::Vector3d computeGravitationalAccelerationZonalSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const std::map< int, double > zonalCoefficientsOfGravityField,
        const double effectiveRadiusOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    // Check that only coefficients for the gravity field are given up to J2 (i.e., size of
    // vector is 3 at max), else throw an error.
    if ( zonalCoefficientsOfGravityField.size( ) > 3 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Currently, accelerations can only be computed up to J4." ) ) );
    }

    // Check if position of body subject to acceleration falls within the effective radius of
    // the central body. If so, throw an error.
    if ( ( positionOfBodySubjectToAcceleration - positionOfBodyExertingAcceleration ).norm( )
         < effectiveRadiusOfBodyExertingAcceleration )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Position of body subject to acceleration is within effective radius."
                            ) ) );
    }

    // Set gravitational acceleration sum equal to central term contribution.
    Eigen::Vector3d gravitationalAccelerationSum
            = computeGravitationalAcceleration(
                positionOfBodySubjectToAcceleration,
                gravitationalParameterOfBodyExertingAcceleration,
                positionOfBodyExertingAcceleration );

    // Add contributions from zonal terms. The switch statement checks
    for ( std::map< int, double >::const_iterator mapIterator
          = zonalCoefficientsOfGravityField.begin( );
          mapIterator != zonalCoefficientsOfGravityField.end( ); mapIterator++ )
    {
        switch ( mapIterator->first )
        {
        case 2:

            gravitationalAccelerationSum += computeGravitationalAccelerationDueToJ2(
                        positionOfBodySubjectToAcceleration,
                        gravitationalParameterOfBodyExertingAcceleration,
                        mapIterator->second,
                        effectiveRadiusOfBodyExertingAcceleration,
                        positionOfBodyExertingAcceleration );

            break;

        case 3:

            gravitationalAccelerationSum += computeGravitationalAccelerationDueToJ3(
                        positionOfBodySubjectToAcceleration,
                        gravitationalParameterOfBodyExertingAcceleration,
                        mapIterator->second,
                        effectiveRadiusOfBodyExertingAcceleration,
                        positionOfBodyExertingAcceleration );

            break;

        case 4:

            gravitationalAccelerationSum += computeGravitationalAccelerationDueToJ4(
                        positionOfBodySubjectToAcceleration,
                        gravitationalParameterOfBodyExertingAcceleration,
                        mapIterator->second,
                        effectiveRadiusOfBodyExertingAcceleration,
                        positionOfBodyExertingAcceleration );

            break;

        default:

            boost::throw_exception(
                        boost::enable_error_info(
                            std::runtime_error(
                                "Degree must be 2, 3, or 4 in current implementation." ) ) );
        };
    }

    // Return total gravitational acceleration computed.
    return gravitationalAccelerationSum;
}

} // namespace acceleration_models
} // namespace astrodynamics
} // namespace tudat
