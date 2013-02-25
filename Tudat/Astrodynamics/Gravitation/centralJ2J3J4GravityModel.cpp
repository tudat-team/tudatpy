/*    Copyright (c) 2010-2013, Delft University of Technology
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
#include <stdexcept>

#include <boost/exception/all.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2J3J4GravityModel.h"

namespace tudat
{
namespace gravitation
{

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

    gravitationalAccelerationDueToJ4( basic_astrodynamics::xCartesianPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.x( )
                 - positionOfBodyExertingAcceleration.x( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ4( basic_astrodynamics::yCartesianPositionIndex )
            *= ( positionOfBodySubjectToAcceleration.y( )
                 - positionOfBodyExertingAcceleration.y( ) ) * factorForXAndYDirections;

    gravitationalAccelerationDueToJ4( basic_astrodynamics::zCartesianPositionIndex )
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

//! Constructor taking position-functions for bodies, and constant parameters of spherical
//! harmonics expansion.
CentralJ2J3J4GravitationalAccelerationModel::CentralJ2J3J4GravitationalAccelerationModel( 
        const StateFunction positionOfBodySubjectToAccelerationFunction,
        const double aGravitationalParameter,
        const double anEquatorialRadius,
        const double aJ2GravityCoefficient,
        const double aJ3GravityCoefficient,
        const double aJ4GravityCoefficient,
        const StateFunction positionOfBodyExertingAccelerationFunction )
    : Base( positionOfBodySubjectToAccelerationFunction,
            aGravitationalParameter,
            positionOfBodyExertingAccelerationFunction ),
      equatorialRadius( anEquatorialRadius ),
      j2GravityCoefficient( aJ2GravityCoefficient ),
      j3GravityCoefficient( aJ3GravityCoefficient ),
      j4GravityCoefficient( aJ4GravityCoefficient )
{
    Base::updateMembers( );
}

//! Get gravitational acceleration.
Eigen::Vector3d CentralJ2J3J4GravitationalAccelerationModel::getAcceleration( )
{
    // Sum and return constituent acceleration terms.
    return computeGravitationalAcceleration(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->positionOfBodyExertingAcceleration )
            + computeGravitationalAccelerationDueToJ2(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->j2GravityCoefficient,
                this->equatorialRadius,
                this->positionOfBodyExertingAcceleration )
            + computeGravitationalAccelerationDueToJ3(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->j3GravityCoefficient,
                this->equatorialRadius,
                this->positionOfBodyExertingAcceleration )
            + computeGravitationalAccelerationDueToJ4(
                this->positionOfBodySubjectToAcceleration,
                this->gravitationalParameter,
                this->j4GravityCoefficient,
                this->equatorialRadius,
                this->positionOfBodyExertingAcceleration );
}

} // namespace gravitation
} // namespace tudat
