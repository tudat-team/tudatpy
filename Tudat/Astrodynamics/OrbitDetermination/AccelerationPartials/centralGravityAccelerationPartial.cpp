/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/centralGravityAccelerationPartial.h"


namespace tudat
{

namespace acceleration_partials
{

//! Calculates partial derivative of point mass gravitational acceleration wrt the position of body undergoing acceleration.
Eigen::Matrix3d calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
        const Eigen::Vector3d& acceleratedBodyPosition,
        const Eigen::Vector3d& acceleratingBodyPositions,
        double gravitationalParameter )
{
    // Calculate relative position
    Eigen::Vector3d relativePosition = acceleratedBodyPosition - acceleratingBodyPositions;

    // Calculate partial (Montenbruck & Gill, Eq. 7.56)
    double relativePositionNorm = relativePosition.norm( );
    double invSquareOfPositionNorm = 1.0 / ( relativePositionNorm * relativePositionNorm );
    double invCubeOfPositionNorm = invSquareOfPositionNorm / relativePositionNorm;
    Eigen::Matrix3d partialMatrix = -gravitationalParameter *
            ( Eigen::Matrix3d::Identity( ) * invCubeOfPositionNorm -
              ( 3.0 * invSquareOfPositionNorm * invCubeOfPositionNorm ) * relativePosition * relativePosition.transpose( ) );

    return partialMatrix;
}

//! Calculates partial derivative of point mass gravitational acceleration wrt gravitational parameter of the central body.
Eigen::Vector3d computePartialOfCentralGravityWrtGravitationalParameter( const Eigen::Vector3d& acceleratedBodyPosition,
                                                                         const Eigen::Vector3d& acceleratingBodyPosition)
{
    // Calculate relative position
    Eigen::Vector3d relativePosition = acceleratingBodyPosition - acceleratedBodyPosition;

    // Calculate partial (Montenbruck & Gill, Eq. 7.76)
    double positionNorm = relativePosition.norm( );
    Eigen::Vector3d partialMatrix = relativePosition / ( positionNorm * positionNorm * positionNorm );
    return partialMatrix;
}

//! Calculates partial derivative of point mass gravitational acceleration wrt gravitational parameter of the central body.
Eigen::Vector3d computePartialOfCentralGravityWrtGravitationalParameter( const Eigen::Vector3d& gravitationalAcceleration,
                                                                         const double gravitationalParameter )
{
    return gravitationalAcceleration / gravitationalParameter;
}

//! Constructor
CentralGravitationPartial::CentralGravitationPartial(
        const boost::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > gravitationalAcceleration,
        const std::string acceleratedBody,
        const std::string acceleratingBody ):
    AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::central_gravity )
{
    accelerationUpdateFunction_ =
            boost::bind( &basic_astrodynamics::AccelerationModel< Eigen::Vector3d>::updateMembers, gravitationalAcceleration, _1 );

    gravitationalParameterFunction_ = gravitationalAcceleration->getGravitationalParameterFunction( );
    centralBodyState_ = gravitationalAcceleration->getStateFunctionOfBodyExertingAcceleration( );
    acceleratedBodyState_ = gravitationalAcceleration->getStateFunctionOfBodyUndergoingAcceleration( );
    accelerationUsesMutualAttraction_ = gravitationalAcceleration->getIsMutualAttractionUsed( );
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
CentralGravitationPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )

{
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;

    // Check dependencies.
    if( parameter->getParameterName( ).first ==  estimatable_parameters::gravitational_parameter )
    {
        // If parameter is gravitational parameter, check and create dependency function .
        partialFunctionPair = getGravitationalParameterPartialFunction( parameter->getParameterName( ) );
    }
    else
    {
        partialFunctionPair = std::make_pair( boost::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }

    return partialFunctionPair;
}

//! Function to create a function returning the current partial w.r.t. a gravitational parameter.
std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
CentralGravitationPartial::getGravitationalParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterId )
{
    boost::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    // Check if parameter is gravitational parameter.
    if( parameterId.first ==  estimatable_parameters::gravitational_parameter )
    {
        // Check if parameter body is central body.
        if( parameterId.second.first == acceleratingBody_ )
        {
            partialFunction = boost::bind( &CentralGravitationPartial::wrtGravitationalParameterOfCentralBody,
                                           this, _1 );
            numberOfColumns = 1;

        }

        // Check if parameter body is accelerated body, and if the mutual acceleration is used.
        if( parameterId.second.first == acceleratedBody_ )
        {
            if( accelerationUsesMutualAttraction_ )
            {
                partialFunction = boost::bind( &CentralGravitationPartial::wrtGravitationalParameterOfCentralBody,
                                               this, _1 );
                numberOfColumns = 1;
            }
        }
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

//! Function to calculate central gravity partial w.r.t. central body gravitational parameter
void CentralGravitationPartial::wrtGravitationalParameterOfCentralBody( Eigen::MatrixXd& gravitationalParameterPartial )
{
    gravitationalParameterPartial = computePartialOfCentralGravityWrtGravitationalParameter(
                currentAcceleratedBodyState_, currentCentralBodyState_ );
}



}

}
