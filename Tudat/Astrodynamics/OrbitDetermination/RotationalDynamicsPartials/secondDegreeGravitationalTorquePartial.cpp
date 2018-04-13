/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/secondDegreeGravitationalTorquePartial.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"


namespace tudat
{

namespace acceleration_partials
{

Eigen::Matrix< double, 3, 4 > getPartialDerivativeOfSecondDegreeGravitationalTorqueWrtQuaternion(
        const double premultiplier,
        const Eigen::Matrix3d& inertiaTensor,
        const Eigen::Vector3d& bodyFixedRelativePosition,
        const Eigen::Vector3d& inertialRelativePosition,
        const std::vector< Eigen::Matrix3d > derivativeOfRotationMatrixWrtQuaternions )
{
    Eigen::Matrix3d scalingMatrix = linear_algebra::getCrossProductMatrix(
                bodyFixedRelativePosition ) *  inertiaTensor - linear_algebra::getCrossProductMatrix(
                inertiaTensor * bodyFixedRelativePosition );
    Eigen::Matrix< double, 3, 4 > partialOfBodyFixedPositionWrtQuaternion = Eigen::Matrix< double, 3, 4 >::Zero( );
    for( unsigned int i = 0; i < derivativeOfRotationMatrixWrtQuaternions.size( ); i++ )
    {
        partialOfBodyFixedPositionWrtQuaternion.block( 0, i, 3, 1 ) =
             derivativeOfRotationMatrixWrtQuaternions.at( i ).transpose( ) * inertialRelativePosition;
    }
    return premultiplier * scalingMatrix * partialOfBodyFixedPositionWrtQuaternion;
}

std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
SecondDegreeGravitationalTorquePartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;

    // Check dependencies.
    if( parameter->getParameterName( ).first ==  estimatable_parameters::gravitational_parameter && (
                parameter->getParameterName( ).second.first == bodyExertingTorque_ ) )
    {
        // If parameter is gravitational parameter, check and create dependency function .
        partialFunctionPair = std::make_pair(
                    boost::bind( &SecondDegreeGravitationalTorquePartial::wrtGravitationalParameterOfCentralBody,
                                       this, _1 ), 1 );
    }
    else
    {
        partialFunctionPair = std::make_pair( boost::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }

    return partialFunctionPair;
}

//! Function to calculate central gravity partial w.r.t. central body gravitational parameter
void SecondDegreeGravitationalTorquePartial::wrtGravitationalParameterOfCentralBody(
        Eigen::MatrixXd& gravitationalParameterPartial )
{
    if( torqueModel_->getCurrentGravitationalParameterOfAttractingBody( ) == 0.0 )
    {
        throw std::runtime_error( "Error when calculating partial of SecondDegreeGravitationalTorquePartial w.r.t. mu: mu=0." );
    }
    gravitationalParameterPartial =
            torqueModel_->getTorque( ) / torqueModel_->getCurrentGravitationalParameterOfAttractingBody( );
}


}

}
