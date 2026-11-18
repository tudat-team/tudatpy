/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/rotational_dynamics_partials/radiationPressureTorquePartial.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

RadiationPressureTorquePartial::RadiationPressureTorquePartial(
        const std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureTorque > torqueModel,
        const std::string& acceleratedBody,
        const std::string& acceleratingBody ):
    TorquePartial( acceleratedBody, acceleratingBody, basic_astrodynamics::radiation_pressure_torque ), torqueModel_( torqueModel ),
    targetModel_( torqueModel->getRadiationPressureAcceleration( )->getTargetModel( ) ),
    paneledTargetModel_( std::dynamic_pointer_cast< electromagnetism::PaneledRadiationPressureTargetModel >( targetModel_ ) ),
    cannonballTargetModel_( std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >( targetModel_ ) )
{
    if( paneledTargetModel_ != nullptr )
    {
        const std::string sourceName = torqueModel->getRadiationPressureAcceleration( )->getSourceName( );
        const auto selfShadowingPerSource = paneledTargetModel_->getSelfShadowingPerSources( );
        useIlluminatedFractions_ =
                selfShadowingPerSource.count( sourceName ) > 0 && selfShadowingPerSource.at( sourceName )->getMaximumNumberOfPixels( ) > 0;
        for( const auto& panel : paneledTargetModel_->getAllPanels( ) )
        {
            if( panel->getTrackedBody( ) != "" && panel->getTrackedBody( ) != sourceName )
            {
                throw std::runtime_error(
                        "Error, radiation pressure torque partial not supported for panels tracking a body other than the source, found "
                        "on " +
                        acceleratedBody );
            }
        }
    }
    currentRotationMatrixDerivativesWrtQuaternion_.resize( 4 );
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int > RadiationPressureTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    using namespace estimatable_parameters;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
            std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
    if( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case source_direction_radiation_pressure_scaling_factor:
                if( parameter->getParameterName( ).second.second == bodyExertingTorque_ )
                {
                    computeSourceDirectionTorque_ = true;
                    partialFunction = std::make_pair(
                            std::bind( &RadiationPressureTorquePartial::wrtSourceDirectionScaling, this, std::placeholders::_1 ), 1 );
                }
                break;
            case source_perpendicular_direction_radiation_pressure_scaling_factor:
                if( parameter->getParameterName( ).second.second == bodyExertingTorque_ )
                {
                    computeSourceDirectionTorque_ = true;
                    partialFunction = std::make_pair(
                            std::bind(
                                    &RadiationPressureTorquePartial::wrtPerpendicularSourceDirectionScaling, this, std::placeholders::_1 ),
                            1 );
                }
                break;
            case radiation_pressure_coefficient:
                if( cannonballTargetModel_ != nullptr )
                {
                    partialFunction = std::make_pair(
                            std::bind( &RadiationPressureTorquePartial::wrtRadiationPressureCoefficient, this, std::placeholders::_1 ), 1 );
                }
                break;
            case specular_reflectivity:
            case diffuse_reflectivity:
                if( paneledTargetModel_ != nullptr )
                {
                    partialFunction = std::make_pair( std::bind( &RadiationPressureTorquePartial::wrtReflectivity,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 parameter->getParameterName( ).second.second,
                                                                 parameter->getParameterName( ).first == specular_reflectivity ),
                                                      1 );
                }
                break;
            default:
                break;
        }
    }
    return partialFunction;
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int > RadiationPressureTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    if( parameter->getParameterName( ).first == estimatable_parameters::arc_wise_radiation_pressure_coefficient &&
        parameter->getParameterName( ).second.first == bodyUndergoingTorque_ && cannonballTargetModel_ != nullptr )
    {
        throw std::runtime_error(
                "Error, partial of radiation pressure torque w.r.t. arc-wise radiation pressure coefficient not yet supported" );
    }
    return std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
}

void RadiationPressureTorquePartial::wrtNonRotationalStateOfAdditionalBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    if( integratedStateType == propagators::translational_state )
    {
        if( stateReferencePoint.first == bodyUndergoingTorque_ )
        {
            partialMatrix.block( 0, 0, 3, 3 ) += currentPartialWrtTargetPosition_;
        }
        else if( stateReferencePoint.first == bodyExertingTorque_ )
        {
            partialMatrix.block( 0, 0, 3, 3 ) -= currentPartialWrtTargetPosition_;
        }
    }
}

void RadiationPressureTorquePartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        torqueModel_->updateMembers( currentTime );

        const std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration > acceleration =
                torqueModel_->getRadiationPressureAcceleration( );
        const std::string sourceName = acceleration->getSourceName( );
        const double sourceDirectionScaling = acceleration->getSourceDirectionScaling( );
        const double perpendicularScaling = acceleration->getPerpendicularSourceDirectionScaling( );

        // Source-target direction, in the global and body-fixed frames
        const Eigen::Quaterniond rotationToGlobalFrame = torqueModel_->getTargetRotationFromLocalToGlobalFrameFunction( )( );
        const Eigen::Matrix3d rotationToGlobalFrameMatrix = rotationToGlobalFrame.toRotationMatrix( );
        const Eigen::Vector3d relativePosition = acceleration->getTargetCenterPositionInSourceFrame( );
        const double distance = relativePosition.norm( );
        const Eigen::Vector3d direction = relativePosition / distance;
        const Eigen::Vector3d localDirection = rotationToGlobalFrameMatrix.transpose( ) * direction;
        currentLocalDirection_ = localDirection;

        // Forces, moment arms and partials of the forces w.r.t. the body-fixed source-target direction, per panel (one for a cannonball)
        std::vector< Eigen::Vector3d > forces;
        std::vector< Eigen::Vector3d > momentArms;
        std::vector< Eigen::Matrix3d > forcePartials;
        if( paneledTargetModel_ != nullptr )
        {
            forces = paneledTargetModel_->getPanelForces( sourceName );
            momentArms = paneledTargetModel_->getPanelCentroidMomentArms( );
            forcePartials.resize( forces.size( ), Eigen::Matrix3d::Zero( ) );

            const double radiationPressure = acceleration->getCurrentRadiationPressure( );
            const std::vector< double >& cosines = paneledTargetModel_->getSurfacePanelCosines( sourceName );
            const std::vector< Eigen::Vector3d >& normals = paneledTargetModel_->getSurfaceNormals( );
            const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& panels = paneledTargetModel_->getAllPanels( );

            // Illuminated fractions differ from one only with self-shadowing
            std::vector< double > illuminatedFractions( forces.size( ), 1.0 );
            if( useIlluminatedFractions_ )
            {
                illuminatedFractions = paneledTargetModel_->getIlluminatedPanelFractions( sourceName );
            }

            for( unsigned int i = 0; i < forces.size( ); i++ )
            {
                const double forceScaling = radiationPressure * panels.at( i )->getPanelArea( ) * illuminatedFractions.at( i );
                if( panels.at( i )->getTrackedBody( ) == sourceName )
                {
                    // Panel normal along the source-target direction: the force rotates with the direction, at constant magnitude
                    forcePartials.at( i ) = forces.at( i ).dot( localDirection ) * Eigen::Matrix3d::Identity( );
                }
                else if( cosines.at( i ) > 0.0 && forceScaling > 0.0 )
                {
                    // Derivative w.r.t. the incoming direction: the unit vector to the source is its negative
                    forcePartials.at( i ) = forceScaling *
                            panels.at( i )->getReflectionLaw( )->evaluateReactionVectorDerivativeWrtTargetPosition(
                                    normals.at( i ),
                                    localDirection,
                                    cosines.at( i ),
                                    forces.at( i ) / forceScaling,
                                    -Eigen::Matrix3d::Identity( ),
                                    -normals.at( i ).transpose( ) );
                }
            }
            currentIlluminatedFractions_ = illuminatedFractions;
        }
        else
        {
            // Cannonball force along the source-target direction, acting at the body-fixed origin
            const Eigen::Vector3d force =
                    rotationToGlobalFrameMatrix.transpose( ) * targetModel_->getCurrentRadiationPressureForce( sourceName );
            forces = { force };
            momentArms = { -torqueModel_->getCenterOfMassFunction( )( ) };
            forcePartials = { force.dot( localDirection ) * Eigen::Matrix3d::Identity( ) };
        }

        // Unscaled torque, and its partial w.r.t. the body-fixed source-target direction
        Eigen::Vector3d torque = Eigen::Vector3d::Zero( );
        Eigen::Matrix3d torqueDirectionPartial = Eigen::Matrix3d::Zero( );

        // Torque of the force components along the source-target direction, and its partial (only needed for scaling)
        const bool computeSourceDirectionTorque = computeSourceDirectionTorque_ || sourceDirectionScaling != perpendicularScaling;
        Eigen::Vector3d sourceDirectionTorque = Eigen::Vector3d::Zero( );
        Eigen::Matrix3d sourceDirectionTorqueDirectionPartial = Eigen::Matrix3d::Zero( );
        for( unsigned int i = 0; i < forces.size( ); i++ )
        {
            const Eigen::Matrix3d momentArmCrossProductMatrix = linear_algebra::getCrossProductMatrix( momentArms.at( i ) );
            torque += momentArms.at( i ).cross( forces.at( i ) );
            torqueDirectionPartial += momentArmCrossProductMatrix * forcePartials.at( i );

            if( computeSourceDirectionTorque )
            {
                const Eigen::Vector3d momentArmCrossDirection = momentArms.at( i ).cross( localDirection );
                const double sourceDirectionForce = forces.at( i ).dot( localDirection );
                sourceDirectionTorque += sourceDirectionForce * momentArmCrossDirection;
                sourceDirectionTorqueDirectionPartial +=
                        momentArmCrossDirection * ( localDirection.transpose( ) * forcePartials.at( i ) + forces.at( i ).transpose( ) ) +
                        sourceDirectionForce * momentArmCrossProductMatrix;
            }
        }

        // Scaled torque: perpendicularScaling * torque + ( sourceDirectionScaling - perpendicularScaling ) * sourceDirectionTorque
        const Eigen::Matrix3d directionPartial = perpendicularScaling * torqueDirectionPartial +
                ( sourceDirectionScaling - perpendicularScaling ) * sourceDirectionTorqueDirectionPartial;
        currentPartialWrtSourceDirectionScaling_ = sourceDirectionTorque;
        currentPartialWrtPerpendicularScaling_ = torque - sourceDirectionTorque;

        // The direction is a unit vector: only its change perpendicular to itself contributes
        const Eigen::Matrix3d unitDirectionPartial =
                directionPartial * ( Eigen::Matrix3d::Identity( ) - localDirection * localDirection.transpose( ) );

        // Orientation dependence, through the rotation of the direction to the body-fixed frame
        linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion(
                linear_algebra::convertQuaternionToVectorFormat( rotationToGlobalFrame ), currentRotationMatrixDerivativesWrtQuaternion_ );
        for( unsigned int i = 0; i < 4; i++ )
        {
            currentPartialWrtQuaternion_.block( 0, i, 3, 1 ) =
                    unitDirectionPartial * currentRotationMatrixDerivativesWrtQuaternion_.at( i ).transpose( ) * direction;
        }

        // Position dependence, through the direction and the inverse-square radiation pressure
        currentPartialWrtTargetPosition_ = unitDirectionPartial * rotationToGlobalFrameMatrix.transpose( ) / distance -
                2.0 * torqueModel_->getTorque( ) * direction.transpose( ) / distance;

        currentTime_ = currentTime;
    }
}

void RadiationPressureTorquePartial::wrtReflectivity( Eigen::MatrixXd& partial, const std::string& panelTypeId, const bool isSpecular )
{
    const std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration > acceleration =
            torqueModel_->getRadiationPressureAcceleration( );
    const std::string sourceName = acceleration->getSourceName( );
    const double radiationPressure = acceleration->getCurrentRadiationPressure( );
    const double sourceDirectionScaling = acceleration->getSourceDirectionScaling( );
    const double perpendicularScaling = acceleration->getPerpendicularSourceDirectionScaling( );

    const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& panels = paneledTargetModel_->getAllPanels( );
    const std::vector< double >& cosines = paneledTargetModel_->getSurfacePanelCosines( sourceName );
    const std::vector< Eigen::Vector3d >& normals = paneledTargetModel_->getSurfaceNormals( );
    const std::vector< Eigen::Vector3d >& momentArms = paneledTargetModel_->getPanelCentroidMomentArms( );

    // Unscaled torque of the partials of the forces of the panels of this type w.r.t. the reflectivity, and of their components
    // along the source-target direction
    Eigen::Vector3d torquePartial = Eigen::Vector3d::Zero( );
    Eigen::Vector3d sourceDirectionTorquePartial = Eigen::Vector3d::Zero( );
    for( unsigned int i = 0; i < panels.size( ); i++ )
    {
        if( panels.at( i )->getPanelTypeId( ) == panelTypeId && cosines.at( i ) > 0.0 )
        {
            const Eigen::Vector3d reactionVectorPartial = isSpecular
                    ? panels.at( i )->getReflectionLaw( )->evaluateReactionVectorPartialWrtSpecularReflectivity( normals.at( i ),
                                                                                                                 currentLocalDirection_ )
                    : panels.at( i )->getReflectionLaw( )->evaluateReactionVectorPartialWrtDiffuseReflectivity( normals.at( i ),
                                                                                                                currentLocalDirection_ );
            const Eigen::Vector3d forcePartial = radiationPressure * panels.at( i )->getPanelArea( ) *
                    currentIlluminatedFractions_.at( i ) * cosines.at( i ) * reactionVectorPartial;
            torquePartial += momentArms.at( i ).cross( forcePartial );
            sourceDirectionTorquePartial += forcePartial.dot( currentLocalDirection_ ) * momentArms.at( i ).cross( currentLocalDirection_ );
        }
    }
    partial = perpendicularScaling * torquePartial + ( sourceDirectionScaling - perpendicularScaling ) * sourceDirectionTorquePartial;
}

void RadiationPressureTorquePartial::wrtRadiationPressureCoefficient( Eigen::MatrixXd& partial )
{
    // The cannonball torque is proportional to the radiation pressure coefficient
    partial = torqueModel_->getTorque( ) / cannonballTargetModel_->getCoefficient( );
}

}  // namespace acceleration_partials

}  // namespace tudat
