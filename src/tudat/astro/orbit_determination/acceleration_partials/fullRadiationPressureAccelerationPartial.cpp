/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/fullRadiationPressureAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

RadiationPressureAccelerationPartial::RadiationPressureAccelerationPartial(
        const std::string acceleratedBody,
        const std::string acceleratingBody,
        const std::shared_ptr< electromagnetism::PaneledSourceRadiationPressureAcceleration > radiationPressureAcceleration,
        const std::shared_ptr< estimatable_parameters::CustomSingleAccelerationPartialCalculatorSet > customAccelerationPartialSet ):
    AccelerationPartial( acceleratedBody, acceleratingBody, radiationPressureAcceleration, basic_astrodynamics::radiation_pressure ),
    radiationPressureAcceleration_( radiationPressureAcceleration ), customAccelerationPartialSet_( customAccelerationPartialSet )
{
    radiationPressureAcceleration->enableSaveForcingQuantities( );

    currentPartialWrtUndergoingState_.setZero( );
    currentPartialWrtExertingState_.setZero( );

    estimatable_parameters::EstimatebleParameterIdentifier singleArcUndergoingBodyIdentifier =
            std::make_pair( estimatable_parameters::initial_body_state, std::make_pair( acceleratedBody, "" ) );
    if( customAccelerationPartialSet->customInitialStatePartials_.count( singleArcUndergoingBodyIdentifier ) > 0 )
    {
        bodyUndergoingPositionPartial_ = customAccelerationPartialSet->customInitialStatePartials_.at( singleArcUndergoingBodyIdentifier );
    }

    estimatable_parameters::EstimatebleParameterIdentifier multiArcUndergoingBodyIdentifier =
            std::make_pair( estimatable_parameters::arc_wise_initial_body_state, std::make_pair( acceleratedBody, "" ) );
    if( customAccelerationPartialSet->customInitialStatePartials_.count( multiArcUndergoingBodyIdentifier ) > 0 )
    {
        bodyUndergoingPositionPartial_ = customAccelerationPartialSet->customInitialStatePartials_.at( multiArcUndergoingBodyIdentifier );
    }

    estimatable_parameters::EstimatebleParameterIdentifier singleArcExertingBodyIdentifier =
            std::make_pair( estimatable_parameters::initial_body_state, std::make_pair( acceleratingBody, "" ) );
    if( customAccelerationPartialSet->customInitialStatePartials_.count( singleArcExertingBodyIdentifier ) > 0 &&
        ( singleArcExertingBodyIdentifier != singleArcUndergoingBodyIdentifier ) )
    {
        bodyExertingPositionPartial_ = customAccelerationPartialSet->customInitialStatePartials_.at( singleArcExertingBodyIdentifier );
    }

    estimatable_parameters::EstimatebleParameterIdentifier multiArcExertingBodyIdentifier =
            std::make_pair( estimatable_parameters::arc_wise_initial_body_state, std::make_pair( acceleratingBody, "" ) );
    if( customAccelerationPartialSet->customInitialStatePartials_.count( multiArcExertingBodyIdentifier ) > 0 &&
        ( multiArcExertingBodyIdentifier != multiArcUndergoingBodyIdentifier ) )
    {
        bodyExertingPositionPartial_ = customAccelerationPartialSet->customInitialStatePartials_.at( multiArcExertingBodyIdentifier );
    }
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
RadiationPressureAccelerationPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int parameterSize = 0;
    if( customAccelerationPartialSet_->customDoubleParameterPartials_.count( parameter->getParameterName( ) ) != 0 )
    {
        partialFunction = std::bind( &RadiationPressureAccelerationPartial::createCustomParameterPartialFunction,
                                     this,
                                     std::placeholders::_1,
                                     customAccelerationPartialSet_->customDoubleParameterPartials_.at( parameter->getParameterName( ) ) );
        parameterSize = 1;
    }
    else if( parameter->getParameterName( ).first == estimatable_parameters::radiation_pressure_coefficient &&
             parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        if( std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >(
                    radiationPressureAcceleration_->getTargetModel( ) ) != nullptr )
        {
            partialFunction = std::bind( &RadiationPressureAccelerationPartial::wrtRadiationPressureCoefficient,
                                         this,
                                         std::placeholders::_1,
                                         std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >(
                                                 radiationPressureAcceleration_->getTargetModel( ) ) );
            parameterSize = 1;
        }
    }
    else if( parameter->getParameterName( ).first == estimatable_parameters::specular_reflectivity &&
             parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        if( std::dynamic_pointer_cast< electromagnetism::PaneledRadiationPressureTargetModel >(
                    radiationPressureAcceleration_->getTargetModel( ) ) != nullptr )
        {
            if( parameter->getParameterName( ).second.second == "" )
            {
                throw std::runtime_error( "Error when creating specular reflectivity partial, panel group name not specified" );
            }
            else
            {
                partialFunction = std::bind( &RadiationPressureAccelerationPartial::wrtSpecularReflectivity,
                                             this,
                                             std::placeholders::_1,
                                             std::dynamic_pointer_cast< electromagnetism::PaneledRadiationPressureTargetModel >(
                                                     radiationPressureAcceleration_->getTargetModel( ) ),
                                             parameter->getParameterName( ).second.second );
                parameterSize = 1;
            };
        }
    }
    else if( parameter->getParameterName( ).first == estimatable_parameters::diffuse_reflectivity &&
             parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        if( std::dynamic_pointer_cast< electromagnetism::PaneledRadiationPressureTargetModel >(
                    radiationPressureAcceleration_->getTargetModel( ) ) != nullptr )
        {
            if( parameter->getParameterName( ).second.second == "" )
            {
                throw std::runtime_error( "Error when creating diffuse reflectivity partial, panel group name not specified" );
            }
            else
            {
                partialFunction = std::bind( &RadiationPressureAccelerationPartial::wrtDiffuseReflectivity,
                                             this,
                                             std::placeholders::_1,
                                             std::dynamic_pointer_cast< electromagnetism::PaneledRadiationPressureTargetModel >(
                                                     radiationPressureAcceleration_->getTargetModel( ) ),
                                             parameter->getParameterName( ).second.second );
                parameterSize = 1;
            };
        }
    }
    // Check if parameter dependency exists.
    else if( parameter->getParameterName( ).second.first == acceleratedBody_ &&
             parameter->getParameterName( ).second.second == acceleratingBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case estimatable_parameters::source_direction_radiation_pressure_scaling_factor:

                partialFunction = std::bind( &computeRadiationPressureAccelerationWrtSourceDirectionScaling,
                                             radiationPressureAcceleration_,
                                             std::placeholders::_1 );
                parameterSize = 1;

                break;
            case estimatable_parameters::source_perpendicular_direction_radiation_pressure_scaling_factor:

                partialFunction = std::bind( &computeRadiationPressureAccelerationWrtSourcePerpendicularDirectionScaling,
                                             radiationPressureAcceleration_,
                                             std::placeholders::_1 );
                parameterSize = 1;

                break;
            default:
                break;
        }
    }

    return std::make_pair( partialFunction, parameterSize );
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
RadiationPressureAccelerationPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int parameterSize = 0;
    if( customAccelerationPartialSet_->customVectorParameterPartials_.count( parameter->getParameterName( ) ) != 0 )
    {
        partialFunction = std::bind( &RadiationPressureAccelerationPartial::createCustomParameterPartialFunction,
                                     this,
                                     std::placeholders::_1,
                                     customAccelerationPartialSet_->customVectorParameterPartials_.at( parameter->getParameterName( ) ) );
        parameterSize = parameter->getParameterSize( );
    }
    else if( parameter->getParameterName( ).first == estimatable_parameters::arc_wise_radiation_pressure_coefficient &&
             parameter->getParameterName( ).second.first == acceleratedBody_ )
    {
        if( std::dynamic_pointer_cast< estimatable_parameters::ArcWiseRadiationPressureCoefficient >( parameter ) != nullptr )
        {
            if( std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >(
                        radiationPressureAcceleration_->getTargetModel( ) ) != nullptr )
            {
                partialFunction =
                        std::bind( &RadiationPressureAccelerationPartial::wrtArcWiseRadiationPressureCoefficient,
                                   this,
                                   std::placeholders::_1,
                                   std::dynamic_pointer_cast< estimatable_parameters::ArcWiseRadiationPressureCoefficient >( parameter ),
                                   std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >(
                                           radiationPressureAcceleration_->getTargetModel( ) ) );
                parameterSize = parameter->getParameterSize( );
            }
        }
        else
        {
            throw std::runtime_error( "Error when making radiation pressure partial, arcwise radiation pressure parameter not consistent" );
        }
    }
    return std::make_pair( partialFunction, parameterSize );
}

void RadiationPressureAccelerationPartial::update( const double currentTime )
{
    radiationPressureAcceleration_->updateMembers( currentTime );

    if( !( currentTime_ == currentTime ) )
    {
        if( bodyUndergoingPositionPartial_ != nullptr )
        {
            currentPartialWrtUndergoingState_ = bodyUndergoingPositionPartial_->computePartial(
                    currentTime, radiationPressureAcceleration_->getAcceleration( ), radiationPressureAcceleration_ );
        }

        if( bodyExertingPositionPartial_ != nullptr )
        {
            currentPartialWrtExertingState_ = bodyExertingPositionPartial_->computePartial(
                    currentTime, radiationPressureAcceleration_->getAcceleration( ), radiationPressureAcceleration_ );
        }

        currentTime_ = currentTime;
    }
}

void RadiationPressureAccelerationPartial::wrtRadiationPressureCoefficient(
        Eigen::MatrixXd& partial,
        std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > targetModel )
{
    if( targetModel->getCoefficient( ) == 0.0 )
    {
        throw std::runtime_error(
                "Error in full radiation pressure partial w.r.t. Cr, partial is only implemented for non-zero coefficient" );
    }

    partial = radiationPressureAcceleration_->getAcceleration( ) / targetModel->getCoefficient( );
}

void RadiationPressureAccelerationPartial::wrtArcWiseRadiationPressureCoefficient(
        Eigen::MatrixXd& partial,
        const std::shared_ptr< estimatable_parameters::ArcWiseRadiationPressureCoefficient > parameter,
        std::shared_ptr< electromagnetism::CannonballRadiationPressureTargetModel > targetModel )
{
    // Get partial w.r.t. radiation pressure coefficient
    Eigen::MatrixXd partialWrtSingleParameter = Eigen::Vector3d::Zero( );
    this->wrtRadiationPressureCoefficient( partialWrtSingleParameter, targetModel );

    // Retrieve current arc
    std::shared_ptr< interpolators::LookUpScheme< double > > currentArcIndexLookUp = parameter->getArcTimeLookupScheme( );
    partial.setZero( );
    if( currentArcIndexLookUp->getMinimumValue( ) <= currentTime_ )
    {
        int currentArc = currentArcIndexLookUp->findNearestLowerNeighbour( currentTime_ );

        if( currentArc >= partial.cols( ) )
        {
            throw std::runtime_error( "Error when getting arc-wise radiation pressure coefficient partials, data not consistent" );
        }

        // Set partial
        partial.block( 0, currentArc, 3, 1 ) = partialWrtSingleParameter;
    }
}

void RadiationPressureAccelerationPartial::wrtSpecularReflectivity(
        Eigen::MatrixXd& partial,
        std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > targetModel,
        const std::string& panelTypeId )
{
    partial.setZero( );

    std::function< double( ) > targetMassFunction = radiationPressureAcceleration_->getTargetMassFunction( );
    double spacecraftMass = targetMassFunction( );

    std::vector< double >& savedPanelOccultedIrradiances = radiationPressureAcceleration_->getSavedPanelOccultedIrradiances( );
    std::vector< Eigen::Vector3d >& savedPanelRelativePositions = radiationPressureAcceleration_->getSavedPanelRelativePositions( );

    Eigen::Quaterniond targetRotationFromLocalToGlobalFrame =
            radiationPressureAcceleration_->getTargetRotationFromLocalToGlobalFrameFunction( )( );

    for( unsigned int i = 0; i < savedPanelOccultedIrradiances.size( ); i++ )
    {
        if( savedPanelOccultedIrradiances.at( i ) >= 0 )
        {
            Eigen::Vector3d forcePartialWrtDiffuseReflectivity = targetRotationFromLocalToGlobalFrame *
                    targetModel->evaluateRadiationPressureForcePartialWrtSpecularReflectivity(
                            savedPanelOccultedIrradiances.at( i ), savedPanelRelativePositions.at( i ), panelTypeId );
            partial += forcePartialWrtDiffuseReflectivity / spacecraftMass;
        }
    }
}
void RadiationPressureAccelerationPartial::wrtDiffuseReflectivity(
        Eigen::MatrixXd& partial,
        std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > targetModel,
        const std::string& panelTypeId )
{
    partial.setZero( );

    std::function< double( ) > targetMassFunction = radiationPressureAcceleration_->getTargetMassFunction( );
    double spacecraftMass = targetMassFunction( );

    std::vector< double >& savedPanelOccultedIrradiances = radiationPressureAcceleration_->getSavedPanelOccultedIrradiances( );
    std::vector< Eigen::Vector3d >& savedPanelRelativePositions = radiationPressureAcceleration_->getSavedPanelRelativePositions( );

    Eigen::Quaterniond targetRotationFromLocalToGlobalFrame =
            radiationPressureAcceleration_->getTargetRotationFromLocalToGlobalFrameFunction( )( );

    for( unsigned int i = 0; i < savedPanelOccultedIrradiances.size( ); i++ )
    {
        if( savedPanelOccultedIrradiances.at( i ) >= 0 )
        {
            Eigen::Vector3d forcePartialWrtDiffuseReflectivity = targetRotationFromLocalToGlobalFrame *
                    targetModel->evaluateRadiationPressureForcePartialWrtDiffuseReflectivity(
                            savedPanelOccultedIrradiances.at( i ), savedPanelRelativePositions.at( i ), panelTypeId );
            partial += forcePartialWrtDiffuseReflectivity / spacecraftMass;
        }
    }
}

}  // namespace acceleration_partials

}  // namespace tudat
