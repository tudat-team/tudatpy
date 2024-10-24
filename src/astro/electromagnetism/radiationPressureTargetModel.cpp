/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"


namespace tudat
{
namespace electromagnetism
{

void RadiationPressureTargetModel::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        updateMembers_(currentTime);
    }
}

void CannonballRadiationPressureTargetModel::updateRadiationPressureForcing(
    double sourceIrradiance,
    const Eigen::Vector3d &sourceToTargetDirection,
    const bool resetForces,
    const std::string sourceName )
{
    if( resetForces )
    {
        resetComputations( sourceName );
    }

    // From Montenbruck (2000), Sec. 3.4
    double radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
    this->currentRadiationPressureForce_[ sourceName ] += currentCoefficient_ * area_ * radiationPressure * sourceToTargetDirection;
    if( computeTorques_ )
    {
        this->currentRadiationPressureTorque_[ sourceName ] += -centerOfMassFunction_( ).cross( this->currentRadiationPressureForce_.at( sourceName ) );
    }
}

void PaneledRadiationPressureTargetModel::updateRadiationPressureForcing(
        double sourceIrradiance,
        const Eigen::Vector3d& sourceToTargetDirectionLocalFrame,
        const bool resetForces,
        const std::string sourceName )
{
    double radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;

    if( resetForces )
    {
        resetComputations( sourceName );
    }
    auto segmentFixedPanelsIterator = segmentFixedPanels_.begin( );

    int counter = 0;
    Eigen::Quaterniond currentOrientation;
    Eigen::Vector3d currentCenterOfMass = Eigen::Vector3d::Constant( TUDAT_NAN );
    if( computeTorques_ )
    {
        currentCenterOfMass = centerOfMassFunction_( );
    }
    Eigen::Vector3d currentPanelForce = Eigen::Vector3d::Zero( );
    Eigen::Vector3d currentPanelTorque = Eigen::Vector3d::Zero( );

    for( unsigned int i = 0; i < segmentFixedPanels_.size( ) + 1; i++ )
    {
        currentOrientation = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
        if( i > 0 )
        {
            currentOrientation = segmentFixedToBodyFixedRotations_.at( segmentFixedPanelsIterator->first )( );
            if( computeTorques_ )
            {
                throw std::runtime_error( "Torques not yet supported for moving vehicle parts." );
            }
        }

        const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& currentPanels_ =
            ( i == 0 ) ? bodyFixedPanels_ : segmentFixedPanels_.at( segmentFixedPanelsIterator->first );
        for( unsigned int j = 0; j < currentPanels_.size( ); j++ )
        {
            surfaceNormals_[ counter ] = currentOrientation * currentPanels_.at( j )->getFrameFixedSurfaceNormal( )( );
            surfacePanelCosines_[ counter ] = (-sourceToTargetDirectionLocalFrame).dot(surfaceNormals_[ counter ]);
            if( computeTorques_ )
            {
                panelCentroidMomentArms_[ counter ] = currentOrientation * ( currentPanels_.at( j )->getFrameFixedPositionVector( )( ) - currentCenterOfMass );
            }
            if (surfacePanelCosines_[ counter ] > 0)
            {
                currentPanelForce = radiationPressure * currentPanels_.at( j )->getPanelArea() * surfacePanelCosines_[ counter ] *
                    currentPanels_.at( j )->getReflectionLaw()->evaluateReactionVector(surfaceNormals_[ counter ], sourceToTargetDirectionLocalFrame );
                this->currentRadiationPressureForce_[ sourceName ] += currentPanelForce;
                if( computeTorques_ )
                {
                    currentPanelTorque = panelCentroidMomentArms_[ counter ].cross( currentPanelForce );
                    this->currentRadiationPressureTorque_[ sourceName ]  += currentPanelTorque;
                }

            }
            else
            {
                currentPanelForce.setZero( );
                if( computeTorques_ )
                {
                    currentPanelTorque.setZero( );
                }
            }
            panelForces_[ counter ] += currentPanelForce;

            if( computeTorques_ )
            {
                panelTorques_[ counter ] += currentPanelTorque;
            }

            counter++;
        }
        if( i > 0 )
        {
            segmentFixedPanelsIterator++;
        }
    }
}

void PaneledRadiationPressureTargetModel::saveLocalComputations( const std::string sourceName, const bool saveCosines )
{
    if( saveCosines )
    {
        surfacePanelCosinesPerSource_[ sourceName ] = surfacePanelCosines_;
    }
    panelForcesPerSource_[ sourceName ] = panelForces_;
    panelTorquesPerSource_[ sourceName ] = panelTorques_;
}


Eigen::Vector3d PaneledRadiationPressureTargetModel::evaluateRadiationPressureForcePartialWrtDiffuseReflectivity(
        double sourceIrradiance,
        const Eigen::Vector3d& sourceToTargetDirectionLocalFrame)
{
    Eigen::Vector3d forcePartialWrtDiffuseReflectivity = Eigen::Vector3d::Zero();

    double radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
    auto segmentFixedPanelsIterator = segmentFixedPanels_.begin( );
    int counter = 0;
    Eigen::Quaterniond currentOrientation;

    for( unsigned int i = 0; i < segmentFixedPanels_.size( ) + 1; i++ )
    {
        currentOrientation = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
        if( i > 0 )
        {
            currentOrientation = segmentFixedToBodyFixedRotations_.at( segmentFixedPanelsIterator->first )( );

        }

        const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& currentPanels_ =
            ( i == 0 ) ? bodyFixedPanels_ : segmentFixedPanels_.at( segmentFixedPanelsIterator->first );
        for( unsigned int j = 0; j < currentPanels_.size( ); j++ )
        {
            surfaceNormals_[ counter ] = currentOrientation * currentPanels_.at( j )->getFrameFixedSurfaceNormal( )( );
            surfacePanelCosines_[ counter ] = (-sourceToTargetDirectionLocalFrame).dot(surfaceNormals_[ counter ]);

            if (surfacePanelCosines_[ counter ] > 0)
            {
                Eigen::Vector3d panelForce = radiationPressure * currentPanels_.at( j )->getPanelArea() * surfacePanelCosines_[ counter ] *
                    currentPanels_.at( j )->getReflectionLaw()->evaluateReactionVectorPartialWrtDiffuseReflectivity(surfaceNormals_[ counter ], sourceToTargetDirectionLocalFrame );
                forcePartialWrtDiffuseReflectivity += panelForce;
            }

            counter++;
        }
        if( i > 0 )
        {
            segmentFixedPanelsIterator++;
        }
    }
    return forcePartialWrtDiffuseReflectivity;
}

Eigen::Vector3d PaneledRadiationPressureTargetModel::evaluateRadiationPressureForcePartialWrtSpecularReflectivity(
        double sourceIrradiance,
        const Eigen::Vector3d& sourceToTargetDirectionLocalFrame)
{
    Eigen::Vector3d forcePartialWrtSpecularReflectivity = Eigen::Vector3d::Zero();

    double radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
    auto segmentFixedPanelsIterator = segmentFixedPanels_.begin( );
    int counter = 0;
    Eigen::Quaterniond currentOrientation;

    for( unsigned int i = 0; i < segmentFixedPanels_.size( ) + 1; i++ )
    {
        currentOrientation = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
        if( i > 0 )
        {
            currentOrientation = segmentFixedToBodyFixedRotations_.at( segmentFixedPanelsIterator->first )( );

        }

        const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& currentPanels_ =
            ( i == 0 ) ? bodyFixedPanels_ : segmentFixedPanels_.at( segmentFixedPanelsIterator->first );
        for( unsigned int j = 0; j < currentPanels_.size( ); j++ )
        {
            surfaceNormals_[ counter ] = currentOrientation * currentPanels_.at( j )->getFrameFixedSurfaceNormal( )( );
            surfacePanelCosines_[ counter ] = (-sourceToTargetDirectionLocalFrame).dot(surfaceNormals_[ counter ]);

            if (surfacePanelCosines_[ counter ] > 0)
            {
                Eigen::Vector3d panelForce = radiationPressure * currentPanels_.at( j )->getPanelArea() * surfacePanelCosines_[ counter ] *
                    currentPanels_.at( j )->getReflectionLaw()->evaluateReactionVectorPartialWrtSpecularReflectivity(surfaceNormals_[ counter ], sourceToTargetDirectionLocalFrame );
                forcePartialWrtSpecularReflectivity += panelForce;
            }

            counter++;
        }
        if( i > 0 )
        {
            segmentFixedPanelsIterator++;
        }
    }
    return forcePartialWrtSpecularReflectivity;
}

void PaneledRadiationPressureTargetModel::updateMembers_(double currentTime)
{

}

} // tudat
} // electromagnetism
