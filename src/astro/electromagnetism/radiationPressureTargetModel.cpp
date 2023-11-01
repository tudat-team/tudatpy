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

Eigen::Vector3d CannonballRadiationPressureTargetModel::evaluateRadiationPressureForce(
    const double sourceIrradiance,
    const Eigen::Vector3d& sourceToTargetDirection)
{
    // From Montenbruck (2000), Sec. 3.4
    radiationPressure_ = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
    return currentCoefficient_ * area_ * radiationPressure_ * sourceToTargetDirection;
}

Eigen::Vector3d PaneledRadiationPressureTargetModel::evaluateRadiationPressureForce(
        double sourceIrradiance,
        const Eigen::Vector3d& sourceToTargetDirectionLocalFrame)
{
    radiationPressure_ = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
    Eigen::Vector3d force = Eigen::Vector3d::Zero();
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
                panelForces_[ counter ] = radiationPressure_ * currentPanels_.at( j )->getPanelArea() * surfacePanelCosines_[ counter ] *
                    currentPanels_.at( j )->getReflectionLaw()->evaluateReactionVector(surfaceNormals_[ counter ], sourceToTargetDirectionLocalFrame );
                force += panelForces_[ counter ];
            }
            else
            {
                panelForces_[ counter ].setZero( );
            }
            counter++;
        }
        if( i > 0 )
        {
            segmentFixedPanelsIterator++;
        }
    }
    return force;
}

void PaneledRadiationPressureTargetModel::updateMembers_(double currentTime)
{

}

} // tudat
} // electromagnetism
