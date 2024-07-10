/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Doornbos, E. N. Thermospheric Density and Wind Determination from Satellite Dynamics, 2011. 
 *      (Page 66+)
 */

#include <map>
#include <string>
#include <vector>

#include <array>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"
#include "tudat/astro/aerodynamics/rarefiedFlowInteractionModel.h"
#include "tudat/astro/aerodynamics/rarefiedFlowAerodynamicCoefficientInterface.h"

namespace tudat
{
namespace aerodynamics
{


// RarefiedFlowAerodynamicCoefficientInterface::RarefiedFlowAerodynamicCoefficientInterface(
//     const std::map< std::string, std::vector< std::shared_ptr< tudat::system_models::VehicleExteriorPanel > > > vehicleExteriorPanels,
//     const std::map< std::string, std::shared_ptr< tudat::ephemerides::RotationalEphemeris > > vehiclePartOrientation,
//     const double referenceLength,
//     const double referenceArea,
//     const Eigen::Vector3d& momentReferencePoint,
//     const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames,
//     const AerodynamicCoefficientFrames forceCoefficientsFrame = negative_aerodynamic_frame_coefficients,
//     const AerodynamicCoefficientFrames momentCoefficientsFrame = body_fixed_frame_coefficients,
//     const bool accountForShadedPanels = false,
//     const std::map< int, std::vector< double > > dataPointsOfInclinationsForShading = std::map< int, std::vector< double > >( )
//     ){


//     }




//! Compute the aerodynamic coefficients of the body itself (without control surfaces) at current flight condition.
/*!
*  Computes the current force and moment coefficients of the body itself (without control surfaces) and is to be
*  implemented in derived classes. Input is a set of independent variables
*  (doubles) which represent the variables from which the coefficients are calculated
*  \param independentVariables Independent variables of force and moment coefficient
*  determination implemented by derived class
*  \param currentTime Time to which coefficients are to be updated
*/
void RarefiedFlowAerodynamicCoefficientInterface::updateCurrentCoefficients(
    const std::vector< double >& independentVariables,
    const double currentTime = TUDAT_NAN )
{
    // Update vehicle part orientations
    // vehicle_.updatePartOrientations( currentTime );

    // Determine inclinations
    determineIncinations( currentTime, independentVariables.at(0), independentVariables.at(1) );

    std::vector< double > numberDensities;

    for ( int i = 4; i < independentVariables.size(); i++ )
    {
        numberDensities.push_back( independentVariables.at(i) );
    }

    determinePanelForceCoefficientVectors( independentVariables.at(2), independentVariables.at(3), numberDensities );

    determinePanelMomentCoefficientVectors( currentTime );

    // sum all panel force and moment coefficient vectors to get the total force and moment coefficient vectors
    
    // initialize total aerodynamic coefficient vector with zeros
    totalAerodynamicCoefficients_ = Eigen::Vector6d::Zero();

    for (auto& vehiclePartEntry : vehicleExteriorPanels_) {
        const std::string& vehiclePartName = vehiclePartEntry.first;

        for ( int i = 0; i < vehicleExteriorPanels_.at( vehiclePartName ).size(); i++ )
        {
            totalAerodynamicCoefficients_.head(3) += vehiclePanelForceCoefficientVectors_[ vehiclePartName ].at( i );
            totalAerodynamicCoefficients_.tail(3) += vehiclePanelMomentCoefficientVectors_[ vehiclePartName ].at( i );
        }
    }

}

//! Get the current aerodynamic coefficients of the body itself (without control surfaces).
/*!
*  Returns the current force and moment coefficients of the body itself (without control surfaces).
*  \return Current force and moment coefficients of the body itself (without control surfaces).
*/

Eigen::Vector6d RarefiedFlowAerodynamicCoefficientInterface::getCurrentAerodynamicCoefficients( )
{
    return totalAerodynamicCoefficients_;
}

void RarefiedFlowAerodynamicCoefficientInterface::determineIncinations(
    double secondsSinceEpoch,
    double angleOfAttack,
    double angleOfSideslip
){
    // Declare free-stream velocity unit vector

    Eigen::Vector3d freestreamVelocityDirection = Eigen::Vector3d::Zero();

    freestreamVelocityDirection(0) = std::cos( angleOfAttack ) * std::cos( angleOfSideslip );
    freestreamVelocityDirection(1) = std::sin( angleOfSideslip );
    freestreamVelocityDirection(2) = std::sin( angleOfAttack ) * std::cos( angleOfSideslip );

    // Drag unit vector is parallel to the freestream velocity direction

    dragUnitVector_ = freestreamVelocityDirection.normalized();

    // Loop over all vehicle part names in vehicleExteriorPanels_

    // Initialize vehicle part rotation matrix to base frame
    // Eigen::Matrix3d rotationMatrixToBaseFrame = Eigen::Matrix3d::Zero();
    
    for (auto& vehiclePartEntry : vehicleExteriorPanels_) {
        const std::string& vehiclePartName = vehiclePartEntry.first;
        const std::vector<std::shared_ptr<tudat::system_models::VehicleExteriorPanel>>& exteriorPanels = vehiclePartEntry.second;

        // get rotation matrix for this part
        // rotationMatrixToBaseFrame = vehicle_.getPartRotationToBaseFrame( vehiclePartName ).toRotationMatrix();

        // Loop over all vehicle panels in vehicleExteriorPanels_ for this vehicle part
        for ( std::shared_ptr< tudat::system_models::VehicleExteriorPanel > vehiclePanel: exteriorPanels )
        {
            
            // Determine panel normal vector and rotate to base frame

            // Eigen::Vector3d panelNormalVector =  rotationMatrixToBaseFrame * vehiclePanel->getFrameFixedSurfaceNormal()();

            Eigen::Vector3d panelNormalVector = vehiclePanel->getFrameFixedSurfaceNormal()();

            // Determine cosine of angle between panel normal vector and freestream velocity direction (drag)
            // gamma_i in Doornbos (2011)
            double panelCosineDragAngle = -dragUnitVector_.dot( panelNormalVector );

            // Initialize variables to be defined in the following if-else statement
            double panelCosineLiftAngle = TUDAT_NAN;
            Eigen::Vector3d liftUnitVector = Eigen::Vector3d::Constant( TUDAT_NAN );


            // If panel is "perfectly" perpendicular to velocity vector, cosine of lift angle is set manually to 0.0
            if ( std::abs( panelCosineDragAngle - 1.0 ) < std::numeric_limits< double >::epsilon() )
            {            
                panelCosineLiftAngle = 0.0;

                // lift unit vector is irrelevant in this case but has to be defined for the code to work
                liftUnitVector = Eigen::Vector3d(1.0, 0.0, 0.0);
                vehiclePanelLiftUnitVectors_[ vehiclePartName ].push_back( liftUnitVector );

                panelCosineDragAngle = 1.0; 

            } else {

                // Determine lift unit vector for this panel
                liftUnitVector = -( ( dragUnitVector_.cross( panelNormalVector ) ).cross( dragUnitVector_ ) ).normalized();

                // Store lift unit vector
                vehiclePanelLiftUnitVectors_[ vehiclePartName ].push_back( liftUnitVector );

                // Determine cosine of angle between panel normal vector and lift unit vector
                // l_i in Doornbos (2011)
                panelCosineLiftAngle = -liftUnitVector.dot( panelNormalVector );
            
            }

            // Store panel cosines of lift and drag angles

            vehiclePanelCosinesOfLiftAndDragAngles_[ vehiclePartName ].push_back( std::make_pair( panelCosineDragAngle, panelCosineLiftAngle ) );

        }
    }
}

void RarefiedFlowAerodynamicCoefficientInterface::determinePanelForceCoefficientVectors(
    double freestreamVelocity, double atmosphericTemperature, std::vector< double > numberDensities
){

    double totalNumberDensity = 0.0;
    for (int j_species = 0; j_species < numberDensities.size(); j_species++)
    {
        totalNumberDensity += numberDensities[j_species];
    }

    // Loop over all vehicle part names in vehicleExteriorPanels_

    for (auto& vehiclePartEntry : vehicleExteriorPanels_) {
        const std::string& vehiclePartName = vehiclePartEntry.first;
        const std::vector<std::shared_ptr<tudat::system_models::VehicleExteriorPanel>>& exteriorPanels = vehiclePartEntry.second;

        // Loop over all vehicle panels in vehicleExteriorPanels_ for this vehicle part
        for ( int i = 0; i < vehicleExteriorPanels_.at( vehiclePartName ).size(); i++ )
        {
            std::shared_ptr< tudat::system_models::VehicleExteriorPanel > vehiclePanel = vehicleExteriorPanels_.at( vehiclePartName ).at( i );
            
            vehiclePanelForceCoefficientVectors_[ vehiclePartName ].push_back( 
                vehiclePanel->getRarefiedFlowInteractionModel()->computePanelForceCoefficientVector( 
                    vehiclePanelCosinesOfLiftAndDragAngles_[ vehiclePartName ].at( i ).first,
                    vehiclePanelCosinesOfLiftAndDragAngles_[ vehiclePartName ].at( i ).second,
                    vehiclePanel->getPanelArea(),
                    vehiclePanel->getPanelTemperature()(),
                    vehiclePanelLiftUnitVectors_[ vehiclePartName ].at( i ),
                    dragUnitVector_,
                    freestreamVelocity,
                    atmosphericTemperature,
                    numberDensities,
                    totalNumberDensity,
                    referenceArea_
                )
            );
        }
    }
}   

void RarefiedFlowAerodynamicCoefficientInterface::determinePanelMomentCoefficientVectors(double secondsSinceEpoch){
    
    Eigen::Vector3d partPositionVector = Eigen::Vector3d::Zero();

    // Loop over all vehicle part names in vehicleExteriorPanels_
    for (auto& vehiclePartEntry : vehicleExteriorPanels_) {
        const std::string& vehiclePartName = vehiclePartEntry.first;
        const std::vector<std::shared_ptr<tudat::system_models::VehicleExteriorPanel>>& exteriorPanels = vehiclePartEntry.second;

        // Loop over all vehicle panels in vehicleExteriorPanels_ for this vehicle part
        for ( int i = 0; i < vehicleExteriorPanels_.at( vehiclePartName ).size(); i++ )
        {
            std::shared_ptr< tudat::system_models::VehicleExteriorPanel > vehiclePanel = vehicleExteriorPanels_.at( vehiclePartName ).at( i );
            
            // Calculate panel position vector in base frame
            Eigen::Vector3d panelPositionVector = 
                // vehiclePartOrientation_.at( vehiclePartName )->getRotationMatrixToBaseFrame(secondsSinceEpoch) * vehiclePanel->getFrameFixedPositionVector()() +
                // partPositionVector;
                vehiclePanel->getFrameFixedPositionVector()();
            
            vehiclePanelMomentCoefficientVectors_[ vehiclePartName ].push_back( 
                vehiclePanel->getRarefiedFlowInteractionModel()->computePanelMomentCoefficientVector( 
                    vehiclePanelForceCoefficientVectors_[ vehiclePartName ].at( i ),
                    panelPositionVector,
                    referenceLength_
                )
            );
        }
    }
}


} // namespace aerodynamics
} // namespace tudat

