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

#ifndef TUDAT_RAREFIEDFLOW_AERODYNAMIC_COEFFICIENT_INTERFACE_H
#define TUDAT_RAREFIEDFLOW_AERODYNAMIC_COEFFICIENT_INTERFACE_H

#include <map>
#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/basics/basicTypedefs.h"
// #include "tudat/astro/system_models/vehicleExteriorPanels.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/astro/aerodynamics/rarefiedFlowInteractionModel.h"

namespace tudat
{
namespace aerodynamics
{

// template< unsigned int NumberOfIndependentVariables >
// class RarefiedFlowAerodynamicCoefficientInterface: public AerodynamicCoefficientInterface
class RarefiedFlowAerodynamicCoefficientInterface: public AerodynamicCoefficientInterface
{
public:
    /*!
     * Constructor of rarefied flow aerodynamic coefficient interface.
     * \param vehicleExteriorPanels Vehicle panels
     * \param vehiclePartOrientation Vehicle part orientation
     * \param referenceLength Reference length
     * \param referenceArea Reference area
     * \param momentReferencePoint Moment reference point
     * \param independentVariableNames Independent variable names
     * \param forceCoefficientsFrame Force coefficients frame
     * \param momentCoefficientsFrame Moment coefficients frame
     * \param accountForShadedPanels Account for shaded panels
     * \param dataPointsOfInclinationsForShading Data points of inclinations for shading
     */
    RarefiedFlowAerodynamicCoefficientInterface(
        tudat::system_models::VehicleSystems vehicle,
        const double referenceLength,
        const double referenceArea,
        const Eigen::Vector3d& momentReferencePoint,
        const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames = {
            angle_of_attack_dependent,
            angle_of_sideslip_dependent,
            temperature_dependent,
            velocity_dependent,
            he_number_density_dependent,
            o_number_density_dependent,
            n2_number_density_dependent,
            o2_number_density_dependent,
            ar_number_density_dependent,
            h_number_density_dependent,
            n_number_density_dependent,
            anomalous_o_number_density_dependent},
        const AerodynamicCoefficientFrames forceCoefficientsFrame = negative_aerodynamic_frame_coefficients,
        const AerodynamicCoefficientFrames momentCoefficientsFrame = body_fixed_frame_coefficients,
        const bool accountForShadedPanels = false
        // const std::map< int, std::vector< double > > dataPointsOfInclinationsForShading = std::map< int, std::vector< double > >( ) 
        ): 
        AerodynamicCoefficientInterface(
            referenceLength, referenceLength, momentReferencePoint, independentVariableNames, forceCoefficientsFrame, momentCoefficientsFrame), 
        vehicle_( vehicle ),
        vehicleExteriorPanels_( vehicle.getVehicleExteriorPanels() ), 
        // vehiclePartOrientations_( vehicle.getVehiclePartOrientations() ), 
        referenceLength_( referenceLength ), referenceArea_( referenceArea ),
        momentReferencePoint_( momentReferencePoint ), 
        independentVariableNames_( independentVariableNames ), 
        forceCoefficientsFrame_( forceCoefficientsFrame ),
        momentCoefficientsFrame_( momentCoefficientsFrame ), accountForShadedPanels_( accountForShadedPanels )
        // dataPointsOfInclinationsForShading_( dataPointsOfInclinationsForShading)
        {
            
            // initializing total aerodynamic coefficient vector with nans
            totalAerodynamicCoefficients_ = Eigen::Vector6d::Constant( TUDAT_NAN );

        }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~RarefiedFlowAerodynamicCoefficientInterface( ) = default;



    void updateCurrentCoefficients(
        const std::vector< double >& independentVariables,
        const double currentTime);
    

    Eigen::Vector6d getCurrentAerodynamicCoefficients( );

private:


    void determineIncinations(
        double secondsSinceEpoch,
        double angleOfAttack,
        double angleOfSideslip
    );
    


    void determinePanelForceCoefficientVectors(
        double freestreamVelocity, double atmosphericTemperature, std::vector< double > numberDensities
    );



    void determinePanelMomentCoefficientVectors(double secondsSinceEpoch);

    
    // Declaration of member variables

    //! Vehicle panels
    std::map< std::string, std::vector< std::shared_ptr< tudat::system_models::VehicleExteriorPanel > > > vehicleExteriorPanels_;

    //! Vehicle part orientation
    // std::map< std::string, std::shared_ptr< tudat::ephemerides::RotationalEphemeris > > vehiclePartOrientations_;

    //! Vehicle panel cosines of lift and drag angles
    std::map< std::string, std::vector< std::pair< double, double > > > vehiclePanelCosinesOfLiftAndDragAngles_;

    //! Vehicle panel lift unit vecotrs
    std::map< std::string, std::vector< Eigen::Vector3d > > vehiclePanelLiftUnitVectors_;

    //! drag unit vector
    Eigen::Vector3d dragUnitVector_;

    //! Vehicle panel force coefficient vectors
    std::map< std::string, std::vector< Eigen::Vector3d > > vehiclePanelForceCoefficientVectors_;

    //! Vehicle panel moment coefficient vectors
    std::map< std::string, std::vector< Eigen::Vector3d > > vehiclePanelMomentCoefficientVectors_;

    //! Total aerodynamic coefficient vector
    Eigen::Vector6d totalAerodynamicCoefficients_;

    //! Reference length
    double referenceLength_;
    //! Reference area
    double referenceArea_;
    //! Moment reference point
    Eigen::Vector3d momentReferencePoint_;

    //! Independent variable names
    std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames_;

    //! Force coefficients frame
    AerodynamicCoefficientFrames forceCoefficientsFrame_;
    //! Moment coefficients frame
    AerodynamicCoefficientFrames momentCoefficientsFrame_;

    //! Account for shaded panels
    bool accountForShadedPanels_;

    //! Data points of inclinations for shading
    std::map< int, std::vector< double > > dataPointsOfInclinationsForShading_;

    //! Vehicle
    tudat::system_models::VehicleSystems vehicle_;

};


} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_RAREFIEDFLOW_AERODYNAMIC_COEFFICIENT_INTERFACE_H
