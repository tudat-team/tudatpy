/*   Copyright (c) 2010-2019, Delft University of Technology
 *   All rigths reserved
 *
 *   This file is part of the Tudat. Redistribution and use in source and
 *   binary forms, with or without modification, are permitted exclusively
 *   under the terms of the Modified BSD license. You should have received
 *   a copy of the license with this file. If not, please or visit:
 *   http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONOUTPUTSETTINGS_H
#define TUDAT_PROPAGATIONOUTPUTSETTINGS_H

#include <string>

#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/basic_astro/torqueModelTypes.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
#include "tudat/astro/orbit_determination/stateDerivativePartial.h"
#endif

namespace tudat
{

namespace propagators
{

// Enum listing the variables that can be exported or use in termination conditions during the propagation
enum VariableType
{
    independentVariable,
    cpuTimeVariable,
    stateVariable,
    dependentVariable,  // -> derivedVariable ?
    stateTransitionMatrix,
    sensitivityMatrix
};

// Functional base class for defining settings for variables
/*
 *  Functional base class for defining settings for variables.
 *  Any variable that requires additional information in addition to what can be provided here, should be
 *  defined by a dedicated derived class.
 */
//! @get_docstring(VariableSettings.__docstring__)
class VariableSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param variableType Type of variable.
     */
    VariableSettings( const VariableType variableType ) :
        variableType_( variableType ) { }

    // Destructor.
    virtual ~VariableSettings( ){ }

    // Type of dependent variable that is to be saved.
    VariableType variableType_;

};


// Enum listing the dependent variables that can be saved during the propagation.
//! @get_docstring(PropagationDependentVariables.__docstring__)
enum PropagationDependentVariables
{
    mach_number_dependent_variable = 0,
    altitude_dependent_variable = 1,
    airspeed_dependent_variable = 2,
    local_density_dependent_variable = 3,
    relative_speed_dependent_variable = 4,
    relative_position_dependent_variable = 5,
    relative_distance_dependent_variable = 6,
    relative_velocity_dependent_variable = 7,
    radiation_pressure_dependent_variable = 8,
    total_acceleration_norm_dependent_variable = 9,
    single_acceleration_norm_dependent_variable = 10,
    total_acceleration_dependent_variable = 11,
    single_acceleration_dependent_variable = 12,
    aerodynamic_force_coefficients_dependent_variable = 13,
    aerodynamic_moment_coefficients_dependent_variable = 14,
    inertial_to_body_fixed_rotation_matrix_variable = 15,
    intermediate_aerodynamic_rotation_matrix_variable = 16,
    relative_body_aerodynamic_orientation_angle_variable = 17,
    body_fixed_airspeed_based_velocity_variable = 18,
    total_aerodynamic_g_load_variable = 19,
    stagnation_point_heat_flux_dependent_variable = 20, // no interface function
    local_temperature_dependent_variable = 21,
    geodetic_latitude_dependent_variable = 22,
    control_surface_deflection_dependent_variable = 23,
    total_mass_rate_dependent_variables = 24,
    tnw_to_inertial_frame_rotation_dependent_variable = 25,
    periapsis_altitude_dependent_variable = 26,
    total_torque_norm_dependent_variable = 27,
    single_torque_norm_dependent_variable = 28,
    total_torque_dependent_variable = 29,
    single_torque_dependent_variable = 30,
    body_fixed_groundspeed_based_velocity_variable = 31,
    keplerian_state_dependent_variable = 32,
    modified_equinocial_state_dependent_variable = 33,
    spherical_harmonic_acceleration_terms_dependent_variable = 34,
    spherical_harmonic_acceleration_norm_terms_dependent_variable = 35,
    body_fixed_relative_cartesian_position = 36,
    body_fixed_relative_spherical_position = 37,
    total_gravity_field_variation_acceleration = 38,
    single_gravity_field_variation_acceleration = 39,
    single_gravity_field_variation_acceleration_terms = 40,
    acceleration_partial_wrt_body_translational_state = 41,
    local_dynamic_pressure_dependent_variable = 42,
    local_aerodynamic_heat_rate_dependent_variable = 43, // no interface function
    euler_angles_to_body_fixed_313 = 44,
    current_body_mass_dependent_variable = 45,
    radiation_pressure_coefficient_dependent_variable = 46,
    rsw_to_inertial_frame_rotation_dependent_variable = 47,
    custom_dependent_variable = 48
};

// Functional base class for defining settings for dependent variables that are to be saved during propagation
/*
 *  Functional base class for defining settings for dependent variables that are to be saved during propagation.
 *  Any dependent variable that requires additional information in addition to what can be provided here, should be
 *  defined by a dedicated derived class.
 */
//! @get_docstring(SingleDependentVariableSaveSettings.__docstring__)
class SingleDependentVariableSaveSettings : public VariableSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param dependentVariableType Type of dependent variable that is to be saved.
     *  \param associatedBody Body associated with dependent variable.
     *  \param secondaryBody Secondary body (not necessarilly required) w.r.t. which parameter is defined (e.g. relative
     *  position, velocity etc. is defined of associatedBody w.r.t. secondaryBody).
     *  \param componentIndex Index of the component to be saved. Only applicable to vectorial dependent variables.
     *  By default -1, i.e. all the components are saved.
     */
    SingleDependentVariableSaveSettings(
            const PropagationDependentVariables dependentVariableType,
            const std::string& associatedBody,
            const std::string& secondaryBody = "",
            const int componentIndex = -1 ):
        VariableSettings( dependentVariable ),
        dependentVariableType_( dependentVariableType ),
        associatedBody_( associatedBody ),
        secondaryBody_( secondaryBody ),
        componentIndex_( componentIndex ) { }

    // Type of dependent variable that is to be saved.
    PropagationDependentVariables dependentVariableType_;

    // Body associated with variable.
    std::string associatedBody_;

    // Secondary body (not necessarilly required) w.r.t. which parameter is defined (e.g. relative  position,
    // velocity etc. is defined of associatedBody w.r.t. secondaryBody).
    std::string secondaryBody_;

    // Index of the component to be saved.
    // Only applicable to vectorial dependent variables.
    // If negative, all the components of the vector are saved.
    int componentIndex_;

};

// Class to define settings for saving a single acceleration (norm or vector) during propagation
/*
 *  Class to define settings for saving a single acceleration (norm or vector) during propagation. NOTE: This acceleration
 *  is returned in the inertial frame!
 */
//! @get_docstring(SingleAccelerationDependentVariableSaveSettings.__docstring__)
class SingleAccelerationDependentVariableSaveSettings: public SingleDependentVariableSaveSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param accelerationModelType Type of acceleration that is to be saved.
     *  \param bodyUndergoingAcceleration Name of body undergoing the acceleration.
     *  \param bodyExertingAcceleration Name of body exerting the acceleration.
     *  \param useNorm Boolean denoting whether to use the norm (if true) or the vector (if false) of the acceleration.
     *  \param componentIndex Index of the component to be saved. Only applicable to vectorial dependent variables.
     *  By default -1, i.e. all the components are saved.
     */
    SingleAccelerationDependentVariableSaveSettings(
            const basic_astrodynamics::AvailableAcceleration accelerationModelType,
            const std::string& bodyUndergoingAcceleration,
            const std::string& bodyExertingAcceleration,
            const bool useNorm = 0,
            const int componentIndex = -1 ):
        SingleDependentVariableSaveSettings(
            ( useNorm == 1 ) ? ( single_acceleration_norm_dependent_variable ) : ( single_acceleration_dependent_variable ),
            bodyUndergoingAcceleration, bodyExertingAcceleration, componentIndex ),
        accelerationModelType_( accelerationModelType )
    { }

    // Type of acceleration that is to be saved.
    basic_astrodynamics::AvailableAcceleration accelerationModelType_;

};

// Class to define settings for saving contributions at separate degree/order to spherical harmonic acceleration.
class SphericalHarmonicAccelerationTermsDependentVariableSaveSettings: public SingleDependentVariableSaveSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param bodyUndergoingAcceleration Name of body undergoing the acceleration.
     *  \param bodyExertingAcceleration Name of body exerting the acceleration.
     *  \param componentIndices List of degree/order terms that are to be saved
     *  \param componentIndex Index of the acceleration vectors component to be saved. By default -1, i.e. all the components
     *  are saved.
     */
    SphericalHarmonicAccelerationTermsDependentVariableSaveSettings(
            const std::string& bodyUndergoingAcceleration,
            const std::string& bodyExertingAcceleration,
            const std::vector< std::pair< int, int > >& componentIndices,
            const int componentIndex = -1,
            const bool useAccelerationNorm = false ):
        SingleDependentVariableSaveSettings(
            useAccelerationNorm ? spherical_harmonic_acceleration_norm_terms_dependent_variable :
                                  spherical_harmonic_acceleration_terms_dependent_variable,
            bodyUndergoingAcceleration, bodyExertingAcceleration,
            componentIndex ), componentIndices_( componentIndices ){ }

    // Constructor.
    /*
     *  Constructor. for saving all terms up to a given degree/order.
     *  \param bodyUndergoingAcceleration Name of body undergoing the acceleration.
     *  \param bodyExertingAcceleration Name of body exerting the acceleration.
     *  \param maximumDegree Maximum degree to which terms are to be saved.
     *  \param maximumOrder Maximum order to which terms are to be saved.
     *  \param componentIndex Index of the acceleration vectors component to be saved. By default -1, i.e. all the components
     *  are saved.
     */
    SphericalHarmonicAccelerationTermsDependentVariableSaveSettings(
            const std::string& bodyUndergoingAcceleration,
            const std::string& bodyExertingAcceleration,
            const int maximumDegree,
            const int maximumOrder,
            const int componentIndex = -1,
            const bool useAccelerationNorm = false ):
        SingleDependentVariableSaveSettings(
            useAccelerationNorm ? spherical_harmonic_acceleration_norm_terms_dependent_variable :
                                  spherical_harmonic_acceleration_terms_dependent_variable,
            bodyUndergoingAcceleration, bodyExertingAcceleration,
            componentIndex )
    {
        for( int i = 0; i <= maximumDegree; i++ )
        {
            for( int j = 0; ( j <= i && j <= maximumOrder ); j++ )
            {
                componentIndices_.push_back( std::make_pair( i, j ) );
            }
        }
    }

    // List of degree/order terms that are to be saved
    std::vector< std::pair< int, int > > componentIndices_;
};

// Class to define settings for saving a single torque (norm or vector) during propagation.
class SingleTorqueDependentVariableSaveSettings: public SingleDependentVariableSaveSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param torqueModelType Type of torque that is to be saved.
     *  \param bodyUndergoingTorque Name of body undergoing the torque.
     *  \param bodyExertingTorque Name of body exerting the torque.
     *  \param useNorm Boolean denoting whether to use the norm (if true) or the vector (if false) of the torque.
     *  \param componentIndex Index of the component to be saved. Only applicable to vectorial dependent variables.
     *  By default -1, i.e. all the components are saved.
     */
    SingleTorqueDependentVariableSaveSettings(
            const basic_astrodynamics::AvailableTorque torqueModelType,
            const std::string& bodyUndergoingTorque,
            const std::string& bodyExertingTorque,
            const bool useNorm = 0,
            const int componentIndex = -1 ):
        SingleDependentVariableSaveSettings(
            ( useNorm == 1 ) ? ( single_torque_norm_dependent_variable ) : ( single_torque_dependent_variable ),
            bodyUndergoingTorque, bodyExertingTorque, componentIndex ),
        torqueModelType_( torqueModelType )
    { }

    // Boolean denoting whether to use the norm (if true) or the vector (if false) of the torque.
    basic_astrodynamics::AvailableTorque torqueModelType_;

};

// Class to define settings for saving a rotation matrix between two AerodynamicsReferenceFrames.
class IntermediateAerodynamicRotationVariableSaveSettings: public SingleDependentVariableSaveSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param associatedBody Body for which the rotation matrix is to be saved.
     *  \param baseFrame Frame from which rotation is to take place.
     *  \param targetFrame Frame to which the rotation is to take place.
     *  \param componentIndex Index of the component to be saved. Only applicable to vectorial dependent variables.
     *  By default -1, i.e. all the components are saved.
     */
    IntermediateAerodynamicRotationVariableSaveSettings(
            const std::string& associatedBody,
            const reference_frames::AerodynamicsReferenceFrames baseFrame,
            const reference_frames::AerodynamicsReferenceFrames targetFrame,
            const int componentIndex = -1 ):
        SingleDependentVariableSaveSettings( intermediate_aerodynamic_rotation_matrix_variable, associatedBody,
                                             "", componentIndex ),
        baseFrame_( baseFrame ), targetFrame_( targetFrame ){ }

    // Frame from which rotation is to take place.
    reference_frames::AerodynamicsReferenceFrames baseFrame_;

    // Frame to which the rotation is to take place.
    reference_frames::AerodynamicsReferenceFrames targetFrame_;

};

// Class to define settings for saving an aerodynamics orientation angle from AerodynamicsReferenceFrameAngles list.
class BodyAerodynamicAngleVariableSaveSettings: public SingleDependentVariableSaveSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param associatedBody Body for which the orientation angle is to be saved.
     *  \param angle Orientation angle that is to be saved.
     *  \param centralBody Body w.r.t. which angles are to be defined (only used to create flight conditions object if none exists
     *  yet).
     */
    BodyAerodynamicAngleVariableSaveSettings(
            const std::string& associatedBody,
            const reference_frames::AerodynamicsReferenceFrameAngles angle,
            const std::string& centralBody = "" ):
        SingleDependentVariableSaveSettings( relative_body_aerodynamic_orientation_angle_variable, associatedBody, centralBody ),
        angle_( angle ){ }

    // Orientation angle that is to be saved.
    reference_frames::AerodynamicsReferenceFrameAngles angle_;

};

// Class to define variations in spherical harmonic acceleration due to single gravity field variation.
class SingleVariationSphericalHarmonicAccelerationSaveSettings: public SingleDependentVariableSaveSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param bodyUndergoingAcceleration Body undergoing acceleration.
     *  \param bodyExertingAcceleration Body exerting acceleration.
     *  \param deformationType Type of gravity field variation.
     *  \param identifier Identifier for gravity field variation.
     */
    SingleVariationSphericalHarmonicAccelerationSaveSettings(
            const std::string& bodyUndergoingAcceleration,
            const std::string& bodyExertingAcceleration,
            const gravitation::BodyDeformationTypes deformationType,
            const std::string& identifier = "" ):
        SingleDependentVariableSaveSettings( single_gravity_field_variation_acceleration,
                                             bodyUndergoingAcceleration, bodyExertingAcceleration ),
        deformationType_( deformationType ), identifier_( identifier ){ }

    // Type of gravity field variation.
    gravitation::BodyDeformationTypes deformationType_;

    // Identifier for gravity field variation.
    std::string identifier_;

};

// Class to define variations in spherical harmonic acceleration due to single gravity field variation at separate degrees/orders.
class SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings: public SingleDependentVariableSaveSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param bodyUndergoingAcceleration Body undergoing acceleration.
     *  \param bodyExertingAcceleration Body exerting acceleration.
     *  \param componentIndices Degrees and orders for which to computed contribution.
     *  \param deformationType Type of gravity field variation.
     *  \param identifier Identifier for gravity field variation.
     */
    SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings(
            const std::string& bodyUndergoingAcceleration,
            const std::string& bodyExertingAcceleration,
            const std::vector< std::pair< int, int > >& componentIndices,
            const gravitation::BodyDeformationTypes deformationType,
            const std::string& identifier = "" ):
        SingleDependentVariableSaveSettings( single_gravity_field_variation_acceleration_terms,
                                             bodyUndergoingAcceleration, bodyExertingAcceleration ),
        componentIndices_( componentIndices ), deformationType_( deformationType ), identifier_( identifier ){ }

    // Constructor.
    /*
     *  Constructor.
     *  \param bodyUndergoingAcceleration Body undergoing acceleration.
     *  \param bodyExertingAcceleration Body exerting acceleration.
     *  \param maximumDegree Maximum degree for which to computed contribution.
     *  \param maximumOrder Maximum order for which to computed contribution.
     *  \param deformationType Type of gravity field variation.
     *  \param identifier Identifier for gravity field variation.
     */
    SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings(
            const std::string& bodyUndergoingAcceleration,
            const std::string& bodyExertingAcceleration,
            const int maximumDegree,
            const int maximumOrder,
            const gravitation::BodyDeformationTypes deformationType,
            const std::string& identifier = "" ):
        SingleDependentVariableSaveSettings( single_gravity_field_variation_acceleration_terms,
                                             bodyUndergoingAcceleration, bodyExertingAcceleration ),
        deformationType_( deformationType ), identifier_( identifier )
    {
        for( int i = 0; i <= maximumDegree; i++ )
        {
            for( int j = 0; ( j <= i && j <= maximumOrder ); j++ )
            {
                componentIndices_.push_back( std::make_pair( i, j ) );
            }
        }
    }

    // Degrees and orders for which to computed contribution.
    std::vector< std::pair< int, int > > componentIndices_;

    // Type of gravity field variation.
    gravitation::BodyDeformationTypes deformationType_;

    // Identifier for gravity field variation.
    std::string identifier_;

};

// Class to define .
class AccelerationPartialWrtStateSaveSettings: public SingleDependentVariableSaveSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param bodyUndergoingAcceleration Body undergoing acceleration.
     *  \param bodyExertingAcceleration Body exerting acceleration.
     *  \param accelerationModelType Type of acceleration that is to be saved.
     *  \param derivativeWrtBody String denoting w.r.t. which body the partial needs to be taken.
     *  \param thirdBody String denoting the third body w.r.t. which the partial needs to be taken (in case
     *      of third body acceleration).
     */
    AccelerationPartialWrtStateSaveSettings(
            const std::string& bodyUndergoingAcceleration,
            const std::string& bodyExertingAcceleration,
            const basic_astrodynamics::AvailableAcceleration accelerationModelType,
            const std::string derivativeWrtBody,
            const std::string thirdBody = "" ):
        SingleDependentVariableSaveSettings(
            acceleration_partial_wrt_body_translational_state, bodyUndergoingAcceleration, bodyExertingAcceleration ),
        accelerationModelType_( accelerationModelType ), derivativeWrtBody_( derivativeWrtBody ),
        thirdBody_( thirdBody ){ }

    // Type of acceleration that is to be saved.
    basic_astrodynamics::AvailableAcceleration accelerationModelType_;

    // String denoting w.r.t. which body the derivative needs to be taken.
    std::string derivativeWrtBody_;

    // String denoting the third body w.r.t. which the partial needs to be taken (in case of third body acceleration).
    std::string thirdBody_;

};

class CustomDependentVariableSaveSettings: public SingleDependentVariableSaveSettings
{
public:


    CustomDependentVariableSaveSettings(
            const std::function< Eigen::VectorXd( ) > customDependentVariableFunction,
            const int dependentVariableSize ):
        SingleDependentVariableSaveSettings(
            custom_dependent_variable, "", "" ),
        customDependentVariableFunction_( customDependentVariableFunction ), dependentVariableSize_( dependentVariableSize ){ }

    const std::function< Eigen::VectorXd( ) > customDependentVariableFunction_;

    const int dependentVariableSize_;

};




// Container class for settings of all dependent variables that are to be saved.
class DependentVariableSaveSettings
{
public:

    // Constructor.
    /*
     *  Constructor.
     *  \param dependentVariables List of settings for parameters that are to be saved.
     *  \param printDependentVariableTypes Variable denoting whether to print the list and vector entries of
     *      dependent variables when propagating.
     */
    DependentVariableSaveSettings(
            const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables,
            const bool printDependentVariableTypes = true ):
        dependentVariables_( dependentVariables ), printDependentVariableTypes_( printDependentVariableTypes ){ }

    // List of settings for parameters that are to be saved.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables_;

    // Variable denoting whether to print the list and vector entries of dependent variables when propagating.
    bool printDependentVariableTypes_;

#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
    // Map of state derivative partials, to be used when saving state derivative partials as dependent variables
    std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > stateDerivativePartials_;
#endif

};

//! @get_docstring(createDependentVariableSaveSettings.__docstring__)
inline std::shared_ptr< DependentVariableSaveSettings > createDependentVariableSaveSettings(
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables,
        const bool printDependentVariableTypes = true )
{
    return std::make_shared< DependentVariableSaveSettings >(
                dependentVariables, printDependentVariableTypes );
}

// Function to get a string representing a 'named identification' of a variable type.
/*
 *  Function to get a string representing a 'named identification' of a variable type.
 *  \param variableType Variable type.
 *  \return String with variable type id.
 */
std::string getVariableName( const VariableType variableType );

// Function to get a string representing a 'named identification' of a variable.
/*
 *  Function to get a string representing a 'named identification' of a variable.
 *  \param variableSettings Variable.
 *  \return String with variable id.
 */
std::string getVariableId( const std::shared_ptr< VariableSettings > variableSettings );

// Function to get a string representing a 'named identification' of a dependent variable type.
/*
 *  Function to get a string representing a 'named identification' of a dependent variable type.
 *  \param propagationDependentVariables Dependent variable type.
 *  \return String with dependent variable type id.
 */
std::string getDependentVariableName( const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings );

// Function to get a string representing a 'named identification' of a dependent variable.
/*
 *  Function to get a string representing a 'named identification' of a dependent variable.
 *  \param dependentVariableSettings Dependent variable.
 *  \return String with dependent variable id.
 */
std::string getDependentVariableId(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings );

//! @get_docstring(machNumberDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > machNumberDependentVariable(
        const std::string& associatedBody,
        const std::string& bodyWithAtmosphere )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                mach_number_dependent_variable, associatedBody, bodyWithAtmosphere );
}

//! @get_docstring(altitudeDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > altitudeDependentVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                altitude_dependent_variable, associatedBody, centralBody );
}

//! @get_docstring(airspeedDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > airspeedDependentVariable(
        const std::string& associatedBody,
        const std::string& bodyWithAtmosphere )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                airspeed_dependent_variable, associatedBody, bodyWithAtmosphere );
}

//! @get_docstring(densityDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > densityDependentVariable(
        const std::string& associatedBody,
        const std::string& bodyWithAtmosphere )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                local_density_dependent_variable, associatedBody, bodyWithAtmosphere );
}

//! @get_docstring(relativeSpeedDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > relativeSpeedDependentVariable(
        const std::string& associatedBody,
        const std::string& relativeBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                relative_speed_dependent_variable, associatedBody, relativeBody );
}

//! @get_docstring(relativePositionDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > relativePositionDependentVariable(
        const std::string& associatedBody,
        const std::string& relativeBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                relative_position_dependent_variable, associatedBody, relativeBody );
}

//! @get_docstring(relativeDistanceDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > relativeDistanceDependentVariable(
        const std::string& associatedBody,
        const std::string& relativeBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                relative_distance_dependent_variable, associatedBody, relativeBody );
}

//! @get_docstring(relativeVelocityDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > relativeVelocityDependentVariable(
        const std::string& associatedBody,
        const std::string& bodyWithAtmosphere )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                relative_velocity_dependent_variable, associatedBody, bodyWithAtmosphere );
}

//! @get_docstring(keplerianStateDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > keplerianStateDependentVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                keplerian_state_dependent_variable, associatedBody, centralBody );
}

//! @get_docstring(modifiedEquinoctialStateDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > modifiedEquinoctialStateDependentVariable(
		const std::string& associatedBody,
		const std::string& centralBody
		)
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			modified_equinocial_state_dependent_variable,
			associatedBody, centralBody);
}

//! @get_docstring(singleAccelerationDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > singleAccelerationDependentVariable(
        const basic_astrodynamics::AvailableAcceleration accelerationModelType,
        const std::string& bodyUndergoingAcceleration,
        const std::string& bodyExertingAcceleration )
{
    return std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                accelerationModelType, bodyUndergoingAcceleration, bodyExertingAcceleration,
                false );
}

//! @get_docstring(singleAccelerationNormDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > singleAccelerationNormDependentVariable(
        const basic_astrodynamics::AvailableAcceleration accelerationModelType,
        const std::string& bodyUndergoingAcceleration,
        const std::string& bodyExertingAcceleration )
{
    return std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                accelerationModelType, bodyUndergoingAcceleration, bodyExertingAcceleration,
                true );
}

//! @get_docstring(sphericalHarmonicAccelerationTermsDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > sphericalHarmonicAccelerationTermsDependentVariable(
        const std::string& bodyUndergoingAcceleration,
        const std::string& bodyExertingAcceleration,
        const std::vector< std::pair< int, int > >& componentIndices )
{
    return std::make_shared< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                bodyUndergoingAcceleration, bodyExertingAcceleration, componentIndices,
                -1, false );
}

//! @get_docstring(sphericalHarmonicAccelerationTermsNormDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > sphericalHarmonicAccelerationTermsNormDependentVariable(
        const std::string& bodyUndergoingAcceleration,
        const std::string& bodyExertingAcceleration,
        const std::vector< std::pair< int, int > >& componentIndices )
{
    return std::make_shared< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                bodyUndergoingAcceleration, bodyExertingAcceleration, componentIndices,
                -1, true );
}

//! @get_docstring(totalGravityFieldVariationAccelerationContributionVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > totalGravityFieldVariationAccelerationContributionVariable(
        const std::string& bodyUndergoingAcceleration,
        const std::string& bodyExertingAcceleration )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                total_gravity_field_variation_acceleration, bodyUndergoingAcceleration, bodyExertingAcceleration );
}

//! @get_docstring(singleGravityFieldVariationAccelerationContributionVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > singleGravityFieldVariationAccelerationContributionVariable(
        const std::string& bodyUndergoingAcceleration,
        const std::string& bodyExertingAcceleration,
        const gravitation::BodyDeformationTypes deformationType,
        const std::string& identifier = "" )
{
    return std::make_shared< SingleVariationSphericalHarmonicAccelerationSaveSettings >(
                bodyUndergoingAcceleration, bodyExertingAcceleration, deformationType, identifier );
}

//! @get_docstring(singleGravityFieldVariationSeparateTermsAccelerationContributionVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > singleGravityFieldVariationSeparateTermsAccelerationContributionVariable(
        const std::string& bodyUndergoingAcceleration,
        const std::string& bodyExertingAcceleration,
        const std::vector< std::pair< int, int > >& componentIndices,
        const gravitation::BodyDeformationTypes deformationType,
        const std::string& identifier = "" )
{
    return std::make_shared< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                bodyUndergoingAcceleration, bodyExertingAcceleration, componentIndices, deformationType,  identifier );
}

//! @get_docstring(totalAccelerationDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > totalAccelerationDependentVariable(
        const std::string& associatedBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                total_acceleration_dependent_variable, associatedBody );
}

//! @get_docstring(totalAccelerationNormDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > totalAccelerationNormDependentVariable(
        const std::string& associatedBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                total_acceleration_norm_dependent_variable, associatedBody );
}

//! @get_docstring(aerodynamicForceCoefficientDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > aerodynamicForceCoefficientDependentVariable(
        const std::string& associatedBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                aerodynamic_force_coefficients_dependent_variable, associatedBody );
}

//! @get_docstring(aerodynamicMomentCoefficientDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > aerodynamicMomentCoefficientDependentVariable(
        const std::string& associatedBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                aerodynamic_moment_coefficients_dependent_variable, associatedBody );
}

//! @get_docstring(inertialToBodyFixedRotationMatrixVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > inertialToBodyFixedRotationMatrixVariable(
        const std::string& associatedBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
            inertial_to_body_fixed_rotation_matrix_variable, associatedBody );
}

//! @get_docstring(intermediateAerodynamicRotationMatrixVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > intermediateAerodynamicRotationMatrixVariable(
        const std::string& associatedBody,
        const reference_frames::AerodynamicsReferenceFrames baseFrame,
        const reference_frames::AerodynamicsReferenceFrames targetFrame )
{
    return std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                associatedBody, baseFrame, targetFrame );
}

//! @get_docstring(latitudeDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > latitudeDependentVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                associatedBody, reference_frames::latitude_angle, centralBody );
}

//! @get_docstring(geodeticLatitudeDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > geodeticLatitudeDependentVariable(
		const std::string& associatedBody,
		const std::string& centralBody )
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			geodetic_latitude_dependent_variable, associatedBody, centralBody );
}

//! @get_docstring(longitudeDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > longitudeDependentVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                associatedBody, reference_frames::longitude_angle, centralBody );
}

//! @get_docstring(headingDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > headingDependentVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                associatedBody, reference_frames::heading_angle, centralBody );
}

//! @get_docstring(flightPathAngleDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > flightPathAngleDependentVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                associatedBody, reference_frames::flight_path_angle, centralBody );
}

//! @get_docstring(angleOfAttackDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > angleOfAttackDependentVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                associatedBody, reference_frames::angle_of_attack, centralBody );
}

//! @get_docstring(sideslipAngleDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > sideslipAngleDependentVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                associatedBody, reference_frames::angle_of_sideslip, centralBody );
}

//! @get_docstring(bankAngleDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > bankAngleDependentVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                associatedBody, reference_frames::bank_angle, centralBody );
}

//! @get_docstring(bodyFixedAirspeedBasedVelocityVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > bodyFixedAirspeedBasedVelocityVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                body_fixed_airspeed_based_velocity_variable, associatedBody );
}

//! @get_docstring(bodyFixedGroundspeedBasedVelocityVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > bodyFixedGroundspeedBasedVelocityVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                body_fixed_groundspeed_based_velocity_variable, associatedBody );
}

//! @get_docstring(tnwToInertialFrameRotationMatrixVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > tnwToInertialFrameRotationMatrixVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                tnw_to_inertial_frame_rotation_dependent_variable, associatedBody, centralBody );
}

//! @get_docstring(rswToInertialFrameRotationMatrixVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > rswToInertialFrameRotationMatrixVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                rsw_to_inertial_frame_rotation_dependent_variable, associatedBody, centralBody );
}

//! @get_docstring(periapsisAltitudeVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > periapsisAltitudeVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                periapsis_altitude_dependent_variable, associatedBody, centralBody );
}

//! @get_docstring(singleTorqueNormVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > singleTorqueNormVariable(
        const basic_astrodynamics::AvailableTorque torqueModelType,
        const std::string& bodyUndergoingTorque,
        const std::string& bodyExertingTorque )
{
    return std::make_shared< SingleTorqueDependentVariableSaveSettings >(
                torqueModelType, bodyUndergoingTorque, bodyExertingTorque, true );
}

//! @get_docstring(singleTorqueVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > singleTorqueVariable(
        const basic_astrodynamics::AvailableTorque torqueModelType,
        const std::string& bodyUndergoingTorque,
        const std::string& bodyExertingTorque )
{
    return std::make_shared< SingleTorqueDependentVariableSaveSettings >(
                torqueModelType, bodyUndergoingTorque, bodyExertingTorque, false );
}

//! @get_docstring(controlSurfaceDeflectionDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > controlSurfaceDeflectionDependentVariable(
		const std::string& associatedBody,
		const std::string& controlSurface )
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			control_surface_deflection_dependent_variable, associatedBody, controlSurface );
}

//! @get_docstring(radiationPressureDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > radiationPressureDependentVariable(
        const std::string& associatedBody,
        const std::string& radiatingBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
               radiation_pressure_dependent_variable, associatedBody, radiatingBody );
}

//! @get_docstring(localTemperatureDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > localTemperatureDependentVariable(
        const std::string& associatedBody )
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			local_temperature_dependent_variable, associatedBody );
}

//! @get_docstring(localDynamicPressureDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > localDynamicPressureDependentVariable(
		const std::string& associatedBody
)
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			local_dynamic_pressure_dependent_variable, associatedBody );
}

//! @get_docstring(localAerodynamicHeatRateDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > localAerodynamicHeatRateDependentVariable(
		const std::string& associatedBody
)
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			local_aerodynamic_heat_rate_dependent_variable, associatedBody );
}

//! @get_docstring(totalAerodynamicGLoadDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > totalAerodynamicGLoadDependentVariable(
		const std::string& associatedBody
)
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			total_aerodynamic_g_load_variable, associatedBody );
}

//! @get_docstring(stagnationPointHeatFluxDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > stagnationPointHeatFluxDependentVariable(
		const std::string& associatedBody
)
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			stagnation_point_heat_flux_dependent_variable, associatedBody );
}

//! @get_docstring(totalMassRateDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > totalMassRateDependentVariable(
		const std::string& associatedBody )
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			total_mass_rate_dependent_variables, associatedBody );
}

//! @get_docstring(totalTorqueNormDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > totalTorqueNormDependentVariable(
		const std::string& associatedBody )
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			total_torque_norm_dependent_variable, associatedBody );
}

//! @get_docstring(totalTorqueDependentVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > totalTorqueDependentVariable(
		const std::string& associatedBody )
{
	return std::make_shared< SingleDependentVariableSaveSettings >(
			total_torque_dependent_variable, associatedBody );
}

//! @get_docstring(centralBodyFixedSphericalPositionVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > centralBodyFixedSphericalPositionVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                body_fixed_relative_spherical_position,  associatedBody, centralBody );
}

//! @get_docstring(centralBodyFixedCartesianPositionVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > centralBodyFixedCartesianPositionVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                body_fixed_relative_cartesian_position,  associatedBody, centralBody );
}

//! @get_docstring(eulerAnglesToBodyFixed313Variable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > eulerAnglesToBodyFixed313Variable(
        const std::string& associatedBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                euler_angles_to_body_fixed_313,  associatedBody );
}

//! @get_docstring(bodyMassVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > bodyMassVariable(
        const std::string& associatedBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                current_body_mass_dependent_variable,  associatedBody );
}

//! @get_docstring(radiationPressureCoefficientVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > radiationPressureCoefficientVariable(
        const std::string& associatedBody,
        const std::string& emittingBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                radiation_pressure_coefficient_dependent_variable,  associatedBody, emittingBody );
}

//! @get_docstring(dynamicPressureVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > dynamicPressureVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                local_dynamic_pressure_dependent_variable,  associatedBody, centralBody );
}

//! @get_docstring(aerodynamicGLoadVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > aerodynamicGLoadVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                total_aerodynamic_g_load_variable,  associatedBody, centralBody );
}

//! @get_docstring(atmosphericTemperatureVariable)
inline std::shared_ptr< SingleDependentVariableSaveSettings > atmosphericTemperatureVariable(
        const std::string& associatedBody,
        const std::string& centralBody )
{
    return std::make_shared< SingleDependentVariableSaveSettings >(
                local_temperature_dependent_variable,  associatedBody, centralBody );
}

inline std::shared_ptr< SingleDependentVariableSaveSettings > customDependentVariable(
        const std::function< Eigen::VectorXd( ) > customDependentVariableFunction,
        const int dependentVariableSize  )
{
    return std::make_shared< CustomDependentVariableSaveSettings >(
                customDependentVariableFunction,  dependentVariableSize );
}





} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONOUTPUTSETTINGS_H
