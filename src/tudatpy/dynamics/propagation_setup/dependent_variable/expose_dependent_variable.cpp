/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_dependent_variable.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/basics/deprecationWarnings.h>
#include <tudat/simulation/propagation_setup.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudat
{
namespace propagators
{

std::shared_ptr< SingleDependentVariableSaveSettings > customDependentVariableDeprecated(
        const std::function< Eigen::VectorXd( ) > customDependentVariableFunction,
        const int dependentVariableSize )
{
    static bool isWarningPrinted = false;
    if( isWarningPrinted == false )
    {
        tudat::utilities::printDeprecationWarning(
                "tudatpy.dynamics.propagation_setup.dependent_"
                "variable.custom",
                "tudatpy.dynamics.propagation_setup.dependent_"
                "variable.custom_dependent_variable" );
        isWarningPrinted = true;
    }

    return customDependentVariable( customDependentVariableFunction, dependentVariableSize );
}

}  // namespace propagators
}  // namespace tudat

namespace tudatpy
{
namespace dynamics
{
namespace propagation_setup
{
namespace dependent_variable
{

void expose_dependent_variable_setup( py::module &m )
{
    //////////////////////////////////////////////////////////////////////////////////////
    /// ENUMS
    /// ////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////

    py::enum_< tp::PropagationDependentVariables >( m,
                                                    "PropagationDependentVariables",
                                                    R"doc(

         Enumeration of available propagation dependent variables.

         Enumeration of propagation dependent variables supported by tudat.





      )doc" )
            // C++ legacy variable names.
            .value( "mach_number_type",
                    tp::PropagationDependentVariables::mach_number_dependent_variable,
                    R"doc(
      )doc" )
            .value( "altitude_type",
                    tp::PropagationDependentVariables::altitude_dependent_variable,
                    R"doc(
      )doc" )
            .value( "airspeed_type",
                    tp::PropagationDependentVariables::airspeed_dependent_variable,
                    R"doc(
      )doc" )
            .value( "local_density_type",
                    tp::PropagationDependentVariables::local_density_dependent_variable,
                    R"doc(
      )doc" )
            .value( "relative_speed_type",
                    tp::PropagationDependentVariables::relative_speed_dependent_variable,
                    R"doc(
      )doc" )
            .value( "relative_position_type",
                    tp::PropagationDependentVariables::relative_position_dependent_variable,
                    R"doc(
      )doc" )
            .value( "relative_distance_type",
                    tp::PropagationDependentVariables::relative_distance_dependent_variable,
                    R"doc(
      )doc" )
            .value( "relative_velocity_type",
                    tp::PropagationDependentVariables::relative_velocity_dependent_variable,
                    R"doc(
      )doc" )
            .value( "radiation_pressure_type",
                    tp::PropagationDependentVariables::radiation_pressure_dependent_variable,
                    R"doc(
      )doc" )
            .value( "total_acceleration_norm_type",
                    tp::PropagationDependentVariables::total_acceleration_norm_dependent_variable,
                    R"doc(
      )doc" )
            .value( "single_acceleration_norm_type",
                    tp::PropagationDependentVariables::single_acceleration_norm_dependent_variable,
                    R"doc(
      )doc" )
            .value( "total_acceleration_type",
                    tp::PropagationDependentVariables::total_acceleration_dependent_variable,
                    R"doc(
      )doc" )
            .value( "single_acceleration_type",
                    tp::PropagationDependentVariables::single_acceleration_dependent_variable,
                    R"doc(
      )doc" )
            .value( "aerodynamic_force_coefficients_type",
                    tp::PropagationDependentVariables::aerodynamic_force_coefficients_dependent_variable,
                    R"doc(
      )doc" )
            .value( "aerodynamic_moment_coefficients_type",
                    tp::PropagationDependentVariables::aerodynamic_moment_coefficients_dependent_variable,
                    R"doc(
      )doc" )
            .value( "rotation_matrix_to_body_fixed_frame_type",
                    tp::PropagationDependentVariables::inertial_to_body_fixed_rotation_matrix_variable,
                    R"doc(
      )doc" )
            .value( "intermediate_aerodynamic_rotation_matrix_type",
                    tp::PropagationDependentVariables::intermediate_aerodynamic_rotation_matrix_variable,
                    R"doc(
      )doc" )
            .value( "relative_body_aerodynamic_orientation_angle_type",
                    tp::PropagationDependentVariables::relative_body_aerodynamic_orientation_angle_variable,
                    R"doc(
      )doc" )
            .value( "body_fixed_airspeed_based_velocity_type",
                    tp::PropagationDependentVariables::body_fixed_airspeed_based_velocity_variable,
                    R"doc(
      )doc" )
            .value( "total_aerodynamic_g_load_type",
                    tp::PropagationDependentVariables::total_aerodynamic_g_load_variable,
                    R"doc(
      )doc" )
            .value( "stagnation_point_heat_flux_type",
                    tp::PropagationDependentVariables::stagnation_point_heat_flux_dependent_variable,
                    R"doc(No documentation found.)doc" )
            .value( "local_temperature_type",
                    tp::PropagationDependentVariables::local_temperature_dependent_variable,
                    R"doc(
      )doc" )
            .value( "geodetic_latitude_type",
                    tp::PropagationDependentVariables::geodetic_latitude_dependent_variable,
                    R"doc(
      )doc" )
            .value( "control_surface_deflection_type",
                    tp::PropagationDependentVariables::control_surface_deflection_dependent_variable,
                    R"doc(
      )doc" )
            .value( "total_mass_rate_type",
                    tp::PropagationDependentVariables::total_mass_rate_dependent_variables,
                    R"doc(
      )doc" )
            .value( "tnw_to_inertial_frame_rotation_type",
                    tp::PropagationDependentVariables::tnw_to_inertial_frame_rotation_dependent_variable,
                    R"doc(
      )doc" )
            .value( "rsw_to_inertial_frame_rotation_type",
                    tp::PropagationDependentVariables::rsw_to_inertial_frame_rotation_dependent_variable,
                    R"doc(
      )doc" )
            .value( "periapsis_altitude_type",
                    tp::PropagationDependentVariables::periapsis_altitude_dependent_variable,
                    R"doc(
      )doc" )
            .value( "apoapsis_altitude_type",
                    tp::PropagationDependentVariables::apoapsis_altitude_dependent_variable,
                    R"doc(
      )doc" )
            .value( "total_torque_norm_type",
                    tp::PropagationDependentVariables::total_torque_norm_dependent_variable,
                    R"doc(
      )doc" )
            .value( "single_torque_norm_type",
                    tp::PropagationDependentVariables::single_torque_norm_dependent_variable,
                    R"doc(
      )doc" )
            .value( "total_torque_type",
                    tp::PropagationDependentVariables::total_torque_dependent_variable,
                    R"doc(
      )doc" )
            .value( "single_torque_type",
                    tp::PropagationDependentVariables::single_torque_dependent_variable,
                    R"doc(
      )doc" )
            .value( "body_fixed_groundspeed_based_velocity_type",
                    tp::PropagationDependentVariables::body_fixed_groundspeed_based_velocity_variable,
                    R"doc(
      )doc" )
            .value( "keplerian_state_type",
                    tp::PropagationDependentVariables::keplerian_state_dependent_variable,
                    R"doc(
      )doc" )
            .value( "modified_equinoctial_state_type",
                    tp::PropagationDependentVariables::modified_equinocial_state_dependent_variable,
                    R"doc(
      )doc" )
            .value( "spherical_harmonic_acceleration_terms_type",
                    tp::PropagationDependentVariables::spherical_harmonic_acceleration_terms_dependent_variable,
                    R"doc(
      )doc" )
            .value( "spherical_harmonic_acceleration_norm_terms_type",
                    tp::PropagationDependentVariables::spherical_harmonic_acceleration_norm_terms_dependent_variable,
                    R"doc(
      )doc" )
            .value( "body_fixed_relative_cartesian_position_type",
                    tp::PropagationDependentVariables::body_fixed_relative_cartesian_position,
                    R"doc(
      )doc" )
            .value( "body_fixed_relative_spherical_position_type",
                    tp::PropagationDependentVariables::body_fixed_relative_spherical_position,
                    R"doc(
      )doc" )
            .value( "total_gravity_field_variation_acceleration_type",
                    tp::PropagationDependentVariables::total_gravity_field_variation_acceleration,
                    R"doc(
      )doc" )
            .value( "single_gravity_field_variation_acceleration_type",
                    tp::PropagationDependentVariables::single_gravity_field_variation_acceleration,
                    R"doc(
      )doc" )
            .value( "single_gravity_field_variation_acceleration_terms_"
                    "type",
                    tp::PropagationDependentVariables::single_gravity_field_variation_acceleration_terms,
                    R"doc(
      )doc" )
            .value( "acceleration_partial_wrt_body_translational_state_"
                    "type",
                    tp::PropagationDependentVariables::acceleration_partial_wrt_body_translational_state,
                    R"doc(
      )doc" )
            .value( "local_dynamic_pressure_type",
                    tp::PropagationDependentVariables::local_dynamic_pressure_dependent_variable,
                    R"doc(
      )doc" )
            //                .value("local_aerodynamic_heat_rate_type",
            //                       tp::PropagationDependentVariables::local_aerodynamic_heat_rate_dependent_variable,
            //                       get_docstring("PropagationDependentVariables.local_aerodynamic_heat_rate_type").c_str())
            .value( "euler_angles_to_body_fixed_type",
                    tp::PropagationDependentVariables::euler_angles_to_body_fixed_313,
                    R"doc(
      )doc" )
            .value( "current_body_mass_type",
                    tp::PropagationDependentVariables::current_body_mass_dependent_variable,
                    R"doc(
      )doc" )
            .value( "radiation_pressure_coefficient_type",
                    tp::PropagationDependentVariables::radiation_pressure_coefficient_dependent_variable,
                    R"doc(
      )doc" )
            .value( "custom_type", tp::PropagationDependentVariables::custom_dependent_variable, R"doc(No documentation found.)doc" )
            .value( "gravity_field_potential_type",
                    tp::PropagationDependentVariables::gravity_field_potential_dependent_variable,
                    R"doc(
      )doc" )
            .value( "gravity_field_laplacian_of_potential_type",
                    tp::PropagationDependentVariables::gravity_field_laplacian_of_potential_dependent_variable,
                    R"doc(
      )doc" )
            .export_values( );

    //////////////////////////////////////////////////////////////////////////////////////
    /// CLASSES
    /// //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////

    py::class_< tp::VariableSettings, std::shared_ptr< tp::VariableSettings > >( m,
                                                                                 "VariableSettings",
                                                                                 R"doc(

         Functional base class to define settings for variables.

         This class is a functional base class for defining settings for variables.
         Any variable that requires additional information in addition to what can be provided here, should be defined by a
         dedicated derived class.





      )doc" );

    py::class_< tp::SingleDependentVariableSaveSettings, std::shared_ptr< tp::SingleDependentVariableSaveSettings >, tp::VariableSettings >(
            m,
            "SingleDependentVariableSaveSettings",
            R"doc(

         `VariableSettings`-derived class to define settings for dependent variables that are to be saved during propagation.

         Functional base class for defining settings for dependent variables that are to be computed and saved during propagation.
         Any dependent variable that requires additional information in addition to what can be provided here, should be
         defined by a dedicated derived class.





      )doc" )
            .def_property_readonly( "dependent_variable_type",
                                    &tp::SingleDependentVariableSaveSettings::getDependentVariableType,
                                    R"doc(
                                    
Type of the dependent variable that is being saved.

:type: PropagationDependentVariables

)doc" )
            .def_property_readonly( "associated_body",
                                    &tp::SingleDependentVariableSaveSettings::getAssociatedBody,
                                    R"doc(
Body associated with variable.

:type: str

)doc" )
            .def_property_readonly( "secondary_body",
                                    &tp::SingleDependentVariableSaveSettings::getSecondaryBody,
                                    R"doc(
                                    
Secondary body w.r.t. which parameter is defined, e.g. relative position is defined for the associated body w.r.t. secondary body.

.. note::

    This parameter is not necessarily defined. In this case an empty string is returned.                    
                                    
:type: str

                                    )doc" )
            .def_property_readonly( "component_index",
                                    &tp::SingleDependentVariableSaveSettings::getComponentIndex,
                                    R"doc(
                                    
Index of the component to be saved.

.. note::

    This parameter is only applicable for vectorial dependent variables. If negative, all components of the vector are saved.                    
                                    
:type: int

                                    )doc" );
    //            .def(py::init<
    //                 const tp::PropagationDependentVariables,
    //                 const std::string &,
    //                 const std::string &,
    //                 const int>(),
    //                 py::arg("dependent_variable_type"),
    //                 py::arg("associated_body"),
    //                 py::arg("secondary_body") = "",
    //                 py::arg("component_idx") = -1);

    py::class_< tp::SingleAccelerationDependentVariableSaveSettings,
                std::shared_ptr< tp::SingleAccelerationDependentVariableSaveSettings >,
                tp::SingleDependentVariableSaveSettings >( m,
                                                           "SingleAccelerationDependentVariableSaveSettings",
                                                           R"doc(

         `SingleDependentVariableSaveSettings`-derived class to save a single acceleration (norm or vector) during propagation.

         Class to define settings for saving a single acceleration (norm or vector) during propagation. Note: this acceleration is returned in the inertial frame!




      )doc" )
            .def_property_readonly( "acceleration_model_type",
                                    &tp::SingleAccelerationDependentVariableSaveSettings::getAccelerationModelType,
                                    R"doc(
                                    
The type of the acceleration that is to be saved.

:type: AvailableAcceleration

)doc" );

    //            .def(py::init<
    //                 const
    //                 tudat::basic_astrodynamics::AvailableAcceleration,
    //                 const std::string &,
    //                 const std::string &,
    //                 const bool,
    //                 const int>(),
    //                 py::arg("acceleration_model_type"),
    //                 py::arg("body_undergoing_acceleration"),
    //                 py::arg("body_exerting_acceleration"),
    //                 py::arg("use_norm") = 0,
    //                 py::arg("component_index") = -1);

    m.def( "get_dependent_variable_id",
           &tp::getDependentVariableId,
           py::arg( "dependent_variable_settings" ),
           R"doc(

 Function to return a string representation of the dependent variable.

 Parameters
 ----------
 dependent_variable_settings : SingleDependentVariableSaveSettings
     Dependent variable settings for which the representation is to be returned.
 Returns
 -------
 str
     String representation of the dependent variable.
 Examples
 --------

 .. code-block:: python

     from tudatpy.dynamics import propagation_setup

     dep_var_settings = propagation_setup.dependent_variable.single_acceleration(propagation_setup.acceleration.point_mass_gravity_type, "SC", "Earth")

     print(propagation_setup.dependent_variable.get_dependent_variable_id(dep_var_settings))
     # prints "Single acceleration in inertial frame of type central gravity , acting on SC, exerted by Earth"

     )doc" );

    m.def( "get_dependent_variable_size",
           &tp::getDependentVariableSaveSize,
           py::arg( "dependent_variable_settings" ),
           py::arg( "bodies" ),
           R"doc(

 Function to get the size of a dependent variable.

 Parameters
 ----------
 dependent_variable_settings : SingleDependentVariableSaveSettings
     Dependent variable settings for which the size is to be returned.
 bodies : SystemOfBodies
     Bodies that contain the bodies associated with the dependent variable.
 Returns
 -------
 int
     Size of the dependent variable.
 Examples
 --------

 .. code-block:: python

     from tudatpy.dynamics import propagation_setup

     bodies = ...

     dep_var_settings = propagation_setup.dependent_variable.single_acceleration(propagation_setup.acceleration.point_mass_gravity_type, "SC", "Earth")

     print(propagation_setup.dependent_variable.get_dependent_variable_size(dep_var_settings, bodies))
     # prints 3

     )doc" );

    m.def( "get_dependent_variable_shape",
           &tp::getDependentVariableShape,
           py::arg( "dependent_variable_settings" ),
           py::arg( "bodies" ),
           R"doc(

 Function to get the shape of a dependent variable.

 Parameters
 ----------
 dependent_variable_settings : SingleDependentVariableSaveSettings
     Dependent variable settings for which the shape is to be returned.
 bodies : SystemOfBodies
     Bodies that contain the bodies associated with the dependent variable.
 Returns
 -------
 tuple[int, int]
     Shape of the dependent variable.
 Examples
 --------

 .. code-block:: python

     from tudatpy.dynamics import propagation_setup

     bodies = ...

     dep_var_settings = propagation_setup.dependent_variable.single_acceleration(propagation_setup.acceleration.point_mass_gravity_type, "SC", "Earth")

     print(propagation_setup.dependent_variable.get_dependent_variable_shape(dep_var_settings, bodies))
     # prints (3, 1)

     )doc" );

    //////////////////////////////////////////////////////////////////////////////////////
    /// FREE FUNCTIONS
    /// ///////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////

    m.def( "mach_number",
           &tp::machNumberDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the Mach number to the dependent variables to save.

 Function to add the Mach number to the dependent variables to save. The calculation of the altitude uses the atmosphere model of the central body and the current state of the body for which the Mach number is to be calculated.

 Parameters
 ----------
 body : str
     Body whose Mach number is to be saved.
 central_body : str
     Body with atmosphere with respect to which the Mach number is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1




 Examples
 --------

 To create settings for saving of a Mach number of a body name 'Spacecraft'
 w.r.t. the atmosphere of body 'Earth', use:

 .. code-block:: python

    # Define save settings for Mach number
    propagation_setup.dependent_variable.mach_number( "Spacecraft", "Earth" )


     )doc" );

    m.def( "altitude",
           &tp::altitudeDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the altitude to the dependent variables to save.

 Function to add the altitude to the dependent variables to save. The calculation of the altitude uses the shape model of the central body and the current state of the body for which the altitude is to be calculated.

 Parameters
 ----------
 body : str
     Body whose altitude is to be saved.
 central_body : str
     Body with respect to which the altitude is computed (requires this body to have a shape model defined).
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "airspeed",
           &tp::airspeedDependentVariable,
           py::arg( "body" ),
           py::arg( "body_with_atmosphere" ),
           R"doc(

 Function to add the airspeed to the dependent variables to save.

 Function to add the airspeed to the dependent variables to save. The calculation of the airspeed uses the rotation and wind models of the central body (to determine the motion of the atmosphere in inertial space), and the current state of the body for which the airspeed is to be calculated.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with atmosphere with respect to which the airspeed is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "density",
           &tp::densityDependentVariable,
           py::arg( "body" ),
           py::arg( "body_with_atmosphere" ),
           R"doc(

 Function to add the local freestream density to the dependent variables to save.

 Function to add the freestream density (at a body's position) to the dependent variables to save. The calculation of the density uses the atmosphere model of the central body, and the current state of the body for which the density is to be calculated.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 body_with_atmosphere : str
     Body with atmosphere with respect to which the density is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "temperature",
           &tp::localTemperatureDependentVariable,
           py::arg( "body" ),
           py::arg( "body_with_atmosphere" ),
           R"doc(

 Function to add the local freestream temperature to the dependent variables to save.

 Function to add the freestream temperature (at a body's position) to the dependent variables to save. The calculation of the temperature uses the atmosphere model of the central body, and the current state of the body for which the temperature is to be calculated.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 body_with_atmosphere : str
     Body with atmosphere with respect to which the temperature is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "dynamic_pressure",
           &tp::localDynamicPressureDependentVariable,
           py::arg( "body" ),
           py::arg( "body_with_atmosphere" ),
           R"doc(

 Function to add the local freestream dynamic pressure to the dependent variables to save.

 Function to add the freestream dynamic pressure (at a body's position) to the dependent variables to save. The calculation of the temperature uses the atmosphere model of the central body, and the current state of the body for which the temperature is to be calculated.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 body_with_atmosphere : str
     Body with atmosphere with respect to which the dynamic pressure is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    //        m.def("local_aerodynamic_heat_rate",
    //              &tp::localAerodynamicHeatRateDependentVariable,
    //              py::arg("body"),
    //              get_docstring("local_aerodynamic_heat_rate").c_str());

    m.def( "local_aerodynamic_g_load",
           &tp::totalAerodynamicGLoadDependentVariable,
           py::arg( "body" ),
           py::arg( "body_with_atmosphere" ),
           R"doc(

 Function to add the total aerodynamic G-load to the dependent variables to save.

 Function to add the total aerodynamic G-load of a body to the dependent variables to save. The calculation uses the atmosphere model of the central body, and the current state of the body for which the temperature is to be calculated.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 body_with_atmosphere : str
     Body with atmosphere exerting the aerodynamic acceleration for which the g-load is to be computed
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "relative_speed",
           &tp::relativeSpeedDependentVariable,
           py::arg( "body" ),
           py::arg( "relative_body" ),
           R"doc(

 Function to add the relative speed to the dependent variables to save.

 Function to add a body's relative speed (norm of the relative velocity vector) with respect to a second body to the dependent variables to save. The relative speed is computed between the bodies' centers of mass.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 relative_body : str
     Body with respect to which the relative speed is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "relative_position",
           &tp::relativePositionDependentVariable,
           py::arg( "body" ),
           py::arg( "relative_body" ),
           R"doc(

 Function to add the relative position vector to the dependent variables to save.

 Function to add a body's relative position vector with respect to a second body to the dependent variables to save. The relative position is computed between the bodies' centers of mass.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 relative_body : str
     Body with respect to which the relative position is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "relative_distance",
           &tp::relativeDistanceDependentVariable,
           py::arg( "body" ),
           py::arg( "relative_body" ),
           R"doc(

 Function to add the relative distance to the dependent variables to save.

 Function to add a body's relative distance (norm of the position vector) with respect to a second body to the dependent variables to save. The relative distance is computed between the bodies' centers of mass.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 relative_body : str
     Body with respect to which the relative distance is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "relative_velocity",
           &tp::relativeVelocityDependentVariable,
           py::arg( "body" ),
           py::arg( "relative_body" ),
           R"doc(

 Function to add the relative velocity vector to the dependent variables to save.

 Function to add a body's relative velocity vector with respect to a second body to the dependent variables to save. The relative velocity is computed between the bodies' centers of mass.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 relative_body : str
     Body with respect to which the relative velocity is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "keplerian_state",
           &tp::keplerianStateDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the Keplerian state to the dependent variables to save.

 Function to add the Keplerian state to the dependent variables to save. The Keplerian state is returned in this order: 1: Semi-major Axis. 2: Eccentricity. 3: Inclination. 4: Argument of Periapsis. 5. Right Ascension of the Ascending Node. 6: True Anomaly.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the Keplerian state is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 6





     )doc" );

    m.def( "modified_equinoctial_state",
           &tp::modifiedEquinoctialStateDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the modified equinoctial state to the dependent variables to save.

 Function to add the modified equinoctial state to the dependent variables to save. The value of the parameter I is automatically chosen as +1 or -1, depending on whether the inclination is smaller or larger than 90 degrees. The elements are returned in the order :math:`p`, :math:`f`, :math:`g`, :math:`h`, :math:`k`, :math:`L`

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the modified equinoctial state is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 6





     )doc" );

    m.def( "single_acceleration",
           &tp::singleAccelerationDependentVariable,
           py::arg( "acceleration_type" ),
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           R"doc(

 Function to add a single acceleration to the dependent variables to save.

 Function to add a single acceleration vector to the dependent variables to save. The requested acceleration is defined by its type, and the bodies undergoing and exerting the acceleration. This acceleration vector represents the acceleration in 3D in the inertial reference frame. NOTE: When requesting a third-body perturbation be saved, you may use either the direct acceleration type, or the third body type. For instance, for saving a point-mass third-body perturbation, you may specify either ``point_mass_gravity_type`` or ``third_body_point_mass_gravity_type`` as acceleration type.

 Parameters
 ----------
 acceleration_type : AvailableAcceleration
     Acceleration type to be saved.
 body_undergoing_acceleration : str
     Body undergoing acceleration.
 body_exerting_acceleration : str
     Body exerting acceleration.
 Returns
 -------
 SingleAccelerationDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3




 Examples
 --------

 To create settings for saving a point mass acceleration acting on body called 'Spacecraft',
 exerted by a body named 'Earth', use:

 .. code-block:: python

    # Define save settings for point-mass acceleration on Spacecraft by Earth
    propagation_setup.dependent_variable.single_acceleration(
            propagation_setup.acceleration.point_mass_gravity_type, 'Spacecraft', 'Earth' )


     )doc" );

    m.def( "single_acceleration_norm",
           &tp::singleAccelerationNormDependentVariable,
           py::arg( "acceleration_type" ),
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           R"doc(

 Function to add a single scalar acceleration to the dependent variables to save.

 Function to add a single scalar acceleration (norm of the acceleration vector) to the dependent variables to save. The requested acceleration is defined by its type, and the bodies undergoing and exerting the acceleration. NOTE: When requesting a third-body perturbation be saved, you may use either the direct acceleration type, or the third body type. For instance, for saving a point-mass third-body perturbation, you may specify either ``point_mass_gravity_type`` or ``third_body_point_mass_gravity_type`` as acceleration type.

 Parameters
 ----------
 acceleration_type : AvailableAcceleration
     Acceleration type to be saved
 body_undergoing_acceleration : str
     Body undergoing acceleration.
 body_exerting_acceleration : str
     Body exerting acceleration.
 Returns
 -------
 SingleAccelerationDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1




 Examples
 --------

 To create settings for saving norm of a point mass acceleration acting on body called 'Spacecraft',
 exerted by a body named 'Earth', use:

 .. code-block:: python

    # Define save settings for point-mass acceleration on Spacecraft by Earth
    propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.point_mass_gravity_type, 'Spacecraft', 'Earth' )


     )doc" );

    m.def( "total_acceleration_norm",
           &tp::totalAccelerationNormDependentVariable,
           py::arg( "body" ),
           R"doc(

 Function to add the total scalar acceleration (norm of the vector) acting on a body to the dependent variables to save.


 Parameters
 ----------
 body : str
     Body undergoing acceleration.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "total_acceleration",
           &tp::totalAccelerationDependentVariable,
           py::arg( "body" ),
           R"doc(

 Function to add the total acceleration vector acting on a body to the dependent variables to save.


 Parameters
 ----------
 body : str
     Body undergoing acceleration.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "single_torque_norm",
           &tp::singleTorqueNormVariable,
           py::arg( "torque_type" ),
           py::arg( "body_undergoing_torque" ),
           py::arg( "body_exerting_torque" ),
           R"doc(

 Function to add a single torque (norm of the torque vector) to the dependent variables to save.


 Parameters
 ----------
 torque_type : AvailableTorque
     Torque type to be saved.
 body_undergoing_torque : str
     Body undergoing torque.
 body_exerting_torque : str
     Body exerting torque.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "single_torque",
           &tp::singleTorqueVariable,
           py::arg( "torque_type" ),
           py::arg( "body_undergoing_torque" ),
           py::arg( "body_exerting_torque" ),
           R"doc(

 Function to add a single torque vector to the dependent variables to save.


 Parameters
 ----------
 torque_type : AvailableTorque
     Torque type to be saved.
 body_undergoing_torque : str
     Body undergoing torque.
 body_exerting_torque : str
     Body exerting torque.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "total_torque_norm",
           &tp::totalTorqueNormDependentVariable,
           py::arg( "body" ),
           R"doc(

 Function to add the total torque (norm of the torque vector) to the dependent variables to save.


 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "total_torque",
           &tp::totalTorqueDependentVariable,
           py::arg( "body" ),
           R"doc(

 Function to add the total torque vector to the dependent variables to save.


 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "spherical_harmonic_terms_acceleration",
           &tp::sphericalHarmonicAccelerationTermsDependentVariable,
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           py::arg( "component_indices" ),
           R"doc(

 Function to add single degree/order contributions of a spherical harmonic acceleration vector to the dependent variables to save.

 Function to add single degree/order contributions of a spherical harmonic acceleration vector to the dependent variables to save. The spherical harmonic acceleration consists of a (truncated) summation of contributions at degree :math:`l` and order :math:`m`. Using this function, you can save the contributions of separate :math:`l,m` entries to the total acceleration. For instance, when requesting dependent variables for :math:`l,m=2,2`, the contribution due to the combined influence of :math:`\bar{C}_{22}` and `\bar{S}_{22}` are provided

 Parameters
 ----------
 body_undergoing_acceleration : str
     Body undergoing acceleration.
 body_exerting_acceleration : str
     Body exerting acceleration.
 component_indices : list[tuple]
     Tuples of (degree, order) indicating the terms to save.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3*``len(component_indices)``




 Examples
 --------

 To create settings for saving spherical harmonic acceleration contributions of degree/order
 2/0, 2/1 and 2/2, acting on a body names 'Spacecraft', exerted by a body named 'Earth',
 use the following for the acceleration. The resulting dependent variable will contain
 nine entries (three acceleration components for 2/0, 2/1 and 2/2, respectively).

 .. code-block:: python

    # Define degree/order combinations for which to save acceleration contributions
    spherical_harmonic_terms = [ (2,0), (2,1), (2,2) ]

    # Define save settings for separate spherical harmonic contributions
    propagation_setup.dependent_variable.spherical_harmonic_terms_acceleration( "Spacecraft", "Earth", spherical_harmonic_terms )


     )doc" );

    m.def( "spherical_harmonic_terms_acceleration_norm",
           &tp::sphericalHarmonicAccelerationTermsNormDependentVariable,
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           py::arg( "component_indices" ),
           R"doc(

 Function to add a single term of the spherical harmonic acceleration (norm of the vector) to the dependent variables to save.

 Function to add single term of the spherical harmonic acceleration (norm of the vector) to the dependent variables to save.

 Parameters
 ----------
 body_undergoing_acceleration : str
     Body undergoing acceleration.
 body_exerting_acceleration : str
     Body exerting acceleration.
 component_indices : list[tuple]
     Tuples of (degree, order) indicating the terms to save.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 ``len(component_indices)``




 Examples
 --------

 To create settings for saving spherical harmonic acceleration contributions of degree/order
 2/0, 2/1 and 2/2, acting on a body names 'Spacecraft', exerted by a body named 'Earth',
 use the following for the acceleration. The resulting dependent variable will contain
 three entries (one acceleration norm for 2/0, 2/1 and 2/2, respectively).

 .. code-block:: python


    # Define degree/order combinations for which to save acceleration contributions
    spherical_harmonic_terms = [ (2,0), (2,1), (2,2) ]

    # Define save settings for separate spherical harmonic contributions
    propagation_setup.dependent_variable.spherical_harmonic_terms_acceleration_norm( "Spacecraft", "Earth", spherical_harmonic_terms )


     )doc" );

    m.def( "aerodynamic_force_coefficients",
           &tp::aerodynamicForceCoefficientDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ) = "",
           R"doc(

 Function to add the aerodynamic force coefficients to the dependent variables to save.

 Function to add the aerodynamic force coefficients to the dependent variables to save. It requires an aerodynamic coefficient interface to be defined for the vehicle. The coefficients are returned in the following order: C_D, C_S, C_l (if coefficient interface defined in aerodynamic frame), or C_X, C_Y, C_Z (if coefficient interface defined in body frame).

 Parameters
 ----------
 body : str
     Body undergoing acceleration.
 central_body : str
     Body exerting acceleration (e.g. body with atmosphere).
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "aerodynamic_moment_coefficients",
           &tp::aerodynamicMomentCoefficientDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ) = "",
           R"doc(

 Function to add the aerodynamic moment coefficients to the dependent variables to save.

 Function to add the aerodynamic force coefficients to the dependent variables to save. It requires an aerodynamic coefficient interface to be defined for the vehicle. The coefficients are returned in the following order: C_l, C_m, C_n , respectively about the X, Y, Z axes of the body-fixed frame, see (see :cite:t:`mooij1994`)

 Parameters
 ----------
 body : str
     Body undergoing acceleration.
 central_body : str
     Body exerting acceleration (e.g. body with atmosphere).
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3


     )doc" );

    m.def( "aerodynamic_force_coefficients_control_surface_free",
           &tp::aerodynamicForceCoefficientControlSurfaceFreeDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ) = "",
           R"doc(No documentation found.)doc" );

    m.def( "aerodynamic_moment_coefficients_control_surface_free",
           &tp::aerodynamicMomentCoefficientControlSurfaceFreeDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ) = "",
           R"doc(No documentation found.)doc" );

    m.def( "aerodynamic_force_coefficients_control_surface_"
           "increment",
           &tp::aerodynamicForceCoefficientControlSurfaceIncrementDependentVariable,
           py::arg( "body" ),
           py::arg( "control_surface_name" ),
           py::arg( "central_body" ) = "",
           R"doc(No documentation found.)doc" );

    m.def( "aerodynamic_moment_coefficients_control_surface_"
           "increment",
           &tp::aerodynamicMomentCoefficientControlSurfaceIncrementDependentVariable,
           py::arg( "body" ),
           py::arg( "control_surface_name" ),
           py::arg( "central_body" ) = "",
           R"doc(No documentation found.)doc" );

    m.def( "latitude",
           &tp::latitudeDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the latitude to the dependent variables to save.

 Function to add the latitude of a body, in the body-fixed frame of a central body, to the dependent variables to save.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the latitude is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "geodetic_latitude",
           &tp::geodeticLatitudeDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the geodetic latitude to the dependent variables to save.

 Function to add the geodetic latitude, in the body-fixed frame of a central body, to the dependent variables to save. If the central body has a spherical shape model, this value is identical to the latitude. If the central body has an oblate spheroid shape model, the calculation of the geodetic latitude uses the flattening of the this shape model to determine the geodetic latitude

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the geodetic latitude is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "longitude",
           &tp::longitudeDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the longitude to the dependent variables to save.

 Function to add the longitude of a body, in the body-fixed frame of a central body, to the dependent variables to save.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the longitude is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "heading_angle",
           &tp::headingDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the heading angle to the dependent variables to save.

 Function to add the heading angle to the dependent variables to save, as defined by :cite:t:`mooij1994`.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the heading angle is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "flight_path_angle",
           &tp::flightPathAngleDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the flight path angle to the dependent variables to save.

 Function to add the flight path angle to the dependent variables to save, as defined by :cite:t:`mooij1994`.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the flight path angle is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "angle_of_attack",
           &tp::angleOfAttackDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the angle of attack to the dependent variables to save.

 Function to add the angle of attack angle to the dependent variables to save, as defined by :cite:t:`mooij1994`.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the angle of attack is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "sideslip_angle",
           &tp::sideslipAngleDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the sideslip angle to the dependent variables to save, as defined by :cite:t:`mooij1994`.


 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the sideslip angle is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "bank_angle",
           &tp::bankAngleDependentVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the bank angle to the dependent variables to save, as defined by :cite:t:`mooij1994`.


 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the bank angle is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "radiation_pressure",
           &tp::radiationPressureDependentVariable,
           py::arg( "target_body" ),
           py::arg( "source_body" ),
           R"doc(

 Function to add the radiation pressure to the dependent variables to save.

 Function to add the local radiation pressure, in N/m^2, to the dependent variables to save. It requires a radiation source model to be defined for the radiating body.

 Parameters
 ----------
 target_body : str
     Name of body at the location of which the radiation pressure is to be returned
 source_body : str
     Name of body from which the radiation originates that causes the radiation pressure at the target body
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "total_gravity_field_variation_acceleration",
           &tp::totalGravityFieldVariationAccelerationContributionVariable,
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           R"doc(

 Function to add the acceleration induced by the total time-variability of a gravity field to the dependent variables to save.

 Function to add the acceleration induced by the total time-variability of a gravity field to the dependent variables to save. This function does not distinguish between different sources of variations of the gravity field, and takes the full time-variation when computing the contribution to the acceleration. To select only one contribution, use the :func:`single_gravity_field_variation_acceleration` function.

 Parameters
 ----------
 body_undergoing_acceleration : str
     Body whose dependent variable should be saved.
 body_exerting_acceleration : str
     Body exerting the acceleration.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "single_gravity_field_variation_acceleration",
           &tp::singleGravityFieldVariationAccelerationContributionVariable,
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           py::arg( "deformation_type" ),
           py::arg( "identifier" ) = "",
           R"doc(

 Function to add the acceleration induced by a single time-variability of a gravity field to the dependent variables to save.

 Function to add the acceleration induced by a single time-variability of a gravity field to the dependent variables to save. The user specifies the type of variability for which the induced acceleration is to be saved.

 Parameters
 ----------
 body_undergoing_acceleration : str
     Body whose dependent variable should be saved.
 body_exerting_acceleration : str
     Body exerting the acceleration.
 deformation_type : BodyDeformationTypes
     Type of gravity field variation for which the acceleration contribution is to be saved
 identifier : str, default=""
     Identifier for the deformation type. To be used in case multiple realizations of a single variation type are present in the given body. Otherwise, this entry can be left empty
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "single_per_term_gravity_field_variation_acceleration",
           &tp::singleGravityFieldVariationSeparateTermsAccelerationContributionVariable,
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           py::arg( "component_indices" ),
           py::arg( "deformation_type" ),
           py::arg( "identifier" ) = "",
           R"doc(

 Function to add the acceleration induced by a single time-variability of a gravity field, at a given list of degrees/orders, to the dependent variables to save. This combines the functionality of the :func:`single_gravity_field_variation_acceleration` and :func:`spherical_harmonic_terms_acceleration` variables


 Parameters
 ----------
 body_undergoing_acceleration : str
     Body whose dependent variable should be saved.
 body_exerting_acceleration : str
     Body exerting the acceleration.
 component_indices : list[tuple]
     Tuples of (degree, order) indicating the terms to save.
 deformation_type : BodyDeformationTypes
     Type of gravity field variation for which the acceleration contribution is to be saved
 identifier : str, default=""
     Identifier for the deformation type. To be used in case multiple realizations of a single variation type are present in the given body. Otherwise, this entry can be left empty
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3*``len(component_indices)``





     )doc" );

    m.def( "total_spherical_harmonic_cosine_coefficient_variations",
           &tp::totalSphericalHarmonicCosineCoefficientVariation,
           py::arg( "body" ),
           py::arg( "minimum_degree" ),
           py::arg( "maximum_degree" ),
           py::arg( "minimum_order" ),
           py::arg( "maximum_order" ),
           R"doc(No documentation found.)doc" );

    m.def( "total_spherical_harmonic_sine_coefficient_variations",
           &tp::totalSphericalHarmonicSineCoefficientVariation,
           py::arg( "body" ),
           py::arg( "minimum_degree" ),
           py::arg( "maximum_degree" ),
           py::arg( "minimum_order" ),
           py::arg( "maximum_order" ),
           R"doc(No documentation found.)doc" );

    m.def( "total_spherical_harmonic_cosine_coefficient_variations_"
           "from_indices",
           &tp::totalSphericalHarmonicCosineCoefficientVariationFromIndices,
           py::arg( "body" ),
           py::arg( "component_indices" ),
           R"doc(No documentation found.)doc" );

    m.def( "total_spherical_harmonic_sine_coefficient_variations_"
           "from_indices",
           &tp::totalSphericalHarmonicSineCoefficientVariationFromIndices,
           py::arg( "body" ),
           py::arg( "component_indices" ),
           R"doc(No documentation found.)doc" );

    m.def( "body_fixed_airspeed_velocity",
           &tp::bodyFixedAirspeedBasedVelocityVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the airspeed velocity vector to the dependent variables to save.

 Function to add the airspeed velocity vector to the dependent variables to save. The airspeed velocity vector is *not provided in an inertial frame*, but instead a frame centered on, and fixed to, the central body. It defines the velocity vector of a body w.r.t. the relative atmosphere It requires the central body to have an atmosphere.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the airspeed is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "body_fixed_groundspeed_velocity",
           &tp::bodyFixedGroundspeedBasedVelocityVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the groundspeed velocity vector to the dependent variables to save.

 Function to add the groundspeed velocity vector to the dependent variables to save. The groundspeed velocity vector is *not provided in an inertial frame*, but instead a frame centered on, and fixed to, the central body. It defines the velocity vector of a body w.r.t. 'the ground' or (alternatively and identically) the relative atmosphere in the case the atmosphere would be perfectly co-rotating with the central body.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the groundspeed is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "inertial_to_body_fixed_rotation_frame",
           &tp::inertialToBodyFixedRotationMatrixVariable,
           py::arg( "body" ),
           R"doc(

 Function to add the rotation matrix from inertial to body-fixed frame to the dependent variables to save.

 Function to add the rotation matrix from inertial to body-fixed frame to the dependent variables to save. This requires the rotation of the body to be defined (either in the environment or the state vector). 
 
 .. note:: 
 
     A rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`).

     You can convert this vector to a matrix using numpy's `reshape` function, e.g.:
     
     .. code-block:: python

        import numpy as np

        # Index where the rotation matrix starts in the dependent variable history
        start_index = ...
        rotation_matrix_components = dependent_variable_array[:, start_index : start_index + 9]

        rotation_matrix_history = rotation_matrix_components.reshape(
            (rotation_matrix_components.shape[0], 3, 3)
        )

 Parameters
 ----------
 body : str
     Body for which the rotation matrix is to be saved.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 9





     )doc" );

    m.def( "tnw_to_inertial_rotation_matrix",
           &tp::tnwToInertialFrameRotationMatrixVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the rotation matrix from the TNW to the inertial frame to the dependent variables to save.

 Function to add the rotation matrix from the TNW to the inertial frame to the dependent variables to save. It has the x-axis pointing along the velocity vector, the z-axis along the orbital angular momentum vector, and the y-axis completing the right-handed system.
 
 .. note:: 
 
     A rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`).

     You can convert this vector to a matrix using numpy's `reshape` function, e.g.:
     
     .. code-block:: python

        import numpy as np

        # Index where the rotation matrix starts in the dependent variable history
        start_index = ...
        rotation_matrix_components = dependent_variable_array[:, start_index : start_index + 9]

        rotation_matrix_history = rotation_matrix_components.reshape(
            (rotation_matrix_components.shape[0], 3, 3)
        )

 Parameters
 ----------
 body : str
     Body for which the rotation matrix is to be saved.
 central_body : str
     Body with respect to which the TNW frame is determined.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 9





     )doc" );

    m.def( "rsw_to_inertial_rotation_matrix",
           &tp::rswToInertialFrameRotationMatrixVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the rotation matrix from the RSW to the inertial frame to the dependent variables to save.

 Function to add the rotation matrix from the RSW to the inertial frame to the dependent variables to save. It has the x-axis pointing along the position vector (away from the central body), the z-axis along the orbital angular momentum vector, and the y-axis completing the right-handed system.

 .. note:: 
 
     A rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`).

     You can convert this vector to a matrix using numpy's `reshape` function, e.g.:
     
     .. code-block:: python

        import numpy as np

        # Index where the rotation matrix starts in the dependent variable history
        start_index = ...
        rotation_matrix_components = dependent_variable_array[:, start_index : start_index + 9]

        rotation_matrix_history = rotation_matrix_components.reshape(
            (rotation_matrix_components.shape[0], 3, 3)
        )

 Parameters
 ----------
 body : str
     Body for which the rotation matrix is to be saved.
 central_body : str
     Body with respect to which the TNW frame is determined.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 9





     )doc" );

    m.def( "inertial_to_body_fixed_313_euler_angles",
           &tp::eulerAnglesToBodyFixed313Variable,
           py::arg( "body" ),
           R"doc(

 Function to add the 3-1-3 Euler angles for the rotation from inertial to body-fixed frame to the dependent variables to save.

 Function to add the 3-1-3 Euler angles for the rotation from inertial to body-fixed frame to the dependent variables to save. This requires the rotation of the body to be defined (either in the environment or the state vector).

 Parameters
 ----------
 body : str
     Body for which the rotation angles are to be saved.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "intermediate_aerodynamic_rotation_matrix_variable",
           &tp::intermediateAerodynamicRotationMatrixVariable,
           py::arg( "body" ),
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           py::arg( "central_body" ) = "",
           R"doc(

 Function to add the rotation matrix between any two reference frames used in aerodynamic calculations.

 Function to add the rotation matrix between any two reference frames used in aerodynamic calculations. The list of available frames is defined by the :class:`AerodynamicsReferenceFrames` enum.

 .. note:: 
 
     A rotation matrix is returned as a nine-entry vector in the dependent variable output, where entry :math:`(i,j)` of the matrix is stored in entry :math:`(3i+j)` of the vector (with :math:`i,j=0,1,2`).

     You can convert this vector to a matrix using numpy's `reshape` function, e.g.:
     
     .. code-block:: python

        import numpy as np

        # Index where the rotation matrix starts in the dependent variable history
        start_index = ...
        rotation_matrix_components = dependent_variable_array[:, start_index : start_index + 9]

        rotation_matrix_history = rotation_matrix_components.reshape(
            (rotation_matrix_components.shape[0], 3, 3)
        )

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 base_frame : AerodynamicsReferenceFrames
     Base reference frame for the rotation.
 target_frame : AerodynamicsReferenceFrames
     Target reference frame for the rotation.
 central_body : str
     Central body w.r.t. which the state of the body is considered.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 9





     )doc" );

    m.def( "periapsis_altitude",
           &tp::periapsisAltitudeVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the altitude of periapsis to the dependent variables to save.

 Function to add the periapsis altitude of the current osculating orbit to the dependent variables to save. The altitude depends on the shape of the central body. This function takes the current (osculating) orbit of the body w.r.t. the central body, and uses this Kepler orbit to extract the position/altitude of periapsis.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the altitude of periapsis is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "apoapsis_altitude",
           &tp::apoapsisAltitudeVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the altitude of apoapsis to the dependent variables to save.

 Function to add the apoapsis altitude of the current osculating orbit to the dependent variables to save. The altitude depends on the shape of the central body. This function takes the current (osculating) orbit of the body w.r.t. the central body, and uses this Kepler orbit to extract the position/altitude of apoapsis.

 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 central_body : str
     Body with respect to which the altitude of apoapsis is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "control_surface_deflection",
           &tp::controlSurfaceDeflectionDependentVariable,
           py::arg( "body" ),
           py::arg( "control_surface" ),
           R"doc(No documentation found.)doc" );

    m.def( "central_body_fixed_spherical_position",
           &tp::centralBodyFixedSphericalPositionVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the spherical, body-fixed position to the dependent variables to save.

 Function to add the spherical position to the dependent variables to save. The spherical position is return as the radius, latitude, longitude, defined in the body-fixed frame of the central body

 Parameters
 ----------
 body : str
     Body whose spherical position is to be saved.
 central_body : str
     Body with respect to which the spherical, body-fixed is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "central_body_fixed_cartesian_position",
           &tp::centralBodyFixedCartesianPositionVariable,
           py::arg( "body" ),
           py::arg( "central_body" ),
           R"doc(

 Function to add the relative Cartesian position, in the central body's fixed frame, to the dependent variables to save.


 Parameters
 ----------
 body : str
     Body whose relative cartesian position is to be saved.
 central_body : str
     Body with respect to which the cartesian, body-fixed is computed.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "body_mass",
           &tp::bodyMassVariable,
           py::arg( "body" ),
           R"doc(

 Function to add the current body mass to the dependent variables to save.


 Parameters
 ----------
 body : str
     Body whose mass should be saved.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "radiation_pressure_coefficient",
           &tp::radiationPressureCoefficientVariable,
           py::arg( "body" ),
           py::arg( "emitting_body" ),
           R"doc(

 Function to add the current radiation pressure coefficient to the dependent variables to save.


 Parameters
 ----------
 body : str
     Body whose dependent variable should be saved.
 emitting_body : str
     Emitting body.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "total_mass_rate",
           &tp::totalMassRateDependentVariable,
           py::arg( "body" ),
           R"doc(

 Function to add the total mass rate to the dependent variables to save.

 Function to add the total mass rate to the dependent variables to save. It requires the body mass to be numerically propagated.

 Parameters
 ----------
 body : str
     Body whose mass rate should be saved.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "custom_dependent_variable",
           &tp::customDependentVariable,
           py::arg( "custom_function" ),
           py::arg( "variable_size" ),
           R"doc(

 Function to compute a custom dependent variable.

 Function to compute a custom dependent variable, which can be implemented by the user as a Python function. The custom dependent variable is typically dependent on the current properties of the environment (e.g. bodies in the environment) or a user-defined guidance class (or similar)

 Parameters
 ----------
 custom_function : Callable[[], numpy.ndarray].
     Function taking no input, and returning the custom dependent variable (as a numpy Nx1 array).
 variable_size : int
     Size N of the array returned by the ``custom_function``
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 ``variable_size``





     )doc" );

    m.def( "custom", &tp::customDependentVariableDeprecated, py::arg( "custom_function" ), py::arg( "variable_size" ) );

    m.def( "gravity_field_potential",
           &tp::gravityFieldPotentialDependentVariable,
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           R"doc(

 Function to add the gravitational potential to the dependent variables to save.

 Function to add the gravitational potential to the dependent variables to save. The gravitational potential is defined by the bodies undergoing and exerting the acceleration.

 Parameters
 ----------
 body_undergoing_acceleration : str
     Body whose dependent variable should be saved.
 body_exerting_acceleration : str
     Body exerting acceleration.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "gravity_field_laplacian_of_potential",
           &tp::gravityFieldLaplacianOfPotentialDependentVariable,
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           R"doc(

 Function to add the laplacian of the gravitational potential to the dependent variables to save.

 Function to add the laplacian of the gravitational potential to the dependent variables to save. The laplacian is defined by the bodies undergoing and exerting the acceleration.

 Parameters
 ----------
 body_undergoing_acceleration : str
     Body whose dependent variable should be saved.
 body_exerting_acceleration : str
     Body exerting acceleration.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "minimum_body_distance",
           &tp::minimumConstellationDistanceDependentVariableSaveSettings,
           py::arg( "body_name" ),
           py::arg( "bodies_to_check" ),
           R"doc(

 Function to compute the minimum distance between a given body, and a set of other bodies.

 Function to compute the minimum distance between a given body, and a set of other bodies. This function takes the instantaneous position of body ``body_name``, and each body in the list ``bodies_to_check``, and computes the body from this list closest to ``body_name``. In this calculation, the positions of the bodies are evaluated at the current propagation time, and therefore **light time is ignored**. In addition, this functions does not consider visibility requirements (e.g. is a planet between two bodies). The dependent variable is of size 2, and consists of: (0) The distance between the body, and the closest other body; (1) The index from ``bodies_to_check`` for which the distance (given by the first index) is closest to ``body`` Typically, this function is used to compute the closest body in a constellation of satellites.

 Parameters
 ----------
 body_name : str
     Body for which the distance to other bodies is to be computed.
 bodies_to_check : list[ str ]
     List of bodies for which it is to be checked which of these bodies is closest to ``body_name``.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 2





     )doc" );

    m.def( "minimum_visible_station_body_distances",
           &tp::minimumConstellationStationDistanceDependentVariableSaveSettings,
           py::arg( "body_name" ),
           py::arg( "station_name" ),
           py::arg( "bodies_to_check" ),
           py::arg( "minimum_elevation_angle" ),
           R"doc(

 Function to compute the minimum distance between a ground station, and a set of other bodies visible from that station.

 Function to compute the minimum distance between a ground station, and a set of other bodies visible from that station This function takes the instantaneous position of the ground station ``station_name`` on ``body_name``, and each body in the list ``bodies_to_check``, and computes the body from this list closest to this ground station, only taking into account those bodies from this list which are visible from teh ground station. For this function, visibility is defined by a single elevation angle cutoff (at the ground station) below which a body is deemed to not be visible. In this calculation, the positions of the bodies are evaluated at the current propagation time, and therefore **light time is ignored**. The dependent variable is of size 3, and consists of: (0) The distance between the ground station, and the closest visible body; (1) The index from ``bodies_to_check`` for which the distance (given by the first index) is closest to thee ground station, and the body is visible. (2) Elevation angle for closest body. In case, no body is visible from the station, this function returns [NaN, -1, NaN]. Typically, this function is used to compute the closest body between a ground station and a constellation of satellites.

 Parameters
 ----------
 body_name : str
     Body on which ground station is located, for which the distance to other bodies is to be computed.
 station_name : str
     Name of ground station, for which the distance to other bodies is to be computed.
 bodies_to_check : list[ str ]
     List of bodies for which it is to be checked which of these bodies is closest to ``station_name`` on ``body_name``.
 minimum_elevation_angle : float
     Minimum elevation angle (at ground station) below which the distance to the ``bodies_to_check`` is not considered.
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 3





     )doc" );

    m.def( "center_of_mass", &tp::centerOfMassVariableSaveSettings, py::arg( "body" ), R"doc(No documentation found.)doc" );

    m.def( "inertia_tensor", &tp::inertiaTensorVariableSaveSettings, py::arg( "body" ), R"doc(No documentation found.)doc" );

    m.def( "received_irradiance",
           &tp::receivedIrradianceDependentVariable,
           py::arg( "target_body" ),
           py::arg( "source_body" ),
           R"doc(

 Function to save the received irradiance from a give source.

 Function to save the received irradiance (in W/m^{2}) from a given source, as used in the computation of radiation pressure

 Parameters
 ----------
 target_body : str
     Name of body at the location of which the irradiance is to be returned
 source_body : str
     Name of body from which the radiation originates that causes the irradiance at the target body
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "received_irradiance_shadow_function",
           &tp::receivedFractionDependentVariable,
           py::arg( "target_body" ),
           py::arg( "source_body" ),
           R"doc(

 Function to save the shadow function that reduces the received irradiance from a given source.

 Function to save the shadow function that reduces the received irradiance from a given source, as a result of occultation of radiation between the source and the target. The shadow function is by definition between 0 (no irradiance is received, total occultation) and 1 (all irradiance is received, no occultation)

 Parameters
 ----------
 target_body : str
     Name of body at the location of which the irradiance is to be returned
 source_body : str
     Name of body from which the radiation originates that causes the irradiance at the target body
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "visible_radiation_source_area",
           &tp::visibleSourceAreaDependentVariable,
           py::arg( "target_body" ),
           py::arg( "source_body" ),
           R"doc(

 Function to save the shadow function that reduces the received irradiance from a given source.

 Function to save the shadow function that reduces the received irradiance from a given source, as a result of occultation of radiation between the source and the target. The shadow function is by definition between 0 (no irradiance is received, total occultation) and 1 (all irradiance is received, no occultation)

 Parameters
 ----------
 target_body : str
     Name of body at the location of which the irradiance is to be returned
 source_body : str
     Name of body from which the radiation originates that causes the irradiance at the target body
 Returns
 -------
 SingleDependentVariableSaveSettings
     Dependent variable settings object.
 Variable Size
 -------------
 1





     )doc" );

    m.def( "vehicle_panel_surface_normals_inertial_frame",
           &tp::vehiclePanelInertialSurfaceNormals,
           py::arg( "body_name" ),
           py::arg( "part_name" ) = "",
           R"doc(No documentation found.)doc" );

    m.def( "vehicle_panel_surface_normals_body_fixed_frame",
           &tp::vehiclePanelInertialSurfaceNormals,
           py::arg( "body_name" ),
           py::arg( "part_name" ) = "",
           R"doc(No documentation found.)doc" );

    m.def( "per_target_panel_radiation_pressure_force",
           &tp::vehiclePanelInertialSurfaceNormals,
           py::arg( "target_name" ),
           py::arg( "source_name" ),
           R"doc(No documentation found.)doc" );

    m.def( "radiation_pressure_source_panel_irradiance",
           &tp::paneledRadiationSourcePerPanelIrradiance,
           py::arg( "target_name" ),
           py::arg( "source_name" ),
           R"doc(No documentation found.)doc" );

    m.def( "radiation_pressure_source_panel_geometry",
           &tp::paneledRadiationSourceGeometry,
           py::arg( "target_name" ),
           py::arg( "source_name" ),
           R"doc(No documentation found.)doc" );

    m.def( "illuminated_panel_fraction",
           &tp::illuminatedPanelFractionDependentVariable,
           py::arg( "target_name" ),
           py::arg( "source_name" ),
           py::arg( "panel_type" ) = "",
           R"doc(No documentation found.)doc" );

    m.def( "cross_section_change",
           &tp::crossSectionChangeDependentVariable,
           py::arg( "target_name" ),
           py::arg( "source_name" ),
           py::arg( "acceleration_type" ) = "radiation_pressure",
           R"doc(No documentation found.)doc" );

    m.def( "full_body_paneled_geometry",
           &tp::fullBodyPaneledGeometryDependentVariable,
           py::arg( "target_name" ),
           R"doc(No documentation found.)doc" );

    m.def( "aerodynamic_coefficients",
           &tp::aerodynamicCoefficientsDependentVariable,
           py::arg( "target_name" ),
           py::arg( "central_body_name" ),
           R"doc(No documentation found.)doc" );

    m.def( "actual_cross_section",
           &tp::actualCrossSectionDependentVariable,
           py::arg( "target_name" ),
           py::arg( "central_body_name" ),
           py::arg( "acceleration_type" ) = "radiation_pressure",
           R"doc(No documentation found.)doc" );
}

}  // namespace dependent_variable
}  // namespace propagation_setup
}  // namespace dynamics
}  // namespace tudatpy
