/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_vehicle_systems.h"

#include <tudat/simulation/environment_setup.h>

// #include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
// #include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/math/basic.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace numerical_simulation
{
namespace environment_setup
{
namespace vehicle_systems
{

void expose_vehicle_systems_setup( py::module& m )
{
    py::class_< tss::BodyPanelGeometrySettings, std::shared_ptr< tss::BodyPanelGeometrySettings > >(
            m,
            "BodyPanelGeometrySettings",
            R"doc(

         Base class for defining the geometrical properties of a single panel on the vehicle's exterior.
         A derived class of this can be instantiated through the :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.frame_fixed_panel_geometry`, :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.time_varying_panel_geometry`, or :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.body_tracking_panel_geometry` functions.


      )doc" );


      py::class_< tss::FrameFixedBodyPanelGeometrySettings,
      tss::BodyPanelGeometrySettings,
      std::shared_ptr< tss::FrameFixedBodyPanelGeometrySettings > >( m,

      "FrameFixedBodyPanelGeometrySettings",
        R"doc(

Derived class for defining the geometrical properties of a single panel on the vehicle's exterior, with a fixed orientation in a given frame.
This class is typically instantiated through the :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.frame_fixed_panel_geometry` function.
)doc" )
  .def_readwrite( "surface_normal",
                  &tss::FrameFixedBodyPanelGeometrySettings::surfaceNormal_,
                  R"doc(
Panel outward surface normal vector (in specified frame)

:type: np.array
)doc" )
  .def_readwrite( "area",
                  &tss::FrameFixedBodyPanelGeometrySettings::area_,
                  R"doc(
Panel surface area

:type: float
)doc" );

py::class_< tss::FrameVariableBodyPanelGeometrySettings,
      tss::BodyPanelGeometrySettings,
      std::shared_ptr< tss::FrameVariableBodyPanelGeometrySettings > >( m,

    "FrameVariableBodyPanelGeometrySettings",
    R"doc(
Derived class for defining the geometrical properties of a single panel on the vehicle's exterior, with a time-variable orientation in a given frame.
This class is typically instantiated through the :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.time_varying_panel_geometry` or :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.body_tracking_panel_geometry` functions.
)doc" )
  .def_readwrite( "surface_normal_function",
                  &tss::FrameVariableBodyPanelGeometrySettings::surfaceNormalFunction_,
                  R"doc(
Function which takes the current epoch as input and returns the panel outward surface normal vector (in specified frame).

:type: Callable[[float], np.ndarray]
)doc" )
  .def_readwrite( "area",
                  &tss::FrameVariableBodyPanelGeometrySettings::area_,
                  R"doc(
Panel surface area

:type: float
)doc" );

    m.def( "frame_fixed_panel_geometry",
           py::overload_cast< const Eigen::Vector3d&, const double, const std::string& >(
                   tss::frameFixedPanelGeometry ),
           py::arg( "surface_normal" ),
           py::arg( "area" ),
           py::arg( "frame_orientation" ) = "",
           R"doc(

 Function for creating settings for a vehicle exterior panel that is fixed to a given frame.

 Function for creating settings for a vehicle exterior panel that is fixed to a given frame, meaning
 that the orientation of the panel is fully defined by the rotation model(s) defined in the vehicle.
 The constant surface normal :math:`\hat{\mathbf{n}}^{\mathcal{F}}` in frame :math:`\mathcal{F}` is provided by the user.
 If the ``frame_orientation`` of this function is left empty, the panel is fixed to the body-frame, and
 :math:`\mathcal{F}` is the  body-fixed frame :math:`\mathcal{B}`.

 Alternatively, the ``frame_orientation`` may be defined as the identifier of the frame fixed to one of the
 vehicle parts (solar array, antenna, etc.). See :func:`~full_panelled_body_settings` for the definition
 of rotation models of vehicle parts.

 Note that this panel model does not contain information on panel location or shape, only its area and surface normal,
 and is therefore not suitable for computation of panel shadowing of torque computations.


 Parameters
 ----------
 surface_normal : np.array
     Panel outward surface normal vector (in specified frame)
 area : float
     Panel surface area
 frame_orientation : str, default = ""
     Identifier of the frame to which the panel is fixed (if body-fixed frame, this can be left empty)
 Returns
 -------
 BodyPanelGeometrySettings
     Object defining settings for panel geometry






     )doc" );

    m.def( "body_tracking_panel_geometry",
           py::overload_cast< const std::string&, const bool, const double, const std::string& >(
                   &tss::bodyTrackingPanelGeometry ),
           py::arg( "body_to_track" ),
           py::arg( "towards_tracked_body" ),
           py::arg( "area" ),
           py::arg( "frame_orientation" ) = "",
           R"doc(

 Function for creating settings for a vehicle exterior panel where the surface normal tracks a given body.

 Function for creating settings for a vehicle exterior panel where the surface normal tracks a given body, for instance
 to define the surface normal of a solar array to always point towards the Sun, or an antenna to always point towards the Earth.
 When using this option, the panel surface normal :math:`\hat{\mathbf{n}}` is computed in an inertial frame based on the tracked
 body, and then (if necessary) rotated to the body-fixed frame.
 Note that this panel model does not contain information on panel location or shape, only its area and surface normal,
 and is therefore not suitable for computation of panel shadowing of torque computations.


 Parameters
 ----------
 body_to_track : str
     Name of the body towards (or away from) which the panel surface normal is to point
 towards_tracked_body : bool
     Boolean defining whether the normal vector points towards (if true) or away from (if false) the tracked body
 area : float
     Panel surface area
 frame_orientation : str, default = ""
     Identifier of the frame in which the panel is defined (with time-variable orientation, defined by tracked body). Note that this option is typically only relevant for internal  book-keeping, and can be left empty
 Returns
 -------
 BodyPanelGeometrySettings
     Object defining settings for panel geometry






     )doc" );

    m.def( "time_varying_panel_geometry",
           py::overload_cast< const std::function< Eigen::Vector3d( ) >&,
                              const double,
                              const std::string& >( &tss::timeVaryingPanelGeometry ),
           py::arg( "surface_normal_function" ),
           py::arg( "area" ),
           py::arg( "frame_orientation" ),
           R"doc(

 Function for creating settings for a vehicle exterior panel that has time-variable orientation in a given frame.

 As :func:`~frame_fixed_panel_geometry`, but with a time-variable outward surface normal :math:`\hat{\mathbf{n}}^{\mathcal{F}}(t)`


 Parameters
 ----------
 surface_normal_function : Callable[[], np.ndarray]
    Function which takes the current epoch as input and returns the panel outward surface normal vector (in specified frame).
 area : float
     Panel surface area
 frame_orientation : str, default = ""
     Identifier of the frame in which the panel normal is defined
 Returns
 -------
 BodyPanelGeometrySettings
     Object defining settings for panel geometry






     )doc" );

    py::class_< tss::BodyPanelSettings, std::shared_ptr< tss::BodyPanelSettings > >(
            m,
            "BodyPanelSettings",
            R"doc(

         Class for defining the complete properties of a single panel on the vehicle's exterior.
         This class is typically instantiated through the :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.body_panel_settings` function.





      )doc" )
      .def_readwrite( "panel_geometry", &tss::BodyPanelSettings::panelGeometry_, R"doc(
        
        Geometric properties of the panel, including surface normal vector and area.

        :type: BodyPanelGeometrySettings

        )doc" )
            .def_readwrite( "reflection_law_settings",
                            &tss::BodyPanelSettings::reflectionLawSettings_,
                            R"doc(
        Reflection law settings of the panel.
            
        :type: BodyPanelReflectionLawSettings
            )doc"

                            )
            .def_readwrite( "panel_type_id", &tss::BodyPanelSettings::panelTypeId_, R"doc(
        Optional identifier for panel type.
        This is typically used to identify the type of panel and can be used to assign a rotation model to a specific panel type, see the :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.full_panelled_body_settings` function.

        :type: str
        )doc" );

    m.def( "body_panel_settings",
           py::overload_cast<
                std::shared_ptr< tss::BodyPanelGeometrySettings >,
                std::shared_ptr< tss::BodyPanelReflectionLawSettings >,
                std::string,
                std::shared_ptr< tss::MaterialProperties > >( &tss::bodyPanelSettings ),
           py::arg( "panel_geometry" ),
           py::arg( "panel_reflection_law" ),
           py::arg( "panel_type_id" ) = "",
           py::arg( "material_properties" ) = nullptr,
           R"doc(

 Function for creating settings for a full panel

 Function for creating settings for a full panel (presently only geometry and reflection properties). The :class:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.BodyPanelReflectionLawSettings` can be created using the :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.specular_diffuse_body_panel_reflection` or :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.lambertian_body_panel_reflection` functions.
 The panel
 can also be endowed with an identifier to specify the type of the panel. This has no direct consequences for the model,
 but may be useful in estimation, to for instance estimate the reflection properties of all panels specified with identified "MLI"
 as a single parameter


 Parameters
 ----------
 panel_geometry : BodyPanelGeometrySettings
     Geometric properties of the panel (size and orientation, at least)
 panel_reflection_law : BodyPanelReflectionLawSettings
     Reflection law settings of the panel
 panel_type_id : str, default = ""
     Optional identifier for panel type
 Returns
 -------
 BodyPanelSettings
     Object defining settings for a panel






     )doc" );

    py::class_< tss::FullPanelledBodySettings, std::shared_ptr< tss::FullPanelledBodySettings > >(
            m,
            "FullPanelledBodySettings",
            R"doc(

         Class for providing the complete settings for a panelled body exterior.
        
         This is typically defined through the :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.full_panelled_body_settings` or :func:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.box_wing_panelled_body_settings` functions.
         The class contains a list of panel settings, and (optionally) a list of rotation model settings for vehicle parts.


     )doc" )
            .def_readwrite( "panel_settings_list", &tss::FullPanelledBodySettings::panelSettingsList_, R"doc(
        List of individual body panel settings on the body.
        
        :type: list[BodyPanelSettings]
        )doc" )
            .def_readwrite( "part_rotation_model_settings",
                            &tss::FullPanelledBodySettings::partRotationModelSettings_,
                            R"doc(
        Dictionary of rotation model settings per vehicle parts.
        The rotation model settings are defined per `panel_type_id`, as defined in the :class:`~tudatpy.numerical_simulation.environment_setup.vehicle_systems.BodyPanelSettings`.
    
        :type: dict[str,RotationModelSettings]

        )doc"

            );

    m.def( "full_panelled_body_settings",
           &tss::fullPanelledBodySettings,
           py::arg( "panel_settings" ),
           py::arg( "part_rotation_model_settings" ) =
                   std::map< std::string, std::shared_ptr< tss::RotationModelSettings > >( ),
           R"doc(

 Function for creating settings for a full panelled vehicle exterior.

 Function for creating settings for a full panelled vehicle exterior, taking a list of panel settings,
 and (optionally) a list of rotation model settings for vehicle parts. The identifiers for the rotation models
 are used to specify the names of part-fixed frames, which are used by the ``frame_orientation`` inputs to
 functions creating settings for :class:`~BodyPanelGeometrySettings`. For instance, assigning a rotation model
 to frame ``LRO_SolarArray`` (dict key for ``part_rotation_model_settings``) allows panels defined in the frame
 with this same frame orientation to be defined. The associated rotation model defines rotations from body-fixed
 frame :math:`\mathcal{B}` to part-fixed frame :math:`\mathcal{F}_{j}` (for part :math:`j`). The rotation from part-fixed
 (where the surface normal is defined) to inertial frame is then computed from
 :math:`\mathbf{R}^{I/\mathcal{F}_{j}}=\mathbf{R}^{I/\mathcal{B}}\mathbf{R}^{\mathcal{B}/\mathcal{F}_{j}}`, where :math:`\mathbf{R}^{I/\mathcal{B}}`
 defines the body's orientation, and :math:`\mathbf{R}^{\mathcal{B}/\mathcal{F}_{j}}` the part orientation (w.r.t. a body-fixed frame)


 Parameters
 ----------
 panel_settings : list[BodyPanelSettings]
     List of settings for body panels.
 part_rotation_model_settings : dict[str,RotationModelSettings], default = dict()
     Rotation model settings per vehicle part (default empty, indicating no part-fixed frames are defined)
 Returns
 -------
 FullPanelledBodySettings
     Object defining full panelled vehicle exterior






     )doc" );

    m.def( "box_wing_panelled_body_settings",
           &tss::bodyWingPanelledGeometry,
           py::arg( "length" ),
           py::arg( "width" ),
           py::arg( "height" ),
           py::arg( "solar_array_area" ),
           py::arg( "box_specular_reflectivity" ),
           py::arg( "box_diffuse_reflectivity" ),
           py::arg( "solar_array_specular_reflectivity" ),
           py::arg( "solar_array_diffuse_reflectivity" ),
           py::arg( "box_instantaneous_reradiation " ) = true,
           py::arg( "solar_array_instantaneous_reradiation " ) = true,
           R"doc(

 Function for creating a simple box-wing spacecraft exterior shape with reflection law settings.

 This function creates a :func:`~full_panelled_body_settings` with ``panel_settings`` generated from simple box-wing
 settings. The assumptions behind the box-wing model are:

 * The spacecraft shape is defined by a rectangular box (cuboid) and solar array
 * The box has its faces parallel to the xy-, xz- and yz-planes
 * The solar array surface normal always points towards the Sun
 * Each box face has identical reflection law settings, defined by :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.specular_diffuse_body_panel_reflection` settings.
 * The solar array has reflection law settings, defined by :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.specular_diffuse_body_panel_reflection` settings.


 Parameters
 ----------
 length : float
     Box length (size in body-fixed x-direction).
 width : float
     Box width (size in body-fixed y-direction).
 height : float
     Box height (size in body-fixed z-direction).
 solar_array_area : float
     Surface area of the solar array.
 box_specular_reflectivity : float
     Box secular reflectivity :math:`\rho`.
 box_diffuse_reflectivity : float
     Box secular reflectivity :math:`\delta`.
 solar_array_specular_reflectivity : float
     Solar array secular reflectivity :math:`\rho`.
 solar_array_diffuse_reflectivity : float
     Solar array secular reflectivity :math:`\delta`.
 box_instantaneous_reradiation : bool
     Boolean denoting whether absorbed radiation is instantaneously retransmitted from box (yes, if true).
 solar_array_instantaneous_reradiation : bool
     Boolean denoting whether absorbed radiation is instantaneously retransmitted from solar array (yes, if true).
 Returns
 -------
 FullPanelledBodySettings
     Object defining full panelled vehicle exterior






     )doc" );

    m.def( "body_panel_settings_list_from_dae",
        &tss::bodyPanelSettingsListFromDae,
        py::arg( "file_path" ),
        py::arg( "frame_origin" ),
        py::arg( "material_properties" ),
        py::arg( "reradiation_settings" ),
        py::arg( "frame_orientation" ) = "",
        R"doc(
Function for creating list of panel body settings
        
Function for creating list of panel body settings from a .dae (COLLADA) file containing the 3D geometry and the
material for the paneled surface.
Parameters
----------
file_path : str
    Path to .dae file with geometry data.
frame_origin : np.array
    Frame origin of the .dae part to be loaded.
material_properties : dict[str, MaterialProperties]
    Dictionary of material properties, as they appear in the .dae file provided.
reradiation_settings : dict[str, bool]
    Dictionary of re-radiation settings for materials, as they appear in the .dae file provided.
frame_orientation : str, default = "
    Identifier of the frame to which the panel is fixed (if body-fixed frame, this can be left empty).
    
Returns
-------
list[BodyPanelSettings]
    List of settings for body panels
    )doc" );

    m.def( "merge_body_panel_setting_lists",
       &tss::mergeBodyPanelSettingsLists,
       py::arg( "list_of_lists_of_body_panel_settings" ),
       R"doc(
Function for merging lists of panel body settings
Parameters
----------
list_of_lists_of_body_panel_settings : list[list[BodyPanelSettings], ...]
    List of lists of body panel settings.
Returns
-------
list[BodyPanelSettings]
    List of settings for body panels assembled from different parts, creating a coherent list of body panel settings.
    )doc" );

    py::class_< tss::MaterialProperties, std::shared_ptr< tss::MaterialProperties > >(m, "MaterialProperties")
        .def_readwrite("specular_reflectivity", &tss::MaterialProperties::specularReflectivity_)
        .def_readwrite("diffuse_reflectivity", &tss::MaterialProperties::diffuseReflectivity_)
        .def_readwrite("energy_accomodation_coefficient", &tss::MaterialProperties::energyAccomodationCoefficient_)
        .def_readwrite("normal_accomodation_coefficient", &tss::MaterialProperties::normalAccomodationCoefficient_)
        .def_readwrite("tangential_accomodation_coefficient", &tss::MaterialProperties::tangentialAccomodationCoefficient_)
        .def_readwrite("normal_velocity_at_wall_ratio", &tss::MaterialProperties::normalVelocityAtWallRatio_);

    m.def( "material_properties",
       &tss::materialProperties,
       py::arg( "specular_reflectivity" ) = -1,
       py::arg( "diffuse_reflectivity" ) = -1,
       py::arg( "energy_accomodation_coefficient" ) = -1,
       py::arg( "normal_accomodation_coefficient" ) = -1,
       py::arg( "tangential_accomodation_coefficient" ) = -1,
       py::arg( "normal_velocity_at_wall_ratio" ) = -1,
       R"doc(No Documentation)doc" );
    
}

}  // namespace vehicle_systems
}  // namespace environment_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
