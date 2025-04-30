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
#include "expose_rigid_body.h"

#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/simulation/environment_setup.h>

// #include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tpc = tudat::physical_constants;

namespace tudatpy
{
namespace numerical_simulation
{
namespace environment_setup
{
namespace rigid_body
{

void expose_rigid_body_setup( py::module &m )
{
    //    py::enum_<tss::RigidBodyPropertiesType>(m,
    //    "RigidBodyPropertiesType",
    //                                     get_docstring("RigidBodyPropertiesType").c_str())
    //        .value("from_function_rigid_body_properties",
    //        tss::RigidBodyPropertiesType::from_function_rigid_body_properties,
    //        get_docstring("RigidBodyPropertiesType.from_function_rigid_body_properties").c_str())
    //        .value("constant_rigid_body_properties",
    //        tss::RigidBodyPropertiesType::constant_rigid_body_properties,
    //        get_docstring("RigidBodyPropertiesType.constant_rigid_body_properties").c_str())
    //        .value("from_gravity_field_rigid_body_properties",
    //        tss::RigidBodyPropertiesType::from_gravity_field_rigid_body_properties,
    //        get_docstring("RigidBodyPropertiesType.from_gravity_field_rigid_body_properties").c_str())
    //        .value("mass_dependent_rigid_body_properties",
    //        tss::RigidBodyPropertiesType::mass_dependent_rigid_body_properties,
    //        get_docstring("RigidBodyPropertiesType.mass_dependent_mass_distribution_properties").c_str())
    //        .export_values();
    //
    py::class_< tss::RigidBodyPropertiesSettings,
                std::shared_ptr< tss::RigidBodyPropertiesSettings > >(
            m,
            "RigidBodyPropertiesSettings",
            R"doc(

         Base class for providing settings for rigid body model creation.

         This class is a functional base class for settings of gravity field models that require no information in addition to their type.
         Gravity field model classes requiring additional information must be created using an object derived from this class.





      )doc" )
            .def_property_readonly( "body_mass_property_type",
                                    &tss::RigidBodyPropertiesSettings::getRigidBodyPropertiesType,
                                    R"doc(

         **read-only**

         Type of rigid body model that is to be created.

         :type: RigidBodyPropertiesType
      )doc" );

    m.def( "constant_rigid_body_properties",
           tss::constantRigidBodyPropertiesSettings,
           py::arg( "mass" ),
           py::arg( "center_of_mass" ) = Eigen::Vector3d::Constant( TUDAT_NAN ),
           py::arg( "inertia_tensor" ) = Eigen::Matrix3d::Constant( TUDAT_NAN ),
           R"doc(

 Function for creating constant rigid body properties.

 Function for creating constant rigid body properties (mass, center of mass, inertia tensor). The center of mass and/or inertia tensor can be left empty by setting them
 to NaN (default), in which case no center of mass or inertia tensor are defined


 Parameters
 ----------
 mass : float
     Constant mass of the body
 center_of_mass : np.array, default = np.full([3, 1], np.nan)
     Constant center of mass of the body (in a body-fixed frame)
 inertia_tensor : np.array, default = np.full([3, 3], np.nan)
     Constant inertia tensor of the body (in a body-fixed frame)
 Returns
 -------
 RigidBodyPropertiesSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings






     )doc" );

    m.def( "custom_time_dependent_rigid_body_properties",
           tss::fromFunctionRigidBodyPropertiesSettings,
           py::arg( "mass_function" ),
           py::arg( "center_of_mass_function" ) = nullptr,
           py::arg( "inertia_tensor_function" ) = nullptr,
           R"doc(

 Function for creating custom (time-dependent) rigid body properties.

 Function for creating custom rigid body properties, where the mass, center of mass and inertia tensor are defined by user-defined functions (as a function of time).
 The center of mass and/or inertia tensor functions can be left empty by setting them
 to None (default), in which case no center of mass or inertia tensor are defined


 Parameters
 ----------
 mass_function : Callable[[float], float]
     Function returning the mass as a function of time (to be used during the propagation)
 center_of_mass_function : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]] = None
     Function returning the center of mass as a function of time (to be used during the propagation)
 inertia_tensor_function : Callable[[float], numpy.ndarray[numpy.float64[3, 3]]] = None
     Function returning the inertia tensor as a function of time (to be used during the propagation)
 Returns
 -------
 RigidBodyPropertiesSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings






     )doc" );

    m.def( "custom_mass_dependent_rigid_body_properties",
           tss::massDependentMassDistributionSettings,
           py::arg( "mass" ),
           py::arg( "center_of_mass_function" ) = nullptr,
           py::arg( "inertia_tensor_function" ) = nullptr,
           R"doc(

 Function for creating custom (time-dependent) rigid body properties.

 Function for creating custom rigid body properties, center of mass and inertia tensor are defined by user-defined functions as a function of mass.
 This functionality is typically used for a body under thrust, where the center of mass and inertia tensor are defined as a function of expended mass.


 Parameters
 ----------
 mass : Callable[[float], float]
     Mass of the body (to be overridden during propagation if mass is propagated)
 center_of_mass_function : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]] = None
     Function returning the center of mass as a function of mass (to be used during the propagation)
 inertia_tensor_function : Callable[[float], numpy.ndarray[numpy.float64[3, 3]]] = None
     Function returning the inertia tensor as a function of mass (to be used during the propagation)
 Returns
 -------
 RigidBodyPropertiesSettings
     Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings






     )doc" );
}

}  // namespace rigid_body
}  // namespace environment_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
