/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_rigid_body.h"

#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/simulation/environment_setup/rigidBodyProperties.h>
#include <tudat/simulation/environment_setup/createGravityField.h>

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
namespace dynamics
{
namespace environment_setup
{
namespace rigid_body
{

void expose_rigid_body_setup( py::module& m )
{
    py::enum_< tss::RigidBodyPropertiesType >( m, "RigidBodyPropertiesType" )
            .value( "from_function_rigid_body_properties", tss::RigidBodyPropertiesType::from_function_rigid_body_properties )
            .value( "constant_rigid_body_properties", tss::RigidBodyPropertiesType::constant_rigid_body_properties )
            .value( "from_gravity_field_rigid_body_properties", tss::RigidBodyPropertiesType::from_gravity_field_rigid_body_properties )
            .value( "mass_dependent_rigid_body_properties", tss::RigidBodyPropertiesType::mass_dependent_rigid_body_properties );

    py::class_< tss::RigidBodyPropertiesSettings, std::shared_ptr< tss::RigidBodyPropertiesSettings > >( m,
                                                                                                         "RigidBodyPropertiesSettings",
                                                                                                         R"doc(

         Base class for providing settings for rigid body model creation.

         Derived settings select how mass, center of mass, and inertia are defined for a body.





      )doc" )
            .def_property_readonly( "body_mass_property_type",
                                    &tss::RigidBodyPropertiesSettings::getRigidBodyPropertiesType,
                                    R"doc(

         **read-only**

         Type of rigid body model that is to be created.

         :type: RigidBodyPropertiesType
      )doc" );

    py::class_< tss::FromGravityFieldRigidBodyPropertiesSettings,
                std::shared_ptr< tss::FromGravityFieldRigidBodyPropertiesSettings >,
                tss::RigidBodyPropertiesSettings >(
            m, "FromGravityFieldRigidBodyPropertiesSettings", R"doc(Settings for properties derived from the body's gravity field.)doc" )
            .def_property( "scaled_mean_moment_of_inertia",
                           &tss::FromGravityFieldRigidBodyPropertiesSettings::getScaledMeanMomentOfInertia,
                           &tss::FromGravityFieldRigidBodyPropertiesSettings::setScaledMeanMomentOfInertia,
                           R"doc(Mean principal moment divided by mass times squared gravity reference radius.)doc" );

    m.def( "from_gravity_field",
           tss::fromGravityFieldRigidBodyPropertiesSettings,
           py::arg( "scaled_mean_moment_of_inertia" ) = TUDAT_NAN,
           R"doc(

Create rigid-body properties derived from the body's gravity field.

The gravitational parameter defines mass. For a spherical-harmonic field, degree-one coefficients define the center
of mass, while complete degree-two coefficients and ``scaled_mean_moment_of_inertia`` define the inertia tensor. A
polyhedron field derives its inertia tensor from its homogeneous geometry. Leaving the scaled mean moment unset keeps
a spherical-harmonic gravity field valid without defining an inertia tensor.

Parameters
----------
scaled_mean_moment_of_inertia : float, default = nan
    Mean principal moment divided by :math:`MR^2`, equivalently
    :math:`(I_{xx}+I_{yy}+I_{zz})/(3MR^2)`.

Returns
-------
FromGravityFieldRigidBodyPropertiesSettings
    Canonical settings for gravity-derived mass, center of mass, and inertia.

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
     Instance of the :class:`~tudatpy.dynamics.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings






     )doc" );

    m.def( "custom_time_dependent_rigid_body_properties",
           tss::fromFunctionRigidBodyPropertiesSettings,
           py::arg( "mass_function" ),
           py::arg_v( "center_of_mass_function", std::function< Eigen::Vector3d( const double ) >( ), "None" ),
           py::arg_v( "inertia_tensor_function", std::function< Eigen::Matrix3d( const double ) >( ), "None" ),
           R"doc(

 Function for creating custom (time-dependent) rigid body properties.

 Function for creating custom rigid body properties, where the mass, center of mass and inertia tensor are defined by user-defined functions (as a function of time).
 The center of mass and/or inertia tensor functions can be left empty by setting them
 to None (default), in which case no center of mass or inertia tensor are defined


 Parameters
 ----------
 mass_function : callable[[float], float]
     Function returning the mass as a function of time in seconds since J2000 TDB, to be used during the propagation.
 center_of_mass_function : callable[[float], numpy.ndarray[numpy.float64[3, 1]]] = None
     Function returning the center of mass as a function of time in seconds since J2000 TDB, to be used during the propagation.
 inertia_tensor_function : callable[[float], numpy.ndarray[numpy.float64[3, 3]]] = None
     Function returning the inertia tensor as a function of time in seconds since J2000 TDB, to be used during the propagation.
 Returns
 -------
 RigidBodyPropertiesSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings






     )doc" );

    m.def( "custom_mass_dependent_rigid_body_properties",
           tss::massDependentMassDistributionSettings,
           py::arg( "mass" ),
           py::arg_v( "center_of_mass_function", std::function< Eigen::Vector3d( const double ) >( ), "None" ),
           py::arg_v( "inertia_tensor_function", std::function< Eigen::Matrix3d( const double ) >( ), "None" ),
           R"doc(

 Function for creating custom (time-dependent) rigid body properties.

 Function for creating custom rigid body properties, center of mass and inertia tensor are defined by user-defined functions as a function of mass.
 This functionality is typically used for a body under thrust, where the center of mass and inertia tensor are defined as a function of expended mass.


 Parameters
 ----------
 mass : float
     Mass of the body (to be overridden during propagation if mass is propagated)
 center_of_mass_function : callable[[float], numpy.ndarray[numpy.float64[3, 1]]] = None
     Function returning the center of mass as a function of mass (to be used during the propagation)
 inertia_tensor_function : callable[[float], numpy.ndarray[numpy.float64[3, 3]]] = None
     Function returning the inertia tensor as a function of mass (to be used during the propagation)
 Returns
 -------
 RigidBodyPropertiesSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.rigid_body.RigidBodyPropertiesSettings` object with the given settings






     )doc" );
}

}  // namespace rigid_body
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
