/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createRelativisticTimeConverter.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;

namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
{

using RelativisticTimeConverterSettings = tss::DirectRelativisticTimeConverterSettings< STATE_SCALAR_TYPE, TIME_TYPE >;
using RelativisticTimePropagatorSettings = tp::RelativisticTimeStatePropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE >;

std::shared_ptr< RelativisticTimeConverterSettings > directRelativisticTimeConverterSettings(
        const std::shared_ptr< RelativisticTimePropagatorSettings >& barycentricToBodycentricSettings,
        const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< TIME_TYPE > >& integratorSettings,
        const std::vector< std::shared_ptr< RelativisticTimePropagatorSettings > >& bodycentricToTopocentricSettings )
{
    return std::make_shared< RelativisticTimeConverterSettings >(
            barycentricToBodycentricSettings, integratorSettings, bodycentricToTopocentricSettings );
}

void setRelativisticTimeConverters( const tss::SystemOfBodies& bodies,
                                    const std::map< std::string, std::shared_ptr< RelativisticTimeConverterSettings > >& settings )
{
    tss::setRelativisticTimeConverters< STATE_SCALAR_TYPE, TIME_TYPE >( bodies, settings );
}

void expose_environment_setup_state_scalar( py::module& m )
{
    m.def( "create_system_of_bodies",
           &tss::createSystemOfBodies< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "body_settings" ),
           R"doc(

 Function that creates a System of bodies from associated settings.

 Function that creates a :class:`~tudatpy.dynamics.environment.SystemOfBodies` of bodies from associated settings in a :class:`~tudatpy.dynamics.environment_setup.BodyListSettings` object.
 This function creates the separate :class:`~tudatpy.dynamics.environment.Body`
 objects and stores them in a :class:`~tudatpy.dynamics.environment.SystemOfBodies` object. This :class:`~tudatpy.dynamics.environment.SystemOfBodies` object represents the full
 physical environment in the simulation, and this function is responsible for creating this environment from the user-defined settings.


 Parameters
 ----------
 body_settings : BodyListSettings
     Object defining the physical environment, with all properties of artificial and natural bodies.
 Returns
 -------
 :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Object containing the objects for bodies and environment models constituting the physical environment






     )doc" );

    m.def( "add_empty_tabulated_ephemeris",
           &tp::addEmptyTabulatedEphemeris< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "ephemeris_origin" ) = "",
           py::arg( "is_part_of_multi_arc" ) = false,
           R"doc(No documentation found.)doc" );

    m.def( "create_tabulated_ephemeris_from_spice",
           &tss::createTabulatedEphemerisFromSpice< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "body" ),
           py::arg( "initial_time" ),
           py::arg( "end_time" ),
           py::arg( "time_step" ),
           py::arg( "observer_name" ),
           py::arg( "reference_frame_name" ),
           py::arg_v( "interpolator_settings", std::make_shared< tudat::interpolators::LagrangeInterpolatorSettings >( 8 ), "..." ),
           R"doc(
Create a tabulated ephemeris from SPICE data.

Parameters
----------
body : dynamics.environment.Body
    Body for which the tabulated ephemeris is created.
initial_time : float
    Start epoch of the tabulation.
end_time : float
    End epoch of the tabulation.
time_step : float
    Time interval between consecutive tabulated states.
observer_name : str
    Name of the SPICE observer.
reference_frame_name : str
    Name of the SPICE reference frame.
interpolator_settings : math.interpolators.InterpolatorSettings, default = math.interpolators.lagrange_interpolation(8, boundary_interpolation=math.interpolators.extrapolate_at_boundary)
    Settings used to interpolate the tabulated ephemeris.
)doc" );

    m.def( "create_body_ephemeris",
           &tss::createBodyEphemeris< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ephemeris_settings" ),
           py::arg( "body_name" ),
           R"doc(

 Function that creates an Ephemeris object.

 Function that creates an :class:`~tudatpy.dynamics.environment.Ephemeris` object, but does *not*
 associate it with any specific body (e.g., it does not go into the environment, but can be used independently of it)


 Parameters
 ----------
 ephemeris_settings : EphemerisSettings
     Object defining the ephemeris settings.
 body_name : str
     Name of body for which the ephemeris is created. Note that this input is only relevant for some ephemeris settings (for instance, a spice ephemeris setting), and it does *not* imply that the ephemeris object is associated with a Body object of this name.
 Returns
 -------
 :class:`~tudatpy.dynamics.environment.Ephemeris`
     Ephemeris object, created according to the provided settings






     )doc" );

    m.def( "create_ground_station_ephemeris",
           py::overload_cast< const std::shared_ptr< tss::Body >, const std::string&, const tss::SystemOfBodies& >(
                   &tss::createReferencePointEphemerisFromId< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "body_with_ground_station" ),
           py::arg( "station_name" ),
           py::arg( "bodies" ) );

    py::class_< RelativisticTimeConverterSettings, std::shared_ptr< RelativisticTimeConverterSettings > >(
            m, "DirectRelativisticTimeConverterSettings" );

    m.def( "direct_relativistic_time_converter_settings",
           &directRelativisticTimeConverterSettings,
           py::arg( "barycentric_to_bodycentric_settings" ),
           py::arg( "integrator_settings" ),
           py::arg( "bodycentric_to_topocentric_settings" ) = std::vector< std::shared_ptr< RelativisticTimePropagatorSettings > >( ),
           R"doc(

 Create settings for a direct relativistic time converter.

 This function combines:

 1. One barycentric↔body-centered conversion settings object, and
 2. Zero or more body-centered↔topocentric conversion settings objects

 into a single converter-settings object for one body.

 The ``barycentric_to_bodycentric_settings`` input should be created with:

 - :func:`~tudatpy.dynamics.propagation_setup.propagator.first_order_bodycentric_relativistic_time_settings`.

 Each entry in ``bodycentric_to_topocentric_settings`` should typically be created with:

 - :func:`~tudatpy.dynamics.propagation_setup.propagator.bodycentered_to_topocentric_time_settings`.

 This function only assembles converter settings. Use
 :func:`~set_relativistic_time_converters` to attach them to bodies.

 Parameters
 ----------
 barycentric_to_bodycentric_settings : RelativisticTimePropagatorSettings
     Settings object defining the barycentric↔body-centered leg.
 integrator_settings : IntegratorSettings
     Numerical integrator settings used when creating the direct converter.
 bodycentric_to_topocentric_settings : list[RelativisticTimePropagatorSettings], optional
     Optional list of settings objects defining body-centered↔topocentric legs.
     Each list entry typically corresponds to one reference point/station.

 Returns
 -------
 DirectRelativisticTimeConverterSettings
     Settings object used by :func:`~set_relativistic_time_converters`.

        )doc" );

    m.def( "set_relativistic_time_converters",
           &setRelativisticTimeConverters,
           py::arg( "bodies" ),
           py::arg( "converter_settings" ),
           R"doc(

 Attach relativistic time converters to bodies.

 This function takes the converter settings assembled with
 :func:`~direct_relativistic_time_converter_settings` and instantiates the
 corresponding converter models in the provided system of bodies.

 For each entry in ``converter_settings``, Tudat sets up:

 - one barycentric↔body-centered conversion leg (first- or second-order), and
 - zero or more body-centered↔topocentric conversion legs.

 The key of each dictionary entry is typically the associated body name, while
 the converter content is defined by the corresponding
 :class:`~DirectRelativisticTimeConverterSettings` object.

 The converter settings used here are typically created from:

 - :func:`~tudatpy.dynamics.propagation_setup.propagator.first_order_bodycentric_relativistic_time_settings`
   for the barycentric↔body-centered leg, and
 - :func:`~tudatpy.dynamics.propagation_setup.propagator.bodycentered_to_topocentric_time_settings`
   for optional topocentric legs.

 After this function returns, each configured body can provide time-scale
 differences through
 :func:`~tudatpy.dynamics.environment.Body.get_time_scale_converter`.

 Parameters
 ----------
 bodies : SystemOfBodies
     The system of bodies to which time converters are attached.
 converter_settings : dict[str, DirectRelativisticTimeConverterSettings]
     Mapping from identifiers (typically body names) to direct converter
     settings objects. Each entry creates one relativistic time converter
     configuration.

 Returns
 -------
 None
     This function modifies ``bodies`` in place by attaching converter models.

        )doc" );
}

}  // namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
