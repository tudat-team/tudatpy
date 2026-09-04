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
           R"doc(Create a system of bodies from body settings.)doc" );

    m.def( "add_empty_tabulated_ephemeris",
           &tp::addEmptyTabulatedEphemeris< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "ephemeris_origin" ) = "",
           py::arg( "is_part_of_multi_arc" ) = false );

    m.def( "create_tabulated_ephemeris_from_spice",
           &tss::createTabulatedEphemerisFromSpice< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "body" ),
           py::arg( "initial_time" ),
           py::arg( "end_time" ),
           py::arg( "time_step" ),
           py::arg( "observer_name" ),
           py::arg( "reference_frame_name" ),
           py::arg_v( "interpolator_settings", std::make_shared< tudat::interpolators::LagrangeInterpolatorSettings >( 8 ), "..." ) );

    m.def( "create_body_ephemeris",
           &tss::createBodyEphemeris< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ephemeris_settings" ),
           py::arg( "body_name" ) );

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
           py::arg( "bodycentric_to_topocentric_settings" ) = std::vector< std::shared_ptr< RelativisticTimePropagatorSettings > >( ) );

    m.def( "set_relativistic_time_converters", &setRelativisticTimeConverters, py::arg( "bodies" ), py::arg( "converter_settings" ) );
}

}  // namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
