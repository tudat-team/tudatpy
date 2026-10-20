#ifndef TUDATPY_ESTIMATION_OBSERVATIONS_BINDINGS_H
#define TUDATPY_ESTIMATION_OBSERVATIONS_BINDINGS_H

#include <string>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace observations
{

//! Shared warning for the temporary legacy observation interfaces.
inline void warnLegacyObservationInterface( const std::string& interfaceName, const std::string& replacementApi = "ObservationDataset" )
{
    const std::string apiAnchor = replacementApi.substr( 0, replacementApi.find_first_of( " ," ) );
    const std::string message = interfaceName + " is deprecated and kept only for backwards compatibility. Use " + replacementApi +
            " instead. API reference: https://py.api.tudat.space/en/latest/estimation/observations.html#"
            "tudatpy.estimation.observations." +
            apiAnchor +
            ". Migration guide: https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-dataset-deprecation.html";
    if( PyErr_WarnEx( PyExc_DeprecationWarning, message.c_str( ), 1 ) < 0 )
    {
        throw py::error_already_set( );
    }
}

void expose_observations_io_bindings( py::module& m );
void expose_observations_simulation_bindings( py::module& m );

}  // namespace observations
}  // namespace estimation
}  // namespace tudatpy

#endif
