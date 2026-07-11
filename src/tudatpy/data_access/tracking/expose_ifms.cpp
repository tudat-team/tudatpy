#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_ifms.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <limits>

#include "scalarTypes.h"
#include "tudat/io/preProcessIfmsFile.h"
#include "tudat/io/readTrackingTxtFile.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace ifms
{

void expose_ifms( py::module& m )
{
    m.def( "read_ifms_file",
           &tio::readIfmsFile,
           py::arg( "file_name" ),
           py::arg( "apply_tropospheric_correction" ) = true,
           py::arg( "remove_invalid_lines" ) = true );

    m.def( "read_ifms_files",
           py::overload_cast< const std::vector< std::string >&,
                              const std::string&,
                              const std::vector< std::string >&,
                              const std::string&,
                              bool,
                              bool,
                              const std::vector< double >&,
                              double,
                              double >( &tio::readIfmsFiles< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "ifms_file_names" ),
           py::arg( "spacecraft_name" ),
           py::arg( "ground_station_names" ),
           py::arg( "earth_name" ) = "Earth",
           py::arg( "apply_tropospheric_correction" ) = true,
           py::arg( "remove_invalid_lines" ) = true,
           py::arg( "frequency_bands" ) = std::vector< double >( ),
           py::arg( "reception_reference_frequency_band" ) = std::numeric_limits< double >::quiet_NaN( ),
           py::arg( "doppler_reference_frequency" ) = std::numeric_limits< double >::quiet_NaN( ) );
}

}  // namespace ifms

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy
