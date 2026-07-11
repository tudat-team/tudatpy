#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_fdets.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/io/preProcessFdetsFile.h"
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingSupplementaryData.h"

namespace py = pybind11;
namespace tdat = tudat::data;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace fdets
{

void expose_fdets( py::module& m )
{
    py::enum_< tio::FdetDateFormat >( m, "FdetDateFormat" )
            .value( "datetime_string", tio::FdetDateFormat::datetime_string )
            .value( "pair_of_numbers", tio::FdetDateFormat::pair_of_numbers )
            .export_values( );

    m.def( "read_fdets_file",
           py::overload_cast< const std::string&, tio::FdetDateFormat >( &tio::readFdetsFile ),
           py::arg( "file_name" ),
           py::arg( "date_format" ) = tio::FdetDateFormat::datetime_string );

    m.def( "read_fdets_files",
           static_cast< std::pair< std::vector< std::shared_ptr< tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE > > >,
                                   std::vector< std::shared_ptr< tdat::TrackingSupplementaryData > > > ( * )(
                   const std::vector< std::string >&,
                   const std::vector< double >&,
                   const tio::FdetDateFormat,
                   const std::string&,
                   const std::vector< std::string >&,
                   const std::vector< std::string >&,
                   const std::string& ) >( &tio::readFdetsFiles< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "fdets_file_names" ),
           py::arg( "base_frequencies" ),
           py::arg( "date_format" ),
           py::arg( "spacecraft_name" ),
           py::arg( "transmitting_station_names" ),
           py::arg( "receiving_station_names" ),
           py::arg( "earth_name" ) = "Earth" );
}

}  // namespace fdets

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy
