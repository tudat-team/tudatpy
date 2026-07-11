#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_generic_text_file.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "tudat/io/readTrackingTxtFile.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace generic_text_file
{

void expose_generic_text_file( py::module& m )
{
    py::enum_< tio::TrackingTxtFileReadFilterType >( m, "TrackingTxtFileReadFilterType" )
            .value( "no_tracking_txt_file_filter", tio::TrackingTxtFileReadFilterType::no_tracking_txt_file_filter )
            .value( "ifms_tracking_txt_file_filter", tio::TrackingTxtFileReadFilterType::ifms_tracking_txt_file_filter )
            .export_values( );

    py::class_< tio::TrackingTxtFileContents, std::shared_ptr< tio::TrackingTxtFileContents > >( m, "TrackingTxtFileContents" )
            .def( py::init< const std::string, const std::vector< std::string >, const char, const std::string >( ),
                  py::arg( "file_name" ),
                  py::arg( "column_types" ),
                  py::arg( "comment_symbol" ) = '#',
                  py::arg( "value_separators" ) = ",:\t " )
            .def_property_readonly( "column_field_types", &tio::TrackingTxtFileContents::getRawColumnTypes )
            .def_property_readonly( "double_datamap", &tio::TrackingTxtFileContents::getDoubleDataMap )
            .def_property_readonly( "raw_datamap", &tio::TrackingTxtFileContents::getRawDataMap )
            .def_property_readonly( "num_rows", &tio::TrackingTxtFileContents::getNumRows );

    m.def( "read_tracking_txt_file",
           &tio::createTrackingTxtFileContents,
           py::arg( "file_name" ),
           py::arg( "column_types" ),
           py::arg( "comment_symbol" ) = '#',
           py::arg( "value_separators" ) = ",:\t ",
           py::arg( "ignore_omitted_columns" ) = false,
           py::arg( "data_filter_method" ) = tio::no_tracking_txt_file_filter );
}

}  // namespace generic_text_file

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy
