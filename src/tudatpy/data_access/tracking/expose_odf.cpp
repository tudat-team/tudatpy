#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_odf.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/io/preProcessOdfFile.h"
#include "tudat/io/readOdfFile.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace odf
{

void expose_odf( py::module& m )
{
    py::class_< tio::OdfRawFileContents, std::shared_ptr< tio::OdfRawFileContents > >( m, "RawOdfFileContents" )
            .def_readonly( "file_reference_date", &tio::OdfRawFileContents::fileReferenceDate_ )
            .def_readonly( "file_reference_time", &tio::OdfRawFileContents::fileReferenceTime_ )
            .def( "write_to_text_file", &tio::OdfRawFileContents::writeOdfToTextFile, py::arg( "output_file" ) );

    m.def( "read_raw_odf_file_contents", &tio::readOdfFile, py::arg( "file_name" ) );
    m.def( "read_odf_files",
           &tio::loadOdfFile< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "odf_file_names" ),
           py::arg( "spacecraft_name" ),
           py::arg( "earth_name" ) = "Earth",
           py::arg( "verbose_output" ) = true );
}

}  // namespace odf

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy
