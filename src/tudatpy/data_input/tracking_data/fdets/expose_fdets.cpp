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

namespace data_input
{

namespace tracking_data
{

namespace fdets
{

void expose_fdets( py::module& m )
{
    py::enum_< tio::FdetDateFormat >( m, "FdetDateFormat", R"doc(Date format used in an Fdets file.)doc" )
            .value( "datetime_string",
                    tio::FdetDateFormat::datetime_string,
                    R"doc(Date is provided in a single UTC datetime string column.)doc" )
            .value( "pair_of_numbers",
                    tio::FdetDateFormat::pair_of_numbers,
                    R"doc(Date is provided as a pair of numeric columns. This format is not currently supported.)doc" )
            .export_values( );

    auto readFdetsData = static_cast< std::pair< std::vector< std::shared_ptr< tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE > > >,
                                                 std::vector< std::shared_ptr< tdat::TrackingSupplementaryData > > > ( * )(
            const std::vector< std::string >&,
            const std::vector< double >&,
            const tio::FdetDateFormat,
            const std::string&,
            const std::vector< std::string >&,
            const std::vector< std::string >&,
            const std::string& ) >( &tio::readFdetsFiles< STATE_SCALAR_TYPE, TIME_TYPE > );

    m.def( "read_fdets_data",
           readFdetsData,
           py::arg( "fdets_file_names" ),
           py::arg( "base_frequencies" ),
           py::arg( "date_format" ),
           py::arg( "spacecraft_name" ),
           py::arg( "transmitting_station_names" ),
           py::arg( "receiving_station_names" ),
           py::arg( "earth_name" ) = "Earth",
           R"doc(
         Load Fdets PRIDE Doppler data files into tracking data and supplementary data objects.

         Each data record is expected to contain a time tag (typically  a UTC datetime string), signal-to-noise
         ratio, normalized spectral maximum, Doppler measured frequency, and Doppler
         noise. If the file contains an additional leading scan-number column, this
         column is detected automatically.
         The reader parses the Fdets text files, assigns the supplied base
         frequencies and station names per file, and converts the resulting
         records to Tudat ``TrackingData`` and ``TrackingSupplementaryData``
         objects.

         Parameters
         ----------
         fdets_file_names : list[str]
             Paths to Fdets files.
         base_frequencies : list[float]
             Base frequencies associated with the input files, per file. This base frequency is added to the frequency column in the file to obtain the observed sky frequencies.
         date_format : FdetDateFormat
             Date format used in the Fdets files (to be used for all files).
         spacecraft_name : str
             Name assigned to the spacecraft link end (to be used for all files).
         transmitting_station_names : list[str]
             Names assigned to transmitting stations, per file.
         receiving_station_names : list[str]
             Names assigned to receiving stations, per file.
         earth_name : str, default="Earth"
             Name assigned to the Earth body (to be used only when it is needed to place the ground stations on a body not named 'Earth').

         Returns
         -------
         tuple[list[TrackingData], list[TrackingSupplementaryData]]
             Tracking data objects and supplementary data objects.
      )doc" );
}

}  // namespace fdets

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy
