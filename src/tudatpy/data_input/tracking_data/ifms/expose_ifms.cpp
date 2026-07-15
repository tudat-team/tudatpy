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

namespace data_input
{

namespace tracking_data
{

namespace ifms
{

void expose_ifms( py::module& m )
{
    auto readIfmsData = py::overload_cast< const std::vector< std::string >&,
                                           const std::string&,
                                           const std::vector< std::string >&,
                                           const std::string&,
                                           bool,
                                           bool,
                                           const std::vector< double >&,
                                           double,
                                           double >( &tio::readIfmsFiles< STATE_SCALAR_TYPE, TIME_TYPE > );

    m.def( "read_ifms_data",
           readIfmsData,
           py::arg( "ifms_file_names" ),
           py::arg( "spacecraft_name" ),
           py::arg( "ground_station_names" ),
           py::arg( "earth_name" ) = "Earth",
           py::arg( "apply_tropospheric_correction" ) = true,
           py::arg( "remove_invalid_lines" ) = true,
           py::arg( "frequency_bands" ) = std::vector< double >( ),
           py::arg( "reception_reference_frequency_band" ) = std::numeric_limits< double >::quiet_NaN( ),
           py::arg( "doppler_reference_frequency" ) = std::numeric_limits< double >::quiet_NaN( ),
           R"doc(
         Load IFMS files into tracking data and supplementary data objects.

         IFMS records contain radio science data, typically closed-loop Doppler
         data, and optionally a station tropospheric-correction column. When
         ``apply_tropospheric_correction`` is ``True``, the tropospheric
         correction is subtracted from the averaged-frequency observable before
         the tracking-data objects are created.
         The file format is described in
         :cite:t:`ifmsOccFtp2006`.
         The reader parses each IFMS file, attaches the corresponding ground
         station name and optional frequency-band metadata, and converts the
         resulting records to Tudat ``TrackingData`` and
         ``TrackingSupplementaryData`` objects.

         Parameters
         ----------
         ifms_file_names : list[str]
             Paths to IFMS files.
         spacecraft_name : str
             Name assigned to the spacecraft link end (to be used for all files).
         ground_station_names : list[str]
             Names assigned to receiving ground stations, per file.
         earth_name : str, default="Earth"
             Name assigned to the Earth body (to be used only when it is needed
             to place the ground stations on a body not named 'Earth').
         apply_tropospheric_correction : bool, default=True
             Whether the troposphere correction column is applied (to be used
             for all files).
         remove_invalid_lines : bool, default=True
             Whether records with invalid frequencies or corrections are skipped
             (to be used for all files).
         frequency_bands : list[float], default=[]
             Frequency-band identifiers used for supplementary data (to be used
             for all files).
         reception_reference_frequency_band : float, optional
             Reference reception frequency band (to be used for all files).
         doppler_reference_frequency : float, optional
             Reference Doppler frequency (to be used for all files).

         Returns
         -------
         tuple[list[TrackingData], list[TrackingSupplementaryData]]
             Tracking data objects and supplementary data objects.
      )doc" );
}

}  // namespace ifms

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy
