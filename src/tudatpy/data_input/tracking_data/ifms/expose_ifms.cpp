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
    m.def( "read_ifms_file",
           &tio::readIfmsFile,
           py::arg( "file_name" ),
           py::arg( "apply_tropospheric_correction" ) = true,
           py::arg( "remove_invalid_lines" ) = true,
           R"doc(Load contents of IFMS file into object

           The keys of the dictionary represent the different columns of the IFMS file, and their values are lists with all the values in the associated column as strings.

           Two of the columns of an IFMS file contain, respectively, the Doppler averaged frequency and a tropospheric correction for the station. When the `apply_tropospheric_correction` option is set to true, the content of the first column is modified by subtracting the values in the second.

           Parameters
           ----------
           file_name : str
               String representing the path to the file to be loaded
           apply_tropospheric_correction : bool
               Whether to modify the averaged Doppler frequency as described above (Default: True)
           remove_invalid_lines : bool
               Boolean defining whether a line is skipped if the transmit frequency, observed frequency, or troposphere correction is undefined (Default: True)

           Returns
           -------
           ifms_contents : TrackingTxtFileContents
               Dictionary with contents of the IFMS file as lists of strings
           )doc" );

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
           py::arg( "doppler_reference_frequency" ) = std::numeric_limits< double >::quiet_NaN( ),
           R"doc(
Load IFMS files into tracking data and supplementary data objects.

Parameters
----------
ifms_file_names : list[str]
    Paths to IFMS files.
spacecraft_name : str
    Name assigned to the spacecraft link end.
ground_station_names : list[str]
    Names assigned to receiving ground stations.
earth_name : str, default="Earth"
    Name assigned to the Earth body.
apply_tropospheric_correction : bool, default=True
    Whether the troposphere correction column is applied.
remove_invalid_lines : bool, default=True
    Whether records with invalid frequencies or corrections are skipped.
frequency_bands : list[float], default=[]
    Frequency-band identifiers used for supplementary data.
reception_reference_frequency_band : float, optional
    Reference reception frequency band.
doppler_reference_frequency : float, optional
    Reference Doppler frequency.

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
