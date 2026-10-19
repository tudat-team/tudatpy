#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_ifms.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <limits>

#include "scalarTypes.h"
#include "tudat/io/preProcessIfmsFile.h"

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
                                           const std::vector< std::string >&,
                                           const std::string&,
                                           double >( &tio::readIfmsFiles< STATE_SCALAR_TYPE, TIME_TYPE > );

    m.def( "read_ifms_data",
           readIfmsData,
           py::arg( "ifms_file_names" ),
           py::arg( "spacecraft_name" ),
           py::arg( "ground_station_names" ),
           py::arg( "earth_name" ) = "Earth",
           py::arg( "apply_tropospheric_correction" ) = true,
           py::arg( "remove_invalid_lines" ) = true,
           py::arg( "frequency_bands" ) = std::vector< std::string >( ),
           py::arg( "reception_reference_frequency_band" ) = std::string( "" ),
           py::arg( "doppler_reference_frequency" ) = std::numeric_limits< double >::quiet_NaN( ),
           R"doc(
         Load Level 2 IFMS Doppler files into tracking data and supplementary data objects.

         The file format is described in `IFMS Doppler Processing Software : Level 1a to Level 2, Table 3-2 and 3-3 <https://archives.esac.esa.int/psa/ftp/MARS-EXPRESS/MRS/MEX-M-MRS-1-2-3-EXT9-4441-V1.0/DOCUMENT/MRS_DOC/MEX_MRS_IGM_DS_3035.PDF>`_.
         The reader parses each IFMS file, attaches the corresponding ground
         station name and optional frequency-band metadata, and converts the
         resulting records to Tudat ``TrackingData`` and
         ``TrackingSupplementaryData`` objects.
         When ``apply_tropospheric_correction`` is ``True``, the tropospheric
         correction is subtracted from the averaged-frequency observable before
         the tracking-data objects are created.

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
         frequency_bands : list[str], default=[]
             Frequency-band identifiers used for supplementary data (to be used
             for all files). Allowed values are ``"S-band"``, ``"X-band"``,
             ``"Ka-band"``, and ``"Ku-band"``.
         reception_reference_frequency_band : str, optional
             Reference reception frequency band (to be used for all files),
             using the same string convention as ``frequency_bands``.
         doppler_reference_frequency : float, optional
             Reference Doppler frequency (to be used for all files).

         Returns
         -------
         tuple[list[TrackingData], list[TrackingSupplementaryData]]
             Tracking data objects and supplementary data objects.

         Example
         -------

         In this example, we load two IFMS files for Mars Express (MEX), tracked by New Norcia (NNOR) in X-band.
         The IFMS files can be downloaded from the `ESA PSA archive <https://archives.esac.esa.int/psa/ftp/MARS-EXPRESS/MRS/MEX-M-MRS-1-2-3-EXT5-3682-V1.0/DATA/LEVEL02/CLOSED_LOOP/IFMS/DP2/>`_.


         .. code-block:: python

            from tudatpy.data_input.tracking_data.ifms import read_ifms_data

            sc_name = "MEX"

            ifms_file_names = [
                "M32ICL2L02_D2X_150050821_00.TAB",
                "M32ICL2L02_D2X_150080523_00.TAB",
            ]

            station_mapping = {"32": "NNOR"}
            ground_station_names_per_file = [station_mapping[f[1:3]] for f in ifms_file_names]

            tracking_data, tracking_supplementary_data = read_ifms_data(
                ifms_file_names,
                spacecraft_name=sc_name,
                ground_station_names=ground_station_names_per_file,
                frequency_bands=["X-band", "X-band"],
                reception_reference_frequency_band="X-band",
                apply_tropospheric_correction=False,
            )


      )doc" );
}

}  // namespace ifms

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy
