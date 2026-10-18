#include "expose_station_loading.h"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "tudat/io/readStationLoadingData.h"

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace station_loading
{
void expose_station_loading( pybind11::module& module )
{
    namespace py = pybind11;
    namespace io = tudat::input_output;

    py::class_< io::StationLoadingSite >( module, "StationLoadingSite", R"doc(

Site of a station loading series, with its position and displacements.

      )doc" )
            .def_readonly( "cartesian_position", &io::StationLoadingSite::cartesianPosition, R"doc(

         **read-only**

         Cartesian position of the site, in metres.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readonly( "validity_radius_metres", &io::StationLoadingSite::validityRadiusMetres, R"doc(

         **read-only**

         Distance from the site within which a station may use its data, in metres.

         :type: float
      )doc" )
            .def_readonly( "displacements", &io::StationLoadingSite::displacements, R"doc(

         **read-only**

         Geocentric radial, east and north displacements in metres, with the index of their epoch as key.

         :type: dict[int, numpy.ndarray[numpy.float64[3, 1]]]
      )doc" );

    py::class_< io::StationLoadingData >( module, "StationLoadingData", R"doc(

Station loading series: displacements of a set of sites at a set of epochs. Objects of this class are created by
:func:`read_ephedisp_station_loading_data`.

      )doc" )
            .def_readonly( "source", &io::StationLoadingData::source, R"doc(

         **read-only**

         Name of the file of the series.

         :type: str
      )doc" )
            .def_readonly( "epochs", &io::StationLoadingData::epochs, R"doc(

         **read-only**

         Epochs of the displacements, in TAI.

         :type: list[astro.time_representation.Time]
      )doc" )
            .def_readonly( "sites", &io::StationLoadingData::sites, R"doc(

         **read-only**

         Sites, with their identifier in the file as key.

         :type: dict[str, StationLoadingSite]
      )doc" );

    module.def( "read_ephedisp_station_loading_data",
                py::overload_cast< const std::string& >( &io::readEphedispStationLoadingData ),
                py::arg( "file_path" ),
                R"doc(

Function to read a station loading series from an EPHEDISP file.

Function to read a station loading series from a file in the `EPHEDISP <https://astrogeo.org/malo/ephedisp_format.txt>`_
format (2005.06.30 version), as distributed for instance by the International Mass Loading Service. The epochs are in
TAI, and are defined by the first epoch and the sampling interval in the file. The displacements are in metres, as
geocentric radial, east and north components at the site locations given in the file.

The physical content of the series (loading components, removed signals, reference origin) is not given in the file.
Avoiding signals that are modelled twice is the responsibility of the user.


Parameters
----------
file_path : str
    Path of the file that is to be read.

Returns
-------
StationLoadingData
    Series read from the file.

      )doc" );
}
}  // namespace station_loading
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
