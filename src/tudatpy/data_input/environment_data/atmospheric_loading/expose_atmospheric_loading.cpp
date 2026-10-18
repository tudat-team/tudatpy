#include "expose_atmospheric_loading.h"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "tudat/io/readAtmosphericTideData.h"
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace atmospheric_loading
{
void expose_atmospheric_loading( pybind11::module& module )
{
    namespace py = pybind11;
    namespace io = tudat::input_output;

    py::class_< io::AtmosphericTideGrid >( module, "AtmosphericTideGrid", R"doc(

Regular latitude/longitude grid of atmospheric S1/S2 loading coefficients. Objects of this class are created by
:func:`read_atmospheric_tide_grid`.

      )doc" )
            .def_readonly( "longitudes_degrees", &io::AtmosphericTideGrid::longitudesDegrees, R"doc(

         **read-only**

         Longitudes of the grid points, in degrees.

         :type: list[float]
      )doc" )
            .def_readonly( "latitudes_degrees", &io::AtmosphericTideGrid::latitudesDegrees, R"doc(

         **read-only**

         Latitudes of the grid points, in degrees.

         :type: list[float]
      )doc" )
            .def_readonly( "coefficients_metres", &io::AtmosphericTideGrid::coefficientsMetres, R"doc(

         **read-only**

         Coefficients at the grid points, in metres, longitude-major. Rows are the radial, north and east components,
         columns the :math:`\cos S_{1}`, :math:`\sin S_{1}`, :math:`\cos S_{2}` and :math:`\sin S_{2}` terms.

         :type: list[numpy.ndarray[numpy.float64[3, 4]]]
      )doc" );

    module.def( "read_atmospheric_tide_grid",
                py::overload_cast< const std::string& >( &io::readAtmosphericTideGrid ),
                py::arg( "file_path" ),
                R"doc(

Function to read an atmospheric S1/S2 loading grid.

Function to read an atmospheric S1/S2 loading grid in the 14-column Bernese format (as distributed for the Ray and
Ponte, 2003, model): longitude, latitude and twelve coefficients in millimetres, being the radial, north and east
components of the :math:`\cos S_{1}`, :math:`\sin S_{1}`, :math:`\cos S_{2}` and :math:`\sin S_{2}` terms. Header lines
start with ``$``. The coefficients are converted to metres.


Parameters
----------
file_path : str
    Path of the file that is to be read.

Returns
-------
AtmosphericTideGrid
    Grid read from the file.

      )doc" );
    module.def( "read_atmospheric_tide_station_data",
                py::overload_cast< const std::string& >( &io::readAtmosphericTideStationData ),
                py::arg( "file_path" ),
                R"doc(

Function to read station-specific atmospheric S1/S2 loading coefficients.

Function to read station-specific atmospheric S1/S2 loading coefficients in the Bernese ATL format: for each station,
a line with its identifier, followed by three lines (radial, north, east) of four coefficients in millimetres. Header
lines start with ``$``. The coefficients are converted to metres.


Parameters
----------
file_path : str
    Path of the file that is to be read.

Returns
-------
dict[str, numpy.ndarray[numpy.float64[3, 4]]]
    Coefficients (metres) per station name; rows radial, north and east, columns :math:`\cos S_{1}`, :math:`\sin S_{1}`,
    :math:`\cos S_{2}` and :math:`\sin S_{2}`.

      )doc" );
}
}  // namespace atmospheric_loading
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
