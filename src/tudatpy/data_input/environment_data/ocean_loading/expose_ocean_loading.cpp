#include "expose_ocean_loading.h"

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "tudat/io/readOceanPoleTideData.h"

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace ocean_loading
{

void expose_ocean_loading( pybind11::module& module )
{
    namespace py = pybind11;
    namespace io = tudat::input_output;
    py::class_< io::OceanPoleTideGrid >( module, "OceanPoleTideGrid", R"doc(

Ocean pole tide loading coefficients on a regular longitude/latitude grid (IERS Conventions 2010, Section 7.1.5).

The coefficients are dimensionless. The grid is stored as read, without interpolation. Objects of this class are
created by :func:`read_ocean_pole_tide_grid`, and used by
:func:`~tudatpy.dynamics.environment_setup.shape_deformation.ocean_pole_tidal_from_data`.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "longitude_count", &io::OceanPoleTideGrid::longitudeCount, R"doc(

         Number of grid longitudes.

         :type: int
      )doc" )
            .def_readwrite( "latitude_count", &io::OceanPoleTideGrid::latitudeCount, R"doc(

         Number of grid latitudes.

         :type: int
      )doc" )
            .def_readwrite( "first_longitude_degrees", &io::OceanPoleTideGrid::firstLongitudeDegrees, R"doc(

         First grid longitude in degrees.

         :type: float
      )doc" )
            .def_readwrite( "first_latitude_degrees", &io::OceanPoleTideGrid::firstLatitudeDegrees, R"doc(

         First grid latitude in degrees.

         :type: float
      )doc" )
            .def_readwrite( "longitude_step_degrees", &io::OceanPoleTideGrid::longitudeStepDegrees, R"doc(

         Longitude spacing of the grid in degrees.

         :type: float
      )doc" )
            .def_readwrite( "latitude_step_degrees", &io::OceanPoleTideGrid::latitudeStepDegrees, R"doc(

         Latitude spacing of the grid in degrees.

         :type: float
      )doc" )
            .def_readwrite( "coefficients", &io::OceanPoleTideGrid::coefficients, R"doc(

         Dimensionless coefficients at each grid point, ordered by latitude and then by longitude (longitude varies
         fastest). Each matrix has the radial, north and east components as rows, and the real and imaginary parts as
         columns.

         :type: list[numpy.ndarray[numpy.float64[3, 2]]]
      )doc" );

    module.def( "read_ocean_pole_tide_grid",
                py::overload_cast< const std::string& >( &io::readOceanPoleTideGrid ),
                py::arg( "file_path" ),
                R"doc(

Function to read an ocean pole tide loading grid.

Function to read an ocean pole tide loading grid in the plain-text format of the IERS conventional grid, which is
available from the `IERS Conventions Centre
<https://iers-conventions.obspm.fr/content/chapter7/additional_info/opoleloadcoefcmcor.txt.gz>`_. Compressed files must
be decompressed before reading. No interpolation or origin conversion is performed (the official grid is given with
respect to the centre of mass).


Parameters
----------
file_path : str
    Path of the (decompressed) grid file that is to be read.

Returns
-------
OceanPoleTideGrid
    Grid read from the file.

      )doc" );
}

}  // namespace ocean_loading
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
