#include "expose_ocean_loading.h"

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "tudat/io/readOceanLoadingData.h"
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
    py::enum_< io::OceanLoadingReferenceOrigin >( module, "OceanLoadingReferenceOrigin", R"doc(

Enumeration of reference origins with respect to which ocean loading displacements are given. In BLQ files, ``CMC: NO``
denotes the centre of the solid Earth, and ``CMC: YES`` the centre of mass of the Earth and oceans.

      )doc" )
            .value( "unspecified", io::OceanLoadingReferenceOrigin::unspecified, "Origin not declared." )
            .value( "solid_earth_centre", io::OceanLoadingReferenceOrigin::solidEarthCentre, "Centre of the solid Earth (CE)." )
            .value( "whole_earth_centre_of_mass",
                    io::OceanLoadingReferenceOrigin::wholeEarthCentreOfMass,
                    "Centre of mass of the Earth and oceans (CM)." );

    py::class_< io::OceanLoadingStation >( module, "OceanLoadingStation", R"doc(

Ocean tidal loading coefficients of a single station, as given in a BLQ file.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "station_identifier", &io::OceanLoadingStation::stationIdentifier, R"doc(

         Identifier of the station, exactly as given in the file.

         :type: str
      )doc" )
            .def_readwrite( "amplitudes_metres", &io::OceanLoadingStation::amplitudesMetres, R"doc(

         Amplitudes in metres. The rows are the up (radial), west and south components; the columns are the eleven
         constituents in the order given by :func:`constituent_names`.

         :type: numpy.ndarray[numpy.float64[3, 11]]
      )doc" )
            .def_readwrite( "phase_lags_degrees", &io::OceanLoadingStation::phaseLagsDegrees, R"doc(

         Greenwich phase lags in degrees (positive for a lag), with the same layout as :attr:`amplitudes_metres`.

         :type: numpy.ndarray[numpy.float64[3, 11]]
      )doc" )
            .def_readwrite( "source", &io::OceanLoadingStation::source, R"doc(

         Name of the file from which the station was read.

         :type: str
      )doc" )
            .def_readwrite( "source_headers", &io::OceanLoadingStation::sourceHeaders, R"doc(

         Header lines of the file, unmodified.

         :type: list[str]
      )doc" )
            .def_readwrite( "ocean_tide_model", &io::OceanLoadingStation::oceanTideModel, R"doc(

         Ocean tide model used to compute the coefficients, as given in the file.

         :type: str
      )doc" )
            .def_readwrite( "greens_function", &io::OceanLoadingStation::greensFunction, R"doc(

         Green functions used to compute the coefficients, as given in the file.

         :type: str
      )doc" )
            .def_readwrite( "reference_origin", &io::OceanLoadingStation::referenceOrigin, R"doc(

         Reference origin of the displacements, as declared in the file.

         :type: OceanLoadingReferenceOrigin
      )doc" )
            .def_readwrite( "provider_longitude_latitude_height", &io::OceanLoadingStation::providerLongitudeLatitudeHeight, R"doc(

         Longitude and latitude (degrees) and height (metres) of the station, as given in the file (only meaningful if
         :attr:`has_provider_coordinates` is true).

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_readwrite( "has_provider_coordinates", &io::OceanLoadingStation::hasProviderCoordinates, R"doc(

         Whether the file gives the coordinates of the station.

         :type: bool
      )doc" );

    py::class_< io::OceanLoadingData >( module, "OceanLoadingData", R"doc(

Ocean tidal loading coefficients of a set of stations, as read from a BLQ file.

Objects of this class are created by :func:`read_ocean_loading_data`, and used by
:func:`~tudatpy.dynamics.environment_setup.shape_deformation.ocean_tidal_from_data`.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::OceanLoadingData::source, R"doc(

         Name of the file from which the data was read.

         :type: str
      )doc" )
            .def_readwrite( "source_headers", &io::OceanLoadingData::sourceHeaders, R"doc(

         Header lines of the file, unmodified.

         :type: list[str]
      )doc" )
            .def_readwrite( "stations", &io::OceanLoadingData::stations, R"doc(

         Coefficients of each station, with the station identifier as key.

         :type: dict[str, OceanLoadingStation]
      )doc" )
            .def( "get_station",
                  &io::OceanLoadingData::getStation,
                  py::arg( "station_identifier" ),
                  py::return_value_policy::reference_internal,
                  R"doc(

         Function to retrieve the coefficients of a single station.

         Parameters
         ----------
         station_identifier : str
             Identifier of the station, exactly as given in the file. An unknown identifier raises an exception.

         Returns
         -------
         OceanLoadingStation
             Coefficients of the station.

      )doc" )
            .def(
                    "validate",
                    []( const io::OceanLoadingData& data ) { io::validateOceanLoadingData( data ); },
                    R"doc(

         Check the header conventions, the coefficients (finite, with non-negative amplitudes) and the station
         identifiers, and raise an exception if they are invalid.

      )doc" );

    module.def( "read_ocean_loading_data",
                py::overload_cast< const std::string& >( &io::readOceanLoadingData ),
                py::arg( "file_path" ),
                R"doc(

Function to read ocean tidal loading coefficients from a BLQ file.

Function to read ocean tidal loading coefficients from a file in the BLQ format (as produced, for instance, by the
Onsala ocean loading service). For each station, the file gives the amplitudes (metres) and Greenwich phase lags
(degrees) of the up, west and south displacements for eleven tidal constituents. The header lines, ocean tide model,
Green functions and reference origin (``CMC`` declaration) are retained; the reference origin is not converted.

Duplicated station identifiers, incomplete station blocks, non-finite values and incompatible column orders or units
raise an exception.


Parameters
----------
file_path : str
    Path of the BLQ file that is to be read.

Returns
-------
OceanLoadingData
    Coefficients of all stations in the file.

      )doc" );

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
      )doc" )
            .def_readwrite( "source", &io::OceanPoleTideGrid::source, R"doc(

         Name of the file from which the grid was read.

         :type: str
      )doc" )
            .def_readwrite( "source_headers", &io::OceanPoleTideGrid::sourceHeaders, R"doc(

         Header lines of the file, unmodified.

         :type: list[str]
      )doc" )
            .def_readwrite( "reference_origin", &io::OceanPoleTideGrid::referenceOrigin, R"doc(

         Reference origin of the displacements, as declared in the file. The official IERS grid is given with respect to
         the centre of mass, and already includes the corresponding correction.

         :type: OceanLoadingReferenceOrigin
      )doc" )
            .def( "validate", &io::validateOceanPoleTideGrid, R"doc(

         Check that the grid is complete, regular and finite, and raise an exception if not.

      )doc" );

    module.def( "read_ocean_pole_tide_grid",
                py::overload_cast< const std::string& >( &io::readOceanPoleTideGrid ),
                py::arg( "file_path" ),
                R"doc(

Function to read an ocean pole tide loading grid.

Function to read an ocean pole tide loading grid in the plain-text format of the IERS conventional grid, which is
available from the `IERS Conventions Centre
<https://iers-conventions.obspm.fr/content/chapter7/additional_info/opoleloadcoefcmcor.txt.gz>`_. Compressed files must
be decompressed before reading. The header lines and the declared reference origin are retained; no interpolation or
origin conversion is performed.


Parameters
----------
file_path : str
    Path of the (decompressed) grid file that is to be read.

Returns
-------
OceanPoleTideGrid
    Grid read from the file.

      )doc" );

    module.def( "constituent_names", &io::getOceanLoadingConstituentNames, R"doc(

Function to retrieve the names of the eleven tidal constituents of a BLQ file, in the order of the file columns.

Returns
-------
list[str]
    Names of the constituents (M2, S2, N2, K2, K1, O1, P1, Q1, Mf, Mm, Ssa).

      )doc" );
}

}  // namespace ocean_loading
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
