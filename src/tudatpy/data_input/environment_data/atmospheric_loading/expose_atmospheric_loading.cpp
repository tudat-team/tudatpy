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
    py::enum_< io::AtmosphericLoadingReferenceOrigin >( module, "AtmosphericLoadingReferenceOrigin", R"doc(

Enumeration of reference origins with respect to which atmospheric loading displacements are given.

      )doc" )
            .value( "unspecified", io::AtmosphericLoadingReferenceOrigin::unspecified, "Origin not declared." )
            .value( "solid_earth_centre", io::AtmosphericLoadingReferenceOrigin::solidEarthCentre, "Centre of the solid Earth (CE)." )
            .value( "whole_earth_centre_of_mass",
                    io::AtmosphericLoadingReferenceOrigin::wholeEarthCentreOfMass,
                    "Centre of mass of the whole Earth, including the atmosphere (CM)." );
    py::enum_< io::AtmosphericTideLatitudeConvention >( module, "AtmosphericTideLatitudeConvention", R"doc(

Enumeration of latitude conventions with which a station is located in an atmospheric loading grid.

      )doc" )
            .value( "unspecified", io::AtmosphericTideLatitudeConvention::unspecified, "Convention not declared." )
            .value( "geocentric", io::AtmosphericTideLatitudeConvention::geocentric, "Geocentric latitude." )
            .value( "geodetic", io::AtmosphericTideLatitudeConvention::geodetic, "Geodetic latitude (requires the body shape)." );
    py::enum_< io::AtmosphericTideDisplacementBasis >( module, "AtmosphericTideDisplacementBasis", R"doc(

Enumeration of local bases in which the three displacement components of atmospheric loading coefficients are given.

      )doc" )
            .value( "unspecified", io::AtmosphericTideDisplacementBasis::unspecified, "Basis not declared." )
            .value( "geocentric_radial_north_east",
                    io::AtmosphericTideDisplacementBasis::geocentricRadialNorthEast,
                    "Geocentric radial, north and east directions." )
            .value( "geodetic_up_north_east",
                    io::AtmosphericTideDisplacementBasis::geodeticUpNorthEast,
                    "Ellipsoidal-normal up, north and east directions (requires the body shape)." );
    py::enum_< io::AtmosphericTideSourceProfile >( module, "AtmosphericTideSourceProfile", R"doc(

Enumeration of identified atmospheric loading products. The model name or file name alone does not identify a product.

      )doc" )
            .value( "unspecified", io::AtmosphericTideSourceProfile::unspecified, "Product not identified." )
            .value( "iers2010_rp03_no_ib_farrell",
                    io::AtmosphericTideSourceProfile::iers2010Rp03NoIbFarrell,
                    "Ray and Ponte (2003) S1/S2 loading without inverted-barometer response, with Farrell Green functions, as "
                    "used in the IERS Conventions 2010 (Section 7.1.3)." );

    py::class_< io::AtmosphericTideMetadata >( module, "AtmosphericTideMetadata", R"doc(

Description of an atmospheric S1/S2 loading data source.

An object of this class can be passed to the readers in this module to declare properties that are not given in the
file (for instance the reference origin, or the expected grid dimensions). The readers complete it with the
information they find in the file header; a declaration that conflicts with the header raises an exception.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::AtmosphericTideMetadata::source, R"doc(

         Name of the file (or other source) of the data.

         :type: str
      )doc" )
            .def_readwrite( "source_headers", &io::AtmosphericTideMetadata::sourceHeaders, R"doc(

         Header lines of the source, unmodified.

         :type: list[str]
      )doc" )
            .def_readwrite( "model_name", &io::AtmosphericTideMetadata::modelName, R"doc(

         Name of the atmospheric tide model (``Ray-Ponte 2003`` when recognized in the header).

         :type: str
      )doc" )
            .def_readwrite( "greens_function", &io::AtmosphericTideMetadata::greensFunction, R"doc(

         Green functions used to compute the loading (``Farrell elastic Green functions`` when recognized in the header).

         :type: str
      )doc" )
            .def_readwrite( "reference_origin", &io::AtmosphericTideMetadata::referenceOrigin, R"doc(

         Reference origin of the displacements.

         :type: AtmosphericLoadingReferenceOrigin
      )doc" )
            .def_readwrite( "latitude_convention", &io::AtmosphericTideMetadata::latitudeConvention, R"doc(

         Latitude convention of the grid.

         :type: AtmosphericTideLatitudeConvention
      )doc" )
            .def_readwrite( "displacement_basis", &io::AtmosphericTideMetadata::displacementBasis, R"doc(

         Local basis of the three displacement components.

         :type: AtmosphericTideDisplacementBasis
      )doc" )
            .def_readwrite( "source_profile", &io::AtmosphericTideMetadata::sourceProfile, R"doc(

         Identified product. Required for the conversion to the centre of mass in the loading model.

         :type: AtmosphericTideSourceProfile
      )doc" )
            .def_readwrite( "expected_longitude_count", &io::AtmosphericTideMetadata::expectedLongitudeCount, R"doc(

         Expected number of grid longitudes (0 if not checked). Since the grid file has no header giving its dimensions,
         setting this value detects a truncated file (361 for the full RP03 grid).

         :type: int
      )doc" )
            .def_readwrite( "expected_latitude_count", &io::AtmosphericTideMetadata::expectedLatitudeCount, R"doc(

         Expected number of grid latitudes (0 if not checked; 181 for the full RP03 grid).

         :type: int
      )doc" )
            .def_readwrite( "original_coefficient_units", &io::AtmosphericTideMetadata::originalCoefficientUnits, R"doc(

         Units of the coefficients in the file (``mm``). Coefficients are stored in metres.

         :type: str
      )doc" )
            .def_readwrite( "original_coefficient_order", &io::AtmosphericTideMetadata::originalCoefficientOrder, R"doc(

         Order of the coefficients in the file.

         :type: str
      )doc" );

    py::class_< io::AtmosphericTideGrid >( module, "AtmosphericTideGrid", R"doc(

Atmospheric S1/S2 loading coefficients on a regular longitude/latitude grid.

The grid is stored as read: both longitude seam columns (0 and 360 degrees) and all polar values are retained, and no
interpolation is performed. Objects of this class are created by :func:`read_atmospheric_tide_grid`, and used by
:func:`~tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_grid`.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "metadata", &io::AtmosphericTideGrid::metadata, R"doc(

         Description of the data source.

         :type: AtmosphericTideMetadata
      )doc" )
            .def_readwrite( "longitudes_degrees", &io::AtmosphericTideGrid::longitudesDegrees, R"doc(

         Grid longitudes in degrees.

         :type: list[float]
      )doc" )
            .def_readwrite( "latitudes_degrees", &io::AtmosphericTideGrid::latitudesDegrees, R"doc(

         Grid latitudes in degrees.

         :type: list[float]
      )doc" )
            .def_readwrite( "coefficients_metres", &io::AtmosphericTideGrid::coefficientsMetres, R"doc(

         Coefficients at each grid point in metres, ordered by longitude and then by latitude. Each matrix has the
         radial, north and east components as rows, and the :math:`\cos S_{1}`, :math:`\sin S_{1}`, :math:`\cos S_{2}`
         and :math:`\sin S_{2}` coefficients as columns.

         :type: list[numpy.ndarray[numpy.float64[3, 4]]]
      )doc" )
            .def( "validate", &io::validateAtmosphericTideGrid, R"doc(

         Check that the grid is complete, regular and finite, and raise an exception if not.

      )doc" );

    py::class_< io::AtmosphericTideStation >( module, "AtmosphericTideStation", R"doc(

Atmospheric S1/S2 loading coefficients of a single station.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "station_identifier", &io::AtmosphericTideStation::stationIdentifier, R"doc(

         Identifier of the station, exactly as given in the file.

         :type: str
      )doc" )
            .def_readwrite( "provider_station_identifier_line", &io::AtmosphericTideStation::providerStationIdentifierLine, R"doc(

         Line of the file that identifies the station, unmodified.

         :type: str
      )doc" )
            .def_readwrite( "coefficients_metres", &io::AtmosphericTideStation::coefficientsMetres, R"doc(

         Coefficients in metres, with the same layout as :attr:`AtmosphericTideGrid.coefficients_metres`.

         :type: numpy.ndarray[numpy.float64[3, 4]]
      )doc" )
            .def_readwrite( "metadata", &io::AtmosphericTideStation::metadata, R"doc(

         Station-specific description of the data source; undeclared entries are taken from the file-level metadata.

         :type: AtmosphericTideMetadata
      )doc" )
            .def_readwrite( "has_provider_coordinates", &io::AtmosphericTideStation::hasProviderCoordinates, R"doc(

         Whether the file gives the coordinates of the station.

         :type: bool
      )doc" )
            .def_readwrite( "provider_longitude_latitude_degrees", &io::AtmosphericTideStation::providerLongitudeLatitudeDegrees, R"doc(

         Longitude and latitude of the station in degrees, as given in the file (only meaningful if
         :attr:`has_provider_coordinates` is true).

         :type: numpy.ndarray[numpy.float64[2, 1]]
      )doc" );

    py::class_< io::AtmosphericTideStationData >( module, "AtmosphericTideStationData", R"doc(

Atmospheric S1/S2 loading coefficients of a set of stations.

Objects of this class are created by :func:`read_atmospheric_tide_station_data`, and used by
:func:`~tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_stations`.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "metadata", &io::AtmosphericTideStationData::metadata, R"doc(

         Description of the data source, common to all stations.

         :type: AtmosphericTideMetadata
      )doc" )
            .def_readwrite( "stations", &io::AtmosphericTideStationData::stations, R"doc(

         Coefficients of each station, with the station identifier as key.

         :type: dict[str, AtmosphericTideStation]
      )doc" )
            .def( "get_station",
                  &io::AtmosphericTideStationData::getStation,
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
         AtmosphericTideStation
             Coefficients of the station.

      )doc" )
            .def( "validate", &io::validateAtmosphericTideStationData, R"doc(

         Check that the data of all stations is complete and finite, and raise an exception if not.

      )doc" );

    module.def( "read_atmospheric_tide_grid",
                py::overload_cast< const std::string&, const io::AtmosphericTideMetadata& >( &io::readAtmosphericTideGrid ),
                py::arg( "file_path" ),
                py::arg( "metadata" ) = io::AtmosphericTideMetadata( ),
                R"doc(

Function to read an atmospheric S1/S2 loading grid.

Function to read an atmospheric S1/S2 loading grid in the 14-column Bernese format (as distributed for the Ray and
Ponte, 2003, model): longitude, latitude and twelve coefficients in millimetres, being the radial, north and east
components of the :math:`\cos S_{1}`, :math:`\sin S_{1}`, :math:`\cos S_{2}` and :math:`\sin S_{2}` terms. The
coefficients are converted to metres; the header lines, model name, Green functions and reference origin found in the
header are stored in the metadata. No interpolation is performed.

Malformed, non-finite, incomplete, duplicated or inconsistent data raise an exception that gives the file and line.


Parameters
----------
file_path : str
    Path of the file that is to be read.
metadata : AtmosphericTideMetadata, default = AtmosphericTideMetadata()
    Declared properties of the data, completed with those found in the file header.

Returns
-------
AtmosphericTideGrid
    Grid read from the file.


Examples
--------
.. code-block:: python

   from tudatpy.data_input.environment_data import atmospheric_loading

   metadata = atmospheric_loading.AtmosphericTideMetadata( )
   metadata.expected_longitude_count = 361
   metadata.expected_latitude_count = 181
   grid = atmospheric_loading.read_atmospheric_tide_grid( "s1_s2_def_ce.dat", metadata )

      )doc" );
    module.def( "read_atmospheric_tide_station_data",
                py::overload_cast< const std::string&, const io::AtmosphericTideMetadata& >( &io::readAtmosphericTideStationData ),
                py::arg( "file_path" ),
                py::arg( "metadata" ) = io::AtmosphericTideMetadata( ),
                R"doc(

Function to read station-specific atmospheric S1/S2 loading coefficients.

Function to read station-specific atmospheric S1/S2 loading coefficients in the Bernese ATL format: for each station,
a line with its identifier, followed by three lines (radial, north, east) of four coefficients in millimetres. The
coefficients are converted to metres.


Parameters
----------
file_path : str
    Path of the file that is to be read.
metadata : AtmosphericTideMetadata, default = AtmosphericTideMetadata()
    Declared properties of the data, completed with those found in the file header.

Returns
-------
AtmosphericTideStationData
    Coefficients of all stations in the file.

      )doc" );
}
}  // namespace atmospheric_loading
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
