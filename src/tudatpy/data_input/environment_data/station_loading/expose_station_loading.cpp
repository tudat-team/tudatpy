#include "expose_station_loading.h"
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
    py::enum_< io::StationLoadingTimeScale >( module, "StationLoadingTimeScale", R"doc(

Enumeration of time scales in which the epochs of a station loading series can be given.

      )doc" )
            .value( "tai", io::StationLoadingTimeScale::tai, "International Atomic Time (used by EPHEDISP files)." )
            .value( "tt", io::StationLoadingTimeScale::tt, "Terrestrial Time." )
            .value( "utc", io::StationLoadingTimeScale::utc, "Coordinated Universal Time." )
            .value( "tdb", io::StationLoadingTimeScale::tdb, "Barycentric Dynamical Time." );
    py::enum_< io::StationLoadingBasis >( module, "StationLoadingBasis", R"doc(

Enumeration of bases in which the displacements of a station loading series can be given.

      )doc" )
            .value( "unspecified", io::StationLoadingBasis::unspecified, "Basis not declared." )
            .value( "geocentric_radial_east_north",
                    io::StationLoadingBasis::geocentricRadialEastNorth,
                    "Geocentric radial, east and north components at the provider's site location (used by EPHEDISP files)." )
            .value( "cartesian", io::StationLoadingBasis::cartesian, "Body-fixed Cartesian components." );
    py::enum_< io::StationLoadingComponent >( module, "StationLoadingComponent", R"doc(

Enumeration of physical loading components that a station loading series can contain.

      )doc" )
            .value( "atmosphere", io::StationLoadingComponent::atmosphere, "Atmospheric pressure loading." )
            .value( "hydrology", io::StationLoadingComponent::hydrology, "Continental water storage loading." )
            .value( "ocean", io::StationLoadingComponent::ocean, "Ocean bottom pressure loading." );
    py::enum_< io::StationLoadingSignalPresence >( module, "StationLoadingSignalPresence", R"doc(

Enumeration describing whether a signal (for instance the diurnal tide) is contained in a station loading series.

      )doc" )
            .value( "unknown", io::StationLoadingSignalPresence::unknown, "Not documented by the provider." )
            .value( "absent", io::StationLoadingSignalPresence::absent, "Removed from, or not modelled in, the series." )
            .value( "present", io::StationLoadingSignalPresence::present, "Contained in the series." );
    py::enum_< io::StationLoadingOrigin >( module, "StationLoadingOrigin", R"doc(

Enumeration of reference origins with respect to which station loading displacements can be given.

      )doc" )
            .value( "unspecified", io::StationLoadingOrigin::unspecified, "Origin not declared." )
            .value( "whole_earth_centre_of_mass",
                    io::StationLoadingOrigin::wholeEarthCentreOfMass,
                    "Centre of mass of the whole Earth, including the loading masses (CM)." )
            .value( "solid_earth_centre_of_mass", io::StationLoadingOrigin::solidEarthCentreOfMass, "Centre of the solid Earth (CE)." )
            .value( "centre_of_figure", io::StationLoadingOrigin::centreOfFigure, "Centre of figure of the Earth's surface (CF)." )
            .value( "custom", io::StationLoadingOrigin::custom, "Other origin, described in the origin definition." );

    py::class_< io::StationLoadingEpoch >( module, "StationLoadingEpoch", R"doc(

Epoch of a station loading series, stored as an integer modified Julian day and the seconds of that day (in the time
scale of the series) to retain the full precision of the source.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "modified_julian_day", &io::StationLoadingEpoch::modifiedJulianDay, R"doc(

         Integer modified Julian day.

         :type: int
      )doc" )
            .def_readwrite( "seconds_of_day", &io::StationLoadingEpoch::secondsOfDay, R"doc(

         Seconds since the start of the day.

         :type: float
      )doc" );

    py::class_< io::StationLoadingMetadata >( module, "StationLoadingMetadata", R"doc(

Description of a station loading series: its product, physical content and reference system.

The physical content cannot be read from EPHEDISP files, and must be supplied by the user (or by a provider-specific
function such as :func:`imls_merra2_atmospheric_loading_metadata`). Unknown entries remain unspecified; they are never
inferred from file names or labels.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::StationLoadingMetadata::source, R"doc(

         Name of the file (or other source) of the data.

         :type: str
      )doc" )
            .def_readwrite( "source_headers", &io::StationLoadingMetadata::sourceHeaders, R"doc(

         Header lines of the source, unmodified.

         :type: list[str]
      )doc" )
            .def_readwrite( "product", &io::StationLoadingMetadata::product, R"doc(

         Name of the loading product.

         :type: str
      )doc" )
            .def_readwrite( "release", &io::StationLoadingMetadata::release, R"doc(

         Release (version) of the loading product.

         :type: str
      )doc" )
            .def_readwrite( "source_profile", &io::StationLoadingMetadata::sourceProfile, R"doc(

         Identifier of the documented provider profile from which the metadata was set (empty if supplied by the user).

         :type: str
      )doc" )
            .def_readwrite( "reference_frame", &io::StationLoadingMetadata::referenceFrame, R"doc(

         Terrestrial reference frame of the provider's site locations.

         :type: str
      )doc" )
            .def_readwrite( "reference_origin", &io::StationLoadingMetadata::referenceOrigin, R"doc(

         Reference origin of the displacements.

         :type: StationLoadingOrigin
      )doc" )
            .def_readwrite( "origin_definition", &io::StationLoadingMetadata::originDefinition, R"doc(

         Description of the reference origin, as documented by the provider.

         :type: str
      )doc" )
            .def_readwrite( "ocean_response", &io::StationLoadingMetadata::oceanResponse, R"doc(

         Assumed ocean response to atmospheric pressure (for instance inverted barometer), as documented by the provider.

         :type: str
      )doc" )
            .def_readwrite( "components", &io::StationLoadingMetadata::components, R"doc(

         Physical loading components contained in the series.

         :type: list[StationLoadingComponent]
      )doc" )
            .def_readwrite( "diurnal", &io::StationLoadingMetadata::diurnal, R"doc(

         Whether diurnal tidal signals are contained in the series.

         :type: StationLoadingSignalPresence
      )doc" )
            .def_readwrite( "semidiurnal", &io::StationLoadingMetadata::semidiurnal, R"doc(

         Whether semi-diurnal tidal signals are contained in the series.

         :type: StationLoadingSignalPresence
      )doc" )
            .def_readwrite( "annual", &io::StationLoadingMetadata::annual, R"doc(

         Whether annual signals are contained in the series.

         :type: StationLoadingSignalPresence
      )doc" )
            .def_readwrite( "semiannual", &io::StationLoadingMetadata::semiannual, R"doc(

         Whether semi-annual signals are contained in the series.

         :type: StationLoadingSignalPresence
      )doc" )
            .def_readwrite( "mean", &io::StationLoadingMetadata::mean, R"doc(

         Whether a mean displacement is contained in the series.

         :type: StationLoadingSignalPresence
      )doc" )
            .def_readwrite( "trend", &io::StationLoadingMetadata::trend, R"doc(

         Whether a linear trend is contained in the series.

         :type: StationLoadingSignalPresence
      )doc" )
            .def_readwrite( "removed_harmonics", &io::StationLoadingMetadata::removedHarmonics, R"doc(

         Names of the harmonic signals that the provider removed from the series.

         :type: list[str]
      )doc" )
            .def_readwrite( "original_displacement_units", &io::StationLoadingMetadata::originalDisplacementUnits, R"doc(

         Units of the displacements in the source. Displacements are stored in metres.

         :type: str
      )doc" );

    py::class_< io::StationLoadingSample >( module, "StationLoadingSample", R"doc(

Displacement of a single site at a single epoch of a station loading series.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "displacement", &io::StationLoadingSample::displacement, R"doc(

         Displacement in metres, in the basis of the series (see :attr:`StationLoadingData.basis`).

         :type: list[float]
      )doc" )
            .def_readwrite( "source_line", &io::StationLoadingSample::sourceLine, R"doc(

         Line number of the sample in the source file (0 if not read from a file).

         :type: int
      )doc" )
            .def_readwrite( "source_record", &io::StationLoadingSample::sourceRecord, R"doc(

         Line of the source file that gives the sample, unmodified.

         :type: str
      )doc" )
            .def_readwrite( "informational_epoch", &io::StationLoadingSample::informationalEpoch, R"doc(

         Calendar epoch given in the record, as text. It is not used: the epoch of a sample is defined by its index.

         :type: str
      )doc" );

    py::class_< io::StationLoadingSite >( module, "StationLoadingSite", R"doc(

Site of a station loading series, with its displacement samples.

The site identifier is the provider's label; it is not used to match the site to a ground station. The location and
validity radius are used for that instead (see
:class:`~tudatpy.dynamics.environment_setup.shape_deformation.StationLoadingEvaluationSettings`).

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "identifier", &io::StationLoadingSite::identifier, R"doc(

         Provider's label of the site.

         :type: str
      )doc" )
            .def_readwrite( "cartesian_position", &io::StationLoadingSite::cartesianPosition, R"doc(

         Body-fixed Cartesian position of the site in metres, as given by the provider.

         :type: list[float]
      )doc" )
            .def_readwrite( "validity_radius_metres", &io::StationLoadingSite::validityRadiusMetres, R"doc(

         Radius around the site location within which its displacements may be used, in metres.

         :type: float
      )doc" )
            .def_readwrite( "source_line", &io::StationLoadingSite::sourceLine, R"doc(

         Line number of the site record in the source file (0 if not read from a file).

         :type: int
      )doc" )
            .def_readwrite( "source_record", &io::StationLoadingSite::sourceRecord, R"doc(

         Line of the source file that defines the site, unmodified.

         :type: str
      )doc" )
            .def_readwrite( "informational_coordinates", &io::StationLoadingSite::informationalCoordinates, R"doc(

         Geographic coordinates given in the site record, as text (not used).

         :type: str
      )doc" )
            .def_readwrite( "samples", &io::StationLoadingSite::samples, R"doc(

         Displacement samples, with the (zero-based) index into :attr:`StationLoadingData.epochs` as key. Samples may be
         missing at the start and end of the series, but not in between.

         :type: dict[int, StationLoadingSample]
      )doc" );

    py::class_< io::StationLoadingData >( module, "StationLoadingData", R"doc(

Station loading series: displacements of a set of sites on a common grid of epochs.

Objects of this class are created by :func:`read_ephedisp_station_loading_data`, or can be filled by the user for
other sources (after which :meth:`validate` should be called). They are used by
:func:`~tudatpy.dynamics.environment_setup.shape_deformation.station_loading`.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "metadata", &io::StationLoadingData::metadata, R"doc(

         Description of the series.

         :type: StationLoadingMetadata
      )doc" )
            .def_readwrite( "time_scale", &io::StationLoadingData::timeScale, R"doc(

         Time scale of the epochs.

         :type: StationLoadingTimeScale
      )doc" )
            .def_readwrite( "basis", &io::StationLoadingData::basis, R"doc(

         Basis of the displacements.

         :type: StationLoadingBasis
      )doc" )
            .def_readwrite( "epochs", &io::StationLoadingData::epochs, R"doc(

         Epochs of the series, strictly increasing.

         :type: list[StationLoadingEpoch]
      )doc" )
            .def_readwrite( "declared_begin_epoch", &io::StationLoadingData::declaredBeginEpoch, R"doc(

         First epoch as declared in the file header, or ``None``.

         :type: StationLoadingEpoch | None
      )doc" )
            .def_readwrite( "declared_end_epoch", &io::StationLoadingData::declaredEndEpoch, R"doc(

         Last epoch as declared in the file header, or ``None``.

         :type: StationLoadingEpoch | None
      )doc" )
            .def_readwrite( "nominal_sample_interval_seconds", &io::StationLoadingData::nominalSampleIntervalSeconds, R"doc(

         Nominal interval between epochs in seconds, or ``None`` for an irregular grid.

         :type: float | None
      )doc" )
            .def_readwrite( "sites", &io::StationLoadingData::sites, R"doc(

         Sites of the series, with the site identifier as key.

         :type: dict[str, StationLoadingSite]
      )doc" )
            .def_readwrite( "informational_begin_epoch", &io::StationLoadingData::informationalBeginEpoch, R"doc(

         Calendar begin epoch given in the file header, as text (not used).

         :type: str
      )doc" )
            .def_readwrite( "informational_end_epoch", &io::StationLoadingData::informationalEndEpoch, R"doc(

         Calendar end epoch given in the file header, as text (not used).

         :type: str
      )doc" )
            .def( "validate", &io::validateStationLoadingData, R"doc(

         Check that the epochs are strictly increasing, that the samples are finite and contiguous, and that the sites
         are consistent, and raise an exception if not.

      )doc" );

    module.def( "read_ephedisp_station_loading_data",
                py::overload_cast< const std::string&, const io::StationLoadingMetadata& >( &io::readEphedispStationLoadingData ),
                py::arg( "file_path" ),
                py::arg( "metadata" ) = io::StationLoadingMetadata( ),
                R"doc(

Function to read a station loading series from an EPHEDISP file.

Function to read a station loading series from a file in the `EPHEDISP <https://astrogeo.org/malo/ephedisp_format.txt>`_
format (2005.06.30 version), as distributed for instance by the International Mass Loading Service. The epochs are in
TAI, and are defined by the first epoch and the sampling interval in the file. The displacements are in metres, as
geocentric radial, east and north components at the site locations given in the file.

The physical content of the series (loading components, removed signals, reference origin) is not given in the file,
and should be supplied through the ``metadata`` argument.


Parameters
----------
file_path : str
    Path of the file that is to be read.
metadata : StationLoadingMetadata, default = StationLoadingMetadata()
    Description of the series, completed with the information found in the file.

Returns
-------
StationLoadingData
    Series read from the file.


Examples
--------
.. code-block:: python

   from tudatpy.data_input.environment_data import station_loading

   metadata = station_loading.imls_merra2_atmospheric_loading_metadata( release="<provider release>" )
   data = station_loading.read_ephedisp_station_loading_data( "provider.eph", metadata=metadata )

      )doc" );
    module.def( "imls_merra2_atmospheric_loading_metadata",
                &io::imlsMerra2AtmosphericLoadingMetadata,
                py::arg( "release" ),
                R"doc(

Function to create the description of the IMLS MERRA2 atmospheric loading product.

Function to create the description (loading components, signal content and reference origin) of the atmospheric
pressure loading series based on MERRA2 distributed by the International Mass Loading Service, as documented by the
provider. Its displacements are given with respect to the centre of mass of the whole Earth, including the atmosphere.


Parameters
----------
release : str
    Release identifier of the product, which must be supplied by the user.

Returns
-------
StationLoadingMetadata
    Description of the product.

      )doc" );
}
}  // namespace station_loading
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
