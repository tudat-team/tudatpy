#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_ilrs.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "tudat/io/readSinexFile.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace ilrs
{

void expose_ilrs( py::module& m )
{
    py::class_< tio::SinexStationState >( m, "SinexStationState", R"doc(
Container for station state data parsed from a SINEX file.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "site_code", &tio::SinexStationState::siteCode_ )
            .def_readwrite( "domes_id", &tio::SinexStationState::domesId_ )
            .def_readwrite( "position", &tio::SinexStationState::position_ )
            .def_readwrite( "velocity", &tio::SinexStationState::velocity_ )
            .def_readwrite( "reference_epoch", &tio::SinexStationState::referenceEpoch_ );

    py::class_< tio::SinexStationEccentricity >( m, "SinexStationEccentricity", R"doc(
Container for SINEX station eccentricity data.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "domes_id", &tio::SinexStationEccentricity::domesId_ )
            .def_readwrite( "station_code", &tio::SinexStationEccentricity::stationCode_ )
            .def_readwrite( "eccentricity", &tio::SinexStationEccentricity::eccentricity_ )
            .def_readwrite( "start_epoch", &tio::SinexStationEccentricity::startEpoch_ )
            .def_readwrite( "end_epoch", &tio::SinexStationEccentricity::endEpoch_ )
            .def_readwrite( "has_open_end", &tio::SinexStationEccentricity::hasOpenEnd_ );

    py::class_< tio::IlrsStationRegistryEntry >( m, "IlrsStationRegistryEntry", R"doc(
Container for one ILRS station registry entry parsed from SINEX ``SITE/ID``.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "station_code", &tio::IlrsStationRegistryEntry::stationCode_ )
            .def_readwrite( "station_name", &tio::IlrsStationRegistryEntry::stationName_ )
            .def_readwrite( "domes_id", &tio::IlrsStationRegistryEntry::domesId_ )
            .def_readwrite( "approximate_longitude", &tio::IlrsStationRegistryEntry::approximateLongitude_ )
            .def_readwrite( "approximate_latitude", &tio::IlrsStationRegistryEntry::approximateLatitude_ )
            .def_readwrite( "approximate_height", &tio::IlrsStationRegistryEntry::approximateHeight_ );

    m.def( "convert_sinex_datetime_to_seconds_since_epoch",
           &tio::convertSinexDateTimeToSecondsSinceEpoch,
           py::arg( "date_time" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Convert a SINEX epoch string (``YY:DDD:SSSSS``) to seconds since a reference Julian day.
           )doc" );

    m.def( "read_sinex_station_data",
           &tio::readSinexStationData,
           py::arg( "file_name" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Read station positions, velocities and reference epochs from a SINEX file.
           )doc" );

    m.def( "read_sinex_station_eccentricities",
           &tio::readSinexStationEccentricities,
           py::arg( "file_name" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Read station eccentricity history from a SINEX file.
           )doc" );

    m.def( "read_ilrs_station_registry_from_sinex_site_id",
           &tio::readIlrsStationRegistryFromSinexSiteId,
           py::arg( "file_name" ),
           R"doc(
Read the ILRS station registry from a SINEX ``SITE/ID`` block as ``{station_code: entry}``.
           )doc" );

    m.def( "read_domes_id_numbers",
           &tio::readDomesIdNumbers,
           py::arg( "file_name" ),
           R"doc(
Read a mapping from station name to DOMES id.
           )doc" );

    m.def( "read_monument_numbers",
           &tio::readMonumentNumbers,
           py::arg( "file_name" ),
           R"doc(
Read a mapping from monument/station code to station name.
           )doc" );

    m.def( "read_ground_station_names",
           &tio::readGroundStationNames,
           py::arg( "file_name" ),
           R"doc(
Read a mapping from DOMES id to station name.
           )doc" );
}

}  // namespace ilrs
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
