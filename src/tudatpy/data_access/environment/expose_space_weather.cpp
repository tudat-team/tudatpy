#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_space_weather.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "tudat/io/solarActivityData.h"

namespace py = pybind11;
namespace tsa = tudat::input_output::solar_activity;

namespace tudatpy
{
namespace data_access
{
namespace environment
{
namespace space_weather
{

void expose_space_weather( py::module& m )
{
    py::class_< tsa::SolarActivityData, std::shared_ptr< tsa::SolarActivityData > >(
            m, "SolarActivityData", R"doc(No documentation available.)doc" )
            .def_readonly( "solar_radio_flux_107_observed", &tsa::SolarActivityData::solarRadioFlux107Observed );

    m.def( "read_solar_activity_data",
           &tsa::readSolarActivityData,
           py::arg( "file_path" ),
           R"doc(
 Reads a space weather data file and produces a dictionary with solar activity data for a range of epochs. Data files can be obtained from http://celestrak.com/SpaceData and should follow the legacy format.

 Parameters
 ----------
 file_path : str
     Path to the space weather data file.
 )doc" );

    py::class_< tsa::SolarActivityContainer, std::shared_ptr< tsa::SolarActivityContainer > >( m, "SolarActivityContainer" )
            .def( py::init< const std::map< double, std::shared_ptr< tsa::SolarActivityData > >& >( ),
                  py::arg( "solar_activity_data_map" ) )
            .def( "get_solar_activity_data",
                  &tsa::SolarActivityContainer::getSolarActivityData,
                  py::arg( "time" ),
                  R"doc(
        Returns the nearest SolarActivityData (in UTC Julian days) for the given time in seconds since J2000.
        )doc" )
            .def( "get_solar_activity_data_map",
                  &tsa::SolarActivityContainer::getSolarActivityDataMap,
                  R"doc(Returns the full map of SolarActivityData.)doc" );
}

}  // namespace space_weather
}  // namespace environment
}  // namespace data_access
}  // namespace tudatpy
