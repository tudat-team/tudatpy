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
namespace data_input
{
namespace environment_data
{
namespace space_weather
{

void expose_space_weather( py::module& m )
{
    py::class_< tsa::SolarActivityData, std::shared_ptr< tsa::SolarActivityData > >( m, "SolarActivityData", R"doc(
         Solar-activity data for a single space-weather epoch.

         Attributes
         ----------
         solar_radio_flux_107_observed : float
             Observed 10.7 cm solar radio flux value.
      )doc" )
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

         Returns
         -------
         dict[float, SolarActivityData]
             Solar activity data indexed by UTC Julian day.
      )doc" );

    py::class_< tsa::SolarActivityContainer, std::shared_ptr< tsa::SolarActivityContainer > >( m,
                                                                                               "SolarActivityContainer",
                                                                                               R"doc(
         Container for time-indexed solar-activity data.

         Parameters
         ----------
         solar_activity_data_map : dict[float, SolarActivityData]
             Solar activity data indexed by UTC Julian day.
      )doc" )
            .def( py::init< const std::map< double, std::shared_ptr< tsa::SolarActivityData > >& >( ),
                  py::arg( "solar_activity_data_map" ),
                  R"doc(
         Create a solar-activity container.

         Parameters
         ----------
         solar_activity_data_map : dict[float, SolarActivityData]
             Solar activity data indexed by UTC Julian day.
      )doc" )
            .def( "get_solar_activity_data",
                  &tsa::SolarActivityContainer::getSolarActivityData,
                  py::arg( "time" ),
                  R"doc(
         Returns the nearest solar-activity data entry.

         Parameters
         ----------
         time : float
             Seconds since J2000.

         Returns
         -------
         SolarActivityData
             Nearest solar-activity data entry.
      )doc" )
            .def( "get_solar_activity_data_map",
                  &tsa::SolarActivityContainer::getSolarActivityDataMap,
                  R"doc(
         Return all solar-activity data.

         Returns
         -------
         dict[float, SolarActivityData]
             Solar activity data indexed by UTC Julian day.
      )doc" );
}

}  // namespace space_weather
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
