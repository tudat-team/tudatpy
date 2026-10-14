#include "expose_earth_orientation.h"
#include "tudat/io/readEarthOrientationData.h"
#include "tudat/io/readHighFrequencyEopData.h"
#include <pybind11/stl.h>
#include <sstream>
namespace py = pybind11;
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace earth_orientation
{
void expose_earth_orientation( py::module& module )
{
    using namespace tudat::input_output;
    py::class_< HighFrequencyEopTerm >( module, "HighFrequencyEopTerm" )
            .def_readonly( "name", &HighFrequencyEopTerm::name )
            .def_readonly( "doodson_number", &HighFrequencyEopTerm::doodsonNumber )
            .def_readonly( "argument_multipliers", &HighFrequencyEopTerm::argumentMultipliers )
            .def_readonly( "period_days", &HighFrequencyEopTerm::periodDays )
            .def_readonly( "sine_amplitudes", &HighFrequencyEopTerm::sineAmplitudes )
            .def_readonly( "cosine_amplitudes", &HighFrequencyEopTerm::cosineAmplitudes );
    py::class_< HighFrequencyEopData >( module,
                                        "HighFrequencyEopData",
                                        "Raw tidal harmonics and provenance. Amplitudes [xp,yp,UT1,excess LOD] use [rad,rad,s,s]; "
                                        "arguments are [GMST+pi,l,l-prime,F,D,Omega]." )
            .def_readonly( "source", &HighFrequencyEopData::source )
            .def_readonly( "header", &HighFrequencyEopData::header )
            .def_readonly( "terms", &HighFrequencyEopData::terms );
    module.def( "read_high_frequency_eop_data",
                py::overload_cast< const std::string& >( &readHighFrequencyEopData ),
                py::arg( "file_name" ),
                "Read an IVS 17-column harmonic table; normalize microarcseconds/microseconds to SI units." );
    module.def(
            "read_high_frequency_eop_data_text",
            []( const std::string& text, const std::string& source ) {
                std::istringstream stream( text );
                return readHighFrequencyEopData( stream, source );
            },
            py::arg( "text" ),
            py::arg( "source" ) = "<text>" );
    module.def( "get_desai_sibois_2016_eop_data",
                &getDesaiSibois2016EopData,
                "Return the pinned 159-term ocean-tide coefficient table, DOI 10.1002/2016JB013125. "
                "This data does not include libration and does not select a rotation model." );
    py::enum_< EarthOrientationFormat >( module, "EarthOrientationFormat", "Supported raw EOP file formats." )
            .value( "c04", EarthOrientationFormat::c04 )
            .value( "finals2000a", EarthOrientationFormat::finals2000a );
    py::class_< EarthOrientationRecord >(
            module,
            "EarthOrientationRecord",
            "A raw UTC epoch. Angles use radians, rates radians/second, and UT1/LOD seconds. Missing optional values are None." )
            .def_readonly( "utc_modified_julian_date", &EarthOrientationRecord::utcModifiedJulianDate )
            .def_readonly( "polar_motion", &EarthOrientationRecord::polarMotion )
            .def_readonly( "ut1_minus_utc", &EarthOrientationRecord::ut1MinusUtc )
            .def_readonly( "length_of_day", &EarthOrientationRecord::lengthOfDay )
            .def_readonly( "celestial_pole_offsets", &EarthOrientationRecord::celestialPoleOffsets )
            .def_readonly( "polar_motion_rates", &EarthOrientationRecord::polarMotionRates )
            .def_readonly( "polar_motion_uncertainty", &EarthOrientationRecord::polarMotionUncertainty )
            .def_readonly( "ut1_minus_utc_uncertainty", &EarthOrientationRecord::ut1MinusUtcUncertainty )
            .def_readonly( "length_of_day_uncertainty", &EarthOrientationRecord::lengthOfDayUncertainty )
            .def_readonly( "celestial_pole_offset_uncertainty", &EarthOrientationRecord::celestialPoleOffsetUncertainty )
            .def_readonly( "polar_motion_rate_uncertainty", &EarthOrientationRecord::polarMotionRateUncertainty )
            .def_readonly( "polar_motion_status", &EarthOrientationRecord::polarMotionStatus )
            .def_readonly( "ut1_status", &EarthOrientationRecord::ut1Status )
            .def_readonly( "celestial_pole_status", &EarthOrientationRecord::celestialPoleStatus );
    py::class_< EarthOrientationData >(
            module,
            "EarthOrientationData",
            "Un-interpolated EOP records indexed by UTC MJD, with source/header/frame provenance and exact raw coverage." )
            .def_readonly( "format", &EarthOrientationData::format )
            .def_readonly( "source", &EarthOrientationData::source )
            .def_readonly( "format_description", &EarthOrientationData::formatDescription )
            .def_readonly( "header", &EarthOrientationData::header )
            .def_readonly( "terrestrial_frame", &EarthOrientationData::terrestrialFrame )
            .def_readonly( "celestial_frame", &EarthOrientationData::celestialFrame )
            .def_readonly( "celestial_pole_convention", &EarthOrientationData::celestialPoleConvention )
            .def_readonly( "records", &EarthOrientationData::records );
    module.def(
            "read_earth_orientation_data",
            py::overload_cast< const std::string&, EarthOrientationFormat >( &readEarthOrientationData ),
            py::arg( "file_name" ),
            py::arg( "format" ) = EarthOrientationFormat::c04,
            "Read a C04 or finals2000A Bulletin A file. Optional values remain missing; no interpolation or extrapolation is performed." );
    module.def(
            "read_earth_orientation_data_text",
            []( const std::string& text, const EarthOrientationFormat format, const std::string& source ) {
                std::istringstream stream( text );
                return readEarthOrientationData( stream, format, source );
            },
            py::arg( "text" ),
            py::arg( "format" ) = EarthOrientationFormat::c04,
            py::arg( "source" ) = "<text>",
            "Read raw EOP text with the same validation and unit normalization as the file reader." );
}
}  // namespace earth_orientation
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
