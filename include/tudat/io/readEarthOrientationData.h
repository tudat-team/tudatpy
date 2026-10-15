#ifndef TUDAT_READ_EARTH_ORIENTATION_DATA_H
#define TUDAT_READ_EARTH_ORIENTATION_DATA_H

#include <array>
#include <istream>
#include <map>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace tudat
{
namespace input_output
{

//! File formats supported by the independent, un-interpolated EOP reader.
enum class EarthOrientationFormat { c04, finals2000a };

//! Earth orientation parameters at a single epoch (one line of an EOP file), in SI units.
struct EarthOrientationRecord {
    //! Epoch of the record, as UTC modified Julian date.
    double utcModifiedJulianDate = std::numeric_limits< double >::quiet_NaN( );

    //! Polar motion x_p and y_p (radians).
    std::array< double, 2 > polarMotion = { std::numeric_limits< double >::quiet_NaN( ), std::numeric_limits< double >::quiet_NaN( ) };

    //! UT1 - UTC (seconds).
    double ut1MinusUtc = std::numeric_limits< double >::quiet_NaN( );

    //! Excess length of day (seconds); empty if not in the file.
    std::optional< double > lengthOfDay;

    //! Celestial pole offsets dX and dY (radians); empty if not in the file.
    std::optional< std::array< double, 2 > > celestialPoleOffsets;

    //! Rates of x_p and y_p (radians per second); empty if not in the file.
    std::optional< std::array< double, 2 > > polarMotionRates;

    //! Uncertainties of the values above, in the same units; empty if not in the file.
    std::optional< std::array< double, 2 > > polarMotionUncertainty;
    std::optional< double > ut1MinusUtcUncertainty;
    std::optional< double > lengthOfDayUncertainty;
    std::optional< std::array< double, 2 > > celestialPoleOffsetUncertainty;
    std::optional< std::array< double, 2 > > polarMotionRateUncertainty;

    //! Whether each value is final (I) or predicted (P), for Bulletin A files; empty for C04 files.
    std::string polarMotionStatus, ut1Status, celestialPoleStatus;
};

//! Contents of an EOP file: its records and a description of the file, without interpolation or extrapolation.
struct EarthOrientationData {
    //! Format of the file.
    EarthOrientationFormat format = EarthOrientationFormat::c04;

    //! Name of the file (or other source) that was read.
    std::string source;

    //! Description of the file format, as detected by the reader.
    std::string formatDescription;

    //! Header lines of the file.
    std::vector< std::string > header;

    //! Terrestrial and celestial reference frames of the data, and the convention for the celestial pole offsets.
    std::string terrestrialFrame = "ITRS";
    std::string celestialFrame = "GCRS";
    std::string celestialPoleConvention = "dX,dY";

    //! Records, with their UTC modified Julian date as key.
    std::map< double, EarthOrientationRecord > records;
};

//! Validate public C++ data before a model consumes it; does not change coverage.
void validateEarthOrientationData( const EarthOrientationData& data );

//! Read legacy 16-column / modern 21-column C04 or Bulletin A finals2000A.
//! Throws with source and line number on malformed, empty, unsorted or duplicate data.
EarthOrientationData readEarthOrientationData( std::istream& stream,
                                               EarthOrientationFormat format,
                                               const std::string& source = "<stream>" );
EarthOrientationData readEarthOrientationData( const std::string& fileName, EarthOrientationFormat format = EarthOrientationFormat::c04 );

}  // namespace input_output
}  // namespace tudat
#endif
