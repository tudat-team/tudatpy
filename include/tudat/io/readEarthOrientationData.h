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

//! An observed or predicted UTC epoch. All quantities are normalized to SI units.
struct EarthOrientationRecord {
    double utcModifiedJulianDate = std::numeric_limits< double >::quiet_NaN( );
    std::array< double, 2 > polarMotion = { std::numeric_limits< double >::quiet_NaN( ),
                                            std::numeric_limits< double >::quiet_NaN( ) };  // xp, yp [rad]
    double ut1MinusUtc = std::numeric_limits< double >::quiet_NaN( );                       // [s]
    std::optional< double > lengthOfDay;                                                    // excess LOD [s]
    std::optional< std::array< double, 2 > > celestialPoleOffsets;                          // dX, dY [rad]
    std::optional< std::array< double, 2 > > polarMotionRates;                              // [rad/s]
    std::optional< std::array< double, 2 > > polarMotionUncertainty;
    std::optional< double > ut1MinusUtcUncertainty;
    std::optional< double > lengthOfDayUncertainty;
    std::optional< std::array< double, 2 > > celestialPoleOffsetUncertainty;
    std::optional< std::array< double, 2 > > polarMotionRateUncertainty;
    std::string polarMotionStatus, ut1Status, celestialPoleStatus;  // I/P for Bulletin A; empty for C04
};

//! Raw data and provenance only: no model, interpolation, extrapolation or padding.
struct EarthOrientationData {
    EarthOrientationFormat format = EarthOrientationFormat::c04;
    std::string source;
    std::string formatDescription;
    std::vector< std::string > header;
    std::string terrestrialFrame = "ITRS";
    std::string celestialFrame = "GCRS";
    std::string celestialPoleConvention = "dX,dY";
    std::map< double, EarthOrientationRecord > records;  // exact UTC MJD coverage
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
