#ifndef TUDAT_READ_HIGH_FREQUENCY_EOP_DATA_H
#define TUDAT_READ_HIGH_FREQUENCY_EOP_DATA_H

#include <array>
#include <istream>
#include <string>
#include <vector>

namespace tudat
{
namespace input_output
{

//! One tidal line, with phase multipliers for [GMST+pi, l, l', F, D, Omega].
struct HighFrequencyEopTerm {
    std::string name;
    std::string doodsonNumber;
    std::array< int, 6 > argumentMultipliers{};
    double periodDays = 0.0;
    //! [xp, yp, UT1, excess LOD], in [rad, rad, s, s].
    std::array< double, 4 > sineAmplitudes{};
    std::array< double, 4 > cosineAmplitudes{};
};

//! Raw harmonic coefficients and provenance; no interpolation or dynamics dependencies.
struct HighFrequencyEopData {
    std::string source;
    std::vector< std::string > header;
    std::vector< HighFrequencyEopTerm > terms;
};

//! Validate finite coefficients, positive periods and unique argument multipliers.
void validateHighFrequencyEopData( const HighFrequencyEopData& data );

//! Read the IVS 17-column format. Input amplitudes are microarcseconds/microseconds.
//! Comments and Doodson labels (including X/E digits) are retained without phase conversion.
HighFrequencyEopData readHighFrequencyEopData( std::istream& stream, const std::string& source = "<stream>" );
HighFrequencyEopData readHighFrequencyEopData( const std::string& fileName );

//! Pinned 159-term Desai--Sibois (2016) ocean model, DOI 10.1002/2016JB013125.
//! This contains no libration terms. No network or external resource file is required.
HighFrequencyEopData getDesaiSibois2016EopData( );

}  // namespace input_output
}  // namespace tudat
#endif
