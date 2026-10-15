/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_READ_AOD1B_GRAVITY_FIELD_H
#define TUDAT_READ_AOD1B_GRAVITY_FIELD_H

#include <array>
#include <cstddef>
#include <istream>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace tudat
{
namespace input_output
{

//! Published non-tidal mass components; GLO is ATM plus OCN, whereas OBA has a different domain.
enum class Aod1bComponent { atmosphere, ocean, global, oceanBottomPressure };

//! A supplied coefficient and its location in the original source, without duplicating line text.
struct Aod1bGravityCoefficient {
    //! Dimensionless Cnm/Snm in the header's declared normalization.
    std::array< double, 2 > values = { 0.0, 0.0 };

    //! One-based input line for provenance and diagnostics.
    std::size_t sourceLine = 0;

    //! Compare the supplied values and their recorded provenance.
    bool operator==( const Aod1bGravityCoefficient& other ) const
    {
        return values == other.values && sourceLine == other.sourceLine;
    }
};

//! One component at a printed GPS epoch; missing coefficients remain distinguishable from zero.
struct Aod1bGravitySample {
    //! Physical component declared by this DATA SET record.
    Aod1bComponent component = Aod1bComponent::global;

    //! Verbatim date/time token and its continuous GPS seconds from 2000-01-01 12:00:00 GPS.
    std::string gpsDate;

    double gpsSecondsSinceJ2000 = 0.0;

    //! Sparse degree/order records, including the published degree-zero/one values.
    std::map< std::pair< int, int >, Aod1bGravityCoefficient > coefficients;

    //! Compare one source component, epoch and all supplied coefficient records.
    bool operator==( const Aod1bGravitySample& other ) const
    {
        return std::tie( component, gpsDate, gpsSecondsSinceJ2000, coefficients ) ==
                std::tie( other.component, other.gpsDate, other.gpsSecondsSinceJ2000, other.coefficients );
    }
};

//! Neutral RL07 potential anomalies relative to the provider's 2007--2014 mean.
struct Aod1bGravityFieldData {
    //! Source identifier and complete original header, including END OF HEADER.
    std::string source;

    std::vector< std::string > headerLines;

    //! Published release, model constants in SI, normalization and expansion limit.
    int release = 0;

    std::string softwareVersion;

    double gravitationalParameter = 0.0;

    double referenceRadius = 0.0;

    bool fullyNormalized = true;

    int maximumDegree = 0;

    //! Original component/epoch records; no normalization or component sums are applied here.
    std::vector< Aod1bGravitySample > samples;

    //! Compare both the scientific records and the retained source metadata.
    bool operator==( const Aod1bGravityFieldData& other ) const
    {
        return std::tie( source,
                         headerLines,
                         release,
                         softwareVersion,
                         gravitationalParameter,
                         referenceRadius,
                         fullyNormalized,
                         maximumDegree,
                         samples ) ==
                std::tie( other.source,
                          other.headerLines,
                          other.release,
                          other.softwareVersion,
                          other.gravitationalParameter,
                          other.referenceRadius,
                          other.fullyNormalized,
                          other.maximumDegree,
                          other.samples );
    }
};

//! Check supplied RL07 metadata, GPS dates, coefficient presence and finite values.
void validateAod1bGravityFieldData( const Aod1bGravityFieldData& data );

//! Read the ASCII RL07 format; gzip decompression remains the caller's responsibility.
Aod1bGravityFieldData readAod1bGravityField( std::istream& stream, const std::string& source = "<stream>" );

//! Open one ASCII product and retain its filename as the source identifier.
Aod1bGravityFieldData readAod1bGravityField( const std::string& fileName );

//! Read successive daily files without merging or copying their raw coefficient records.
std::vector< Aod1bGravityFieldData > readAod1bGravityFields( const std::vector< std::string >& fileNames );

}  // namespace input_output
}  // namespace tudat
#endif
