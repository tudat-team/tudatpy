/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_READ_ICGEM_GRAVITY_FIELD_H
#define TUDAT_READ_ICGEM_GRAVITY_FIELD_H
#include <array>
#include <cstddef>
#include <istream>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "tudat/astro/basic_astro/tideSystem.h"

namespace tudat
{
namespace input_output
{
enum class IcgemNormalization { fullyNormalized, unnormalized };
enum class IcgemCoefficientErrors { none, calibrated, formal, calibratedAndFormal };
//! Published values and optional sigmaC/sigmaS pairs, in the declared normalization and record units.
struct IcgemGravityCoefficient {
    //! Static/reference values, trends or periodic amplitudes as specified by the containing record.
    double cosine = 0.0;

    double sine = 0.0;

    //! Declared calibrated/formal standard deviations; absence is distinct from a supplied zero.
    std::optional< std::array< double, 2 > > calibratedErrors;

    std::optional< std::array< double, 2 > > formalErrors;

    //! One-based line and original text provide the source provenance already used by static records.
    std::size_t sourceLine = 0;

    std::string originalRecord;

    //! Compare supplied values, uncertainty columns and source provenance exactly.
    bool operator==( const IcgemGravityCoefficient& other ) const
    {
        return std::tie( cosine, sine, calibratedErrors, formalErrors, sourceLine, originalRecord ) ==
                std::tie( other.cosine, other.sine, other.calibratedErrors, other.formalErrors, other.sourceLine, other.originalRecord );
    }
};

//! One format-defined gfct/trnd/dot/asin/acos record, retaining its independent dates and units.
struct IcgemTemporalCoefficient {
    //! Original keyword and spherical-harmonic degree/order.
    std::string type;

    int degree = 0, order = 0;

    //! Supplied values, uncertainty columns and source provenance, without normalization or evaluation.
    IcgemGravityCoefficient coefficient;

    //! Printed dates, left uninterpreted as a physical time scale or year convention.
    //! In 1.0 only gfct supplies referenceDate; the other records inherit it when consumed.
    std::string referenceDate;

    std::string validityEndDate;

    //! Supplied period in years, present only for asin/acos.
    std::optional< double > periodYears;

    //! Compare the independent reference/validity tokens as well as the numerical record.
    bool operator==( const IcgemTemporalCoefficient& other ) const
    {
        return std::tie( type, degree, order, coefficient, referenceDate, validityEndDate, periodYears ) ==
                std::tie( other.type,
                          other.degree,
                          other.order,
                          other.coefficient,
                          other.referenceDate,
                          other.validityEndDate,
                          other.periodYears );
    }
};

//! Neutral sparse static and temporal records. Missing entries remain distinct from published zeros.
struct IcgemGravityFieldData {
    //! Source identity, untouched header/comment text and declared product/model/format names.
    std::string source;

    std::vector< std::string > headerLines;

    std::string productType;

    std::string modelName;

    std::string format = "icgem1.0";

    //! Published SI constants and maximum expansion degree.
    double gravitationalParameter = 0.0;

    double referenceRadius = 0.0;

    int maximumDegree = 0;

    //! Error-column declaration applies equally to static and temporal records.
    IcgemCoefficientErrors errors = IcgemCoefficientErrors::none;

    //! These defaults are defined by the 2023 format, not inferred from coefficients.
    IcgemNormalization normalization = IcgemNormalization::fullyNormalized;

    basic_astrodynamics::TideSystem tideSystem = basic_astrodynamics::TideSystem::unknown;

    //! Static gfc values and ordered temporal records; temporal baselines are never inserted into this static map.
    std::map< std::pair< int, int >, IcgemGravityCoefficient > coefficients;

    std::vector< IcgemTemporalCoefficient > temporalCoefficients;

    //! Compare complete source records so archived consuming settings retain their original data.
    bool operator==( const IcgemGravityFieldData& other ) const
    {
        return std::tie( source,
                         headerLines,
                         productType,
                         modelName,
                         format,
                         gravitationalParameter,
                         referenceRadius,
                         maximumDegree,
                         errors,
                         normalization,
                         tideSystem,
                         coefficients,
                         temporalCoefficients ) ==
                std::tie( other.source,
                          other.headerLines,
                          other.productType,
                          other.modelName,
                          other.format,
                          other.gravitationalParameter,
                          other.referenceRadius,
                          other.maximumDegree,
                          other.errors,
                          other.normalization,
                          other.tideSystem,
                          other.coefficients,
                          other.temporalCoefficients );
    }
};

//! Validate header conventions, coefficient values and independent temporal record coverage.
void validateIcgemGravityFieldData( const IcgemGravityFieldData& data );

//! Read ICGEM 1.0/2.0 records without evaluating temporal terms; unknown keywords and trailing text are comments.
IcgemGravityFieldData readIcgemGravityField( std::istream& stream, const std::string& source = "<stream>" );

//! Open a gravity product and retain its filename for record provenance.
IcgemGravityFieldData readIcgemGravityField( const std::string& fileName );

}  // namespace input_output
}  // namespace tudat
#endif
