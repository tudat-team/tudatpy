/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_CREATE_ICGEM_TEMPORAL_GRAVITY_FIELD_VARIATIONS_H
#define TUDAT_CREATE_ICGEM_TEMPORAL_GRAVITY_FIELD_VARIATIONS_H

#include <cereal/types/array.hpp>
#include <cereal/types/optional.hpp>
#include "tudat/astro/gravitation/icgemTemporalGravityFieldVariations.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"

namespace tudat
{
namespace input_output
{

//! Archive published coefficient values, uncertainty columns and existing record provenance.
template< class Archive >
void serialize( Archive& archive, IcgemGravityCoefficient& coefficient )
{
    archive( coefficient.cosine,
             coefficient.sine,
             coefficient.calibratedErrors,
             coefficient.formalErrors,
             coefficient.sourceLine,
             coefficient.originalRecord );
}

//! Archive independent temporal dates and periods without choosing their consuming interpretation.
template< class Archive >
void serialize( Archive& archive, IcgemTemporalCoefficient& record )
{
    archive( record.type,
             record.degree,
             record.order,
             record.coefficient,
             record.referenceDate,
             record.validityEndDate,
             record.periodYears );
}

//! Archive the raw product once; restored settings do not reopen provider files.
template< class Archive >
void serialize( Archive& archive, IcgemGravityFieldData& data )
{
    archive( data.source,
             data.headerLines,
             data.productType,
             data.modelName,
             data.format,
             data.gravitationalParameter,
             data.referenceRadius,
             data.maximumDegree,
             data.errors,
             data.normalization,
             data.tideSystem,
             data.coefficients,
             data.temporalCoefficients );
}

}  // namespace input_output
namespace simulation_setup
{

//! Select temporal ICGEM coefficients with explicit source time/year conventions and a nominal-field baseline.
//! Evaluation is analytic: format 2.0 uses declared half-open gfct coverage and format 1.0 remains unbounded.
//! A second interpolation stage is not supported.
class IcgemTemporalGravityFieldVariationSettings : public GravityFieldVariationSettings
{
public:
    //! Retain preloaded data and optional exact-token overrides in source-scale seconds since J2000.
    IcgemTemporalGravityFieldVariationSettings( input_output::IcgemGravityFieldData data,
                                                basic_astrodynamics::TimeScales sourceTimeScale,
                                                gravitation::IcgemYearConvention yearConvention,
                                                int maximumDegree,
                                                int maximumOrder,
                                                std::map< std::string, double > dateOverrides = {} );

    //! Inspect the supplied static/temporal records and their declared source conventions.
    const input_output::IcgemGravityFieldData& getData( ) const
    {
        return data_;
    }

    //! Retrieve the explicitly selected time scale for interpreting provider date tokens.
    basic_astrodynamics::TimeScales getSourceTimeScale( ) const
    {
        return sourceTimeScale_;
    }

    //! Retrieve the explicitly selected definition of trend/periodic years.
    gravitation::IcgemYearConvention getYearConvention( ) const
    {
        return yearConvention_;
    }

    //! Retrieve the inclusive highest selected degree; the block begins at degree two.
    int getMaximumDegree( ) const
    {
        return maximumDegree_;
    }

    //! Retrieve the inclusive highest selected order; the block begins at order zero.
    int getMaximumOrder( ) const
    {
        return maximumOrder_;
    }

    //! Inspect the caller's interpretation of exact irregular date tokens, in source-scale seconds since J2000.
    const std::map< std::string, double >& getDateOverrides( ) const
    {
        return dateOverrides_;
    }

    //! Return a stable product label through the existing variation-identifier lookup.
    std::string getIdentifier( ) const;

    //! Compare complete preloaded records and the physical interpretation chosen for their consumption.
    bool equals( const GravityFieldVariationSettings& other ) const override;

private:
    //! Preloaded source product and its explicitly chosen time/year interpretation.
    input_output::IcgemGravityFieldData data_;

    basic_astrodynamics::TimeScales sourceTimeScale_ = basic_astrodynamics::tai_scale;

    gravitation::IcgemYearConvention yearConvention_ = gravitation::IcgemYearConvention::julianYear;

    //! Inclusive upper bounds and caller-supplied interpretations of irregular printed date tokens.
    int maximumDegree_ = 2, maximumOrder_ = 0;

    std::map< std::string, double > dateOverrides_;

    friend class cereal::access;

    //! Provide an empty archive destination; a usable runtime model validates all restored settings.
    IcgemTemporalGravityFieldVariationSettings( ): GravityFieldVariationSettings( gravitation::icgem_temporal ) {}

    //! Save the complete construction inputs without reopening external coefficient files.
    template< class Archive >
    void save( Archive& archive ) const
    {
        archive( cereal::base_class< GravityFieldVariationSettings >( this ),
                 data_,
                 sourceTimeScale_,
                 yearConvention_,
                 maximumDegree_,
                 maximumOrder_,
                 dateOverrides_ );
    }

    //! Restore the complete construction inputs with the same explicit consuming conventions.
    template< class Archive >
    void load( Archive& archive )
    {
        archive( cereal::base_class< GravityFieldVariationSettings >( this ),
                 data_,
                 sourceTimeScale_,
                 yearConvention_,
                 maximumDegree_,
                 maximumOrder_,
                 dateOverrides_ );
    }
};

//! Construct analytic temporal settings; callers must select both the source scale and definition of a year.
std::shared_ptr< GravityFieldVariationSettings > icgemTemporalGravityFieldVariationSettings(
        input_output::IcgemGravityFieldData data,
        basic_astrodynamics::TimeScales sourceTimeScale,
        gravitation::IcgemYearConvention yearConvention,
        int maximumDegree,
        int maximumOrder,
        std::map< std::string, double > dateOverrides = {} );

}  // namespace simulation_setup
}  // namespace tudat
#endif
