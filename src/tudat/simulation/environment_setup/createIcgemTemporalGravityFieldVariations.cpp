/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/simulation/environment_setup/createIcgemTemporalGravityFieldVariations.h"

namespace tudat
{
namespace simulation_setup
{

//! Store the preloaded source once with the caller's explicit date, time-scale and year interpretation.
IcgemTemporalGravityFieldVariationSettings::IcgemTemporalGravityFieldVariationSettings(
        input_output::IcgemGravityFieldData data,
        const basic_astrodynamics::TimeScales sourceTimeScale,
        const gravitation::IcgemYearConvention yearConvention,
        const int maximumDegree,
        const int maximumOrder,
        std::map< std::string, double > dateOverrides ):
    GravityFieldVariationSettings( gravitation::icgem_temporal ), data_( std::move( data ) ), sourceTimeScale_( sourceTimeScale ),
    yearConvention_( yearConvention ), maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder ),
    dateOverrides_( std::move( dateOverrides ) )
{}

//! Identify the source product through the existing per-variation lookup.
std::string IcgemTemporalGravityFieldVariationSettings::getIdentifier( ) const
{
    return "ICGEM_TEMPORAL_" + data_.modelName;
}

//! Compare all source records and their consuming interpretation, including exact-token overrides.
bool IcgemTemporalGravityFieldVariationSettings::equals( const GravityFieldVariationSettings& other ) const
{
    const auto* settings = dynamic_cast< const IcgemTemporalGravityFieldVariationSettings* >( &other );

    return settings != nullptr && GravityFieldVariationSettings::equals( other ) && data_ == settings->data_ &&
            sourceTimeScale_ == settings->sourceTimeScale_ && yearConvention_ == settings->yearConvention_ &&
            maximumDegree_ == settings->maximumDegree_ && maximumOrder_ == settings->maximumOrder_ &&
            dateOverrides_ == settings->dateOverrides_;
}

//! Construct ordinary variation settings without introducing an interpolation across temporal discontinuities.
std::shared_ptr< GravityFieldVariationSettings > icgemTemporalGravityFieldVariationSettings(
        input_output::IcgemGravityFieldData data,
        const basic_astrodynamics::TimeScales sourceTimeScale,
        const gravitation::IcgemYearConvention yearConvention,
        const int maximumDegree,
        const int maximumOrder,
        std::map< std::string, double > dateOverrides )
{
    return std::make_shared< IcgemTemporalGravityFieldVariationSettings >(
            std::move( data ), sourceTimeScale, yearConvention, maximumDegree, maximumOrder, std::move( dateOverrides ) );
}

}  // namespace simulation_setup
}  // namespace tudat
