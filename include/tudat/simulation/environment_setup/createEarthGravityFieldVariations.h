/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_CREATE_EARTH_GRAVITY_FIELD_VARIATIONS_H
#define TUDAT_CREATE_EARTH_GRAVITY_FIELD_VARIATIONS_H

#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"
#include "tudat/simulation/environment_setup/createAod1bGravityFieldVariations.h"
#include "tudat/simulation/environment_setup/createIcgemTemporalGravityFieldVariations.h"
#include "tudat/astro/gravitation/oceanPoleTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/harmonicGravityFieldVariations.h"

#include <optional>
#include <cereal/types/array.hpp>

namespace tudat
{
namespace simulation_setup
{
//! Full IERS 6.23 ocean-pole response, with explicit input integrals, load response and constants.
//! The selected body-centered block starts at degree two/order zero.
class OceanPoleTideGravityFieldVariationSettings : public GravityFieldVariationSettings
{
public:
    //! Retain published integrals, degree-indexed load Love numbers and explicit forcing constants.
    OceanPoleTideGravityFieldVariationSettings( const input_output::OceanPoleGravityData& data,
                                                const std::vector< double >& loadLoveNumbers,
                                                const gravitation::OceanPoleGravityConstants& constants,
                                                int maximumDegree,
                                                int maximumOrder,
                                                const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr ):
        GravityFieldVariationSettings( gravitation::ocean_pole_tide, interpolation ), data_( data ), loadLoveNumbers_( loadLoveNumbers ),
        constants_( constants ), maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder )
    {}

    //! Return the published equilibrium integrals and their provenance.
    const input_output::OceanPoleGravityData& getData( ) const
    {
        return data_;
    }

    //! Return degree-indexed load Love numbers.
    const std::vector< double >& getLoadLoveNumbers( ) const
    {
        return loadLoveNumbers_;
    }

    //! Return the explicit forcing constants and source GM/radius.
    const gravitation::OceanPoleGravityConstants& getConstants( ) const
    {
        return constants_;
    }

    //! Return the largest selected spherical-harmonic degree.
    int getMaximumDegree( ) const
    {
        return maximumDegree_;
    }

    //! Return the largest selected spherical-harmonic order.
    int getMaximumOrder( ) const
    {
        return maximumOrder_;
    }

    //! Compare the exact source integrals, load response, physical constants and consuming settings.
    bool equals( const GravityFieldVariationSettings& other ) const override;

private:
    //! Construct an archive destination; normal construction requires explicit source inputs.
    OceanPoleTideGravityFieldVariationSettings( ): GravityFieldVariationSettings( gravitation::ocean_pole_tide ) {}

    friend class cereal::access;

    //! Save complete input data and output selection, retaining all source provenance.
    template< class Archive >
    void save( Archive& archive ) const
    {
        archive( cereal::base_class< GravityFieldVariationSettings >( this ),
                 data_,
                 loadLoveNumbers_,
                 constants_,
                 maximumDegree_,
                 maximumOrder_ );
    }

    //! Restore complete input data and output selection, retaining all source provenance.
    template< class Archive >
    void load( Archive& archive )
    {
        archive( cereal::base_class< GravityFieldVariationSettings >( this ),
                 data_,
                 loadLoveNumbers_,
                 constants_,
                 maximumDegree_,
                 maximumOrder_ );
    }

    //! Published equilibrium integrals, degree-indexed load Love numbers and explicit forcing constants.
    input_output::OceanPoleGravityData data_;

    std::vector< double > loadLoveNumbers_;

    gravitation::OceanPoleGravityConstants constants_;

    //! Inclusive upper bounds of the body-centered output block.
    int maximumDegree_ = 0, maximumOrder_ = 0;
};

//! Final ocean or direct atmospheric potential constituents, evaluated with the selected Earth rotation clock.
//! The Earth coefficient block starts at degree two/order zero; no extra loading, phase or nodal factors are applied.
class HarmonicGravityFieldVariationSettings : public GravityFieldVariationSettings
{
public:
    //! Retain canonical constituents, explicit product conventions and optional deferred interpolation.
    HarmonicGravityFieldVariationSettings( std::vector< gravitation::HarmonicGravityConstituent > constituents,
                                           const gravitation::HarmonicGravityMetadata& metadata,
                                           int maximumDegree,
                                           int maximumOrder,
                                           const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr );

    //! Return the final ordered inventory, including any explicitly applied admittance targets.
    const std::vector< gravitation::HarmonicGravityConstituent >& getConstituents( ) const
    {
        return constituents_;
    }

    //! Return the source product's physical component and reference conventions.
    const gravitation::HarmonicGravityMetadata& getMetadata( ) const
    {
        return metadata_;
    }

    //! Return a stable component, product revision and final inventory label for model lookup.
    std::string getIdentifier( ) const;

    //! Return the largest selected spherical-harmonic degree.
    int getMaximumDegree( ) const
    {
        return maximumDegree_;
    }

    //! Return the largest selected spherical-harmonic order.
    int getMaximumOrder( ) const
    {
        return maximumOrder_;
    }

    //! Compare the exact final inventory, conventions, selected block and interpolation settings.
    bool equals( const GravityFieldVariationSettings& other ) const override;

private:
    //! Construct an empty archive destination; a usable model still validates the restored data.
    HarmonicGravityFieldVariationSettings( ) = default;

    friend class cereal::access;

    //! Save the final canonical data and consuming settings.
    template< class Archive >
    void save( Archive& archive ) const
    {
        archive( cereal::base_class< GravityFieldVariationSettings >( this ), constituents_, metadata_, maximumDegree_, maximumOrder_ );
    }

    //! Restore the final canonical data and consuming settings.
    template< class Archive >
    void load( Archive& archive )
    {
        archive( cereal::base_class< GravityFieldVariationSettings >( this ), constituents_, metadata_, maximumDegree_, maximumOrder_ );
    }

    //! Final dimensionless coefficient matrices and their source GM/radius, normalization and loading declarations.
    std::vector< gravitation::HarmonicGravityConstituent > constituents_;

    gravitation::HarmonicGravityMetadata metadata_;

    //! Inclusive upper bounds; Earth setup uses the degree-two/order-zero lower bounds.
    int maximumDegree_ = 0, maximumOrder_ = 0;
};

//! Convert explicit potential data and optional ordered admittance once, storing only the final constituents.
//! Missing selected source coefficients reject; a caller needing generic zero-fill can construct canonical settings explicitly.
std::shared_ptr< GravityFieldVariationSettings > harmonicPotentialGravityFieldVariationSettings(
        const input_output::HarmonicPotentialData& data,
        const gravitation::HarmonicGravityMetadata& metadata,
        int maximumDegree,
        int maximumOrder,
        const std::vector< std::string >& selectedDoodsonIds = {},
        const std::optional< gravitation::HarmonicGravityAdmittance >& admittance = std::nullopt,
        const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr );

//! Select the complete anelastic Moon/Sun response and optional deferred interpolation.
//! Requires an ITRS field with an explicit tide-free or zero-tide convention and TDB-input Earth rotation.
inline std::shared_ptr< GravityFieldVariationSettings > iers2010SolidBodyTideGravityFieldVariationSettings(
        const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr )
{
    return std::make_shared< GravityFieldVariationSettings >( gravitation::iers_2010, interpolation );
}

//! Select the conventional solid pole response using the rotation model's pole.
inline std::shared_ptr< GravityFieldVariationSettings > solidPoleTideGravityFieldVariationSettings(
        const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr )
{
    return std::make_shared< GravityFieldVariationSettings >( gravitation::pole_tide, interpolation );
}

//! Select the supplied equilibrium ocean-pole response and body-centered truncation.
inline std::shared_ptr< GravityFieldVariationSettings > oceanPoleTideGravityFieldVariationSettings(
        const input_output::OceanPoleGravityData& data,
        const std::vector< double >& loadLoveNumbers,
        const gravitation::OceanPoleGravityConstants& constants,
        int maximumDegree,
        int maximumOrder,
        const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr )
{
    return std::make_shared< OceanPoleTideGravityFieldVariationSettings >(
            data, loadLoveNumbers, constants, maximumDegree, maximumOrder, interpolation );
}

//! Construct conventional and harmonic Earth variations; interpolation remains deferred.
std::shared_ptr< gravitation::GravityFieldVariations > createEarthGravityFieldVariationsModel(
        const std::shared_ptr< GravityFieldVariationSettings >& settings,
        const std::string& body,
        const SystemOfBodies& bodies );

//! Reject repeated conventional responses and overlapping final harmonic inventories within one physical component.
void validateEarthGravityFieldVariationCollection( const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& settings );

}  // namespace simulation_setup
}  // namespace tudat
#endif
