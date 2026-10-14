/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_CREATE_EARTH_GRAVITY_FIELD_VARIATIONS_H
#define TUDAT_CREATE_EARTH_GRAVITY_FIELD_VARIATIONS_H

#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"
#include "tudat/astro/gravitation/oceanPoleTideGravityFieldVariations.h"

namespace tudat
{
namespace simulation_setup
{
//! Complete conventional anelastic Earth response to Moon and Sun, including frequency terms.
//! Requires an ITRS field with an explicit tide-free or zero-tide convention and TDB-input Earth rotation.
class Iers2010SolidBodyGravityFieldVariationSettings : public GravityFieldVariationSettings
{
public:
    explicit Iers2010SolidBodyGravityFieldVariationSettings( const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr ):
        GravityFieldVariationSettings( gravitation::iers_2010, interpolation )
    {}
};

//! Solid pole response using the selected corrected pole and shared 2018 secular pole.
class SolidPoleTideGravityFieldVariationSettings : public GravityFieldVariationSettings
{
public:
    explicit SolidPoleTideGravityFieldVariationSettings( const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr ):
        GravityFieldVariationSettings( gravitation::pole_tide, interpolation )
    {}
};

//! Full IERS 6.23 ocean-pole response, with explicit input integrals, load response and constants.
//! The selected body-centered block starts at degree two/order zero.
class OceanPoleTideGravityFieldVariationSettings : public GravityFieldVariationSettings
{
public:
    OceanPoleTideGravityFieldVariationSettings( const input_output::OceanPoleGravityData& data,
                                                const std::vector< double >& loadLoveNumbers,
                                                const gravitation::OceanPoleGravityConstants& constants,
                                                int maximumDegree,
                                                int maximumOrder,
                                                const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr ):
        GravityFieldVariationSettings( gravitation::ocean_pole_tide, interpolation ), data_( data ), loadLoveNumbers_( loadLoveNumbers ),
        constants_( constants ), maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder )
    {}

    const input_output::OceanPoleGravityData& getData( ) const
    {
        return data_;
    }
    const std::vector< double >& getLoadLoveNumbers( ) const
    {
        return loadLoveNumbers_;
    }
    const gravitation::OceanPoleGravityConstants& getConstants( ) const
    {
        return constants_;
    }
    int getMaximumDegree( ) const
    {
        return maximumDegree_;
    }
    int getMaximumOrder( ) const
    {
        return maximumOrder_;
    }

private:
    input_output::OceanPoleGravityData data_;
    std::vector< double > loadLoveNumbers_;
    gravitation::OceanPoleGravityConstants constants_;
    int maximumDegree_, maximumOrder_;
};

//! Select the complete Moon/Sun response and optional deferred interpolation.
inline std::shared_ptr< GravityFieldVariationSettings > iers2010SolidBodyTideGravityFieldVariationSettings(
        const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr )
{
    return std::make_shared< Iers2010SolidBodyGravityFieldVariationSettings >( interpolation );
}
//! Select the conventional solid pole response using the rotation model's pole.
inline std::shared_ptr< GravityFieldVariationSettings > solidPoleTideGravityFieldVariationSettings(
        const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr )
{
    return std::make_shared< SolidPoleTideGravityFieldVariationSettings >( interpolation );
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

//! Construction and overlap checks local to the three Earth presets; interpolation remains deferred.
std::shared_ptr< gravitation::GravityFieldVariations > createEarthGravityFieldVariationsModel(
        const std::shared_ptr< GravityFieldVariationSettings >& settings,
        const std::string& body,
        const SystemOfBodies& bodies );
void validateEarthGravityFieldVariationCollection( const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& settings );

}  // namespace simulation_setup
}  // namespace tudat
#endif
