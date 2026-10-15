/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_CREATE_AOD1B_GRAVITY_FIELD_VARIATIONS_H
#define TUDAT_CREATE_AOD1B_GRAVITY_FIELD_VARIATIONS_H

#include <cereal/types/array.hpp>
#include "tudat/io/readAod1bGravityField.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"

namespace tudat
{
namespace input_output
{

//! Archive source values and provenance without making the neutral reader depend on cereal.
template< class Archive >
void serialize( Archive& archive, Aod1bGravityCoefficient& coefficient )
{
    archive( CEREAL_NVP( coefficient.values ), CEREAL_NVP( coefficient.sourceLine ) );
}

//! Archive a component's printed GPS date, continuous epoch and supplied coefficient records.
template< class Archive >
void serialize( Archive& archive, Aod1bGravitySample& sample )
{
    archive( CEREAL_NVP( sample.component ),
             CEREAL_NVP( sample.gpsDate ),
             CEREAL_NVP( sample.gpsSecondsSinceJ2000 ),
             CEREAL_NVP( sample.coefficients ) );
}

//! Archive the complete raw product so settings remain independent of external files.
template< class Archive >
void serialize( Archive& archive, Aod1bGravityFieldData& data )
{
    archive( CEREAL_NVP( data.source ),
             CEREAL_NVP( data.headerLines ),
             CEREAL_NVP( data.release ),
             CEREAL_NVP( data.softwareVersion ),
             CEREAL_NVP( data.gravitationalParameter ),
             CEREAL_NVP( data.referenceRadius ),
             CEREAL_NVP( data.fullyNormalized ),
             CEREAL_NVP( data.maximumDegree ),
             CEREAL_NVP( data.samples ) );
}

}  // namespace input_output
namespace simulation_setup
{

//! Select one RL07 non-tidal component and an expansion beginning at degree two/order zero.
//! Products are moved into the settings; the runtime model retains only scaled interpolation tables.
class Aod1bGravityFieldVariationSettings : public GravityFieldVariationSettings
{
public:
    //! Store preloaded daily products without introducing a second interpolation stage.
    Aod1bGravityFieldVariationSettings( std::vector< input_output::Aod1bGravityFieldData > data,
                                        input_output::Aod1bComponent component,
                                        int maximumDegree,
                                        int maximumOrder );

    //! Inspect the retained source records without copying daily products.
    const std::vector< input_output::Aod1bGravityFieldData >& getData( ) const
    {
        return data_;
    }

    //! Retrieve the selected physical mass component.
    input_output::Aod1bComponent getComponent( ) const
    {
        return component_;
    }

    //! Return the stable RL07 component label used by GravityFieldVariationsSet lookup.
    std::string getIdentifier( ) const;

    //! Retrieve the highest included spherical-harmonic degree.
    int getMaximumDegree( ) const
    {
        return maximumDegree_;
    }

    //! Retrieve the highest included spherical-harmonic order.
    int getMaximumOrder( ) const
    {
        return maximumOrder_;
    }

    //! Compare the selected component, truncation and complete source records.
    bool equals( const GravityFieldVariationSettings& other ) const override;

private:
    //! Raw daily products, selected component and maximum output degree/order.
    std::vector< input_output::Aod1bGravityFieldData > data_;

    input_output::Aod1bComponent component_ = input_output::Aod1bComponent::global;

    int maximumDegree_ = 2, maximumOrder_ = 0;

    //! Cereal reconstructs the same preloaded settings, rather than reopening provider files.
    friend class cereal::access;

    //! Provide an empty object solely for archive reconstruction.
    Aod1bGravityFieldVariationSettings( ): GravityFieldVariationSettings( gravitation::aod1b_rl07 ) {}

    //! Archive the base dispatch type and all AOD1B construction inputs.
    template< class Archive >
    void save( Archive& archive ) const
    {
        archive( cereal::base_class< GravityFieldVariationSettings >( this ),
                 CEREAL_NVP( data_ ),
                 CEREAL_NVP( component_ ),
                 CEREAL_NVP( maximumDegree_ ),
                 CEREAL_NVP( maximumOrder_ ) );
    }

    //! Restore the same preloaded source records and selected coefficient block.
    template< class Archive >
    void load( Archive& archive )
    {
        archive( cereal::base_class< GravityFieldVariationSettings >( this ),
                 CEREAL_NVP( data_ ),
                 CEREAL_NVP( component_ ),
                 CEREAL_NVP( maximumDegree_ ),
                 CEREAL_NVP( maximumOrder_ ) );
    }
};

//! Construct settings from preloaded daily RL07 files and an explicit physical component.
std::shared_ptr< GravityFieldVariationSettings > aod1bGravityFieldVariationSettings(
        std::vector< input_output::Aod1bGravityFieldData > data,
        input_output::Aod1bComponent component,
        int maximumDegree,
        int maximumOrder );

//! Reject duplicate components, GLO plus ATM/OCN, and unsupported OBA decompositions.
void validateAod1bGravityFieldVariationCollection( const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& settings );

}  // namespace simulation_setup
}  // namespace tudat
#endif
