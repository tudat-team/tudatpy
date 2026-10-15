/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/simulation/environment_setup/createAod1bGravityFieldVariations.h"

#include <set>
#include <stdexcept>

namespace tudat
{
namespace simulation_setup
{

//! Move the preloaded products into ordinary variation settings with native GPS interpolation.
Aod1bGravityFieldVariationSettings::Aod1bGravityFieldVariationSettings( std::vector< input_output::Aod1bGravityFieldData > data,
                                                                        const input_output::Aod1bComponent component,
                                                                        const int maximumDegree,
                                                                        const int maximumOrder ):
    GravityFieldVariationSettings( gravitation::aod1b_rl07 ), data_( std::move( data ) ), component_( component ),
    maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder )
{}

//! Identify the selected component independently of daily filenames and truncation.
std::string Aod1bGravityFieldVariationSettings::getIdentifier( ) const
{
    switch( component_ )
    {
        case input_output::Aod1bComponent::atmosphere:
            return "AOD1B_RL07_ATM";
        case input_output::Aod1bComponent::ocean:
            return "AOD1B_RL07_OCN";
        case input_output::Aod1bComponent::global:
            return "AOD1B_RL07_GLO";
        case input_output::Aod1bComponent::oceanBottomPressure:
            return "AOD1B_RL07_OBA";
        default:
            throw std::runtime_error( "Unknown AOD1B mass component" );
    }
}

//! Compare value-owned source products as well as the selected physical interpretation.
bool Aod1bGravityFieldVariationSettings::equals( const GravityFieldVariationSettings& other ) const
{
    // Only this settings type has the component interpretation and raw records being compared.
    const auto* settings = dynamic_cast< const Aod1bGravityFieldVariationSettings* >( &other );

    return settings != nullptr && GravityFieldVariationSettings::equals( other ) && component_ == settings->component_ &&
            maximumDegree_ == settings->maximumDegree_ && maximumOrder_ == settings->maximumOrder_ && data_ == settings->data_;
}

//! Retain preloaded products so construction and archive loading never depend on external paths.
std::shared_ptr< GravityFieldVariationSettings > aod1bGravityFieldVariationSettings(
        std::vector< input_output::Aod1bGravityFieldData > data,
        const input_output::Aod1bComponent component,
        const int maximumDegree,
        const int maximumOrder )
{
    return std::make_shared< Aod1bGravityFieldVariationSettings >( std::move( data ), component, maximumDegree, maximumOrder );
}

//! Validate only the AOD1B subset, leaving independent tidal and generic variations unchanged.
void validateAod1bGravityFieldVariationCollection( const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& settings )
{
    // ATM and OCN form the one supported additive decomposition; GLO and OBA stand alone.
    std::set< input_output::Aod1bComponent > components;

    for( const auto& entry : settings )
    {
        const auto aod = std::dynamic_pointer_cast< Aod1bGravityFieldVariationSettings >( entry );

        if( aod && !components.insert( aod->getComponent( ) ).second )
        {
            throw std::runtime_error( "Duplicate AOD1B mass component" );
        }
    }
    if( components.size( ) > 1 &&
        ( components.count( input_output::Aod1bComponent::global ) ||
          components.count( input_output::Aod1bComponent::oceanBottomPressure ) ) )
    {
        throw std::runtime_error( "AOD1B GLO overlaps ATM/OCN; OBA has no supported additive decomposition with another AOD1B product" );
    }
}

}  // namespace simulation_setup
}  // namespace tudat
