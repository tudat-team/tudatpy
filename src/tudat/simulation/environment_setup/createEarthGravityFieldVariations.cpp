/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/simulation/environment_setup/createEarthGravityFieldVariations.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/gravitation/earthSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/aod1bGravityFieldVariations.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/interface/sofa/fundamentalArguments.h"
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace tudat
{
namespace simulation_setup
{
namespace
{
//! Build Doodson phases from the selected rotation clock, with TT and UT1 kept distinct.
Eigen::Vector6d earthGravityDoodsonArguments( double time,
                                              const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter )
{
    const Time tai = converter->getTaiTime( basic_astrodynamics::tdb_scale, time );
    const Time tt = basic_astrodynamics::convertTAItoTT( tai );
    return sofa_interface::calculateDoodsonFundamentalArguments(
            static_cast< double >( tt ), static_cast< double >( tt ), static_cast< double >( converter->getUt1FromTai( tai ) ) );
}
//! Subtract the 2018 secular pole and convert m1/m2 from radians to the requested units.
Eigen::Vector2d earthGravityWobble( double time,
                                    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter,
                                    const std::shared_ptr< earth_orientation::PolarMotionCalculator >& polar,
                                    double radiansToOutputUnits )
{
    const Time tt = basic_astrodynamics::convertTAItoTT( converter->getTaiTime( basic_astrodynamics::tdb_scale, time ) );
    const double radiansPerMas = mathematical_constants::PI / ( 180.0 * 3600000.0 );
    Eigen::Vector2d wobble = polar->getPositionOfCipInItrsAtTt( tt ) -
            earth_orientation::getSecularPolePositionInMas( static_cast< double >( tt ) ) * radiansPerMas;
    wobble.y( ) = -wobble.y( );
    return wobble * radiansToOutputUnits;
}
}  // namespace

//! Compare complete ocean-pole source records and physical response settings exactly.
bool OceanPoleTideGravityFieldVariationSettings::equals( const GravityFieldVariationSettings& other ) const
{
    const auto* rhs = dynamic_cast< const OceanPoleTideGravityFieldVariationSettings* >( &other );

    if( !rhs || !GravityFieldVariationSettings::equals( other ) || maximumDegree_ != rhs->maximumDegree_ ||
        maximumOrder_ != rhs->maximumOrder_ || loadLoveNumbers_ != rhs->loadLoveNumbers_ || data_.source != rhs->data_.source ||
        data_.headerLines != rhs->data_.headerLines || data_.maximumDegree != rhs->data_.maximumDegree ||
        data_.coefficients.size( ) != rhs->data_.coefficients.size( ) ||
        std::tie( constants_.angularVelocity,
                  constants_.gravitationalConstant,
                  constants_.oceanDensity,
                  constants_.equatorialGravity,
                  constants_.gravitationalParameter,
                  constants_.referenceRadius,
                  constants_.gamma2 ) !=
                std::tie( rhs->constants_.angularVelocity,
                          rhs->constants_.gravitationalConstant,
                          rhs->constants_.oceanDensity,
                          rhs->constants_.equatorialGravity,
                          rhs->constants_.gravitationalParameter,
                          rhs->constants_.referenceRadius,
                          rhs->constants_.gamma2 ) )
    {
        return false;
    }
    auto right = rhs->data_.coefficients.begin( );

    for( const auto& left : data_.coefficients )
    {
        if( left.first != right->first || left.second.a != right->second.a || left.second.b != right->second.b ||
            left.second.sourceLine != right->second.sourceLine )
        {
            return false;
        }
        ++right;
    }
    return true;
}

//! Preserve final coefficients and select the gravity variation type from their declared physical component.
HarmonicGravityFieldVariationSettings::HarmonicGravityFieldVariationSettings(
        std::vector< gravitation::HarmonicGravityConstituent > constituents,
        const gravitation::HarmonicGravityMetadata& metadata,
        int maximumDegree,
        int maximumOrder,
        const std::shared_ptr< ModelInterpolationSettings >& interpolation ):
    GravityFieldVariationSettings( metadata.component == gravitation::HarmonicGravityComponent::oceanMass ? gravitation::ocean_tide
                                                                                                          : gravitation::atmospheric_tide,
                                   interpolation ),
    constituents_( std::move( constituents ) ), metadata_( metadata ), maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder )
{
    if( metadata_.component != gravitation::HarmonicGravityComponent::oceanMass &&
        metadata_.component != gravitation::HarmonicGravityComponent::directAtmosphericMass )
    {
        throw std::invalid_argument( "Harmonic gravity settings require an explicit ocean or direct atmospheric mass component" );
    }
}

//! Label each harmonic model with its declared component, revision and complete ordered target inventory.
std::string HarmonicGravityFieldVariationSettings::getIdentifier( ) const
{
    std::string identifier = metadata_.component == gravitation::HarmonicGravityComponent::oceanMass ? "ocean_mass:" : "direct_air_mass:";

    identifier += metadata_.revision;
    for( const auto& wave : constituents_ )
    {
        identifier += ":" + wave.doodson;
    }
    return identifier;
}

//! Compare every stored product declaration and canonical matrix, including the final constituent order.
bool HarmonicGravityFieldVariationSettings::equals( const GravityFieldVariationSettings& other ) const
{
    const auto* rhs = dynamic_cast< const HarmonicGravityFieldVariationSettings* >( &other );

    if( !rhs || !GravityFieldVariationSettings::equals( other ) || maximumDegree_ != rhs->maximumDegree_ ||
        maximumOrder_ != rhs->maximumOrder_ || metadata_.component != rhs->metadata_.component ||
        metadata_.gravitationalParameter != rhs->metadata_.gravitationalParameter ||
        metadata_.referenceRadius != rhs->metadata_.referenceRadius || metadata_.fullyNormalized != rhs->metadata_.fullyNormalized ||
        metadata_.revision != rhs->metadata_.revision || metadata_.loading != rhs->metadata_.loading ||
        constituents_.size( ) != rhs->constituents_.size( ) )
    {
        return false;
    }
    for( std::size_t i = 0; i < constituents_.size( ); ++i )
    {
        const auto& left = constituents_[ i ];

        const auto& right = rhs->constituents_[ i ];

        if( left.doodson != right.doodson || left.multipliers != right.multipliers || left.cCos.rows( ) != right.cCos.rows( ) ||
            left.cCos.cols( ) != right.cCos.cols( ) || left.cSin.rows( ) != right.cSin.rows( ) || left.cSin.cols( ) != right.cSin.cols( ) ||
            left.sCos.rows( ) != right.sCos.rows( ) || left.sCos.cols( ) != right.sCos.cols( ) || left.sSin.rows( ) != right.sSin.rows( ) ||
            left.sSin.cols( ) != right.sSin.cols( ) || left.cCos != right.cCos || left.cSin != right.cSin || left.sCos != right.sCos ||
            left.sSin != right.sSin )
        {
            return false;
        }
    }
    return true;
}

//! Convert a selected raw potential inventory and replace it by the requested ordered admittance targets once.
std::shared_ptr< GravityFieldVariationSettings > harmonicPotentialGravityFieldVariationSettings(
        const input_output::HarmonicPotentialData& data,
        const gravitation::HarmonicGravityMetadata& metadata,
        int maximumDegree,
        int maximumOrder,
        const std::vector< std::string >& selectedDoodsonIds,
        const std::optional< gravitation::HarmonicGravityAdmittance >& admittance,
        const std::shared_ptr< ModelInterpolationSettings >& interpolation )
{
    auto constituents = gravitation::convertHarmonicPotentialData( data, maximumDegree, maximumOrder, false, selectedDoodsonIds );

    if( admittance )
    {
        constituents = gravitation::applyHarmonicGravityAdmittance( constituents, *admittance );
    }
    return std::make_shared< HarmonicGravityFieldVariationSettings >(
            std::move( constituents ), metadata, maximumDegree, maximumOrder, interpolation );
}

//! Construct the Earth response using the selected rotation and consuming field conventions.
std::shared_ptr< gravitation::GravityFieldVariations > createEarthGravityFieldVariationsModel(
        const std::shared_ptr< GravityFieldVariationSettings >& settings,
        const std::string& body,
        const SystemOfBodies& bodies )
{
    using namespace basic_astrodynamics;
    if( body != "Earth" || !bodies.count( body ) )
    {
        throw std::runtime_error( "Conventional Earth gravity variations require body Earth" );
    }
    const auto earth = bodies.at( body );
    const auto field =
            std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( earth->getGravityFieldModel( ) );
    const auto rotation = std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( earth->getRotationalEphemeris( ) );
    if( !field || !rotation )
    {
        throw std::runtime_error( "Earth gravity tides require a time-dependent spherical field and GCRS/ITRS rotation" );
    }
    if( rotation->getInputTimeScale( ) != tdb_scale )
    {
        throw std::runtime_error( "Earth gravity tides require a TDB-input rotation; select TDB in the Earth rotation settings" );
    }
    if( field->getFixedReferenceFrame( ) != "ITRS" || bodies.getFrameOrientation( ) != rotation->getBaseFrameOrientation( ) )
    {
        throw std::runtime_error( "Earth gravity tides require ITRS coefficients and identical global/rotation base axes labels" );
    }
    const auto angles = rotation->getAnglesCalculator( );
    if( !angles )
    {
        throw std::runtime_error( "Earth gravity tides require the selected Earth orientation calculator" );
    }
    const auto converter = angles->getTerrestrialTimeScaleConverter( );
    const auto polar = angles->getPolarMotionCalculator( );
    if( !converter || !polar )
    {
        throw std::runtime_error( "Earth gravity tides require the selected rotation clock and pole" );
    }
    constexpr double sourceMu = 3.986004415e14;
    constexpr double sourceRadius = 6378136.3;
    switch( settings->getBodyDeformationType( ) )
    {
        case gravitation::icgem_temporal: {
            // Preserve analytic interval ownership; interpolating again could bridge a jump or missing base piece.
            const auto temporal = std::dynamic_pointer_cast< IcgemTemporalGravityFieldVariationSettings >( settings );

            if( !temporal || settings->getInterpolatorSettings( ) != nullptr )
            {
                throw std::runtime_error( "Expected temporal ICGEM settings without a second interpolation stage" );
            }
            if( field->getTideSystem( ) == TideSystem::unknown || temporal->getData( ).tideSystem != field->getTideSystem( ) )
            {
                throw std::runtime_error( "Temporal ICGEM source and nominal field require matching known tide systems" );
            }
            return gravitation::createIcgemTemporalGravityFieldVariations( temporal->getData( ),
                                                                           temporal->getSourceTimeScale( ),
                                                                           temporal->getYearConvention( ),
                                                                           temporal->getMaximumDegree( ),
                                                                           temporal->getMaximumOrder( ),
                                                                           temporal->getDateOverrides( ),
                                                                           field->getNominalCosineCoefficients( ),
                                                                           field->getNominalSineCoefficients( ),
                                                                           field->getGravitationalParameter( ),
                                                                           field->getReferenceRadius( ),
                                                                           converter );
        }
        case gravitation::aod1b_rl07: {
            // RL07 already interpolates in GPS seconds; another TDB interpolation stage would change that model.
            const auto aod = std::dynamic_pointer_cast< Aod1bGravityFieldVariationSettings >( settings );

            if( !aod || settings->getInterpolatorSettings( ) != nullptr )
            {
                throw std::runtime_error( "Expected AOD1B RL07 settings without a second interpolation stage" );
            }
            return gravitation::createAod1bGravityFieldVariations( aod->getData( ),
                                                                   aod->getComponent( ),
                                                                   aod->getMaximumDegree( ),
                                                                   aod->getMaximumOrder( ),
                                                                   field->getGravitationalParameter( ),
                                                                   field->getReferenceRadius( ),
                                                                   converter );
        }
        case gravitation::ocean_tide:
        case gravitation::atmospheric_tide: {
            const auto harmonic = std::dynamic_pointer_cast< HarmonicGravityFieldVariationSettings >( settings );

            if( !harmonic )
            {
                throw std::runtime_error( "Expected canonical harmonic gravity settings" );
            }
            return std::make_shared< gravitation::HarmonicGravityFieldVariations >(
                    harmonic->getConstituents( ),
                    harmonic->getMetadata( ),
                    std::bind( earthGravityDoodsonArguments, std::placeholders::_1, converter ),
                    field->getGravitationalParameter( ),
                    field->getReferenceRadius( ),
                    2,
                    0,
                    harmonic->getMaximumDegree( ),
                    harmonic->getMaximumOrder( ) );
        }
        case gravitation::iers_2010: {
            if( field->getTideSystem( ) != TideSystem::tideFree && field->getTideSystem( ) != TideSystem::zeroTide )
            {
                throw std::runtime_error( "Complete IERS Earth solid tide requires an explicit tide-free or zero-tide nominal field" );
            }
            std::function< Eigen::Vector6d( double ) > earthState;
            std::function< Eigen::Quaterniond( double ) > orientation;
            std::vector< std::function< Eigen::Vector6d( double ) > > raiserStates;
            std::vector< std::function< double( ) > > raiserMus;
            const bool interpolated = settings->getInterpolatorSettings( ) != nullptr;
            if( interpolated )
            {
                if( !earth->getEphemeris( ) )
                {
                    throw std::runtime_error( "Interpolated Earth solid tide requires an Earth ephemeris" );
                }
                earthState = std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >, earth, std::placeholders::_1 );
                orientation = std::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame, rotation, std::placeholders::_1 );
            }
            else
            {
                earthState = std::bind( &Body::getState, earth );
                orientation = std::bind( &Body::getCurrentRotationToLocalFrame, earth );
            }
            const std::vector< std::string > raisers{ "Moon", "Sun" };
            for( const auto& name : raisers )
            {
                if( !bodies.count( name ) || !bodies.at( name )->getGravityFieldModel( ) )
                {
                    throw std::runtime_error( "Complete IERS Earth solid tide requires Moon and Sun gravity models" );
                }
                const auto raiser = bodies.at( name );
                if( interpolated )
                {
                    if( !raiser->getEphemeris( ) )
                    {
                        throw std::runtime_error( "Interpolated Earth solid tide requires raiser ephemerides" );
                    }
                    raiserStates.push_back(
                            std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >, raiser, std::placeholders::_1 ) );
                }
                else
                {
                    raiserStates.push_back( std::bind( &Body::getState, raiser ) );
                }
                raiserMus.push_back( std::bind( &Body::getGravitationalParameter, raiser ) );
            }
            return std::make_shared< gravitation::EarthSolidBodyTideGravityFieldVariations >(
                    earthState,
                    orientation,
                    raiserStates,
                    sourceRadius,
                    []( ) { return sourceMu; },
                    raiserMus,
                    raisers,
                    std::bind( earthGravityDoodsonArguments, std::placeholders::_1, converter ),
                    field->getReferenceRadius( ),
                    field->getGravitationalParameter( ),
                    field->getTideSystem( ) );
        }
        case gravitation::pole_tide: {
            return std::make_shared< gravitation::SolidPoleTideGravityFieldVariations >(
                    std::bind( earthGravityWobble, std::placeholders::_1, converter, polar, 180.0 * 3600.0 / mathematical_constants::PI ),
                    sourceMu,
                    sourceRadius,
                    field->getGravitationalParameter( ),
                    field->getReferenceRadius( ) );
        }
        case gravitation::ocean_pole_tide: {
            const auto ocean = std::dynamic_pointer_cast< OceanPoleTideGravityFieldVariationSettings >( settings );
            if( !ocean )
            {
                throw std::runtime_error( "Expected ocean pole gravity settings" );
            }
            return std::make_shared< gravitation::OceanPoleTideGravityFieldVariations >(
                    ocean->getData( ),
                    std::bind( earthGravityWobble, std::placeholders::_1, converter, polar, 1.0 ),
                    ocean->getLoadLoveNumbers( ),
                    ocean->getConstants( ),
                    field->getGravitationalParameter( ),
                    field->getReferenceRadius( ),
                    2,
                    0,
                    ocean->getMaximumDegree( ),
                    ocean->getMaximumOrder( ) );
        }
        default:
            throw std::runtime_error( "Unsupported Earth gravity preset" );
    }
}

//! Reject repeated conventional responses and duplicate final harmonic waves within one physical component.
void validateEarthGravityFieldVariationCollection( const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& settings )
{
    // Enforce the published non-tidal component decomposition before the conventional-tide early return.
    validateAod1bGravityFieldVariationCollection( settings );

    std::set< gravitation::BodyDeformationTypes > unique;

    std::set< std::pair< gravitation::HarmonicGravityComponent, std::string > > harmonicInventory;

    for( const auto& entry : settings )
    {
        if( !entry )
        {
            throw std::runtime_error( "Null gravity variation settings" );
        }
        const auto type = entry->getBodyDeformationType( );

        if( const auto harmonic = std::dynamic_pointer_cast< HarmonicGravityFieldVariationSettings >( entry ) )
        {
            // Radiational ocean tides remain ocean mass; direct air mass is a separate physical signal.
            for( const auto& wave : harmonic->getConstituents( ) )
            {
                if( !harmonicInventory.emplace( harmonic->getMetadata( ).component, wave.doodson ).second )
                {
                    throw std::runtime_error( "Overlapping final harmonic gravity constituents within the same physical component" );
                }
            }
        }
        if( type == gravitation::iers_2010 || type == gravitation::pole_tide || type == gravitation::ocean_pole_tide ||
            type == gravitation::icgem_temporal )
        {
            if( !unique.insert( type ).second )
            {
                throw std::runtime_error( "Duplicate complete solid tide, pole tide or temporal ICGEM model" );
            }
        }
    }
    if( !unique.count( gravitation::iers_2010 ) )
    {
        return;
    }
    for( const auto& entry : settings )
    {
        std::vector< std::string > raisers;
        if( const auto basic = std::dynamic_pointer_cast< BasicSolidBodyGravityFieldVariationSettings >( entry ) )
        {
            raisers = basic->getDeformingBodies( );
        }
        else if( const auto coupled = std::dynamic_pointer_cast< ModeCoupledSolidBodyGravityFieldVariationSettings >( entry ) )
        {
            raisers = coupled->getDeformingBodies( );
        }
        for( const auto& raiser : raisers )
        {
            if( raiser == "Moon" || raiser == "Sun" )
            {
                throw std::runtime_error( "Complete IERS Earth tide overlaps basic or mode-coupled Moon/Sun forcing" );
            }
        }
    }
}
}  // namespace simulation_setup
}  // namespace tudat
