/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/simulation/environment_setup/createEarthGravityFieldVariations.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/gravitation/earthSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/interface/sofa/fundamentalArguments.h"
#include <set>

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
    const Time tai = converter->getTaiTime( basic_astrodynamics::tdb_scale, Time( time ) );
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
    const Time tt = basic_astrodynamics::convertTAItoTT( converter->getTaiTime( basic_astrodynamics::tdb_scale, Time( time ) ) );
    const double radiansPerMas = mathematical_constants::PI / ( 180.0 * 3600000.0 );
    Eigen::Vector2d wobble = polar->getPositionOfCipInItrsAtTt( tt ) -
            earth_orientation::getSecularPolePositionInMas( static_cast< double >( tt ) ) * radiansPerMas;
    wobble.y( ) = -wobble.y( );
    return wobble * radiansToOutputUnits;
}
}  // namespace
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
        case gravitation::iers_2010: {
            if( !std::dynamic_pointer_cast< Iers2010SolidBodyGravityFieldVariationSettings >( settings ) )
            {
                throw std::runtime_error( "Expected complete IERS Earth solid-tide settings" );
            }
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
            if( !std::dynamic_pointer_cast< SolidPoleTideGravityFieldVariationSettings >( settings ) )
            {
                throw std::runtime_error( "Expected solid pole gravity settings" );
            }
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

//! Reject repeated responses and Moon/Sun forcing already included by the complete IERS model.
void validateEarthGravityFieldVariationCollection( const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& settings )
{
    std::set< gravitation::BodyDeformationTypes > unique;
    for( const auto& entry : settings )
    {
        if( !entry )
        {
            throw std::runtime_error( "Null gravity variation settings" );
        }
        const auto type = entry->getBodyDeformationType( );
        if( type == gravitation::iers_2010 || type == gravitation::pole_tide || type == gravitation::ocean_pole_tide )
        {
            if( !unique.insert( type ).second )
            {
                throw std::runtime_error( "Duplicate complete solid or pole gravity tide" );
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
