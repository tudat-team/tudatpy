/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createGroundStations.h"
#include <set>
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/ground_stations/meteorologicalConditions.h"
#include "tudat/io/readTabulatedWeatherData.h"

namespace tudat
{

namespace simulation_setup
{

std::shared_ptr< ground_stations::StationMeteoData > createGroundStationMeteoData(
        const std::shared_ptr< WeatherDataSettings > weatherDataSettings )
{
    if( weatherDataSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating ground station meteo data: no weather data settings were provided." );
    }

    switch( weatherDataSettings->getWeatherDataType( ) )
    {
        case dsn_weather_data: {
            std::shared_ptr< DsnWeatherDataSettings > dsnWeatherDataSettings =
                    std::dynamic_pointer_cast< DsnWeatherDataSettings >( weatherDataSettings );
            if( dsnWeatherDataSettings == nullptr || dsnWeatherDataSettings->getWeatherData( ) == nullptr )
            {
                throw std::runtime_error( "Error when creating DSN meteo data: incompatible or empty weather data settings." );
            }

            std::map< ground_stations::MeteoDataEntries, int > dsnMeteoEntries = { { ground_stations::temperature_meteo_data, 1 },
                                                                                   { ground_stations::pressure_meteo_data, 2 },
                                                                                   { ground_stations::water_vapor_pressure_meteo_data, 3 },
                                                                                   { ground_stations::relative_humidity_meteo_data, 4 },
                                                                                   { ground_stations::dew_point_meteo_data, 0 } };

            return std::make_shared< ground_stations::ContinuousInterpolatedMeteoData >(
                    interpolators::createOneDimensionalInterpolator( dsnWeatherDataSettings->getWeatherData( )->meteoDataMap_,
                                                                     dsnWeatherDataSettings->getInterpolatorSettings( ) ),
                    dsnMeteoEntries );
        }
        case estrack_weather_data: {
            std::shared_ptr< EstrackWeatherDataSettings > estrackWeatherDataSettings =
                    std::dynamic_pointer_cast< EstrackWeatherDataSettings >( weatherDataSettings );
            if( estrackWeatherDataSettings == nullptr || estrackWeatherDataSettings->getWeatherData( ) == nullptr )
            {
                throw std::runtime_error( "Error when creating ESTRACK meteo data: incompatible or empty weather data settings." );
            }

            std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > > interpolators;
            std::vector< std::map< double, Eigen::VectorXd > > meteoDataList =
                    estrackWeatherDataSettings->getWeatherData( )->getMeteoDataPerFile( );
            for( unsigned int i = 0; i < meteoDataList.size( ); i++ )
            {
                interpolators.push_back( interpolators::createOneDimensionalInterpolator(
                        meteoDataList.at( i ), estrackWeatherDataSettings->getInterpolatorSettings( ) ) );
            }

            std::map< ground_stations::MeteoDataEntries, int > estrackMeteoEntries = { { ground_stations::temperature_meteo_data, 2 },
                                                                                       { ground_stations::pressure_meteo_data, 1 },
                                                                                       { ground_stations::relative_humidity_meteo_data,
                                                                                         0 } };

            return std::make_shared< ground_stations::PiecewiseInterpolatedMeteoData >( interpolators, estrackMeteoEntries );
        }
        default:
            throw std::runtime_error( "Error when creating ground station meteo data: weather data settings type not recognized." );
    }
}

void setWeatherDataInGroundStation( const std::shared_ptr< ground_stations::GroundStation > groundStation,
                                    const std::shared_ptr< WeatherDataSettings > weatherDataSettings )
{
    if( groundStation == nullptr )
    {
        throw std::runtime_error( "Error when setting weather data in ground station: ground station is nullptr." );
    }
    groundStation->setMeteoData( createGroundStationMeteoData( weatherDataSettings ) );
}

std::shared_ptr< GroundStationSettings > getGroundStationSettingsByName(
        const std::vector< std::shared_ptr< GroundStationSettings > >& groundStationSettings,
        const std::string& groundStation )
{
    for( unsigned int i = 0; i < groundStationSettings.size( ); i++ )
    {
        if( groundStationSettings.at( i ) != nullptr && groundStationSettings.at( i )->getStationName( ) == groundStation )
        {
            return groundStationSettings.at( i );
        }
    }

    throw std::runtime_error( "Error when setting weather data in ground station settings: station " + groundStation +
                              " not found in settings list." );
}

void setDsnWeatherDataInGroundStationSettings(
        const std::vector< std::shared_ptr< GroundStationSettings > >& groundStationSettings,
        const std::map< int, std::shared_ptr< input_output::DsnWeatherData > >& weatherDataPerComplex,
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
        const std::map< int, std::vector< std::string > >& groundStationsPerComplex )
{
    for( auto weatherDataIt = weatherDataPerComplex.begin( ); weatherDataIt != weatherDataPerComplex.end( ); ++weatherDataIt )
    {
        int dsnComplex = weatherDataIt->first;
        std::shared_ptr< input_output::DsnWeatherData > weatherData = weatherDataIt->second;

        std::vector< std::string > groundStations;
        if( groundStationsPerComplex.count( dsnComplex ) )
        {
            groundStations = groundStationsPerComplex.at( dsnComplex );
        }
        else
        {
            throw std::runtime_error( "Error when setting weather data in ground station settings: no ground stations in complex." );
        }

        std::shared_ptr< DsnWeatherDataSettings > weatherDataSettings =
                std::make_shared< DsnWeatherDataSettings >( weatherData, interpolatorSettings );
        for( const std::string& groundStation : groundStations )
        {
            getGroundStationSettingsByName( groundStationSettings, groundStation )->setWeatherDataSettings( weatherDataSettings );
        }
    }
}

void setDsnWeatherDataInGroundStationSettings( const std::vector< std::shared_ptr< GroundStationSettings > >& groundStationSettings,
                                               const std::vector< std::string >& weatherFileNames,
                                               std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings,
                                               const std::map< int, std::vector< std::string > >& groundStationsPerComplex )
{
    setDsnWeatherDataInGroundStationSettings( groundStationSettings,
                                              input_output::readDsnWeatherDataFiles( weatherFileNames ),
                                              interpolatorSettings,
                                              groundStationsPerComplex );
}

void setEstrackWeatherDataInGroundStationSettings( const std::vector< std::shared_ptr< GroundStationSettings > >& groundStationSettings,
                                                   const std::shared_ptr< input_output::EstrackWeatherData > weatherData,
                                                   const std::string groundStation,
                                                   std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    getGroundStationSettingsByName( groundStationSettings, groundStation )
            ->setWeatherDataSettings( std::make_shared< EstrackWeatherDataSettings >( weatherData, interpolatorSettings ) );
}

void setEstrackWeatherDataInGroundStationSettings( const std::vector< std::shared_ptr< GroundStationSettings > >& groundStationSettings,
                                                   const std::vector< std::string >& weatherFileNames,
                                                   const std::string groundStation,
                                                   std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    setEstrackWeatherDataInGroundStationSettings( groundStationSettings,
                                                  std::make_shared< input_output::EstrackWeatherData >( weatherFileNames ),
                                                  groundStation,
                                                  interpolatorSettings );
}

//! Function to create a ground station from pre-defined station state object, and add it to a Body object
void createGroundStation( const std::shared_ptr< Body > body,
                          const std::string groundStationName,
                          const std::shared_ptr< ground_stations::GroundStationState > groundStationState )
{
    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            std::make_shared< ground_stations::PointingAnglesCalculator >(
                    std::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame,
                               body->getRotationalEphemeris( ),
                               std::placeholders::_1 ),
                    std::bind( &ground_stations::GroundStationState::getRotationFromBodyFixedToTopocentricFrame,
                               groundStationState,
                               std::placeholders::_1 ) );
    body->addGroundStation(
            groundStationName,
            std::make_shared< ground_stations::GroundStation >( groundStationState, pointingAnglesCalculator, groundStationName ) );
}

std::shared_ptr< ground_stations::StationMotionModel > createGroundStationMotionModels(
        const std::shared_ptr< Body > body,
        const std::vector< std::shared_ptr< GroundStationMotionSettings > > stationMotionSettings =
                std::vector< std::shared_ptr< GroundStationMotionSettings > >( ),
        const simulation_setup::SystemOfBodies& bodies = simulation_setup::SystemOfBodies( ) )
{
    std::vector< std::shared_ptr< ground_stations::StationMotionModel > > stationMotionModelList;

    for( unsigned int i = 0; i < stationMotionSettings.size( ); i++ )
    {
        std::shared_ptr< ground_stations::StationMotionModel > currentStationMotionModel;
        switch( stationMotionSettings.at( i )->getModelType( ) )
        {
            case station_history_motion: {
                const auto settings =
                        std::dynamic_pointer_cast< StationHistoryGroundStationMotionSettings >( stationMotionSettings.at( i ) );
                if( !settings ) throw std::runtime_error( "Incompatible station history motion settings" );
                std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter;
                if( const auto rotation =
                            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( body->getRotationalEphemeris( ) ) )
                    converter = rotation->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( );
                currentStationMotionModel = std::make_shared< ground_stations::StationHistoryMotionModel >(
                        settings->trajectory_, settings->timeFunction_, converter );
                break;
            }
            case body_deformation_station_motion: {
                std::shared_ptr< BodyDeformationStationMotionSettings > bodyDeformationModelSettings =
                        std::dynamic_pointer_cast< BodyDeformationStationMotionSettings >( stationMotionSettings.at( i ) );
                if( bodyDeformationModelSettings == nullptr )
                {
                    throw std::runtime_error( "Error when making body deformation station motion model, settings type is incompatible" );
                }
                currentStationMotionModel = std::make_shared< ground_stations::BodyDeformationStationMotionModel >(
                        std::bind( &Body::getBodyDeformationModelsReference, body ),
                        bodyDeformationModelSettings->throwExceptionWhenNotAvailable_ );
                break;
            }
            case linear_station_motion: {
                std::shared_ptr< LinearGroundStationMotionSettings > linearStationMotionSettings =
                        std::dynamic_pointer_cast< LinearGroundStationMotionSettings >( stationMotionSettings.at( i ) );
                if( linearStationMotionSettings == nullptr )
                {
                    throw std::runtime_error( "Error when making linear ground station motion model, settings type is incompatible" );
                }
                currentStationMotionModel = std::make_shared< ground_stations::LinearStationMotionModel >(
                        linearStationMotionSettings->linearVelocity_, linearStationMotionSettings->referenceEpoch_ );
                break;
            }
            case piecewise_constant_station_motion: {
                std::shared_ptr< PiecewiseConstantGroundStationMotionSettings > piecewiseConstantStationMotionSettings =
                        std::dynamic_pointer_cast< PiecewiseConstantGroundStationMotionSettings >( stationMotionSettings.at( i ) );
                if( piecewiseConstantStationMotionSettings == nullptr )
                {
                    throw std::runtime_error(
                            "Error when making piecewise constant ground station motion model, settings type is incompatible" );
                }
                currentStationMotionModel = std::make_shared< ground_stations::PiecewiseConstantStationMotionModel >(
                        piecewiseConstantStationMotionSettings->displacementList_ );
                break;
            }
            case bodycentric_to_barycentric_station_position_motion: {
                std::shared_ptr< BodyCentricToBarycentricGroundStationMotionSettings > relativisticStationMotionSettings =
                        std::dynamic_pointer_cast< BodyCentricToBarycentricGroundStationMotionSettings >( stationMotionSettings.at( i ) );
                if( relativisticStationMotionSettings == nullptr )
                {
                    throw std::runtime_error(
                            "Error when making relativistic (body-centered to barycentered) ground station motion model, settings type is "
                            "incompatible" );
                }

                std::function< Eigen::Vector6d( const Time& ) > bodyBarycentricStateFunction;
                if( bodies.getFrameOrigin( ) == "SSB" )
                {
                    bodyBarycentricStateFunction =
                            std::bind( &Body::getStateInBaseFrameFromEphemeris< double, Time >, body, std::placeholders::_1 );
                }
                else
                {
                    throw std::runtime_error(
                            "Error when getting body state function for (body-centered to barycentered) ground station motion model, only "
                            "SSB global frame origin is currently supported" );
                }

                auto rotation = body->getRotationalEphemeris( );
                if( !rotation ) throw std::invalid_argument( "Relativistic station correction requires a rotation model." );
                const auto sameAxes = []( const std::string& first, const std::string& second ) {
                    const std::set< std::string > icrfAxes{ "GCRS", "ICRF", "J2000" };
                    return first == second || ( icrfAxes.count( first ) && icrfAxes.count( second ) );
                };
                if( !sameAxes( rotation->getBaseFrameOrientation( ), bodies.getFrameOrientation( ) ) )
                    throw std::invalid_argument( "Relativistic station rotation base axes must match the bodies environment axes." );
                const auto validateStateAxes = [ & ]( const std::shared_ptr< Body >& source ) {
                    if( source->getEphemeris( ) &&
                        !sameAxes( source->getEphemeris( )->getReferenceFrameOrientation( ), bodies.getFrameOrientation( ) ) )
                        throw std::invalid_argument(
                                "Relativistic station barycentric ephemeris axes must match the bodies environment axes." );
                    if( source->getIsBodyGlobalFrameOrigin( ) != 1 && source->getEphemeris( ) &&
                        source->getEphemeris( )->getReferenceFrameOrigin( ) != bodies.getFrameOrigin( ) &&
                        source->getEphemerisFrameToBaseFrame( )->getBaseFrameId( ).empty( ) )
                        throw std::invalid_argument( "Relativistic station ephemeris origin conversion has not been configured." );
                };
                validateStateAxes( body );
                if( const auto earth = std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( rotation ) )
                {
                    rotation = std::make_shared< ephemerides::GcrsToItrsRotationModel >(
                            earth->getAnglesCalculator( ), basic_astrodynamics::tdb_scale, earth->getBaseFrameOrientation( ) );
                }
                std::vector< ground_stations::RelativisticStationExternalBody > externalBodies;
                std::set< std::shared_ptr< Body > > selectedExternalBodies;
                for( const auto& name : relativisticStationMotionSettings->externalBodyNames_ )
                {
                    if( !bodies.count( name ) ) throw std::invalid_argument( "Relativistic station external body not found: " + name );
                    const auto external = bodies.at( name );
                    if( external == body )
                        throw std::invalid_argument( "Relativistic external potential must exclude the station host: " + name );
                    if( !selectedExternalBodies.insert( external ).second )
                        throw std::invalid_argument( "Duplicate external body in relativistic station potential: " + name );
                    if( !external->getGravityFieldModel( ) ) throw std::invalid_argument( "Relativistic external body lacks GM: " + name );
                    validateStateAxes( external );
                    externalBodies.push_back(
                            { name,
                              std::bind( &Body::getStateInBaseFrameFromEphemeris< double, Time >, external, std::placeholders::_1 ),
                              std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                         external->getGravityFieldModel( ) ) } );
                }
                currentStationMotionModel = std::make_shared< ground_stations::BodyCentricToBarycentricRelativisticStationMotion >(
                        bodyBarycentricStateFunction,
                        std::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrameFromExtendedTime,
                                   rotation,
                                   std::placeholders::_1 ),
                        std::bind( &ephemerides::RotationalEphemeris::getDerivativeOfRotationToTargetFrameFromExtendedTime,
                                   rotation,
                                   std::placeholders::_1 ),
                        externalBodies,
                        relativisticStationMotionSettings->coordinateCompatibility_ );
                break;
            }
            case custom_station_motion: {
                std::shared_ptr< CustomGroundStationMotionSettings > customStationMotionSettings =
                        std::dynamic_pointer_cast< CustomGroundStationMotionSettings >( stationMotionSettings.at( i ) );
                if( customStationMotionSettings == nullptr )
                {
                    throw std::runtime_error( "Error when making custom ground station motion model, settings type is incompatible" );
                }
                currentStationMotionModel = std::make_shared< ground_stations::CustomStationMotionModel >(
                        customStationMotionSettings->customDisplacementModel_ );
                break;
            }
            default:
                throw std::runtime_error( "Error when making ground station motion model, settings type not recognized" );
        }
        stationMotionModelList.push_back( currentStationMotionModel );
    }
    return std::make_shared< ground_stations::CombinedStationMotionModel >( stationMotionModelList );
}

std::shared_ptr< ground_stations::GroundStationState > createGroundStationState(
        const std::shared_ptr< Body > body,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType,
        const std::vector< std::shared_ptr< GroundStationMotionSettings > >& stationMotionSettings,
        const simulation_setup::SystemOfBodies& bodies )
{
    std::shared_ptr< ground_stations::StationMotionModel > stationMotionModel =
            createGroundStationMotionModels( body, stationMotionSettings, bodies );
    return std::make_shared< ground_stations::GroundStationState >(
            groundStationPosition, positionElementType, body->getShapeModel( ), stationMotionModel );
}

std::shared_ptr< ground_stations::GroundStationState > createGroundStationState(
        const std::shared_ptr< Body > body,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType,
        const simulation_setup::SystemOfBodies& bodies )
{
    return createGroundStationState(
            body, groundStationPosition, positionElementType, std::vector< std::shared_ptr< GroundStationMotionSettings > >( ), bodies );
}

//! Function to create a ground station and add it to a Body object
void createGroundStation( const std::shared_ptr< Body > body,
                          const std::string groundStationName,
                          const Eigen::Vector3d groundStationPosition,
                          const coordinate_conversions::PositionElementTypes positionElementType,
                          const std::vector< std::shared_ptr< GroundStationMotionSettings > > stationMotionSettings,
                          const SystemOfBodies& bodies )
{
    createGroundStation( body,
                         groundStationName,
                         createGroundStationState( body, groundStationPosition, positionElementType, stationMotionSettings, bodies ) );
}

//! Function to create a set of ground stations and add them to the corresponding Body objects
void createGroundStations( const SystemOfBodies& bodies,
                           const std::map< std::pair< std::string, std::string >, Eigen::Vector3d >& groundStationsWithPosition,
                           const coordinate_conversions::PositionElementTypes positionElementType )
{
    for( std::map< std::pair< std::string, std::string >, Eigen::Vector3d >::const_iterator stationIterator =
                 groundStationsWithPosition.begin( );
         stationIterator != groundStationsWithPosition.end( );
         stationIterator++ )
    {
        if( bodies.count( stationIterator->first.first ) > 0 )
        {
            createGroundStation( bodies.at( stationIterator->first.first ),
                                 stationIterator->first.second,
                                 stationIterator->second,
                                 positionElementType );
        }
    }
}

void createGroundStation( const std::shared_ptr< Body > body,
                          const std::shared_ptr< GroundStationSettings > groundStationSettings,
                          const SystemOfBodies& bodies )
{
    if( body->getGroundStationMap( ).count( groundStationSettings->getStationName( ) ) != 0 )
    {
        throw std::runtime_error( "Error when creating ground station " + groundStationSettings->getStationName( ) + " on body " +
                                  body->getBodyName( ) + ", station already exists." );
    }
    else
    {
        createGroundStation( body,
                             groundStationSettings->getStationName( ),
                             groundStationSettings->getGroundStationPosition( ),
                             groundStationSettings->getPositionElementType( ),
                             groundStationSettings->getStationMotionSettings( ),
                             bodies );
        if( groundStationSettings->getWeatherDataSettings( ) != nullptr )
        {
            setWeatherDataInGroundStation( body->getGroundStation( groundStationSettings->getStationName( ) ),
                                           groundStationSettings->getWeatherDataSettings( ) );
        }
    }
}

std::vector< std::pair< std::string, std::string > > getGroundStationsLinkEndList( const std::shared_ptr< Body > body )
{
    std::vector< std::pair< std::string, std::string > > stationList;

    std::map< std::string, std::shared_ptr< ground_stations::GroundStation > > groundStationsMap = body->getGroundStationMap( );
    for( auto stationIterator : groundStationsMap )
    {
        stationList.push_back( std::make_pair( body->getBodyName( ), stationIterator.first ) );
    }
    return stationList;
}
//
// std::vector< double >  getTargetElevationAngles(
//        const std::shared_ptr< Body > observingBody,
//        const std::shared_ptr< Body > targetBody,
//        const std::string groundStationName,
//        const std::vector< double > times )
//{
//    if( observingBody->getGroundStationMap( ).count( groundStationName ) == 0 )
//    {
//        throw std::runtime_error( "Error when computing elevating angle, station " + groundStationName +
//                                  " not found on body " + observingBody->getBodyName( ) );
//    }
//
//
//    std::function< Eigen::Vector6d( const double& ) > groundStationStateFunction =
//            getLinkEndCompleteEphemerisFunction(
//                observingBody, std::make_pair( observingBody->getBodyName( ), groundStationName ) );
//
//    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
//            observingBody->getGroundStationMap( ).at( groundStationName )->getPointingAnglesCalculator( );
//    Eigen::Vector3d relativePosition;
//    std::vector< double > elevationAngles;
//    for( unsigned int i = 0; i < times.size( ); i++ )
//    {
//        targetBody->getStateInBaseFrameFromEphemeris( times.at( i ) );
//        groundStationStateFunction( times.at( i ) );
//        relativePosition = ( targetBody->getStateInBaseFrameFromEphemeris( times.at( i ) ) -
//                groundStationStateFunction( times.at( i ) ) ).segment( 0, 3 );
//        elevationAngles.push_back( pointingAnglesCalculator->calculateElevationAngleFromInertialVector( relativePosition, times.at( i ) )
//        );
//    }
//    return elevationAngles;
//}

//! Function to determine whether the Sun is below a given elevation angle as seen from a ground station.
// Note: this function is just a convenient wrapper that allows the user to check whether it is night or day at a given ground station.
// It uses isTargetVisibleFromGroundStation and allows the user to avoid extracting the states at a given time for both Sun and Ground
// Station. The same result could be achieved by extracting the states and using those directly inside isTargetVisibleFromGroundStation
// (this would also be faster).
bool isItDarkAtGroundStation( const std::string& groundStationName,
                              const double timeTdb,
                              const double maxSunElevation,
                              const SystemOfBodies& bodies,
                              const std::string& centralBodyName )
{
    std::shared_ptr< Body > centralBody = bodies.getBody( centralBodyName );

    if( centralBody->getGroundStationMap( ).count( groundStationName ) == 0 )
    {
        throw std::runtime_error( "Error when checking darkness at ground station " + groundStationName + " on body " + centralBodyName +
                                  ": station not found." );
    }

    // Inertial state of the ground station, obtained from the central body's ephemeris and rotation model, combined with the
    // station's body-fixed position.
    std::function< Eigen::Vector6d( const double ) > groundStationStateFunction = getLinkEndCompleteEphemerisFunction< double, double >(
            observation_models::LinkEndId( centralBodyName, groundStationName ), bodies );
    const Eigen::Vector3d groundStationPosition = groundStationStateFunction( timeTdb ).segment( 0, 3 );

    // Inertial position of the Sun.
    const Eigen::Vector3d sunPosition = bodies.getBody( "Sun" )->getPositionInBaseFrameFromEphemeris< double, double >( timeTdb );

    const std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            centralBody->getGroundStation( groundStationName )->getPointingAnglesCalculator( );

    return !ground_stations::isTargetVisibleFromGroundStation(
            groundStationPosition, sunPosition, timeTdb, maxSunElevation, pointingAnglesCalculator );
}

}  // namespace simulation_setup

}  // namespace tudat
