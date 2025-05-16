/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readViennaMappingFunctionData.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"
#include "tudat/astro/observation_models/corrections/firstOrderRelativisticCorrection.h"
#include "tudat/astro/observation_models/corrections/solarCoronaCorrection.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/math/basic/coordinateConversions.h"

namespace tudat
{

namespace observation_models
{

//! Function to create object that computes a single (type of) correction to the light-time
std::shared_ptr< LightTimeCorrection > createLightTimeCorrections( const std::shared_ptr< LightTimeCorrectionSettings > correctionSettings,
                                                                   const simulation_setup::SystemOfBodies& bodies,
                                                                   const LinkEnds& linkEnds,
                                                                   const LinkEndType& transmittingLinkEndType,
                                                                   const LinkEndType& receivingLinkEndType,
                                                                   const ObservableType observableType )
{
    using namespace tudat::ephemerides;
    using namespace tudat::gravitation;

    std::shared_ptr< LightTimeCorrection > lightTimeCorrection;

    LinkEndId transmitter = linkEnds.at( transmittingLinkEndType );
    LinkEndId receiver = linkEnds.at( receivingLinkEndType );

    // Identify type of light time correction to be created.
    switch( correctionSettings->getCorrectionType( ) )
    {
        case first_order_relativistic: {
            // Check input consistency
            if( std::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings ) != nullptr )
            {
                // Retrieve list of bodies causing light time perturbation
                std::vector< std::string > perturbingBodies =
                        std::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings )
                                ->getPerturbingBodies( );

                std::vector< std::function< Eigen::Vector6d( const double ) > > perturbingBodyStateFunctions;
                std::vector< std::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

                // Retrieve mass and state functions for each perturbing body.
                for( unsigned int i = 0; i < perturbingBodies.size( ); i++ )
                {
                    if( bodies.count( perturbingBodies[ i ] ) == 0 )
                    {
                        throw std::runtime_error( "Error when making 1st order relativistic light time correction, could not find body " +
                                                  perturbingBodies.at( i ) );
                    }
                    else
                    {
                        // Set state function.
                        perturbingBodyStateFunctions.push_back(
                                std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, double >,
                                           bodies.at( perturbingBodies[ i ] ),
                                           std::placeholders::_1 ) );

                        // Set gravitational parameter function.
                        perturbingBodyGravitationalParameterFunctions.push_back(
                                std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                           bodies.at( perturbingBodies[ i ] )->getGravityFieldModel( ) ) );
                    }
                }

                // Create light-time correction function
                lightTimeCorrection = std::make_shared< FirstOrderLightTimeCorrectionCalculator >(
                        perturbingBodyStateFunctions,
                        perturbingBodyGravitationalParameterFunctions,
                        perturbingBodies,
                        transmitter.bodyName_,
                        receiver.bodyName_,
                        std::bind( &relativity::PPNParameterSet::getParameterGamma, relativity::ppnParameterSet ),
                        std::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings )
                                ->getBendingFlag() );
            }
            else
            {
                throw std::runtime_error( "Error, correction settings type (1st order relativistic) does not coincide with data type." );
            }

            break;
        }
        case tabulated_tropospheric: {
            std::shared_ptr< TabulatedTroposphericCorrectionSettings > troposphericCorrectionSettings =
                    std::dynamic_pointer_cast< TabulatedTroposphericCorrectionSettings >( correctionSettings );
            if( correctionSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating tabulated tropospheric correction: incompatible settings type." );
            }

            if( !isRadiometricObservableType( observableType ) )
            {
                throw std::runtime_error(
                        "Error when creating tabulated tropospheric correction: selected observable type is not radiometric." );
            }

            // If one of the link ends is the body with the atmosphere then create the tropospheric correction
            if( transmitter.bodyName_ != receiver.bodyName_ &&
                ( transmitter.bodyName_ == troposphericCorrectionSettings->getBodyWithAtmosphere( ) ||
                  receiver.bodyName_ == troposphericCorrectionSettings->getBodyWithAtmosphere( ) ) )
            {
                bool isUplinkCorrection;
                LinkEndId groundStation;
                if( transmitter.bodyName_ == troposphericCorrectionSettings->getBodyWithAtmosphere( ) )
                {
                    isUplinkCorrection = true;
                    groundStation = transmitter;
                }
                else
                {
                    isUplinkCorrection = false;
                    groundStation = receiver;
                }

                std::shared_ptr< TroposhericElevationMapping > troposphericElevationMapping =
                        createTroposphericElevationMapping( troposphericCorrectionSettings->getTroposphericMappingModelType( ),
                                                            bodies,
                                                            transmitter,
                                                            receiver,
                                                            isUplinkCorrection );

                ObservableType baseObservableType = getBaseObservableType( observableType );

                std::shared_ptr< TabulatedMediaReferenceCorrectionManager > dryCorrectionCalculator, wetCorrectionCalculator;
                std::pair< std::string, std::string > stationSpacecraftPair = std::make_pair( groundStation.stationName_, "" );

                if( troposphericCorrectionSettings->getTroposphericDryCorrection( ).count( stationSpacecraftPair ) &&
                    troposphericCorrectionSettings->getTroposphericDryCorrection( )
                            .at( stationSpacecraftPair )
                            .count( baseObservableType ) &&
                    troposphericCorrectionSettings->getTroposphericWetCorrection( ).count( stationSpacecraftPair ) &&
                    troposphericCorrectionSettings->getTroposphericWetCorrection( )
                            .at( stationSpacecraftPair )
                            .count( baseObservableType ) &&
                    troposphericCorrectionSettings->getTroposphericDryCorrectionAdjustment( ).count( stationSpacecraftPair ) &&
                    troposphericCorrectionSettings->getTroposphericDryCorrectionAdjustment( )
                            .at( stationSpacecraftPair )
                            .count( baseObservableType ) &&
                    troposphericCorrectionSettings->getTroposphericWetCorrectionAdjustment( ).count( stationSpacecraftPair ) &&
                    troposphericCorrectionSettings->getTroposphericWetCorrectionAdjustment( )
                            .at( stationSpacecraftPair )
                            .count( baseObservableType ) )
                {
                    lightTimeCorrection = std::make_shared< TabulatedTroposphericCorrection >(
                            troposphericCorrectionSettings->getTroposphericDryCorrection( )
                                    .at( stationSpacecraftPair )
                                    .at( baseObservableType ),
                            troposphericCorrectionSettings->getTroposphericWetCorrection( )
                                    .at( stationSpacecraftPair )
                                    .at( baseObservableType ),
                            troposphericCorrectionSettings->getTroposphericDryCorrectionAdjustment( )
                                    .at( stationSpacecraftPair )
                                    .at( baseObservableType ),
                            troposphericCorrectionSettings->getTroposphericWetCorrectionAdjustment( )
                                    .at( stationSpacecraftPair )
                                    .at( baseObservableType ),
                            troposphericElevationMapping,
                            isUplinkCorrection );
                }
                else
                {
                    throw std::runtime_error( "Error when creating tabulated tropospheric corrections for " + groundStation.stationName_ +
                                              "ground station and " + getObservableName( observableType ) +
                                              " observable: tabulated data not available. " );
                }
            }
            // Set correction to nullptr if correction isn't valid for selected link ends
            else
            {
                lightTimeCorrection = nullptr;
            }

            break;
        }
        case saastamoinen_tropospheric: {
            std::shared_ptr< SaastamoinenTroposphericCorrectionSettings > troposphericCorrectionSettings =
                    std::dynamic_pointer_cast< SaastamoinenTroposphericCorrectionSettings >( correctionSettings );
            if( correctionSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating Saastamoinen tropospheric correction: incompatible settings type." );
            }

            if( !isRadiometricObservableType( observableType ) )
            {
                throw std::runtime_error(
                        "Error when creating saastamoinen tropospheric correction: selected observable type is not radiometric." );
            }

            // If one of the link ends is the body with the atmosphere then create the tropospheric correction
            if( transmitter.bodyName_ != receiver.bodyName_ &&
                ( transmitter.bodyName_ == troposphericCorrectionSettings->getBodyWithAtmosphere( ) ||
                  receiver.bodyName_ == troposphericCorrectionSettings->getBodyWithAtmosphere( ) ) )
            {
                bool isUplinkCorrection;
                LinkEndId groundStation;
                if( transmitter.bodyName_ == troposphericCorrectionSettings->getBodyWithAtmosphere( ) )
                {
                    isUplinkCorrection = true;
                    groundStation = transmitter;
                }
                else
                {
                    isUplinkCorrection = false;
                    groundStation = receiver;
                }

                std::shared_ptr< TroposhericElevationMapping > troposphericElevationMapping =
                        createTroposphericElevationMapping( troposphericCorrectionSettings->getTroposphericMappingModelType( ),
                                                            bodies,
                                                            transmitter,
                                                            receiver,
                                                            isUplinkCorrection );

                // Using nominal geodetic position, i.e. ignoring station motion
                std::function< Eigen::Vector3d( double ) > groundStationGeodeticPositionFunction =
                        std::bind( &ground_stations::GroundStationState::getNominalGeodeticPosition,
                                   bodies.getBody( groundStation.bodyName_ )
                                           ->getGroundStation( groundStation.stationName_ )
                                           ->getNominalStationState( ) );

                std::function< double( const double ) > waterVaporPartialPressureFunction;
                if( troposphericCorrectionSettings->getWaterVaporPartialPressureModelType( ) == tabulated )
                {
                    waterVaporPartialPressureFunction = bodies.getBody( groundStation.bodyName_ )
                                                                ->getGroundStation( groundStation.stationName_ )
                                                                ->getWaterVaporPartialPressureFunction( );
                }
                else if( troposphericCorrectionSettings->getWaterVaporPartialPressureModelType( ) == bean_and_dutton )
                {
                    waterVaporPartialPressureFunction = ground_stations::getBeanAndDuttonWaterVaporPartialPressureFunction(
                            bodies.getBody( groundStation.bodyName_ )
                                    ->getGroundStation( groundStation.stationName_ )
                                    ->getRelativeHumidityFunction( ),
                            bodies.getBody( groundStation.bodyName_ )
                                    ->getGroundStation( groundStation.stationName_ )
                                    ->getTemperatureFunction( ) );
                }
                else
                {
                    throw std::runtime_error(
                            "Error when creating Saastamoinen tropospheric correction: water vapor partial "
                            "pressure model type not recognized." );
                }

                lightTimeCorrection = std::make_shared< SaastamoinenTroposphericCorrection >(
                        groundStationGeodeticPositionFunction,
                        bodies.getBody( groundStation.bodyName_ )->getGroundStation( groundStation.stationName_ )->getPressureFunction( ),
                        bodies.getBody( groundStation.bodyName_ )
                                ->getGroundStation( groundStation.stationName_ )
                                ->getTemperatureFunction( ),
                        waterVaporPartialPressureFunction,
                        troposphericElevationMapping,
                        isUplinkCorrection );
            }
            // Set correction to nullptr if correction isn't valid for selected link ends
            else
            {
                lightTimeCorrection = nullptr;
            }

            break;
        }
        case tabulated_ionospheric: {
            std::shared_ptr< TabulatedIonosphericCorrectionSettings > ionosphericCorrectionSettings =
                    std::dynamic_pointer_cast< TabulatedIonosphericCorrectionSettings >( correctionSettings );
            if( correctionSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating tabulated ionospheric correction: incompatible settings type." );
            }

            // If one of the link ends is the body with the atmosphere then create the tropospheric correction
            if( transmitter.bodyName_ != receiver.bodyName_ &&
                ( transmitter.bodyName_ == ionosphericCorrectionSettings->getBodyWithAtmosphere( ) ||
                  receiver.bodyName_ == ionosphericCorrectionSettings->getBodyWithAtmosphere( ) ) )
            {
                bool isUplinkCorrection;
                LinkEndId groundStation, spacecraft;
                if( transmitter.bodyName_ == ionosphericCorrectionSettings->getBodyWithAtmosphere( ) )
                {
                    isUplinkCorrection = true;
                    groundStation = transmitter;
                    spacecraft = receiver;
                }
                else
                {
                    isUplinkCorrection = false;
                    groundStation = receiver;
                    spacecraft = transmitter;
                }

                ObservableType baseObservableType = getBaseObservableType( observableType );

                std::shared_ptr< TabulatedMediaReferenceCorrectionManager > correctionCalculator;
                std::pair< std::string, std::string > stationSpacecraftPair =
                        std::make_pair( groundStation.stationName_, spacecraft.bodyName_ );

                if( ionosphericCorrectionSettings->getReferenceRangeCorrection( ).count( stationSpacecraftPair ) &&
                    ionosphericCorrectionSettings->getReferenceRangeCorrection( ).at( stationSpacecraftPair ).count( baseObservableType ) )
                {
                    lightTimeCorrection = std::make_shared< TabulatedIonosphericCorrection >(
                            ionosphericCorrectionSettings->getReferenceRangeCorrection( )
                                    .at( stationSpacecraftPair )
                                    .at( baseObservableType ),
                            createLinkFrequencyFunction( bodies, linkEnds, transmittingLinkEndType, receivingLinkEndType ),
                            baseObservableType,
                            isUplinkCorrection,
                            ionosphericCorrectionSettings->getReferenceFrequency( ) );
                }
                else
                {
                    throw std::runtime_error( "Error when creating tabulated ionospheric corrections for " + groundStation.stationName_ +
                                              " ground station, " + spacecraft.bodyName_ + " spacecraft, and " +
                                              getObservableName( observableType ) + " observable: tabulated data not available. " );
                }
            }
            // Set correction to nullptr if correction isn't valid for selected link ends
            else
            {
                lightTimeCorrection = nullptr;
            }

            break;
        }
        case jakowski_vtec_ionospheric: {
            std::shared_ptr< JakowskiIonosphericCorrectionSettings > ionosphericCorrectionSettings =
                    std::dynamic_pointer_cast< JakowskiIonosphericCorrectionSettings >( correctionSettings );
            if( correctionSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating Jakowski VTEC ionospheric correction: incompatible settings type." );
            }

            // If one of the link ends is the body with the atmosphere then create the tropospheric correction
            if( transmitter.bodyName_ != receiver.bodyName_ &&
                ( transmitter.bodyName_ == ionosphericCorrectionSettings->getBodyWithAtmosphere( ) ||
                  receiver.bodyName_ == ionosphericCorrectionSettings->getBodyWithAtmosphere( ) ) )
            {
                bool isUplinkCorrection;
                LinkEndId groundStation, spacecraft;
                if( transmitter.bodyName_ == ionosphericCorrectionSettings->getBodyWithAtmosphere( ) )
                {
                    isUplinkCorrection = true;
                    groundStation = transmitter;
                    spacecraft = receiver;
                }
                else
                {
                    isUplinkCorrection = false;
                    groundStation = receiver;
                    spacecraft = transmitter;
                }

                ObservableType baseObservableType = getBaseObservableType( observableType );

                std::shared_ptr< TabulatedMediaReferenceCorrectionManager > correctionCalculator;
                std::pair< std::string, std::string > stationSpacecraftPair =
                        std::make_pair( groundStation.stationName_, spacecraft.bodyName_ );

                // Create Jakowski VTEC calculator

                std::shared_ptr< input_output::solar_activity::SolarActivityContainer > solarActivityContainer =
                        std::make_shared< input_output::solar_activity::SolarActivityContainer >(
                                ionosphericCorrectionSettings->getSolarActivityData( ) );
                std::function< double( double ) > flux10p7Function = [ = ]( double time ) {
                    return solarActivityContainer->getSolarActivityData( time )->solarRadioFlux107Observed;
                };

                // Elevation computed wrt to body-fixed frame of body with ground station
                // In principle the elevation is defined wrt a body-centered inertial frame, but since this frame has the same
                // z axis as the body-fixed frame we can use either of them to compute the elevation. This wouldn't be valid
                // for computing the right ascension though.
                std::function< double( double ) > sunElevationFunction = [ = ]( double time ) {
                    Eigen::Vector3d bodyFixedCartesianRelativePosition =
                            bodies.getBody( groundStation.bodyName_ )->getRotationalEphemeris( )->getRotationMatrixToTargetFrame( time ) *
                            ( bodies.getBody( "Sun" )->getStateInBaseFrameFromEphemeris( time ) -
                              bodies.getBody( groundStation.bodyName_ )->getStateInBaseFrameFromEphemeris( time ) )
                                    .segment( 0, 3 );
                    Eigen::Vector3d bodyFixedSphericalRelativePosition =
                            coordinate_conversions::convertCartesianToSpherical( bodyFixedCartesianRelativePosition );
                    return mathematical_constants::PI / 2.0 - bodyFixedSphericalRelativePosition.y( );
                };

                std::shared_ptr< JakowskiVtecCalculator > vtecCalculator =
                        std::make_shared< JakowskiVtecCalculator >( sunElevationFunction,
                                                                    flux10p7Function,
                                                                    ionosphericCorrectionSettings->getUseUtcTimeForLocalTime( ),
                                                                    ionosphericCorrectionSettings->getGeomagneticPoleLatitude( ),
                                                                    ionosphericCorrectionSettings->getGeomagneticPoleLongitude( ),
                                                                    ionosphericCorrectionSettings->getIonosphereHeight( ) );

                // Create ionospheric correction

                std::function< double( Eigen::Vector3d, double ) > elevationFunction =
                        std::bind( &ground_stations::PointingAnglesCalculator::calculateElevationAngleFromInertialVector,
                                   bodies.getBody( groundStation.bodyName_ )
                                           ->getGroundStation( groundStation.stationName_ )
                                           ->getPointingAnglesCalculator( ),
                                   std::placeholders::_1,
                                   std::placeholders::_2 );
                std::function< double( Eigen::Vector3d, double ) > azimuthFunction =
                        std::bind( &ground_stations::PointingAnglesCalculator::calculateAzimuthAngleFromInertialVector,
                                   bodies.getBody( groundStation.bodyName_ )
                                           ->getGroundStation( groundStation.stationName_ )
                                           ->getPointingAnglesCalculator( ),
                                   std::placeholders::_1,
                                   std::placeholders::_2 );

                // Nominal geodetic position, i.e. ignoring station motion
                std::function< Eigen::Vector3d( double ) > groundStationGeodeticPositionFunction =
                        std::bind( &ground_stations::GroundStationState::getNominalGeodeticPosition,
                                   bodies.getBody( groundStation.bodyName_ )
                                           ->getGroundStation( groundStation.stationName_ )
                                           ->getNominalStationState( ) );

                // Get equatorial radius
                double equatorialRadius;
                std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel =
                        bodies.getBody( groundStation.bodyName_ )->getShapeModel( );
                if( std::dynamic_pointer_cast< basic_astrodynamics::SphericalBodyShapeModel >( shapeModel ) != nullptr )
                {
                    equatorialRadius =
                            std::dynamic_pointer_cast< basic_astrodynamics::SphericalBodyShapeModel >( shapeModel )->getAverageRadius( );
                }
                else if( std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( shapeModel ) != nullptr )
                {
                    equatorialRadius = std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( shapeModel )
                                               ->getEquatorialRadius( );
                }
                else
                {
                    throw std::runtime_error( "Error when creating Jakowski ionospheric corrections for body " + groundStation.bodyName_ +
                                              ": shape model not recognized" );
                }

                lightTimeCorrection = std::make_shared< MappedVtecIonosphericCorrection >(
                        vtecCalculator,
                        createLinkFrequencyFunction( bodies, linkEnds, transmittingLinkEndType, receivingLinkEndType ),
                        elevationFunction,
                        azimuthFunction,
                        groundStationGeodeticPositionFunction,
                        baseObservableType,
                        isUplinkCorrection,
                        equatorialRadius,
                        ionosphericCorrectionSettings->getFirstOrderDelayCoefficient( ) );
            }
            // Set correction to nullptr if correction isn't valid for selected link ends
            else
            {
                lightTimeCorrection = nullptr;
            }

            break;
        }
        case inverse_power_series_solar_corona: {
            std::shared_ptr< InversePowerSeriesSolarCoronaCorrectionSettings > coronaCorrectionSettings =
                    std::dynamic_pointer_cast< InversePowerSeriesSolarCoronaCorrectionSettings >( correctionSettings );
            if( correctionSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating inverse power series solar corona correction: incompatible settings type." );
            }

            std::string sunBodyName = "Sun";  // coronaCorrectionSettings->getSunBodyName( );

            std::function< Eigen::Vector6d( const double ) > sunStateFunction =
                    std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, double >,
                               bodies.at( sunBodyName ),
                               std::placeholders::_1 );

            std::shared_ptr< basic_astrodynamics::BodyShapeModel > sunShapeModel = bodies.at( sunBodyName )->getShapeModel( );
            if( sunShapeModel == nullptr )
            {
                throw std::runtime_error(
                        "Error when creating inverse power series solar corona correction: no shape model "
                        "was found for the Sun" );
            }

            lightTimeCorrection = std::make_shared< InversePowerSeriesSolarCoronaCorrection >(
                    observableType,
                    sunStateFunction,
                    coronaCorrectionSettings->getCoefficients( ),
                    coronaCorrectionSettings->getPositiveExponents( ),
                    coronaCorrectionSettings->getCriticalPlasmaDensityDelayCoefficient( ),
                    sunShapeModel->getAverageRadius( ) );

            break;
        }
        default: {
            std::string errorMessage =
                    "Error, light time correction type " + std::to_string( correctionSettings->getCorrectionType( ) ) + " not recognized.";
            throw std::runtime_error( errorMessage );

            break;
        }
    }
    return lightTimeCorrection;
}

std::shared_ptr< TroposhericElevationMapping > createTroposphericElevationMapping(
        const TroposphericMappingModel troposphericMappingModelType,
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEndId& transmitter,
        const LinkEndId& receiver,
        const bool isUplinkCorrection )
{
    std::shared_ptr< TroposhericElevationMapping > troposphericMappingModel;

    LinkEndId groundStation;
    if( isUplinkCorrection )
    {
        groundStation = transmitter;
    }
    else
    {
        groundStation = receiver;
    }

    switch( troposphericMappingModelType )
    {
        case simplified_chao: {
            std::function< double( Eigen::Vector3d, double ) > elevationFunction =
                    std::bind( &ground_stations::PointingAnglesCalculator::calculateElevationAngleFromInertialVector,
                               bodies.getBody( groundStation.bodyName_ )
                                       ->getGroundStation( groundStation.stationName_ )
                                       ->getPointingAnglesCalculator( ),
                               std::placeholders::_1,
                               std::placeholders::_2 );

            troposphericMappingModel = std::make_shared< SimplifiedChaoTroposphericMapping >( elevationFunction, isUplinkCorrection );

            break;
        }
        case niell: {
            std::function< double( Eigen::Vector3d, double ) > elevationFunction =
                    std::bind( &ground_stations::PointingAnglesCalculator::calculateElevationAngleFromInertialVector,
                               bodies.getBody( groundStation.bodyName_ )
                                       ->getGroundStation( groundStation.stationName_ )
                                       ->getPointingAnglesCalculator( ),
                               std::placeholders::_1,
                               std::placeholders::_2 );

            // Using nominal geodetic position, i.e. ignoring station motion
            std::function< Eigen::Vector3d( double ) > groundStationGeodeticPositionFunction = std::bind(
                    &ground_stations::GroundStationState::getNominalGeodeticPosition,
                    bodies.getBody( groundStation.bodyName_ )->getGroundStation( groundStation.stationName_ )->getNominalStationState( ) );

            troposphericMappingModel = std::make_shared< NiellTroposphericMapping >(
                    elevationFunction, groundStationGeodeticPositionFunction, isUplinkCorrection );

            break;
        }
        default:
            throw std::runtime_error( "Error when creating tropospheric elevation mapping model: model type " +
                                      std::to_string( troposphericMappingModelType ) + " not recognized." );

            break;
    }

    return troposphericMappingModel;
}

void setVmfTroposphereCorrections( const std::vector< std::string >& dataFiles,
                                   const bool fileHasMeteo,
                                   const bool fileHasGradient,
                                   const simulation_setup::SystemOfBodies& bodies,
                                   const bool setTropospherData,
                                   const bool setMeteoData,
                                   const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    std::map< std::string, input_output::VMFData > vmfData;
    input_output::readVMFFiles( dataFiles, vmfData, fileHasMeteo, fileHasGradient );

    for( auto it: vmfData )
    {
        std::string currentStationName = it.first;
        if( bodies.at( "Earth" )->getGroundStationMap( ).count( currentStationName ) > 0 )
        {
            std::shared_ptr< ground_stations::GroundStation > currentStation =
                    bodies.at( "Earth" )->getGroundStationMap( ).at( currentStationName );

            std::map< double, Eigen::VectorXd > processedTroposphereData;
            std::map< double, Eigen::VectorXd > processedMeteoData;
            it.second.getFullDataSet( processedTroposphereData, processedMeteoData );

            if( setMeteoData )
            {
                std::map< ground_stations::MeteoDataEntries, int > vmfMeteoEntries = { { ground_stations::temperature_meteo_data, 1 },
                                                                                       { ground_stations::pressure_meteo_data, 0 },
                                                                                       { ground_stations::water_vapor_pressure_meteo_data,
                                                                                         2 } };

                std::shared_ptr< ground_stations::StationMeteoData > meteoData =
                        std::make_shared< ground_stations::ContinuousInterpolatedMeteoData >(
                                interpolators::createOneDimensionalInterpolator( processedMeteoData, interpolatorSettings ),
                                vmfMeteoEntries );

                currentStation->setMeteoData( meteoData );
            }

            if( setTropospherData )
            {
                currentStation->setTroposphereData( std::make_shared< ground_stations::InterpolatedStationTroposphereData >(
                        interpolators::createOneDimensionalInterpolator( processedTroposphereData, interpolatorSettings ),
                        fileHasMeteo,
                        fileHasGradient ) );
            }
        }
    }
}

std::function< double( std::vector< FrequencyBands >, double ) > createLinkFrequencyFunction(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds& linkEnds,
        const LinkEndType& transmittingLinkEndType,
        const LinkEndType& receivingLinkEndType )
{
    std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittedFrequencyCalculator =
            bodies.getBody( linkEnds.at( transmitter ).bodyName_ )
                    ->getGroundStation( linkEnds.at( transmitter ).stationName_ )
                    ->getTransmittingFrequencyCalculator( );

    std::vector< std::function< double( FrequencyBands, FrequencyBands ) > > turnaroundRatioFunctions;

    for( auto retransmitterLinkEndsIt = ++linkEnds.begin( ); retransmitterLinkEndsIt->first != receivingLinkEndType;
         ++retransmitterLinkEndsIt )
    {
        std::shared_ptr< system_models::VehicleSystems > vehicleSystems;
        // Check if retransmitter is a body
        if( retransmitterLinkEndsIt->second.stationName_ == "" ||
            !simulation_setup::isReferencePointGroundStation(
                    bodies, retransmitterLinkEndsIt->second.bodyName_, retransmitterLinkEndsIt->second.stationName_ ) )
        {
            vehicleSystems = bodies.getBody( retransmitterLinkEndsIt->second.bodyName_ )->getVehicleSystems( );
        }
        // If retransmitter is a ground station of the body
        else
        {
            vehicleSystems = bodies.getBody( retransmitterLinkEndsIt->second.bodyName_ )
                                     ->getGroundStation( retransmitterLinkEndsIt->second.stationName_ )
                                     ->getVehicleSystems( );
        }

        if( vehicleSystems == nullptr )
        {
            throw std::runtime_error( "Error when creating link frequency function: vehicle systems are not defined for " +
                                      retransmitterLinkEndsIt->second.bodyName_ + "/" + retransmitterLinkEndsIt->second.stationName_ +
                                      " link end." );
        }

        turnaroundRatioFunctions.push_back( vehicleSystems->getTransponderTurnaroundRatio( ) );
    }

    std::function< double( std::vector< FrequencyBands >, double ) > linkFrequencyFunction =
            [ = ]( const std::vector< FrequencyBands >& frequencyBands, const double time ) {
                double frequency = transmittedFrequencyCalculator->getTemplatedCurrentFrequency< double, double >( time );

                if( frequencyBands.size( ) + 1 < turnaroundRatioFunctions.size( ) )
                {
                    throw std::runtime_error( "Error when computing link frequency, found incompatible number of frequency bands: " +
                                              std::to_string( frequencyBands.size( ) ) + " and turnaround ratios " +
                                              std::to_string( turnaroundRatioFunctions.size( ) ) );
                }
                for( unsigned int i = 0; i < turnaroundRatioFunctions.size( ); ++i )
                {
                    frequency *= turnaroundRatioFunctions.at( i )( frequencyBands.at( i ), frequencyBands.at( i + 1 ) );
                }

                return frequency;
            };

    return linkFrequencyFunction;
}

}  // namespace observation_models

}  // namespace tudat
