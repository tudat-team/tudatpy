/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"
#include "tudat/astro/observation_models/corrections/firstOrderRelativisticCorrection.h"
#include "tudat/astro/relativity/metric.h"

namespace tudat
{

namespace observation_models
{

//! Function to create object that computes a single (type of) correction to the light-time
std::shared_ptr< LightTimeCorrection > createLightTimeCorrections(
        const std::shared_ptr< LightTimeCorrectionSettings > correctionSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds& linkEnds,
        const LinkEndType& transmittingLinkEndType,
        const LinkEndType& receivingLinkEndType,
        const ObservableType observableType )
{

    using namespace tudat::ephemerides;
    using namespace tudat::gravitation;

    std::shared_ptr< LightTimeCorrection > lightTimeCorrection;

    LinkEndId transmitter= linkEnds.at( transmittingLinkEndType );
    LinkEndId receiver = linkEnds.at( receivingLinkEndType );

    if ( observableType == two_way_doppler && requiresMultiLegIterations( correctionSettings->getCorrectionType( ) ) )
    {
        throw std::runtime_error(
                "Error when creating 2-way Doppler light time corrections: selected correction (" +
                getLightTimeCorrectionName( correctionSettings->getCorrectionType( ) ) +
                ") requires multi-leg light iterations, which aren't implemented for 2-way Doppler." );
    }

    // Identify type of light time correction to be created.
    switch( correctionSettings->getCorrectionType( ) )
    {
    case first_order_relativistic:
    {
        // Check input consistency
        if( std::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings ) != nullptr )
        {
            // Retrieve list of bodies causing light time perturbation
            std::vector< std::string > perturbingBodies =
                    std::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings )->
                    getPerturbingBodies( );

            std::vector< std::function< Eigen::Vector6d( const double ) > > perturbingBodyStateFunctions;
            std::vector< std::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

            // Retrieve mass and state functions for each perturbing body.
            for( unsigned int i = 0; i < perturbingBodies.size( ); i++ )
            {
                if( bodies.count( perturbingBodies[ i ] ) == 0 )
                {
                    throw std::runtime_error(
                                "Error when making 1st order relativistic light time correction, could not find body " +
                                perturbingBodies.at( i ) );
                }
                else
                {
                    // Set state function.
                    perturbingBodyStateFunctions.push_back(
                                std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, double >,
                                                                         bodies.at( perturbingBodies[ i ] ), std::placeholders::_1 ) );

                    // Set gravitational parameter function.
                    perturbingBodyGravitationalParameterFunctions.push_back(
                                std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                             bodies.at( perturbingBodies[ i ] )->
                                             getGravityFieldModel( ) ) );
                }
            }

            // Create light-time correction function
            lightTimeCorrection = std::make_shared< FirstOrderLightTimeCorrectionCalculator >(
                        perturbingBodyStateFunctions, perturbingBodyGravitationalParameterFunctions, perturbingBodies,
                        transmitter.bodyName_, receiver.bodyName_,
                        std::bind( &relativity::PPNParameterSet::getParameterGamma, relativity::ppnParameterSet ) );

        }
        else
        {
            throw std::runtime_error(
                        "Error, correction settings type (1st order relativistic) does not coincide with data type." );
        }

        break;
    }
    case tabulated_tropospheric:
    {
        std::shared_ptr< TabulatedTroposphericCorrectionSettings > troposphericCorrectionSettings =
            std::dynamic_pointer_cast< TabulatedTroposphericCorrectionSettings >( correctionSettings );
        if ( correctionSettings == nullptr )
        {
            throw std::runtime_error(
                    "Error when creating tabulated tropospheric correction: incompatible settings type." );
        }

        if ( !isRadiometricObservableType( observableType ) )
        {
            throw std::runtime_error(
                    "Error when creating tabulated tropospheric correction: selected observable type is not radiometric." );
        }

        // If one of the link ends is the body with the atmosphere then create the tropospheric correction
        if ( transmitter.bodyName_ != receiver.bodyName_ && (
                transmitter.bodyName_ == troposphericCorrectionSettings->getBodyWithAtmosphere( ) ||
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
                                                        bodies, transmitter, receiver, isUplinkCorrection );

            ObservableType baseObservableType = getBaseObservableType( observableType );

            std::shared_ptr< TabulatedMediaReferenceCorrectionManager > dryCorrectionCalculator, wetCorrectionCalculator;
            std::pair< std::string, std::string > stationSpacecraftPair = std::make_pair( groundStation.stationName_, "" );

            if ( troposphericCorrectionSettings->getTroposphericDryCorrection( ).count( stationSpacecraftPair ) &&
                troposphericCorrectionSettings->getTroposphericDryCorrection( ).at( stationSpacecraftPair ).count( baseObservableType ) &&
                troposphericCorrectionSettings->getTroposphericWetCorrection( ).count( stationSpacecraftPair ) &&
                troposphericCorrectionSettings->getTroposphericWetCorrection( ).at( stationSpacecraftPair ).count( baseObservableType ) )
            {
                lightTimeCorrection = std::make_shared< TabulatedTroposphericCorrection >(
                    troposphericCorrectionSettings->getTroposphericDryCorrection( ).at(
                            stationSpacecraftPair ).at( baseObservableType ),
                    troposphericCorrectionSettings->getTroposphericWetCorrection( ).at(
                            stationSpacecraftPair ).at( baseObservableType ),
                    troposphericElevationMapping,
                    isUplinkCorrection );
            }
            else
            {
                throw std::runtime_error(
                        "Error when creating tabulated tropospheric corrections for " + groundStation.stationName_ +
                        "ground station and " + getObservableName( observableType ) + " observable: tabulated data not available. " );
            }
        }
        // Set correction to nullptr if correction isn't valid for selected link ends
        else
        {
            lightTimeCorrection = nullptr;
        }

        break;
    }
    case saastamoinen_tropospheric:
    {
        std::shared_ptr< SaastamoinenTroposphericCorrectionSettings > troposphericCorrectionSettings =
            std::dynamic_pointer_cast< SaastamoinenTroposphericCorrectionSettings >( correctionSettings );
        if ( correctionSettings == nullptr )
        {
            throw std::runtime_error(
                    "Error when creating Saastamoinen tropospheric correction: incompatible settings type." );
        }

        if ( !isRadiometricObservableType( observableType ) )
        {
            throw std::runtime_error(
                    "Error when creating tabulated tropospheric correction: selected observable type is not radiometric." );
        }

        // If one of the link ends is the body with the atmosphere then create the tropospheric correction
        if ( transmitter.bodyName_ != receiver.bodyName_ && (
                transmitter.bodyName_ == troposphericCorrectionSettings->getBodyWithAtmosphere( ) ||
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
                                                        bodies, transmitter, receiver, isUplinkCorrection );

            // Using nominal geodetic position, i.e. ignoring station motion
            std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction = std::bind(
                    &ground_stations::GroundStationState::getNominalGeodeticPosition,
                    bodies.getBody( groundStation.bodyName_ )->getGroundStation(
                            groundStation.stationName_ )->getNominalStationState( ) );

            std::function< double ( const double ) > waterVaporPartialPressureFunction;
            if ( troposphericCorrectionSettings->getWaterVaporPartialPressureModelType( ) == tabulated )
            {
                waterVaporPartialPressureFunction = bodies.getBody( groundStation.bodyName_ )->getGroundStation(
                        groundStation.stationName_ )->getWaterVaporPartialPressureFunction( );
            }
            else if ( troposphericCorrectionSettings->getWaterVaporPartialPressureModelType( ) == bean_and_dutton )
            {
                waterVaporPartialPressureFunction = getBeanAndDuttonWaterVaporPartialPressureFunction(
                        bodies.getBody( groundStation.bodyName_ )->getGroundStation(
                            groundStation.stationName_ )->getRelativeHumidityFunction( ),
                        bodies.getBody( groundStation.bodyName_ )->getGroundStation(
                            groundStation.stationName_ )->getTemperatureFunction( ) );
            }
            else
            {
                throw std::runtime_error( "Error when creating Saastamoinen tropospheric correction: water vapor partial "
                                          "pressure model type not recognized." );
            }

            lightTimeCorrection = std::make_shared< SaastamoinenTroposphericCorrection >(
                    groundStationGeodeticPositionFunction,
                    bodies.getBody( groundStation.bodyName_ )->getGroundStation( groundStation.stationName_ )->getPressureFunction( ),
                    bodies.getBody( groundStation.bodyName_ )->getGroundStation( groundStation.stationName_ )->getTemperatureFunction( ),
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
    case tabulated_ionospheric:
    {
        std::shared_ptr< TabulatedIonosphericCorrectionSettings > ionosphericCorrectionSettings =
            std::dynamic_pointer_cast< TabulatedIonosphericCorrectionSettings >( correctionSettings );
        if ( correctionSettings == nullptr )
        {
            throw std::runtime_error(
                    "Error when creating tabulated ionospheric correction: incompatible settings type." );
        }

        // If one of the link ends is the body with the atmosphere then create the tropospheric correction
        if ( transmitter.bodyName_ != receiver.bodyName_ && (
                transmitter.bodyName_ == ionosphericCorrectionSettings->getBodyWithAtmosphere( ) ||
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
            std::pair< std::string, std::string > stationSpacecraftPair = std::make_pair(
                    groundStation.stationName_, spacecraft.bodyName_ );

            if ( ionosphericCorrectionSettings->getReferenceRangeCorrection( ).count( stationSpacecraftPair ) &&
                ionosphericCorrectionSettings->getReferenceRangeCorrection( ).at( stationSpacecraftPair ).count( baseObservableType ) )
            {
                lightTimeCorrection = std::make_shared< TabulatedIonosphericCorrection >(
                    ionosphericCorrectionSettings->getReferenceRangeCorrection( ).at( stationSpacecraftPair ).at(
                            baseObservableType ),
                    createLinkFrequencyFunction( bodies, linkEnds, transmittingLinkEndType, receivingLinkEndType ),
                    baseObservableType,
                    isUplinkCorrection,
                    ionosphericCorrectionSettings->getReferenceFrequency( ) );
            }
            else
            {
                throw std::runtime_error(
                        "Error when creating tabulated tropospheric corrections for " + groundStation.stationName_ +
                        "ground station and " + getObservableName( observableType ) + " observable: tabulated data not available. " );
            }

        }
        // Set correction to nullptr if correction isn't valid for selected link ends
        else
        {
            lightTimeCorrection = nullptr;
        }

        break;
    }
    default:
    {
        std::string errorMessage = "Error, light time correction type " +
                std::to_string( correctionSettings->getCorrectionType( ) ) + " not recognized.";
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
    if ( isUplinkCorrection )
    {
        groundStation = transmitter;
    }
    else
    {
        groundStation = receiver;
    }

    switch( troposphericMappingModelType )
    {
        case simplified_chao:
        {
            std::function< double ( Eigen::Vector3d, double ) > elevationFunction = std::bind(
                    &ground_stations::PointingAnglesCalculator::calculateElevationAngle,
                    bodies.getBody( groundStation.bodyName_ )->getGroundStation(
                            groundStation.stationName_ )->getPointingAnglesCalculator( ),
                            std::placeholders::_1, std::placeholders::_2 );

            troposphericMappingModel = std::make_shared< SimplifiedChaoTroposphericMapping >(
                    elevationFunction, isUplinkCorrection );

            break;
        }
        case niell:
        {
            std::function< double ( Eigen::Vector3d, double ) > elevationFunction = std::bind(
                    &ground_stations::PointingAnglesCalculator::calculateElevationAngle,
                    bodies.getBody( groundStation.bodyName_ )->getGroundStation(
                            groundStation.stationName_ )->getPointingAnglesCalculator( ),
                            std::placeholders::_1, std::placeholders::_2 );

            // Using nominal geodetic position, i.e. ignoring station motion
            std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction = std::bind(
                    &ground_stations::GroundStationState::getNominalGeodeticPosition,
                    bodies.getBody( groundStation.bodyName_ )->getGroundStation(
                            groundStation.stationName_ )->getNominalStationState( ) );

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

std::function< double ( std::vector< FrequencyBands >, double ) > createLinkFrequencyFunction(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds& linkEnds,
        const LinkEndType& transmittingLinkEndType,
        const LinkEndType& receivingLinkEndType )
{
    std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittedFrequencyCalculator = bodies.getBody(
            linkEnds.at( transmitter ).bodyName_ )->getGroundStation( linkEnds.at( transmitter ).stationName_
                    )->getTransmittingFrequencyCalculator( );

    std::vector< std::function< double ( FrequencyBands, FrequencyBands ) > > turnaroundRatioFunctions;

    for ( auto retransmitterLinkEndsIt = linkEnds.begin( )++; retransmitterLinkEndsIt->first != receivingLinkEndType;
            ++retransmitterLinkEndsIt )
    {
        turnaroundRatioFunctions.push_back( bodies.getBody( retransmitterLinkEndsIt->second.bodyName_
            )->getVehicleSystems( )->getTransponderTurnaroundRatio( ) );
    }

    std::function< double ( std::vector< FrequencyBands >, double ) > linkFrequencyFunction = [=] (
            const std::vector< FrequencyBands >& frequencyBands, const double time )
    {
        double frequency = transmittedFrequencyCalculator->getTemplatedCurrentFrequency< double, double >( time );

        for ( unsigned int i = 0; i < turnaroundRatioFunctions.size( ); ++i )
        {
            frequency *= turnaroundRatioFunctions.at( i )( frequencyBands.at( i ), frequencyBands.at( i + 1 ) );
        }

        return frequency;
    };

    return linkFrequencyFunction;
}

} // namespace observation_models

} // namespace tudat
