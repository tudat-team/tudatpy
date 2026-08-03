/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/observationViabilityCalculator.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace simulation_setup
{

namespace
{

struct LightTimeCorrectionSourceIndex {
    unsigned int calculatorIndex_;

    int componentIndex_;
};

//! Resolve, for the given `LightTimeCorrectionComponentsDependentVariableSettings`, which entries
//! of the LightTimeCalculator per-correction caches should populate the dependent-variable vector.
//! If the settings object carries no filter, the order and size matches the calculator order and
//! each calculator's registered-correction order; otherwise all corrections of each requested type
//! are selected, grouped in filter order (missing types produce a clear error).
std::vector< LightTimeCorrectionSourceIndex > resolveLightTimeCorrectionSourceIndices(
        const std::shared_ptr< LightTimeCorrectionComponentsDependentVariableSettings >& lightTimeSettings,
        const std::vector< std::shared_ptr< observation_models::LightTimeCalculatorBase > >& lightTimeCalculators )
{
    const auto& filter = lightTimeSettings->correctionTypeFilter_;

    std::vector< LightTimeCorrectionSourceIndex > sourceIndices;
    if( filter.empty( ) )
    {
        for( unsigned int calculatorIndex = 0; calculatorIndex < lightTimeCalculators.size( ); calculatorIndex++ )
        {
            if( lightTimeCalculators.at( calculatorIndex ) == nullptr )
            {
                throw std::runtime_error(
                        "Error when adding light_time_correction_components dependent variable: null LightTimeCalculator "
                        "registered on this leg." );
            }
            const auto registeredCorrections = lightTimeCalculators.at( calculatorIndex )->getLightTimeCorrectionList( );
            for( unsigned int i = 0; i < registeredCorrections.size( ); i++ )
            {
                sourceIndices.push_back( { calculatorIndex, static_cast< int >( i ) } );
            }
        }
        return sourceIndices;
    }

    std::map< observation_models::LightTimeCorrectionType, std::vector< LightTimeCorrectionSourceIndex > > indicesByType;
    for( unsigned int calculatorIndex = 0; calculatorIndex < lightTimeCalculators.size( ); calculatorIndex++ )
    {
        if( lightTimeCalculators.at( calculatorIndex ) == nullptr )
        {
            throw std::runtime_error(
                    "Error when adding light_time_correction_components dependent variable: null LightTimeCalculator "
                    "registered on this leg." );
        }
        const auto registeredCorrections = lightTimeCalculators.at( calculatorIndex )->getLightTimeCorrectionList( );
        for( unsigned int i = 0; i < registeredCorrections.size( ); i++ )
        {
            indicesByType[ registeredCorrections[ i ]->getLightTimeCorrectionType( ) ].push_back(
                    { calculatorIndex, static_cast< int >( i ) } );
        }
    }

    std::map< observation_models::LightTimeCorrectionType, bool > processedRequestedTypes;
    for( const observation_models::LightTimeCorrectionType requestedType : filter )
    {
        if( processedRequestedTypes[ requestedType ] )
        {
            continue;
        }

        auto matchedIndices = indicesByType.find( requestedType );
        if( matchedIndices == indicesByType.end( ) || matchedIndices->second.empty( ) )
        {
            throw std::runtime_error( "Error when adding light_time_correction_components dependent variable: requested correction type '" +
                                      observation_models::getLightTimeCorrectionName( requestedType ) +
                                      "' is not registered on this leg." );
        }

        sourceIndices.insert( sourceIndices.end( ), matchedIndices->second.begin( ), matchedIndices->second.end( ) );
        processedRequestedTypes[ requestedType ] = true;
    }
    return sourceIndices;
}

}  // namespace

void checkObservationDependentVariableEnvironment(
        const SystemOfBodies& bodies,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings )
{
    if( isObservationDependentVariableGroundStationProperty( variableSettings->variableType_ ) && variableSettings->isLinkEndDefined_ )
    {
        std::string bodyName = variableSettings->linkEndId_.bodyName_;
        std::string stationName = variableSettings->linkEndId_.getReferencePointName( );

        if( bodies.count( bodyName ) == 0 )
        {
            throw std::runtime_error( "Error in observation dependent variables when creating function for " +
                                      getObservationDependentVariableId( variableSettings ) + ", could not find body " + bodyName );
        }

        if( bodies.at( bodyName )->getGroundStationMap( ).count( stationName ) == 0 )
        {
            throw std::runtime_error( "Error in observation dependent variables when creating function for " +
                                      getObservationDependentVariableId( variableSettings ) + ", could not find station " + stationName +
                                      " on body " + bodyName );
        }
    }
}

std::pair< int, int > getLinkEndStateTimeIndices( const observation_models::ObservableType observableType,
                                                  const observation_models::LinkDefinition linkEnds,
                                                  const observation_models::LinkEndId linkEndId,
                                                  const observation_models::LinkEndType linkEndRole,
                                                  const observation_models::LinkEndType originatingLinkEndRole,
                                                  const IntegratedObservationPropertyHandling integratedObservableHandling )
{
    // Get link-end indices consistent with settings
    std::vector< std::pair< int, int > > currentStateTimeIndex =
            observation_models::getLinkStateAndTimeIndicesForLinkEnd( linkEnds.linkEnds_, observableType, linkEndId );

    // If no indices are found, throw exception
    if( currentStateTimeIndex.size( ) == 0 )
    {
        throw std::runtime_error( "Error in getting link and state ID for " + observation_models::getObservableName( observableType ) +
                                  " Could not find link end: (" + linkEndId.bodyName_ + ", " + linkEndId.getReferencePointName( ) + ")" );
    }

    // Filter list of indices by link end role
    if( linkEndRole != observation_models::unidentified_link_end )
    {
        std::vector< std::pair< int, int > > stateTimeIndex;
        std::vector< int > linkEndIndices =
                observation_models::getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndRole, linkEnds.size( ) );
        for( unsigned int i = 0; i < currentStateTimeIndex.size( ); i++ )
        {
            for( unsigned int j = 0; j < linkEndIndices.size( ); j++ )
            {
                // If link end index and transmission index are equal, retain current link
                if( currentStateTimeIndex.at( i ).first == linkEndIndices.at( j ) )
                {
                    stateTimeIndex.push_back( currentStateTimeIndex.at( i ) );
                }
            }
        }
        currentStateTimeIndex = stateTimeIndex;
    }

    // Filter list of indices by originating link end role
    if( originatingLinkEndRole != observation_models::unidentified_link_end )
    {
        std::vector< std::pair< int, int > > stateTimeIndex;
        std::vector< int > linkEndIndices =
                observation_models::getLinkEndIndicesForLinkEndTypeAtObservable( observableType, originatingLinkEndRole, linkEnds.size( ) );
        for( unsigned int i = 0; i < currentStateTimeIndex.size( ); i++ )
        {
            for( unsigned int j = 0; j < linkEndIndices.size( ); j++ )
            {
                // If originating link end index and transmission index are equal, retain current link
                if( currentStateTimeIndex.at( i ).second == linkEndIndices.at( j ) )
                {
                    stateTimeIndex.push_back( currentStateTimeIndex.at( i ) );
                }
            }
        }
        currentStateTimeIndex = stateTimeIndex;
    }

    if( currentStateTimeIndex.size( ) == 0 )
    {
        throw std::runtime_error(
                "Error when getting link end state and time indices for observation dependent variable calculation, "
                "no indices could be found for required relevant and originating link ends (must be inconsistent)." );
    }

    // Check if observation is integrated
    if( observation_models::isObservableOfIntegratedType( observableType ) )
    {
        std::vector< std::pair< int, int > > stateTimeIndex;
        if( currentStateTimeIndex.size( ) != 2 )
        {
            throw std::runtime_error( "Error when getting integrated observable state and time indices; 2 remaining indices required" );
        }
        else if( integratedObservableHandling == interval_undefined )
        {
            throw std::runtime_error(
                    "Error when getting integrated observable state and time indices; choice of start or end of arc undefined" );
        }
        else if( integratedObservableHandling == interval_start )
        {
            stateTimeIndex.push_back( currentStateTimeIndex.at( 0 ) );
        }
        else if( integratedObservableHandling == interval_end )
        {
            stateTimeIndex.push_back( currentStateTimeIndex.at( 1 ) );
        }
        currentStateTimeIndex = stateTimeIndex;
    }
    else
    {
        if( currentStateTimeIndex.size( ) != 1 )
        {
            throw std::runtime_error( "Error when getting observable state and time indices; 1 index required" );
        }
    }

    return currentStateTimeIndex.at( 0 );
}

ObservationDependentVariableFunction getStationObservationAngleFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition linkEnds )
{
    // Check if relevant properties of environment exist
    checkObservationDependentVariableEnvironment( bodies, variableSettings );

    // Retrieve link-end ID for station

    observation_models::LinkEndId linkEndIdToUse;
    if( variableSettings->isLinkEndDefined_ )
    {
        linkEndIdToUse = variableSettings->linkEndId_;
    }
    else
    {
        linkEndIdToUse = linkEnds.at( variableSettings->linkEndType_ );
    }

    std::string bodyName = linkEndIdToUse.bodyName_;
    std::string stationName = linkEndIdToUse.getReferencePointName( );

    if( bodies.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when getting station angle function, body " + bodyName + " does not exist " );
    }

    if( bodies.at( bodyName )->getGroundStationMap( ).count( stationName ) == 0 )
    {
        throw std::runtime_error( "Error when getting station angle function, station " + stationName + " does not exist on body " +
                                  bodyName );
    }

    std::pair< int, int > linkEndIndicesToUse = getLinkEndStateTimeIndices( observableType,
                                                                            linkEnds,
                                                                            linkEndIdToUse,
                                                                            variableSettings->linkEndType_,
                                                                            variableSettings->originatingLinkEndType_,
                                                                            variableSettings->integratedObservableHandling_ );

    // Retrieve pointing angles calculator for station, and setup output function
    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            bodies.at( bodyName )->getGroundStationMap( ).at( stationName )->getPointingAnglesCalculator( );
    ObservationDependentVariableFunction outputFunction;

    if( variableSettings->variableType_ == station_elevation_angle )
    {
        outputFunction =
                [ = ]( const std::vector< double >& linkEndTimes,
                       const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                       const Eigen::VectorXd& observationValue,
                       const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySimulationSettings ) {
                    return ( Eigen::VectorXd( 1 ) << ground_stations::calculateGroundStationElevationAngle(
                                     pointingAnglesCalculator, linkEndStates, linkEndTimes, linkEndIndicesToUse ) )
                            .finished( );
                };
    }
    else if( variableSettings->variableType_ == station_azimuth_angle )
    {
        outputFunction =
                [ = ]( const std::vector< double >& linkEndTimes,
                       const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                       const Eigen::VectorXd& observationValue,
                       const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySimulationSettings ) {
                    return ( Eigen::VectorXd( 1 ) << ground_stations::calculateGroundStationAzimuthAngle(
                                     pointingAnglesCalculator, linkEndStates, linkEndTimes, linkEndIndicesToUse ) )
                            .finished( );
                };
    }
    return outputFunction;
}

ObservationDependentVariableFunction getInterlinkObservationVariableFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< InterlinkObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition linkEnds )
{
    std::pair< int, int > linkEndIndicesToUse = getLinkEndStateTimeIndices( observableType,
                                                                            linkEnds,
                                                                            linkEnds.at( variableSettings->linkEndType_ ),
                                                                            variableSettings->linkEndType_,
                                                                            variableSettings->originatingLinkEndType_,
                                                                            variableSettings->integratedObservableHandling_ );

    ObservationDependentVariableFunction outputFunction;

    switch( variableSettings->variableType_ )
    {
        case target_range: {
            if( variableSettings->relativeBody_ != "" )
            {
                throw std::runtime_error(
                        "Error when parsing target range observation dependent variable, relative body must not be defined" );
            }
            outputFunction = [ = ]( const std::vector< double >& linkEndTimes,
                                    const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                                    const Eigen::VectorXd& observationValue,
                                    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >
                                            ancillarySimulationSettings ) {
                return ( Eigen::VectorXd( 1 ) << ( linkEndStates.at( linkEndIndicesToUse.first ).segment( 0, 3 ) -
                                                   linkEndStates.at( linkEndIndicesToUse.second ).segment( 0, 3 ) )
                                                         .norm( ) )
                        .finished( );
            };
            break;
        }
        case body_avoidance_angle_variable: {
            if( bodies.count( variableSettings->relativeBody_ ) == 0 )
            {
                throw std::runtime_error( "Error when parsing body avoidance observation dependent variable w.r.t. " +
                                          variableSettings->relativeBody_ + ", body is not defined" );
            }

            if( bodies.at( variableSettings->relativeBody_ )->getEphemeris( ) == nullptr )
            {
                throw std::runtime_error( "Error when parsing body avoidance observation dependent variable w.r.t. " +
                                          variableSettings->relativeBody_ + ", body has no ephemeris" );
            }

            outputFunction = [ = ]( const std::vector< double >& linkEndTimes,
                                    const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                                    const Eigen::VectorXd& observationValue,
                                    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >
                                            ancillarySimulationSettings ) {
                double cosineOfAvoidanceAngle = observation_models::computeCosineBodyAvoidanceAngle(
                        linkEndStates.at( linkEndIndicesToUse.first ).segment( 0, 3 ),
                        linkEndStates.at( linkEndIndicesToUse.second ).segment( 0, 3 ),
                        bodies.at( variableSettings->relativeBody_ )
                                ->getStateInBaseFrameFromEphemeris< double, double >(
                                        ( linkEndTimes.at( linkEndIndicesToUse.first ) + linkEndTimes.at( linkEndIndicesToUse.second ) ) /
                                        2.0 )
                                .segment( 0, 3 ) );
                if( std::fabs( cosineOfAvoidanceAngle ) >= 1.0 )
                {
                    return ( Eigen::Vector1d( ) << utilities::sgn( cosineOfAvoidanceAngle ) * mathematical_constants::PI / 2.0 )
                            .finished( );
                }
                else
                {
                    return ( Eigen::Vector1d( ) << std::acos( cosineOfAvoidanceAngle ) ).finished( );
                }
            };
            break;
        }
        case link_body_center_distance: {
            if( bodies.count( variableSettings->relativeBody_ ) == 0 )
            {
                throw std::runtime_error( "Error when parsing link-body distance observation dependent variable w.r.t. " +
                                          variableSettings->relativeBody_ + ", body is not defined" );
            }

            if( bodies.at( variableSettings->relativeBody_ )->getEphemeris( ) == nullptr )
            {
                throw std::runtime_error( "Error when parsing link-body distance observation dependent variable w.r.t. " +
                                          variableSettings->relativeBody_ + ", body has no ephemeris" );
            }

            outputFunction = [ = ]( const std::vector< double >& linkEndTimes,
                                    const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                                    const Eigen::VectorXd& observationValue,
                                    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >
                                            ancillarySimulationSettings ) {
                double minimumDistance = observation_models::computeMinimumLinkDistanceToPoint(
                        linkEndStates.at( linkEndIndicesToUse.first ).segment( 0, 3 ),
                        linkEndStates.at( linkEndIndicesToUse.second ).segment( 0, 3 ),
                        bodies.at( variableSettings->relativeBody_ )
                                ->getStateInBaseFrameFromEphemeris< double, double >(
                                        ( linkEndTimes.at( linkEndIndicesToUse.first ) + linkEndTimes.at( linkEndIndicesToUse.second ) ) /
                                        2.0 )
                                .segment( 0, 3 ) );
                return ( Eigen::Vector1d( ) << minimumDistance ).finished( );
            };
            break;
        }
        case link_limb_distance: {
            if( bodies.count( variableSettings->relativeBody_ ) == 0 )
            {
                throw std::runtime_error( "Error when parsing link-body distance observation dependent variable w.r.t. " +
                                          variableSettings->relativeBody_ + ", body is not defined" );
            }

            if( bodies.at( variableSettings->relativeBody_ )->getEphemeris( ) == nullptr )
            {
                throw std::runtime_error( "Error when parsing link-body distance observation dependent variable w.r.t. " +
                                          variableSettings->relativeBody_ + ", body has no ephemeris" );
            }

            if( bodies.at( variableSettings->relativeBody_ )->getShapeModel( ) == nullptr )
            {
                throw std::runtime_error( "Error when parsing body link limb distance observation dependent variable w.r.t. " +
                                          variableSettings->relativeBody_ + ", body has no shape model" );
            }
            auto shapeModel = bodies.at( variableSettings->relativeBody_ )->getShapeModel( );
            outputFunction = [ = ]( const std::vector< double >& linkEndTimes,
                                    const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                                    const Eigen::VectorXd& observationValue,
                                    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >
                                            ancillarySimulationSettings ) {
                double minimumDistance = observation_models::computeMinimumLinkDistanceToPoint(
                        linkEndStates.at( linkEndIndicesToUse.first ).segment( 0, 3 ),
                        linkEndStates.at( linkEndIndicesToUse.second ).segment( 0, 3 ),
                        bodies.at( variableSettings->relativeBody_ )
                                ->getStateInBaseFrameFromEphemeris< double, double >(
                                        ( linkEndTimes.at( linkEndIndicesToUse.first ) + linkEndTimes.at( linkEndIndicesToUse.second ) ) /
                                        2.0 )
                                .segment( 0, 3 ) );
                return ( Eigen::Vector1d( ) << minimumDistance - shapeModel->getAverageRadius( ) ).finished( );
            };
            break;
        }
        case link_angle_with_orbital_plane: {
            if( bodies.count( variableSettings->relativeBody_ ) == 0 )
            {
                throw std::runtime_error( "Error when parsing link-orbital plane angle observation dependent variable w.r.t. " +
                                          variableSettings->relativeBody_ + ", body is not defined" );
            }

            if( bodies.at( variableSettings->relativeBody_ )->getEphemeris( ) == nullptr )
            {
                throw std::runtime_error( "Error when parsing link-orbital plane angle observation dependent variable w.r.t. " +
                                          variableSettings->relativeBody_ + ", body has no ephemeris" );
            }

            outputFunction = [ = ]( const std::vector< double >& linkEndTimes,
                                    const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                                    const Eigen::VectorXd& observationValue,
                                    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >
                                            ancillarySimulationSettings ) {
                Eigen::Vector3d vectorToTarget = linkEndStates.at( linkEndIndicesToUse.second ).segment( 0, 3 ) -
                        linkEndStates.at( linkEndIndicesToUse.first ).segment( 0, 3 );
                Eigen::Vector6d targetStateWrtCentralBody = linkEndStates.at( linkEndIndicesToUse.second ) -
                        bodies.at( variableSettings->relativeBody_ )
                                ->getStateInBaseFrameFromEphemeris< double, double >( linkEndTimes.at( linkEndIndicesToUse.second ) );
                Eigen::Vector3d orbitNormal =
                        targetStateWrtCentralBody.segment< 3 >( 0 ).cross( targetStateWrtCentralBody.segment< 3 >( 3 ) );
                return ( Eigen::Vector1d( ) << linear_algebra::computeAngleBetweenVectors( orbitNormal, vectorToTarget ) -
                                 mathematical_constants::PI / 2.0 )
                        .finished( );
            };
            break;
        }
        default:
            throw std::runtime_error( "Error when parsing interlink observation dependent variable, did not recognize variable" +
                                      getObservationDependentVariableId( variableSettings ) );
    }
    return outputFunction;
}

ObservationDependentVariableFunction getObservationVectorDependentVariableFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition linkEnds )
{
    ObservationDependentVariableFunction outputFunction;
    switch( variableSettings->variableType_ )
    {
        case station_elevation_angle: {
            std::shared_ptr< StationAngleObservationDependentVariableSettings > angleSettings =
                    std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >( variableSettings );
            if( angleSettings == nullptr )
            {
                throw std::runtime_error( "Error in observation dependent variables, incorrect type found for station_elevation_angle" );
            }
            outputFunction = getStationObservationAngleFunction( bodies, angleSettings, observableType, linkEnds );
            break;
        }
        case station_azimuth_angle: {
            std::shared_ptr< StationAngleObservationDependentVariableSettings > angleSettings =
                    std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >( variableSettings );
            if( angleSettings == nullptr )
            {
                throw std::runtime_error( "Error in observation dependent variables, incorrect type found for station_azimuth_angle" );
            }

            outputFunction = getStationObservationAngleFunction( bodies, angleSettings, observableType, linkEnds );
            break;
        }
        case target_range: {
            std::shared_ptr< InterlinkObservationDependentVariableSettings > linkSettings =
                    std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >( variableSettings );
            if( linkSettings == nullptr )
            {
                throw std::runtime_error( "Error in observation dependent variables, incorrect type found for target_range" );
            }

            outputFunction = getInterlinkObservationVariableFunction( bodies, linkSettings, observableType, linkEnds );
            break;
        }
        case body_avoidance_angle_variable: {
            std::shared_ptr< InterlinkObservationDependentVariableSettings > linkSettings =
                    std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >( variableSettings );
            if( linkSettings == nullptr )
            {
                throw std::runtime_error(
                        "Error in observation dependent variables, incorrect type found for body_avoidance_angle_variable" );
            }

            outputFunction = getInterlinkObservationVariableFunction( bodies, linkSettings, observableType, linkEnds );
            break;
        }
        case link_limb_distance: {
            std::shared_ptr< InterlinkObservationDependentVariableSettings > linkSettings =
                    std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >( variableSettings );
            if( linkSettings == nullptr )
            {
                throw std::runtime_error( "Error in observation dependent variables, incorrect type found for link_limb_distance" );
            }

            outputFunction = getInterlinkObservationVariableFunction( bodies, linkSettings, observableType, linkEnds );
            break;
        }
        case link_body_center_distance: {
            std::shared_ptr< InterlinkObservationDependentVariableSettings > linkSettings =
                    std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >( variableSettings );
            if( linkSettings == nullptr )
            {
                throw std::runtime_error( "Error in observation dependent variables, incorrect type found for link_body_center_distance" );
            }

            outputFunction = getInterlinkObservationVariableFunction( bodies, linkSettings, observableType, linkEnds );
            break;
        }
        case link_angle_with_orbital_plane: {
            std::shared_ptr< InterlinkObservationDependentVariableSettings > linkSettings =
                    std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >( variableSettings );
            if( linkSettings == nullptr )
            {
                throw std::runtime_error(
                        "Error in observation dependent variables, incorrect type found for link_angle_with_orbital_plane" );
            }

            outputFunction = getInterlinkObservationVariableFunction( bodies, linkSettings, observableType, linkEnds );
            break;
        }
        case integration_time_dependent_variable: {
            if( !observation_models::isObservableOfIntegratedType( observableType ) )
            {
                throw std::runtime_error( "Error in observation dependent variables, requested integration time for observable " +
                                          observation_models::getObservableName( observableType, linkEnds.size( ) ) );
            }

            outputFunction = [ = ]( const std::vector< double >& linkEndTimes,
                                    const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                                    const Eigen::VectorXd& observationValue,
                                    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >
                                            ancillarySimulationSettings ) {
                return ( Eigen::VectorXd( 1 ) << ancillarySimulationSettings->getAncillaryDoubleData(
                                 observation_models::doppler_integration_time, true ) )
                        .finished( );
            };
            break;
        }
        case retransmission_delays_dependent_variable: {
            if( !observation_models::observableCanHaveRetransmissionDelay( observableType ) )
            {
                throw std::runtime_error( "Error in observation dependent variables, requested retransmission time for observable " +
                                          observation_models::getObservableName( observableType, linkEnds.size( ) ) );
            }

            outputFunction = [ = ]( const std::vector< double >& linkEndTimes,
                                    const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                                    const Eigen::VectorXd& observationValue,
                                    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >
                                            ancillarySimulationSettings ) {
                int numberOfLinkEnds = linkEnds.size( );
                Eigen::VectorXd zeroDelay = Eigen::VectorXd::Zero( numberOfLinkEnds - 2 );

                std::vector< double > retransmissionDelays =
                        ancillarySimulationSettings->getAncillaryDoubleVectorData( observation_models::link_ends_delays, false );
                if( retransmissionDelays.size( ) > 0 )
                {
                    if( static_cast< int >( retransmissionDelays.size( ) ) != numberOfLinkEnds - 2 )
                    {
                        throw std::runtime_error( "Error in observation dependent variables, retransmission time size for observable " +
                                                  observation_models::getObservableName( observableType, linkEnds.size( ) ) +
                                                  " is of inconsistent size. Should be " + std::to_string( numberOfLinkEnds - 2 ) +
                                                  " but is " + std::to_string( retransmissionDelays.size( ) ) );
                    }

                    return utilities::convertStlVectorToEigenVector( retransmissionDelays );
                }
                else
                {
                    return zeroDelay;
                }
            };
            break;
        }
        case link_end_epochs_dependent_variable: {
            if( observableType != n_way_range && observableType != dsn_n_way_range )
            {
                throw std::runtime_error( "Error in observation dependent variables, requested link-end epochs for observable " +
                                          observation_models::getObservableName( observableType, linkEnds.size( ) ) +
                                          ". This setting is currently only available for n-way range observables." );
            }

            outputFunction = [ = ]( const std::vector< double >& linkEndTimes,
                                    const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                                    const Eigen::VectorXd& observationValue,
                                    const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >
                                            ancillarySimulationSettings ) {
                const int expectedNumberOfLinkEndTimes = 2 * ( static_cast< int >( linkEnds.size( ) ) - 1 );
                if( static_cast< int >( linkEndTimes.size( ) ) != expectedNumberOfLinkEndTimes )
                {
                    throw std::runtime_error( "Error in observation dependent variables, link-end epochs size for observable " +
                                              observation_models::getObservableName( observableType, linkEnds.size( ) ) +
                                              " is of inconsistent size. Should be " + std::to_string( expectedNumberOfLinkEndTimes ) +
                                              " but is " + std::to_string( linkEndTimes.size( ) ) + "." );
                }
                return utilities::convertStlVectorToEigenVector( linkEndTimes );
            };
            break;
        }
        case light_time_correction_components: {
            throw std::runtime_error(
                    "Error in observation dependent variables: light_time_correction_components is handled by "
                    "ObservationDependentVariableCalculator::addDependentVariable, not by the generic vector factory." );
            break;
        }
        default:
            throw std::runtime_error( "Error when parsing vector observation dependent variable, did not recognize variable" +
                                      getObservationDependentVariableId( variableSettings ) );
    }
    return outputFunction;
}

std::map< std::pair< int, int >, std::shared_ptr< ObservationDependentVariableSettings > >
ObservationDependentVariableBookkeeping::getSettingsIndicesAndSizes( ) const
{
    std::map< std::pair< int, int >, std::shared_ptr< ObservationDependentVariableSettings > > settingsStartIndices;
    for( unsigned int i = 0; i < dependentVariableStartIndices_.size( ); i++ )
    {
        settingsStartIndices[ std::make_pair( dependentVariableStartIndices_[ i ], dependentVariableSizes_[ i ] ) ] = settingsList_[ i ];
    }
    return settingsStartIndices;
}

std::pair< int, int > ObservationDependentVariableBookkeeping::addDependentVariable(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const int sizeOverride )
{
    // Check if the requested dependent variable can be used for given link
    if( !doesObservationDependentVariableExistForGivenLink( observableType_, linkEnds_.linkEnds_, variableSettings ) )
    {
        return std::make_pair( 0, 0 );
    }

    // Defer `light_time_correction_components` whose size cannot be resolved robustly without the
    // observation model's LightTimeCalculator (i.e. no caller-supplied size override). Even with
    // a non-empty type filter, the final size depends on how many registered corrections match
    // each requested type on the selected leg. The setting is kept in `deferredSettings_` and
    // turned into a real entry by `ObservationDependentVariableCalculator` once the leg map is
    // known.
    if( sizeOverride < 0 && variableSettings->variableType_ == light_time_correction_components )
    {
        auto lightTimeSettings = std::dynamic_pointer_cast< LightTimeCorrectionComponentsDependentVariableSettings >( variableSettings );
        if( lightTimeSettings == nullptr )
        {
            throw std::runtime_error(
                    "Error when adding light_time_correction_components dependent variable to bookkeeping: "
                    "settings object must be a LightTimeCorrectionComponentsDependentVariableSettings." );
        }
        addDeferredSetting( variableSettings );
        return std::make_pair( 0, 0 );
    }

    // Retrieve the current index in list of dependent variables and size of new parameter
    int currentIndex = totalDependentVariableSize_;
    int parameterSize = ( sizeOverride >= 0 ) ? sizeOverride : getObservationDependentVariableSize( variableSettings, linkEnds_.linkEnds_ );

    dependentVariableStartIndices_.push_back( totalDependentVariableSize_ );
    dependentVariableSizes_.push_back( parameterSize );
    settingsList_.push_back( variableSettings );
    totalDependentVariableSize_ += parameterSize;

    return std::make_pair( currentIndex, parameterSize );
}

void ObservationDependentVariableBookkeeping::addDependentVariables(
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList )
{
    // Parse all settings to be added and check if they are already included
    for( auto settingsToAdd : settingsList )
    {
        // Check if settings already exist for given observation set
        bool settingsDetected = false;
        for( auto existingSettings : settingsList_ )
        {
            if( existingSettings->areSettingsCompatible( settingsToAdd ) )
            {
                settingsDetected = true;
            }
        }

        // Add required dependent variable if not yet included
        if( !settingsDetected )
        {
            addDependentVariable( settingsToAdd );
        }
    }
}

std::pair< int, int > ObservationDependentVariableBookkeeping::getDependentVariableIndices(
        const std::shared_ptr< ObservationDependentVariableSettings > dependentVariables )
{
    std::pair< int, int > startAndSizePair = std::make_pair( 0, 0 );
    for( unsigned int i = 0; i < settingsList_.size( ); i++ )
    {
        if( settingsList_.at( i )->areSettingsCompatible( dependentVariables ) )
        {
            if( startAndSizePair.second == 0 )
            {
                startAndSizePair = std::make_pair( dependentVariableStartIndices_.at( i ), dependentVariableSizes_.at( i ) );
            }
            else
            {
                throw std::runtime_error( "Error when finding observation dependent variable; multiple candidates found" );
            }
        }
    }
    return startAndSizePair;
}

namespace
{

//! Build the add-function for a `light_time_correction_components` dependent variable given the
//! resolved slot (currentIndex, parameterSize), LightTimeCalculator, and source-index map.
ObservationDependentVariableAddFunction makeLightTimeCorrectionComponentsAddFunction(
        const int currentIndex,
        const int parameterSize,
        const std::vector< std::shared_ptr< observation_models::LightTimeCalculatorBase > >& lightTimeCalculators,
        const std::vector< LightTimeCorrectionSourceIndex >& sourceIndices )
{
    if( parameterSize != static_cast< int >( sourceIndices.size( ) ) )
    {
        throw std::runtime_error(
                "Error when saving light_time_correction_components: resolved source count does not match "
                "the dependent-variable bookkeeping size." );
    }

    return [ currentIndex, parameterSize, lightTimeCalculators, sourceIndices ](
                   Eigen::VectorXd& dependentVariables,
                   const std::vector< double >& /*linkEndTimes*/,
                   const std::vector< Eigen::Matrix< double, 6, 1 > >& /*linkEndStates*/,
                   const Eigen::VectorXd& /*observable*/,
                   const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > /*ancillary*/ ) {
        for( int i = 0; i < parameterSize; i++ )
        {
            if( dependentVariables( currentIndex + i ) == dependentVariables( currentIndex + i ) )
            {
                throw std::runtime_error( "Error when saving observation dependent variables; overriding existing value" );
            }
        }
        for( int i = 0; i < parameterSize; i++ )
        {
            const LightTimeCorrectionSourceIndex srcIdx = sourceIndices[ static_cast< size_t >( i ) ];
            if( srcIdx.calculatorIndex_ >= lightTimeCalculators.size( ) || lightTimeCalculators.at( srcIdx.calculatorIndex_ ) == nullptr )
            {
                throw std::runtime_error( "Error when saving light_time_correction_components: cached calculator index out of range." );
            }
            const std::vector< double >& components =
                    lightTimeCalculators.at( srcIdx.calculatorIndex_ )->getCurrentLightTimeCorrectionComponents( );
            if( srcIdx.componentIndex_ < 0 || srcIdx.componentIndex_ >= static_cast< int >( components.size( ) ) )
            {
                throw std::runtime_error(
                        "Error when saving light_time_correction_components: cached component index out of range. "
                        "Did the observation model evaluate before the dependent variable was computed?" );
            }
            dependentVariables( currentIndex + i ) = components[ static_cast< size_t >( srcIdx.componentIndex_ ) ];
        }
    };
}

}  // namespace

void ObservationDependentVariableCalculator::registerLightTimeCorrectionComponents(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    auto lightTimeSettings = std::dynamic_pointer_cast< LightTimeCorrectionComponentsDependentVariableSettings >( variableSettings );
    if( lightTimeSettings == nullptr )
    {
        throw std::runtime_error(
                "Error when adding light_time_correction_components dependent variable: settings object must be a "
                "LightTimeCorrectionComponentsDependentVariableSettings." );
    }

    const auto legKey = std::make_pair( lightTimeSettings->originatingLinkEndType_, lightTimeSettings->linkEndType_ );
    auto calculatorIt = legLightTimeCalculators_.find( legKey );
    if( calculatorIt == legLightTimeCalculators_.end( ) || calculatorIt->second.empty( ) )
    {
        throw std::runtime_error(
                "Error when adding light_time_correction_components dependent variable: no LightTimeCalculator registered for leg (" +
                observation_models::getLinkEndTypeString( lightTimeSettings->originatingLinkEndType_ ) + " -> " +
                observation_models::getLinkEndTypeString( lightTimeSettings->linkEndType_ ) +
                "). The observation model for this observable may not expose a leg matching these link-end types." );
    }

    const std::vector< std::shared_ptr< observation_models::LightTimeCalculatorBase > > lightTimeCalculators = calculatorIt->second;
    const std::vector< LightTimeCorrectionSourceIndex > sourceIndices =
            resolveLightTimeCorrectionSourceIndices( lightTimeSettings, lightTimeCalculators );
    lightTimeSettings->resolvedSize_ = static_cast< int >( sourceIndices.size( ) );

    const std::pair< int, int > indices =
            dependentVariableBookkeeping_->addDependentVariable( variableSettings, lightTimeSettings->resolvedSize_ );
    if( indices.second == 0 )
    {
        // Leg does not apply to this observable (shouldn't normally happen for a correctly-typed leg).
        return;
    }

    dependentVariableAddFunctions_.push_back(
            makeLightTimeCorrectionComponentsAddFunction( indices.first, indices.second, lightTimeCalculators, sourceIndices ) );
}

void ObservationDependentVariableCalculator::addDependentVariable(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const SystemOfBodies& bodies )
{
    // `light_time_correction_components` needs the per-leg LightTimeCalculator map to compute its
    // size and lambda.
    if( variableSettings->variableType_ == light_time_correction_components )
    {
        if( !doesObservationDependentVariableExistForGivenLink( dependentVariableBookkeeping_->getObservableType( ),
                                                                dependentVariableBookkeeping_->getLinkEnds( ).linkEnds_,
                                                                variableSettings ) )
        {
            return;
        }
        if( legLightTimeCalculators_.empty( ) )
        {
            throw std::runtime_error(
                    "Error when adding light_time_correction_components dependent variable: no leg "
                    "LightTimeCalculator map is available for observable '" +
                    observation_models::getObservableName( dependentVariableBookkeeping_->getObservableType( ) ) +
                    "'. Construct ObservationDependentVariableCalculator with a populated leg map." );
        }
        registerLightTimeCorrectionComponents( variableSettings );
        return;
    }

    const std::pair< int, int > currentIndexAndSize = dependentVariableBookkeeping_->addDependentVariable( variableSettings );
    if( currentIndexAndSize.second > 0 )
    {
        addDependentVariableFunction( variableSettings, bodies, currentIndexAndSize.first, currentIndexAndSize.second );
    }
}

void ObservationDependentVariableCalculator::addDependentVariables(
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList,
        const SystemOfBodies& bodies )
{
    // Parse all settings to be added and check if they are already included
    for( auto settingsToAdd : settingsList )
    {
        // Check if settings already exist for given observation set
        bool settingsDetected = false;
        for( auto existingSettings : dependentVariableBookkeeping_->getDependentVariableSettings( ) )
        {
            if( existingSettings->areSettingsCompatible( settingsToAdd ) )
            {
                settingsDetected = true;
            }
        }

        // Add required dependent variable if not yet included
        if( !settingsDetected )
        {
            addDependentVariable( settingsToAdd, bodies );
        }
    }
}

Eigen::VectorXd ObservationDependentVariableCalculator::calculateDependentVariables(
        const std::vector< double >& linkEndTimes,
        const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
        const Eigen::VectorXd& observation,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > observationAncillarySimulationSettings )
{
    Eigen::VectorXd dependentVariables =
            Eigen::VectorXd::Constant( dependentVariableBookkeeping_->getTotalDependentVariableSize( ), TUDAT_NAN );

    for( unsigned int i = 0; i < dependentVariableAddFunctions_.size( ); i++ )
    {
        dependentVariableAddFunctions_.at( i )(
                dependentVariables, linkEndTimes, linkEndStates, observation, observationAncillarySimulationSettings );
    }
    return dependentVariables;
}

void ObservationDependentVariableCalculator::addDependentVariableFunction(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const SystemOfBodies& bodies,
        const int currentIndex,
        const int parameterSize )
{
    // `light_time_correction_components` is not produced by the generic vector factory; it is
    // backed by the LightTimeCalculator's per-correction cache. In the constructor path, the
    // bookkeeping already holds the slot — we only need to build the add-function lambda.
    if( variableSettings->variableType_ == light_time_correction_components )
    {
        auto lightTimeSettings = std::dynamic_pointer_cast< LightTimeCorrectionComponentsDependentVariableSettings >( variableSettings );
        if( lightTimeSettings == nullptr )
        {
            throw std::runtime_error(
                    "Error when building light_time_correction_components add-function: settings object must be a "
                    "LightTimeCorrectionComponentsDependentVariableSettings." );
        }
        const auto legKey = std::make_pair( lightTimeSettings->originatingLinkEndType_, lightTimeSettings->linkEndType_ );
        auto calculatorIt = legLightTimeCalculators_.find( legKey );
        if( calculatorIt == legLightTimeCalculators_.end( ) || calculatorIt->second.empty( ) )
        {
            throw std::runtime_error(
                    "Error when building light_time_correction_components add-function: no LightTimeCalculator "
                    "registered for leg (" +
                    observation_models::getLinkEndTypeString( lightTimeSettings->originatingLinkEndType_ ) + " -> " +
                    observation_models::getLinkEndTypeString( lightTimeSettings->linkEndType_ ) +
                    "). Construct the calculator with a populated leg map before rehydrating." );
        }
        const std::vector< std::shared_ptr< observation_models::LightTimeCalculatorBase > > lightTimeCalculators = calculatorIt->second;
        const std::vector< LightTimeCorrectionSourceIndex > sourceIndices =
                resolveLightTimeCorrectionSourceIndices( lightTimeSettings, lightTimeCalculators );
        dependentVariableAddFunctions_.push_back(
                makeLightTimeCorrectionComponentsAddFunction( currentIndex, parameterSize, lightTimeCalculators, sourceIndices ) );
        return;
    }

    // Create function to compute dependent variable
    ObservationDependentVariableFunction observationDependentVariableFunction =
            getObservationVectorDependentVariableFunction( bodies,
                                                           variableSettings,
                                                           dependentVariableBookkeeping_->getObservableType( ),
                                                           dependentVariableBookkeeping_->getLinkEnds( ).linkEnds_ );

    // Create function to compute dependent variable and add to existing list
    ObservationDependentVariableAddFunction dependentVariableAddFunction =
            [ = ]( Eigen::VectorXd& dependentVariables,
                   const std::vector< double >& linkEndTimes,
                   const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                   const Eigen::VectorXd& observable,
                   const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySimulationSettings ) {
                // Check if computation is not overriding existing values
                for( int i = 0; i < parameterSize; i++ )
                {
                    if( dependentVariables( currentIndex + i ) == dependentVariables( currentIndex + i ) )
                    {
                        throw std::runtime_error( "Error when saving observation dependent variables; overriding existing value" );
                    }
                }
                dependentVariables.segment( currentIndex, parameterSize ) =
                        observationDependentVariableFunction( linkEndTimes, linkEndStates, observable, ancillarySimulationSettings );
            };

    // Add new dependent variable function and settings to list
    dependentVariableAddFunctions_.push_back( dependentVariableAddFunction );
}

}  // namespace simulation_setup

}  // namespace tudat
