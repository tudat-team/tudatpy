/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace observation_models
{

//std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::VectorXd,
//std::pair< std::vector< double >, LinkEndType > > > > getTudatCompatibleObservationsAndTimes(
//        const std::vector< std::tuple< ObservableType, LinkEnds, Eigen::VectorXd,
//        std::vector< double >, LinkEndType > >& tudatpyObservationsAndTimes )
//{
//    std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::VectorXd,
//    std::pair< std::vector< double >, LinkEndType > > > > tudatCompatibleObservationsAndTimes ;

//    for( unsigned int i = 0; i < tudatpyObservationsAndTimes.size( ); i++ )
//    {
//        auto currentTuple = tudatpyObservationsAndTimes.at( i );
//        tudatCompatibleObservationsAndTimes[ std::get< 0 >( currentTuple ) ][ std::get< 1 >( currentTuple ) ] =
//                std::make_pair(  std::get< 2 >( currentTuple ), std::make_pair(
//                                     std::get< 3 >( currentTuple ), std::get< 4 >( currentTuple ) ) );
//    }
//    return tudatCompatibleObservationsAndTimes;
//}

bool doesLinkEndTypeDefineId( const ObservableType observableType )
{
    bool linkEndTypeDefinesId = true;
    switch( observableType )
    {
    case one_way_range:
        break;
    case angular_position:
        break;
    case position_observable:
        break;
    case one_way_doppler:
        break;
    case one_way_differenced_range:
        break;
    case n_way_range:
        linkEndTypeDefinesId = false;
        break;
    case two_way_doppler:
    case doppler_measured_frequency:
        linkEndTypeDefinesId = false;
        break;
    case euler_angle_313_observable:
        break;
    case velocity_observable:
        break;
    case relative_angular_position:
        break;
    case relative_position_observable:
        break;
    case n_way_differenced_range:
        linkEndTypeDefinesId = false;
        break;
    default:
        throw std::runtime_error( "Error when determining if link end type defines id; observable " +
                                  getObservableName( observableType ) + " not found" );
    }
    return linkEndTypeDefinesId;
}

bool isObservableTypeMultiLink( const ObservableType observableType )
{
    bool isIntegratedTypeisMultiLink = false;
    switch( observableType )
    {
    case one_way_range:
        break;
    case angular_position:
        break;
    case position_observable:
        break;
    case one_way_doppler:
        break;
    case one_way_differenced_range:
        break;
    case n_way_range:
        isIntegratedTypeisMultiLink = true;
        break;
    case two_way_doppler:
    case doppler_measured_frequency:
        isIntegratedTypeisMultiLink = true;
        break;
    case euler_angle_313_observable:
        break;
    case velocity_observable:
        break;
    case relative_angular_position:
        isIntegratedTypeisMultiLink = true;
        break;
    case relative_position_observable:
        break;
    case n_way_differenced_range:
        isIntegratedTypeisMultiLink = true;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type is multilink; observable " +
                                  getObservableName( observableType ) + " not found" );
    }
    return isIntegratedTypeisMultiLink;
}

bool isObservableOfIntegratedType( const ObservableType observableType )
{
    bool isIntegratedType = true;
    switch( observableType )
    {
    case one_way_range:
        isIntegratedType = false;
        break;
    case angular_position:
        isIntegratedType = false;
        break;
    case position_observable:
        isIntegratedType = false;
        break;
    case one_way_doppler:
        isIntegratedType = false;
        break;
    case one_way_differenced_range:
        isIntegratedType = true;
        break;
    case n_way_range:
        isIntegratedType = false;
        break;
    case two_way_doppler:
    case doppler_measured_frequency:
        isIntegratedType = false;
        break;
    case euler_angle_313_observable:
        isIntegratedType = false;
        break;
    case velocity_observable:
        isIntegratedType = false;
        break;
    case relative_angular_position:
        isIntegratedType = false;
        break;
    case relative_position_observable:
        isIntegratedType = false;
        break;
    case n_way_differenced_range:
        isIntegratedType = true;
        break;
    case dsn_one_way_averaged_doppler:
        isIntegratedType = true;
        break;
    case dsn_n_way_averaged_doppler:
        isIntegratedType = true;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type is integrated; observable " +
                                  getObservableName( observableType ) + " not found" );
    }
    return isIntegratedType;
}

bool requiresTransmittingStation( const ObservableType observableType )
{
    bool requiresTransmittingStation = true;
    switch( observableType )
    {
    case one_way_range:
    case angular_position:
    case position_observable:
    case one_way_doppler:
    case one_way_differenced_range:
    case n_way_range:
    case two_way_doppler:
    case euler_angle_313_observable:
    case velocity_observable:
    case relative_angular_position:
    case n_way_differenced_range:
    case dsn_one_way_averaged_doppler:
        requiresTransmittingStation = false;
        break;
    case dsn_n_way_averaged_doppler:
    case doppler_measured_frequency:
        requiresTransmittingStation = true;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type requires transmitting station: observable " +
                                  getObservableName( observableType ) + " not found." );
    }
    return requiresTransmittingStation;
}

bool requiresFirstReceivingStation( const ObservableType observableType )
{
    bool requiresFirstReceivingStation = true;
    switch( observableType )
    {
    case one_way_range:
    case angular_position:
    case position_observable:
    case one_way_doppler:
    case one_way_differenced_range:
    case n_way_range:
    case two_way_doppler:
    case doppler_measured_frequency:
    case euler_angle_313_observable:
    case velocity_observable:
    case relative_angular_position:
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
        requiresFirstReceivingStation = false;
        break;
    case dsn_one_way_averaged_doppler:
        requiresFirstReceivingStation = true;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type requires first receiving station: observable " +
                                  getObservableName( observableType ) + " not found." );
    }
    return requiresFirstReceivingStation;
}

bool requiresSecondReceivingStation( const ObservableType observableType )
{
    bool requiresSecondReceivingStation = true;
    switch( observableType )
    {
    case one_way_range:
    case angular_position:
    case position_observable:
    case one_way_doppler:
    case one_way_differenced_range:
    case n_way_range:
    case two_way_doppler:
    case doppler_measured_frequency:
    case euler_angle_313_observable:
    case velocity_observable:
    case relative_angular_position:
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
    case dsn_one_way_averaged_doppler:
        requiresSecondReceivingStation = false;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type requires second receiving station: observable " +
                                  getObservableName( observableType ) + " not found." );
    }
    return requiresSecondReceivingStation;
}

bool isRadiometricObservableType( const ObservableType observableType )
{
    bool isRadiometric;
    switch( observableType )
    {
    case one_way_range:
    case one_way_doppler:
    case one_way_differenced_range:
    case n_way_range:
    case two_way_doppler:
    case doppler_measured_frequency:
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
    case dsn_one_way_averaged_doppler:
        isRadiometric = true;
        break;
    case angular_position:
    case position_observable:
    case euler_angle_313_observable:
    case velocity_observable:
    case relative_angular_position:
        isRadiometric = false;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type is radiometric: observable " +
                                  getObservableName( observableType ) + " not found." );
    }
    return isRadiometric;
}

bool isPhaseVelocityBasedObservableType( const ObservableType observableType )
{
    bool isPhaseVelocityBased;
    switch( observableType )
    {
    case one_way_doppler:
    case one_way_differenced_range:
    case two_way_doppler:
    case doppler_measured_frequency:
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
    case dsn_one_way_averaged_doppler:
        isPhaseVelocityBased = true;
        break;
    case one_way_range:
    case n_way_range:
    case angular_position:
    case position_observable:
    case euler_angle_313_observable:
    case velocity_observable:
    case relative_angular_position:
        isPhaseVelocityBased = false;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type is radiometric: observable " +
                                  getObservableName( observableType ) + " not found." );
    }
    return isPhaseVelocityBased;
}

bool isGroupVelocityBasedObservableType( const ObservableType observableType )
{
    bool isGroupVelocityBased;
    switch( observableType )
    {
    case one_way_range:
    case n_way_range:
        isGroupVelocityBased = true;
        break;
    case two_way_doppler:
    case doppler_measured_frequency:
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
    case dsn_one_way_averaged_doppler:
    case one_way_doppler:
    case one_way_differenced_range:
    case angular_position:
    case position_observable:
    case euler_angle_313_observable:
    case velocity_observable:
    case relative_angular_position:
        isGroupVelocityBased = false;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type is radiometric: observable " +
                                  getObservableName( observableType ) + " not found." );
    }
    return isGroupVelocityBased;
}

bool observableCanHaveRetransmissionDelay( const ObservableType observableType )
{

    bool canHaveDelay = false;
    switch( observableType )
    {
    case one_way_range:
        break;
    case angular_position:
        break;
    case position_observable:
        break;
    case one_way_doppler:
        break;
    case one_way_differenced_range:
        break;
    case n_way_range:
        canHaveDelay = true;
        break;
    case two_way_doppler:
        canHaveDelay = true;
        break;
    case euler_angle_313_observable:
        break;
    case velocity_observable:
        break;
    case relative_angular_position:
        break;
    case relative_position_observable:
        break;
    case n_way_differenced_range:
        canHaveDelay = true;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type can have retransmission delay; observable " +
                                  getObservableName( observableType ) + " not found" );
    }
    return canHaveDelay;
}

bool linkEndIdDefinesSingleLink( const ObservableType observableType )
{
    if( isObservableOfIntegratedType( observableType ) || isObservableTypeMultiLink( observableType ) || !doesLinkEndTypeDefineId( observableType ) )
    {
        return false;
    }
    else
    {
        return true;
    }
}


//bool areObservableLinksContinuous( const ObservableType observableType )
//{
//    bool isTypeContinuous = true;
//    switch( observableType )
//    {
//    case one_way_range:
//        isTypeContinuous = true;
//        break;
//    case angular_position:
//        isTypeContinuous = true;
//        break;
//    case position_observable:
//        isTypeContinuous = true;
//        break;
//    case one_way_doppler:
//        isTypeContinuous = true;
//        break;
//    case one_way_differenced_range:
//        isTypeContinuous = true;
//        break;
//    case n_way_range:
//        isTypeContinuous = true;
//        break;
//    case two_way_doppler:
//        isTypeContinuous = true;
//        break;
//    case euler_angle_313_observable:
//        isTypeContinuous = true;
//        break;
//    case velocity_observable:
//        isTypeContinuous = true;
//        break;
//    case relative_angular_position:
//        isTypeContinuous = true;
//        break;
//    case n_way_differenced_range:
//        isIntegratedType = true;
//        break;
//    default:
//        throw std::runtime_error( "Error when determining if observable type is continuous; observable " +
//                                  getObservableName( observableType ) + " not found" );
//    }
//    return isTypeContinuous;
//}

std::string getNWayString( const int numberOfLinkEnds )
{
    std::string numberOfWays = "";
    switch( numberOfLinkEnds )
    {
    case 2:
        numberOfWays = "One";
        break;
    case 3:
        numberOfWays = "Two";
        break;
    case 4:
        numberOfWays = "Three";
        break;
    case 5:
        numberOfWays = "Four";
        break;
    case 6:
        numberOfWays = "Five";
        break;
    case 7:
        numberOfWays = "Six";
        break;
    default:
        numberOfWays = "N";
    }
    return numberOfWays;
}

//! Function to get the name (string) associated with a given observable type.
std::string getObservableName( const ObservableType observableType, const int numberOfLinkEnds )
{
    std::string observableName;
    switch( observableType )
    {
    case one_way_range:
        observableName = "OneWayRange";
        break;
    case angular_position:
        observableName = "AngularPosition";
        break;
    case position_observable:
        observableName = "CartesianPosition";
        break;
    case velocity_observable:
        observableName = "CartesianVelocity";
        break;
    case one_way_doppler:
        observableName = "OneWayDoppler";
        break;
    case one_way_differenced_range:
        observableName = "OneWayDifferencedRange";
        break;
    case two_way_doppler:
        observableName = "TwoWayDoppler";
        break;
    case n_way_range:
    {
        observableName = getNWayString( numberOfLinkEnds ) + "WayRange";
        break;
    }
    case euler_angle_313_observable:
        observableName = "EulerAngle313";
        break;
    case relative_angular_position:
        observableName = "RelativeAngularPosition";
        break;
    case relative_position_observable:
        observableName = "RelativeCartesianPosition";
        break;
    case n_way_differenced_range:
        observableName = getNWayString( numberOfLinkEnds ) + "WayDifferencedRange";
        break;
    case dsn_one_way_averaged_doppler:
        observableName = "DsnOneWayAveragedDoppler";
        break;
    case dsn_n_way_averaged_doppler:
        observableName = "Dsn" + getNWayString( numberOfLinkEnds ) + "WayAveragedDoppler";
        break;
    case doppler_measured_frequency:
        observableName = getNWayString( numberOfLinkEnds ) + "WayDopplerMeasuredFrequency";
        break;
    default:
        std::string errorMessage =
                "Error, could not find observable type " + std::to_string( observableType ) +
                " when getting name from type";
        std::cerr<<errorMessage<<std::endl;

//        throw std::runtime_error( errorMessage );
    }

    return observableName;
}

//! Function to get the observable type associated with the name (string) of observable.
ObservableType getObservableType( const std::string& observableName )
{
    ObservableType observableType;

    if( observableName == "OneWayRange" )
    {
        observableType = one_way_range;
    }
    else if( observableName == "AngularPosition" )
    {
        observableType = angular_position;
    }
    else if( observableName == "CartesianPosition" )
    {
        observableType = position_observable;
    }
    else if( observableName == "CartesianVelocity" )
    {
        observableType = velocity_observable;
    }
    else if( observableName ==  "OneWayDoppler" )
    {
        observableType = one_way_doppler;
    }
    else if( observableName ==  "TwoWayDoppler" )
    {
        observableType = two_way_doppler;
    }
    else if( observableName ==  "OneWayDifferencedRange" )
    {
        observableType = one_way_differenced_range;
    }
    else if( observableName == "EulerAngle313" )
    {
        observableType = euler_angle_313_observable;
    }
    else if( observableName == "RelativeAngularPosition" )
    {
        observableType = relative_angular_position;
    }
    else if ( observableName == "RelativeCartesianPosition" )
    {
        observableType = relative_position_observable;
    }
    else
    {
        std::string errorMessage =
                "Error, could not find observable name " + observableName +
                " when getting type from name";
        throw std::runtime_error( errorMessage );
    }

    return observableType;
}

ObservableType getUndifferencedObservableType( const ObservableType differencedObservableType )
{
    ObservableType undifferencedObservableType = undefined_observation_model;
    switch( differencedObservableType )
    {
    case one_way_differenced_range:
        undifferencedObservableType = one_way_range;
        break;
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
        undifferencedObservableType = n_way_range;
        break;
    case relative_angular_position:
        undifferencedObservableType = angular_position;
        break;
    default:
        throw std::runtime_error( "Error when getting undifferenced observable type for " + getObservableName(
                                      differencedObservableType ) + ", no such type exists" );

    }
    return undifferencedObservableType;
}

ObservableType getDifferencedObservableType( const ObservableType undifferencedObservableType )
{
    ObservableType differencedObservableType = undefined_observation_model;
    switch( undifferencedObservableType )
    {
    case one_way_range:
        differencedObservableType = one_way_differenced_range;
        break;
    case n_way_range:
        differencedObservableType = n_way_differenced_range;
        break;
    case angular_position:
        differencedObservableType = relative_angular_position;
        break;
    default:
        throw std::runtime_error( "Error when getting differenced observable type for " + getObservableName(
                                      undifferencedObservableType ) + ", no such type exists" );

    }
    return differencedObservableType;
}

ObservableType getUnconcatenatedObservableType( const ObservableType observableType )
{
    ObservableType unconcatenatedObservableType = undefined_observation_model;
    switch( observableType )
    {
    case n_way_differenced_range:
    case n_way_range:
    case dsn_n_way_averaged_doppler:
        unconcatenatedObservableType = one_way_range;
        break;
    default:
        throw std::runtime_error( "Error when getting unconcatenated observable type for " + getObservableName(
            observableType ) + ", no such type exists" );

    }
    return unconcatenatedObservableType;
}

ObservableType getBaseObservableType( const ObservableType observableType )
{
    ObservableType baseObservableType = undefined_observation_model;
    switch( observableType )
    {
    case one_way_range:
    case n_way_range:
        baseObservableType = one_way_range;
        break;
    case one_way_doppler:
    case one_way_differenced_range:
    case two_way_doppler:
    case n_way_differenced_range:
    case dsn_one_way_averaged_doppler:
    case dsn_n_way_averaged_doppler:
    case doppler_measured_frequency:
        baseObservableType = one_way_doppler;
        break;
    case angular_position:
    case relative_angular_position:
        baseObservableType = angular_position;
        break;
    default:
        throw std::runtime_error( "Error when getting base observable type for " + getObservableName( observableType ) +
            ", no such type exists" );

    }
    return baseObservableType;
}

std::pair< std::vector< int >, std::vector< int > > getUndifferencedTimeAndStateIndices(
        const ObservableType differencedObservableType,
        const int numberOfLinkEnds )
{
    std::vector< int > firstIndices;
    std::vector< int > secondIndices;

    switch( differencedObservableType )
    {
    case one_way_differenced_range:
        firstIndices = { 0, 1 };
        secondIndices = { 2, 3 };
        break;
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
    {
        int numberOfLinkEndTimesStates = 2 + ( numberOfLinkEnds - 2 ) * 2;
        for( int i = 0; i < numberOfLinkEndTimesStates; i++ )
        {
           firstIndices.push_back( i );
           secondIndices.push_back( i + numberOfLinkEndTimesStates );
        }
        break;
    }
    case relative_angular_position:
        firstIndices = { 0, 2 };
        secondIndices = { 1, 2 };
        break;
    default:
        throw std::runtime_error( "Error when getting undifferenced time and state entries for: " + getObservableName(
                                      differencedObservableType ) + ", no such type found" );
    }
    return std::make_pair( firstIndices, secondIndices );
}


std::pair< LinkEnds, LinkEnds > getUndifferencedLinkEnds( const ObservableType differencedObservableType, const LinkEnds& differencedLinkEnds )
{
    std::pair< LinkEnds, LinkEnds > linkEndsPair;
    switch( differencedObservableType )
    {
    case one_way_differenced_range:
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
        linkEndsPair = { differencedLinkEnds, differencedLinkEnds };
        break;
    case relative_angular_position:
    {
        LinkEnds firstLinkEnds;
        firstLinkEnds[ transmitter ] = differencedLinkEnds.at( transmitter );
        firstLinkEnds[ receiver ] = differencedLinkEnds.at( receiver );

        LinkEnds secondLinkEnds;
        secondLinkEnds[ transmitter ] = differencedLinkEnds.at( transmitter2 );
        secondLinkEnds[ receiver ] = differencedLinkEnds.at( receiver );
        linkEndsPair = { firstLinkEnds, secondLinkEnds };
        break;
    }
    default:
        throw std::runtime_error( "Error when getting undifferenced link ends for " + getObservableName(
            differencedObservableType ) + ", no such type exists" );

    }
    return linkEndsPair;
}


std::vector< LinkEnds > getUnconcatenatedLinkEnds( const ObservableType concatenatedObservableType, const LinkEnds& concatenatedLinkEnds )
{
    std::vector< LinkEnds >  linkEndsList;
    switch( concatenatedObservableType )
    {
    case n_way_differenced_range:
    case n_way_range:
    case dsn_n_way_averaged_doppler:
    {
        auto linkEndIterator = concatenatedLinkEnds.begin( );
        for( unsigned int i = 0; i < concatenatedLinkEnds.size( ) - 1; i++ )
        {
            LinkEnds currentOneWayLinkEnds;
            currentOneWayLinkEnds[ transmitter ] = linkEndIterator->second;
            linkEndIterator++;
            currentOneWayLinkEnds[ receiver ] = linkEndIterator->second;
            linkEndsList.push_back( currentOneWayLinkEnds );

        }
        break;
    }
    default:
        throw std::runtime_error( "Error when getting unconcatenated link ends for " + getObservableName(
            concatenatedObservableType ) + ", no such type exists" );

    }
    return linkEndsList;
}





//! Function to get the size of an observable of a given type.
int getObservableSize( const ObservableType observableType )
{
    int observableSize = -1;
    switch( observableType )
    {
    case one_way_range:
        observableSize = 1;
        break;
    case angular_position:
        observableSize = 2;
        break;
    case position_observable:
        observableSize = 3;
        break;
    case velocity_observable:
        observableSize = 3;
        break;
    case one_way_doppler:
        observableSize = 1;
        break;
    case two_way_doppler:
    case doppler_measured_frequency:
        observableSize = 1;
        break;
    case one_way_differenced_range:
        observableSize = 1;
        break;
    case n_way_range:
        observableSize = 1;
        break;
    case euler_angle_313_observable:
        observableSize = 3;
        break;
    case relative_angular_position:
        observableSize = 2;
        break;
    case n_way_differenced_range:
        observableSize = 1;
        break;
    case dsn_one_way_averaged_doppler:
        observableSize = 1;
        break;
    case dsn_n_way_averaged_doppler:
        observableSize = 1;
        break;
    case relative_position_observable:
        observableSize = 3;
        break;
    default:
       std::string errorMessage = "Error, did not recognize observable " + std::to_string( observableType )
               + ", when getting observable size";
       throw std::runtime_error( errorMessage );
    }
    return observableSize;
}

//! Function to get the indices in link end times/states for a given link end type and observable type
std::vector< int > getLinkEndIndicesForLinkEndTypeAtObservable(
        const ObservableType observableType, const LinkEndType linkEndType, const int numberOfLinkEnds )
{
    std::vector< int > linkEndIndices;

    switch( observableType )
    {
    case one_way_range:
        switch( linkEndType )
        {
        case transmitter:
            linkEndIndices.push_back( 0 );
            break;
        case receiver:
            linkEndIndices.push_back( 1 );
            break;
        default:
            std::string errorMessage =
                    "Error, could not find link end type index for link end " +
                    std::to_string( linkEndType ) + " of observable " +
                    std::to_string( observableType );
            throw std::runtime_error( errorMessage );
        }
        break;
    case one_way_doppler:
        switch( linkEndType )
        {
        case transmitter:
            linkEndIndices.push_back( 0 );
            break;
        case receiver:
            linkEndIndices.push_back( 1 );
            break;
        default:
            std::string errorMessage =
                    "Error, could not find link end type index for link end " +
                    std::to_string( linkEndType ) + " of observable " +
                    std::to_string( observableType );
            throw std::runtime_error( errorMessage );
        }
        break;
    case two_way_doppler:
    case doppler_measured_frequency:
        switch( linkEndType )
        {
        case transmitter:
            linkEndIndices.push_back( 0 );
            break;
        case reflector1:
            linkEndIndices.push_back( 1 );
            linkEndIndices.push_back( 2 );
            break;
        case receiver:
            linkEndIndices.push_back( 3 );
            break;
        default:
            std::string errorMessage =
                    "Error, could not find link end type index for link end " +
                    std::to_string( linkEndType ) + " of observable " +
                    std::to_string( observableType );
            throw std::runtime_error( errorMessage );
        }
        break;
    case one_way_differenced_range:
    case dsn_one_way_averaged_doppler:
        switch( linkEndType )
        {
        case transmitter:
            linkEndIndices.push_back( 0 );
            linkEndIndices.push_back( 2 );
            break;
        case receiver:
            linkEndIndices.push_back( 1 );
            linkEndIndices.push_back( 3 );
            break;
        default:
            std::string errorMessage =
                    "Error, could not find link end type index for link end " +
                    std::to_string( linkEndType ) + " of observable " +
                    std::to_string( observableType );
            throw std::runtime_error( errorMessage );
        }
        break;
    case angular_position:
        switch( linkEndType )
        {
        case transmitter:
            linkEndIndices.push_back( 0 );
            break;
        case receiver:
            linkEndIndices.push_back( 1 );
            break;
        default:
            std::string errorMessage =
                    "Error, could not find link end type index for link end " +
                    std::to_string( linkEndType ) + " of observable " +
                    std::to_string( observableType );
            throw std::runtime_error( errorMessage );
        }
        break;
    case position_observable:
        if( linkEndType == observed_body )
        {
            linkEndIndices.push_back( 0 );
        }
        else
        {
            std::string errorMessage =
                    "Error, could not find link end type index for link end " +
                    std::to_string( linkEndType ) + " of observable " +
                    std::to_string( observableType );
            throw std::runtime_error( errorMessage );
        }
        break;
    case euler_angle_313_observable:
        if( linkEndType == observed_body )
        {
            linkEndIndices.push_back( 0 );
        }
        else
        {
            std::string errorMessage =
                    "Error, could not find link end type index for link end " +
                    std::to_string( linkEndType ) + " of observable " +
                    std::to_string( observableType );
            throw std::runtime_error( errorMessage );
        }
        break;
    case velocity_observable:
        if( linkEndType == observed_body )
        {
            linkEndIndices.push_back( 0 );
        }
        else
        {
            std::string errorMessage =
                    "Error, could not find link end type index for link end " +
                    std::to_string( linkEndType ) + " of observable " +
                    std::to_string( observableType );
            throw std::runtime_error( errorMessage );
        }
        break;
    case n_way_range:
        if( numberOfLinkEnds < 2 )
        {
            throw std::runtime_error( "Error when getting n way differenced range link end indices, not enough link ends" );
        }
        if( linkEndType == transmitter )
        {
            linkEndIndices.push_back( 0 );
        }
        else if( linkEndType == receiver )
        {
            linkEndIndices.push_back( 2 * ( numberOfLinkEnds - 1 ) - 1 );
        }
        else
        {
            int linkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndType, numberOfLinkEnds );
            linkEndIndices.push_back( 2 * linkEndIndex - 1 );
            linkEndIndices.push_back( 2 * linkEndIndex );
        }
        break;
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
        if( numberOfLinkEnds < 2 )
        {
            throw std::runtime_error( "Error when getting " + std::to_string( observableType ) +
                " link end indices, not enough link ends" );
        }
        if( linkEndType == transmitter )
        {
            linkEndIndices.push_back( 0 );
            linkEndIndices.push_back( 2 * ( numberOfLinkEnds - 1 ) );

        }
        else if( linkEndType == receiver )
        {
            linkEndIndices.push_back( 2 * ( numberOfLinkEnds - 1 ) - 1 );
            linkEndIndices.push_back( 4 * ( numberOfLinkEnds - 1 ) - 1 );
        }
        else
        {
            int linkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndType, numberOfLinkEnds );
            linkEndIndices.push_back( 2 * linkEndIndex - 1 );
            linkEndIndices.push_back( 2 * linkEndIndex );
            linkEndIndices.push_back( 2 * ( ( numberOfLinkEnds - 1 ) + linkEndIndex ) - 1 );
            linkEndIndices.push_back( 2 * ( ( numberOfLinkEnds - 1 ) + linkEndIndex ) );
        }
        break;
    case relative_angular_position:
        switch( linkEndType )
        {
            case transmitter:
                linkEndIndices.push_back( 0 );
                break;
            case transmitter2:
                linkEndIndices.push_back( 1 );
                break;
            case receiver:
                linkEndIndices.push_back( 2 );
                break;
            default:
                std::string errorMessage =
                        "Error, could not find link end type index for link end " +
                        std::to_string( linkEndType ) + " of observable " +
                        std::to_string( observableType );
                throw std::runtime_error( errorMessage );
        }
        break;
    case relative_position_observable:
        switch( linkEndType )
        {
            case observed_body:
                linkEndIndices.push_back( 0 );
                break;
            case observer:
                linkEndIndices.push_back( 1 );
                break;
            default:
                std::string errorMessage = "Error, could not find link end type index for link end " + std::to_string( linkEndType ) + " of observable " +
                        std::to_string( observableType );
                throw std::runtime_error( errorMessage );
        }
        break;
    default:
        std::string errorMessage =
                "Error, could not find link end type index for link end types of observable " +
                std::to_string( observableType );
        throw std::runtime_error( errorMessage );
    }

    return linkEndIndices;
}


LinkEndType getDefaultReferenceLinkEndType(
        const ObservableType observableType )
{
    LinkEndType referenceLinkEndType;
    switch( observableType )
    {
    case one_way_range:
        referenceLinkEndType = receiver;
        break;
    case angular_position:
        referenceLinkEndType = receiver;
        break;
    case one_way_doppler:
    case doppler_measured_frequency:
        referenceLinkEndType = receiver;
        break;
    case one_way_differenced_range:
        referenceLinkEndType = receiver;
        break;
    case n_way_range:
        referenceLinkEndType = receiver;
        break;
    case two_way_doppler:
        referenceLinkEndType = receiver;
        break;
    case position_observable:
        referenceLinkEndType = observed_body;
        break;
    case velocity_observable:
        referenceLinkEndType = observed_body;
        break;
    case euler_angle_313_observable:
        referenceLinkEndType = observed_body;
        break;
    case relative_angular_position:
        referenceLinkEndType = receiver;
        break;
    case n_way_differenced_range:
        referenceLinkEndType = receiver;
        break;
    case dsn_one_way_averaged_doppler:
        referenceLinkEndType = receiver;
        break;
    case dsn_n_way_averaged_doppler:
        referenceLinkEndType = receiver;
        break;
    case relative_position_observable:
        referenceLinkEndType = observed_body;
        break;
    default:
        throw std::runtime_error( "Error, default reference link end not defined for observable " +
                                  std::to_string( observableType ) );
    }
    return referenceLinkEndType;
}

int getNumberOfLinksInObservable(
        const ObservableType observableType, const int numberOfLinkEnds )
{
    int numberOfLinks = -1;
    switch( observableType )
    {
    case one_way_range:
        numberOfLinks = 1;
        break;
    case angular_position:
        numberOfLinks = 1;
        break;
    case one_way_doppler:
        numberOfLinks = 1;
        break;
    case one_way_differenced_range:
        numberOfLinks = 1;
        break;
    case dsn_one_way_averaged_doppler:
        numberOfLinks = 1;
        break;
    case n_way_range:
        if( numberOfLinkEnds < 0 )
        {
            throw std::runtime_error( "Error when determining number of links for n-way range: number of link ends not provided" );
        }
        numberOfLinks = numberOfLinkEnds - 1;
        break;
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
        if( numberOfLinkEnds < 0 )
        {
            throw std::runtime_error( "Error when determining number of links for " + std::to_string( observableType ) +
                ": number of link ends not provided." );
        }
        numberOfLinks = numberOfLinkEnds - 1;
        break;
    case two_way_doppler:
    case doppler_measured_frequency:
        numberOfLinks = 2;
        break;
    case position_observable:
        numberOfLinks = 0;
        break;
    case velocity_observable:
        numberOfLinks = 0;
        break;
    case euler_angle_313_observable:
        numberOfLinks = 0;
        break;
    case relative_angular_position:
        numberOfLinks = 3;
        break;
    case relative_position_observable:
        numberOfLinks = 0;
        break;
    default:
        throw std::runtime_error( "Error, number of links not defined for observable " +
                                  std::to_string( observableType ) );
    }
    return numberOfLinks;
}

std::vector< LinkEndType > getLinkEndTypesForGivenLinkEndId(
        const LinkEnds& linkEnds,
        const LinkEndId linkEndToCheck )
{
    std::vector< LinkEndType > linkEndTypeList;
    for( auto linkEndIterator : linkEnds )
    {
        if( linkEndToCheck == linkEndIterator.second )
        {
            linkEndTypeList.push_back( linkEndIterator.first );
        }
    }
    return linkEndTypeList;
}

std::vector< int > getLinkEndIndicesForLinkEndIdAtObservable(
        const ObservableType observableType, const LinkEnds& linkEnds, const LinkEndId linkEndToCheck )
{
    std::vector< LinkEndType > linkEndTypes = getLinkEndTypesForGivenLinkEndId(
            linkEnds, linkEndToCheck );
    std::vector< int > totalLinkEndIndices;
    for( unsigned int i = 0; i < linkEndTypes.size( ); i++ )
    {
        std::vector< int > currentLinkEndIndices = getLinkEndIndicesForLinkEndTypeAtObservable(
                observableType, linkEndTypes.at( i ), linkEnds.size( ) );
        totalLinkEndIndices.insert(totalLinkEndIndices.end(), currentLinkEndIndices.begin(), currentLinkEndIndices.end());
    }
    return totalLinkEndIndices;
}


//! Function to retrieve the link end indices in link end states/times that are to be used in viability calculation
std::vector< std::pair< int, int > > getLinkStateAndTimeIndicesForLinkEnd(
        const LinkEnds& linkEnds, const ObservableType observableType,  const LinkEndId linkEndToCheck )
{
    std::vector< std::pair< int, int > > linkEndIndices;

    switch( observableType )
    {
    case one_way_range:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) ||
                ( ( linkEnds.at( transmitter ).bodyName_ == linkEndToCheck.bodyName_ ) && ( linkEndToCheck.stationName_ == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                linkEndToCheck.stationName_ == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way range link end types for link end indices" );
        }
        break;
    case one_way_doppler:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                  ( linkEndToCheck.stationName_ == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                linkEndToCheck.stationName_ == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way doppler link end types for link end indices" );
        }
        break;
    case two_way_doppler:
    case doppler_measured_frequency:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                  ( linkEndToCheck.stationName_ == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }

        if( linkEnds.at( reflector1 ) == linkEndToCheck || ( ( linkEnds.at( reflector1 ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                linkEndToCheck.stationName_ == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 2, 3 ) );
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }

        if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                linkEndToCheck.stationName_ == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 3, 2 ) );
        }

        if( linkEndIndices.size( ) == 0 )
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way doppler link end types for link end indices" );
        }
        break;
    case one_way_differenced_range:
        if( linkEnds.at( transmitter ) == linkEndToCheck || ( ( linkEnds.at( transmitter ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                              linkEndToCheck.stationName_ == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
            linkEndIndices.push_back( std::make_pair( 2, 3 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                linkEndToCheck.stationName_ == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
            linkEndIndices.push_back( std::make_pair( 3, 2 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant " + std::to_string( observableType ) + " link end types for link end indices" );
        }
        break;
    case dsn_one_way_averaged_doppler:
        if( linkEnds.at( transmitter ) == linkEndToCheck || ( ( linkEnds.at( transmitter ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                              linkEndToCheck.stationName_ == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
            linkEndIndices.push_back( std::make_pair( 2, 3 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                linkEndToCheck.stationName_ == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
            linkEndIndices.push_back( std::make_pair( 3, 2 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant " + std::to_string( observableType ) + " link end types for link end indices" );
        }
        break;
    case n_way_range:
    {
        std::vector< int > matchingLinkEndIndices = getNWayLinkEndIndicesFromLinkEndId( linkEndToCheck, linkEnds );
        if( matchingLinkEndIndices.size( ) > 0 )
        {
            for( unsigned int i = 0; i < matchingLinkEndIndices.size( ); i++ )
            {
                if( matchingLinkEndIndices.at( i ) == 0 )
                {
                    linkEndIndices.push_back( std::make_pair( 0, 1 ) );
                }
                else if( matchingLinkEndIndices.at( i ) == static_cast< int >( linkEnds.size( ) )  - 1 )
                {
                    linkEndIndices.push_back( std::make_pair( 2 * ( linkEnds.size( ) - 1 ) - 1,
                                                              2 * ( linkEnds.size( ) - 1 ) - 2 ) );
                }
                else
                {
                    linkEndIndices.push_back(
                                std::make_pair( 2 * matchingLinkEndIndices.at( i ),
                                                2 * matchingLinkEndIndices.at( i ) + 1 ) );
                    linkEndIndices.push_back(
                                std::make_pair( 2 * matchingLinkEndIndices.at( i ) - 1,
                                                2 * matchingLinkEndIndices.at( i ) - 2 ) );
                }
            }
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant n-way range link end types for link end indices" );
        }
        break;
    }
    case n_way_differenced_range:
    case dsn_n_way_averaged_doppler:
    {
        int undifferenceNumberOfEntries = 2 * ( linkEnds.size( ) - 1 );
        std::vector< int > matchingLinkEndIndices = getNWayLinkEndIndicesFromLinkEndId( linkEndToCheck, linkEnds );
        if( matchingLinkEndIndices.size( ) > 0 )
        {
            for( unsigned int i = 0; i < matchingLinkEndIndices.size( ); i++ )
            {
                if( matchingLinkEndIndices.at( i ) == 0 )
                {
                    linkEndIndices.push_back( std::make_pair( 0, 1 ) );
                    linkEndIndices.push_back( std::make_pair( undifferenceNumberOfEntries, undifferenceNumberOfEntries + 1 ) );

                }
                else if( matchingLinkEndIndices.at( i ) == static_cast< int >( linkEnds.size( ) )  - 1 )
                {
                    linkEndIndices.push_back( std::make_pair( 2 * ( linkEnds.size( ) - 1 ) - 1,
                                                              2 * ( linkEnds.size( ) - 1 ) - 2 ) );
                    linkEndIndices.push_back( std::make_pair( undifferenceNumberOfEntries + 2 * ( linkEnds.size( ) - 1 ) - 1,
                                                              undifferenceNumberOfEntries + 2 * ( linkEnds.size( ) - 1 ) - 2 ) );
                }
                else
                {
                    linkEndIndices.push_back( std::make_pair(
                                                  2 * matchingLinkEndIndices.at( i ),
                                                  2 * matchingLinkEndIndices.at( i ) + 1 ) );
                    linkEndIndices.push_back( std::make_pair(
                                                  2 * matchingLinkEndIndices.at( i ) - 1,
                                                  2 * matchingLinkEndIndices.at( i ) - 2 ) );
                    linkEndIndices.push_back( std::make_pair(
                                                  undifferenceNumberOfEntries + 2 * matchingLinkEndIndices.at( i ),
                                                  undifferenceNumberOfEntries + 2 * matchingLinkEndIndices.at( i ) + 1 ) );
                    linkEndIndices.push_back( std::make_pair(
                                                  undifferenceNumberOfEntries + 2 * matchingLinkEndIndices.at( i ) - 1,
                                                  undifferenceNumberOfEntries + 2 * matchingLinkEndIndices.at( i ) - 2 ) );
                }
            }
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant " + std::to_string( observableType ) +
                " link end types for link end indices" );
        }
        break;
    }
    case angular_position:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                  ( linkEndToCheck.stationName_ == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                linkEndToCheck.stationName_ == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant angular position link end types for link end indices" );
        }
        break;
    case position_observable:

        throw std::runtime_error( "Error, parsed irrelevant position observable link end types for link end indices" );
        break;
    case euler_angle_313_observable:

        throw std::runtime_error( "Error, parsed irrelevant euler angle observable link end types for link end indices" );
        break;
    case velocity_observable:

        throw std::runtime_error( "Error, parsed irrelevant position observable link end types for link end indices" );
        break;
    case relative_angular_position:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                  ( linkEndToCheck.stationName_ == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 2 ) );
        }
        else if( ( linkEnds.at( transmitter2 ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter2 ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                  ( linkEndToCheck.stationName_ == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 2 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).bodyName_ == linkEndToCheck.bodyName_ ) &&
                                                                linkEndToCheck.stationName_ == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 2, 0 ) );
            linkEndIndices.push_back( std::make_pair( 2, 1 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant angular position link end types for link end indices" );
        }
        break;
    case relative_position_observable:

        throw std::runtime_error( "Error, parsed irrelevant relative position observable link end types for link end indices" );
        break;
    default:

        throw std::runtime_error( "Error, observable type " + std::to_string(
                                      observableType ) + " not recognized when making viability link ends" );

    }

    return linkEndIndices;
}

std::map< LinkEndType, int > getSingleLinkStateEntryIndices( const ObservableType observableType )
{
    std::map< LinkEndType, int > singleLinkStateEntries;
    if( observableType == one_way_range ||
            observableType == angular_position ||
            observableType == one_way_doppler )
    {
        singleLinkStateEntries = oneWayLinkStateEntries;
    }
    else if( observableType == position_observable ||
             observableType == velocity_observable ||
             observableType == euler_angle_313_observable )
    {
        singleLinkStateEntries = observedBodyLinkStateEntries;
    }
    else if ( observableType == relative_position_observable )
    {
        singleLinkStateEntries = observedObserverBodiesLinkStateEntries;
    }
    else
    {
        throw std::runtime_error( "Error when getting single link state entries, observable type " +
                                  getObservableName( observableType ) + " not recognized." );
    }
    return singleLinkStateEntries;
}


}

}
