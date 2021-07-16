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
        isIntegratedType = false;
        break;
    case euler_angle_313_observable:
        isIntegratedType = false;
        break;
    case velocity_observable:
        isIntegratedType = false;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type is integrated; observable " +
                                  getObservableName( observableType ) + " not found" );
    }
    return isIntegratedType;
}

bool areObservableLinksContinuous( const ObservableType observableType )
{
    bool isTypeContinuous = true;
    switch( observableType )
    {
    case one_way_range:
        isTypeContinuous = true;
        break;
    case angular_position:
        isTypeContinuous = true;
        break;
    case position_observable:
        isTypeContinuous = true;
        break;
    case one_way_doppler:
        isTypeContinuous = true;
        break;
    case one_way_differenced_range:
        isTypeContinuous = true;
        break;
    case n_way_range:
        isTypeContinuous = true;
        break;
    case two_way_doppler:
        isTypeContinuous = true;
        break;
    case euler_angle_313_observable:
        isTypeContinuous = true;
        break;
    case velocity_observable:
        isTypeContinuous = true;
        break;
    default:
        throw std::runtime_error( "Error when determining if observable type is continuous; observable " +
                                  getObservableName( observableType ) + " not found" );
    }
    return isTypeContinuous;
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
        std::string numberOfWays = "N";
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

        observableName = numberOfWays + "WayRange";
        break;
    }
    case euler_angle_313_observable:
        observableName = "EulerAngle313";
        break;
    default:
        std::string errorMessage =
                "Error, could not find observable type " + std::to_string( observableType ) +
                " when getting name from type";
        throw std::runtime_error( errorMessage );
    }

    return observableName;
}

//! Function to get the observable type.ssociated with the name (string) of observable.
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
    else
    {
        std::string errorMessage =
                "Error, could not find observable name " + observableName +
                " when getting type from name";
        throw std::runtime_error( errorMessage );
    }

    return observableType;
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
            throw std::runtime_error( errorMessage );        }
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
            throw std::runtime_error( "Error when getting n way range link end indices, not enough link ends" );
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
    case n_way_range:
        if( numberOfLinkEnds < 0 )
        {
            throw std::runtime_error( "Error when determining number of links for n-way range: number of link ends not provided" );
        }
        numberOfLinks = numberOfLinkEnds - 1;
        break;
    case two_way_doppler:
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

void checkObservationResidualDiscontinuities(
        Eigen::Block< Eigen::VectorXd > observationResidualBlock,
        const ObservableType observableType )
{
    if( observableType == angular_position || observableType == euler_angle_313_observable )
    {
        for( int i = 1; i < observationResidualBlock.rows( ); i++ )
        {
            if( std::fabs( observationResidualBlock( i, 0 ) - observationResidualBlock( i - 1, 0 ) ) > 6.0 )
            {
                if( observationResidualBlock( i, 0 ) > 0 )
                {
                    observationResidualBlock( i, 0 ) = observationResidualBlock( i, 0 ) - 2.0 * mathematical_constants::PI;
                }
                else
                {
                    observationResidualBlock( i, 0 ) = observationResidualBlock( i, 0 ) + 2.0 * mathematical_constants::PI;
                }
            }
            else if( std::fabs( observationResidualBlock( i, 0 ) - observationResidualBlock( i - 1, 0 ) ) > 3.0 )
            {
                std::cerr<<"Warning, detected jump in observation residual of size "<<std::fabs( observationResidualBlock( i, 0 ) - observationResidualBlock( i - 1, 0 ) )<<
                           " for observable type "<<observableType<<std::endl;
            }
        }
    }
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
                ( ( linkEnds.at( transmitter ).first == linkEndToCheck.first ) && ( linkEndToCheck.second == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way range link end types for link end indices" );
        }
        break;
    case one_way_doppler:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter ).first == linkEndToCheck.first ) &&
                                                                  ( linkEndToCheck.second == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way doppler link end types for link end indices" );
        }
        break;
    case two_way_doppler:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter ).first == linkEndToCheck.first ) &&
                                                                  ( linkEndToCheck.second == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }

        if( linkEnds.at( reflector1 ) == linkEndToCheck || ( ( linkEnds.at( reflector1 ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 2, 3 ) );
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }

        if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 3, 2 ) );
        }

        if( linkEndIndices.size( ) == 0 )
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way doppler link end types for link end indices" );
        }
        break;
    case one_way_differenced_range:
        if( linkEnds.at( transmitter ) == linkEndToCheck || ( ( linkEnds.at( transmitter ).first == linkEndToCheck.first ) &&
                                                              linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
            linkEndIndices.push_back( std::make_pair( 2, 3 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
            linkEndIndices.push_back( std::make_pair( 3, 2 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way differenced link end types for link end indices" );
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
    case angular_position:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter ).first == linkEndToCheck.first ) &&
                                                                  ( linkEndToCheck.second == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
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
    default:
        throw std::runtime_error( "Error, observable type " + std::to_string(
                                      observableType ) + " not recognized when making viability link ends" );

    }

    return linkEndIndices;
}


}

}
