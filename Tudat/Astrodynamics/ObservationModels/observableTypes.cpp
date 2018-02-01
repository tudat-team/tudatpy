/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace observation_models
{

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

}

}
