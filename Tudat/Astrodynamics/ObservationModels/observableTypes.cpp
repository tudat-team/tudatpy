#include <iostream>
#include <algorithm>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"

namespace tudat
{

std::string getObservableName( const ObservableType observableType )
{
    std::string observableName;
    switch( observableType )
    {
    case oneWayRange:
        observableName = "OneWayRange";
        break;
    case twoWayRange:
        observableName = "TwoWayRange";
        break;
    case threeWayRange:
        observableName = "ThreeWayRange";
        break;
    case oneWayDiffencedRangeRate:
        observableName = "OneWayRangeDifferencedRangeRate";
        break;
    case angular_position:
        observableName = "AngularPosition";
        break;
    case stateObservable:
        observableName = "State";
        break;
    default:
        std::cerr<<"Error, could not find observable type "<<observableType<<" when getting name of type"<<std::endl;
    }
    return observableName;
}

ObservableType getObservableType( const std::string& observableName )
{
    ObservableType observableType;

    if( observableName == "OneWayRange" )
    {
        observableType = oneWayRange;
    }
    else if( observableName == "TwoWayRange" )
    {
        observableType = twoWayRange;
    }
    else if( observableName == "ThreeWayRange" )
    {
        observableType = threeWayRange;
    }
    else if( observableName == "OneWayRangeDifferencedRangeRate" )
    {
        observableType = oneWayDiffencedRangeRate;
    }
    else if( observableName == "AngularPosition" )
    {
        observableType = angular_position;
    }
    else if( observableName == "State" )
    {
        observableType = stateObservable;
    }
    else
    {
        std::cerr<<"Error, could not find observable name "<<observableName<<" when getting type of name"<<std::endl;
    }

    return observableType;
}

//std::vector< LinkEndId > getActiveLinkEnds(
//        const std::map< ObservableType, std::vector< LinkEnds > >& linkEndMap )
//{
//    std::vector< LinkEndId > activeLinkEndList;
//    std::vector< LinkEndId >::iterator findIterator;

//    for( std::map< ObservableType, std::vector< LinkEnds > >::const_iterator linkEndIterator = linkEndMap.begin( );
//         linkEndIterator != linkEndMap.end( ); linkEndIterator++ )
//    {
//        for( unsigned int i = 0; i < linkEndIterator->second.size( ); i++ )
//        {
//            if( linkEndIterator->first == oneWayRange || linkEndIterator->first == twoWayRange ||
//                    linkEndIterator->first == threeWayRange || linkEndIterator->first == oneWayDiffencedRangeRate ||
//                    linkEndIterator->first == angular_position )
//            {
//                findIterator = std::find( activeLinkEndList.begin( ), activeLinkEndList.end( ), linkEndIterator->second[ i ].at( receiver ) );
//                if( findIterator == activeLinkEndList.end( ) )
//                {
//                    activeLinkEndList.push_back( linkEndIterator->second[ i ].at( receiver ) );
//                }

//                findIterator = std::find( activeLinkEndList.begin( ), activeLinkEndList.end( ), linkEndIterator->second[ i ].at( transmitter ) );
//                if( findIterator == activeLinkEndList.end( ) )
//                {
//                    activeLinkEndList.push_back( linkEndIterator->second[ i ].at( transmitter ) );
//                }
//            }
//            else if( linkEndIterator->first == stateObservable )
//            {
//                findIterator = std::find( activeLinkEndList.begin( ), activeLinkEndList.end( ), linkEndIterator->second[ i ].at( observed_body ) );
//                if( findIterator == activeLinkEndList.end( ) )
//                {
//                    activeLinkEndList.push_back( linkEndIterator->second[ i ].at( observed_body ) );
//                }
//            }
//            else if( linkEndIterator->first == positionObservable )
//            {
//                findIterator = std::find( activeLinkEndList.begin( ), activeLinkEndList.end( ), linkEndIterator->second[ i ].at( observed_body ) );
//                if( findIterator == activeLinkEndList.end( ) )
//                {
//                    activeLinkEndList.push_back( linkEndIterator->second[ i ].at( observed_body ) );
//                }
//            }
//            else
//            {
//                std::cerr<<"Error when identifying active link ends, observable "<<linkEndIterator->first<<" not found "<<std::endl;
//            }
//        }
//    }

//    return activeLinkEndList;
//}

//int getNumberOfLinkEndTimesForBaseObservable( const ObservableType observableType, const int numberOfLinkEnds )
//{
//    int numberOfLinkEndTimesForObservable_;

//    switch( observableType )
//    {
//    case positionObservable:
//        numberOfLinkEndTimesForObservable_ = 1;
//        break;
//    case oneWayRange:
//        numberOfLinkEndTimesForObservable_ = 2;
//        break;
//    case twoWayRange:
//        numberOfLinkEndTimesForObservable_ = 3;
//        break;
//    case threeWayRange:
//        numberOfLinkEndTimesForObservable_ = 3;
//        break;
//    case angular_position:
//        numberOfLinkEndTimesForObservable_ = 2;
//        break;
//    case oneWayDiffencedRangeRate:
//        numberOfLinkEndTimesForObservable_ = 4;
//        break;
//    case twoWayDifferencedRangeRate:
//        numberOfLinkEndTimesForObservable_ = 6;
//        break;
//    case oneWayTimeTransfer:
//        numberOfLinkEndTimesForObservable_ = 2;
//        break;
//    case nWayRange:
//        if( numberOfLinkEnds < 2 )
//        {
//            std::cerr<<"Error when getting number of link end times for n way range,found "<<numberOfLinkEnds<<" link ends"<<std::endl;
//        }
//        else
//        {
//            numberOfLinkEndTimesForObservable_ = 2 * numberOfLinkEnds;
//        }
//        break;
//    case nWayRangeRate:
//        if( numberOfLinkEnds < 2 )
//        {
//            std::cerr<<"Error when getting number of link end times for n way range,found "<<numberOfLinkEnds<<" link ends"<<std::endl;
//        }
//        else
//        {
//            numberOfLinkEndTimesForObservable_ = 4 * numberOfLinkEnds;
//        }
//        break;
//    default:
//        std::cerr<<"Error, could not find number of link end times for obsdervation of type "<<observableType<<std::endl;
//    }

//    return numberOfLinkEndTimesForObservable_;
//}

//std::vector< int > getLinkEndIndicesForLinkEndTypeAtObservable(
//        const ObservableType observableType, const LinkEndType linkEndType, const int numberOfLinkEnds )
//{
//    std::vector< int > linkEndIndices;

//    switch( observableType )
//    {
//    case oneWayRange:
//        switch( linkEndType )
//        {
//        case transmitter:
//            linkEndIndices.push_back( 0 );
//            break;
//        case receiver:
//            linkEndIndices.push_back( 1 );
//            break;
//        default:
//            std::cerr<<"Error, could not find link end type index for link end "<<linkEndType<<" of observable "<<observableType<<std::endl;
//        }
//        break;
//    case oneWayTimeTransfer:
//        switch( linkEndType )
//        {
//        case transmitter:
//            linkEndIndices.push_back( 0 );
//            break;
//        case receiver:
//            linkEndIndices.push_back( 1 );
//            break;
//        default:
//            std::cerr<<"Error, could not find link end type index for link end "<<linkEndType<<" of observable "<<observableType<<std::endl;
//        }
//        break;
//    case twoWayRange:
//        switch( linkEndType )
//        {
//        case transmitter:
//            linkEndIndices.push_back( 0 );
//            break;
//        case reflector:
//            linkEndIndices.push_back( 1 );
//            break;
//        case receiver:
//            linkEndIndices.push_back( 2 );
//            break;
//        default:
//            std::cerr<<"Error, could not find link end type index for link end "<<linkEndType<<" of observable "<<observableType<<std::endl;
//        }
//        break;
//    case threeWayRange:
//        switch( linkEndType )
//        {
//        case transmitter:
//            linkEndIndices.push_back( 0 );
//            break;
//        case reflector:
//            linkEndIndices.push_back( 1 );
//            break;
//        case receiver:
//            linkEndIndices.push_back( 2 );
//            break;
//        default:
//            std::cerr<<"Error, could not find link end type index for link end "<<linkEndType<<" of observable "<<observableType<<std::endl;
//        }
//        break;
//    case angular_position:
//        switch( linkEndType )
//        {
//        case transmitter:
//            linkEndIndices.push_back( 0 );
//            break;
//        case receiver:
//            linkEndIndices.push_back( 1 );
//            break;
//        default:
//            std::cerr<<"Error, could not find link end type index for link end "<<linkEndType<<" of observable "<<observableType<<std::endl;
//        }
//        break;
//    case oneWayDiffencedRangeRate:
//        switch( linkEndType )
//        {
//        case transmitter:
//            linkEndIndices.push_back( 0 );
//            linkEndIndices.push_back( 2 );
//            break;
//        case receiver:
//            linkEndIndices.push_back( 1 );
//            linkEndIndices.push_back( 3 );
//            break;
//        default:
//            std::cerr<<"Error, could not find link end type index for link end "<<linkEndType<<" of observable "<<observableType<<std::endl;
//        }
//        break;
//    case twoWayDifferencedRangeRate:
//        switch( linkEndType )
//        {
//        case transmitter:
//            linkEndIndices.push_back( 0 );
//            linkEndIndices.push_back( 3 );
//            break;
//        case reflector:
//            linkEndIndices.push_back( 1 );
//            linkEndIndices.push_back( 4 );
//            break;
//        case receiver:
//            linkEndIndices.push_back( 2 );
//            linkEndIndices.push_back( 5 );
//            break;
//        default:
//            std::cerr<<"Error, could not find link end type index for link end "<<linkEndType<<" of observable "<<observableType<<std::endl;
//        }
//        break;
//    case nWayRange:
//        if( numberOfLinkEnds < 2 )
//        {
//            std::cerr<<"Error when getting n way range link end indices, not enough link ends"<<std::endl;
//        }
//        if( linkEndType == transmitter )
//        {
//            linkEndIndices.push_back( 0 );
//        }
//        else if( linkEndType == receiver )
//        {
//            linkEndIndices.push_back( 2 * ( numberOfLinkEnds - 1 ) - 1 );
//        }
//        else
//        {
//            int linkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndType, numberOfLinkEnds );
//            linkEndIndices.push_back( 2 * linkEndIndex - 1 );
//            linkEndIndices.push_back( 2 * linkEndIndex );
//        }
//        break;
//    case nWayRangeRate:
//        if( numberOfLinkEnds < 2 )
//        {
//            std::cerr<<"Error when getting n way range link end indices, not enough link ends"<<std::endl;
//        }
//        if( linkEndType == transmitter )
//        {
//            linkEndIndices.push_back( 0 );
//            linkEndIndices.push_back( 2 * ( numberOfLinkEnds - 1 ) );
//        }
//        else if( linkEndType == receiver )
//        {
//            linkEndIndices.push_back( 2 * ( numberOfLinkEnds - 1 ) - 1 );
//            linkEndIndices.push_back( 4 * ( numberOfLinkEnds - 1 ) - 1 );
//        }
//        else
//        {
//            int linkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndType, numberOfLinkEnds );
//            linkEndIndices.push_back( 2 * linkEndIndex - 1 );
//            linkEndIndices.push_back( 2 * linkEndIndex );
//            linkEndIndices.push_back( 2 * ( ( numberOfLinkEnds - 1 ) + linkEndIndex ) - 1 );
//            linkEndIndices.push_back( 2 * ( ( numberOfLinkEnds - 1 ) + linkEndIndex ) );
//        }
//        break;
//    case positionObservable:
//        if( linkEndType == observed_body )
//        {
//            linkEndIndices.push_back( 0 );
//        }
//        else
//        {
//            std::cerr<<"Error, could not find link end type index for link end "<<linkEndType<<" of observable "<<observableType<<std::endl;
//        }

//    default:
//        std::cerr<<"Error, could not find link end type index of observable "<<observableType<<std::endl;
//    }

//    return linkEndIndices;
//}

int getObservableSize( const ObservableType observableType )
{
    int observableSize = TUDAT_NAN;
    switch( observableType )
    {
    case oneWayRange:
        observableSize = 1;
        break;
    case twoWayRange:
        observableSize = 1;
        break;
    case threeWayRange:
        observableSize = 1;
        break;
    case oneWayDiffencedRangeRate:
        observableSize = 1;
        break;
    case angular_position:
        observableSize = 2;
        break;
    case stateObservable:
        observableSize = 6;
        break;
    case multiBaselineRangeRate:
        observableSize = 1;
        break;
    case twoWayDifferencedRangeRate:
        observableSize = 1;
        break;
    case  nWayRange:
        observableSize = 1;
        break;
    case nWayRangeRate:
        observableSize = 1;
        break;
    case oneWayTimeTransfer:
        observableSize = 1;
        break;
    case differencedBaselineObservation:
        observableSize = -1;
        break;
    case positionObservable:
        observableSize = 3;
        break;
    case oneWayDoppler:
        observableSize = 1;
        break;
    default:
        std::cerr<<"Error, did not recognize observable "<<observableType<<", when getting observable size"<<std::endl;
    }
    return observableSize;
}

}
