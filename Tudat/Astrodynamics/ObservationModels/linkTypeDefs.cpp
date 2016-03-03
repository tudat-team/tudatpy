//#include <iostream>
//#include <algorithm>

//#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

//namespace tudat
//{

//std::pair< std::map< std::string, std::vector< std::string > >, std::vector< std::string > > getBodiesWithAndWithoutGroundStations(
//        const std::vector< LinkEndId >& linkEndList )
//{
//    std::map< std::string, std::vector< std::string > > celestialBodiesWithGroundStations;
//    std::vector< std::string > bodiesWithoutGroundStation;

//    // Get identifiers of bodies without ground stations (vehicles) and bodies with grounds stations (celestial bodies)
//    for( unsigned int i = 0; i < linkEndList.size( ); i++ )
//    {
//        if( linkEndList.at( i ).second == "" )
//        {
//            if( std::find( bodiesWithoutGroundStation.begin( ), bodiesWithoutGroundStation.end( ), linkEndList.at( i ).first ) ==
//                    bodiesWithoutGroundStation.end( ) )
//            {
//                bodiesWithoutGroundStation.push_back( linkEndList.at( i ).first );
//            }
//        }
//        else
//        {
//            bool addGroundStation = 0;

//            if( celestialBodiesWithGroundStations.count( linkEndList.at( i ).first ) == 0 )
//            {
//                addGroundStation = 1;
//            }
//            else if( std::find( celestialBodiesWithGroundStations[ linkEndList.at( i ).first ].begin( ),
//                                celestialBodiesWithGroundStations[ linkEndList.at( i ).first ].end( ),
//                                linkEndList.at( i ).first ) == bodiesWithoutGroundStation.end( ) )
//            {
//                addGroundStation = 1;
//            }

//            if( addGroundStation )
//            {
//                celestialBodiesWithGroundStations[ linkEndList.at( i ).first ].push_back( linkEndList.at( i ).second );
//            }
//        }
//    }

//    return std::make_pair( celestialBodiesWithGroundStations, bodiesWithoutGroundStation );
//}


//std::string getLinkEndNames( const LinkEnds linkEndsToPrint )
//{
//    std::string linkEndName;
//    for( LinkEnds::const_iterator linkEndIterator = linkEndsToPrint.begin( ); linkEndIterator != linkEndsToPrint.end( ); linkEndIterator++ )
//    {
//        linkEndName = linkEndName + linkEndIterator->second.first + "," + linkEndIterator->second.second + " ";
//    }

//    return linkEndName;
//}

//void printLinkEnds( const LinkEnds linkEndsToPrint )
//{
//    for( LinkEnds::const_iterator linkEndIterator = linkEndsToPrint.begin( ); linkEndIterator != linkEndsToPrint.end( ); linkEndIterator++ )
//    {
//        printLinkEndType( linkEndIterator->first );
//        std::cout<<": ("<<linkEndIterator->second.first<<", "<<linkEndIterator->second.second<<"), ";
//    }
//}

//std::string getLinkEndTypeName( const LinkEndType linkEndType )
//{
//    std::string name;
//    switch( linkEndType )
//    {
//    case transmitter:
//    {
//        name = "Transmitter";
//        break;
//    }
//    case reflector:
//    {
//        name = "Reflector";
//        break;
//    }
//    case receiver:
//    {
//        name = "Receiver";
//        break;
//    }
//    case observed_body:
//    {
//        name = "Observed body";
//        break;
//    }
//    default:
//        std::cerr<<"Error when getting link end type name, "<<linkEndType<<" not recognized"<<std::endl;
//    }
//    return name;
//}

//LinkEndType getLinkEndType( const std::string linkEndTypeName )
//{
//    LinkEndType type;
//    if( linkEndTypeName == "Transmitter" )
//    {
//        type = transmitter;
//    }
//    else if( linkEndTypeName == "Reflector" )
//    {
//        type = reflector;
//    }
//    else if( linkEndTypeName == "Receiver" )
//    {
//        type = receiver;
//    }
//    else if( linkEndTypeName == "Observed body" )
//    {
//        type = observed_body;
//    }
//    else
//    {
//        std::cerr<<"Error when getting link end type, "<<linkEndTypeName<<" not recognized"<<std::endl;
//    }
//    return type;
//}

//void printLinkEndType( const LinkEndType linkEndType )
//{
//    std::cout<<getLinkEndTypeName( linkEndType );
//}

//int convertReceiverEnumToIndex( const LinkEndType linkEndType )
//{
//    if( linkEndType < receiver1 || linkEndType > receiver10 )
//    {
//        std::cerr<<"Error when converting receiver enum to index, enum is inconsistent"<<std::endl;
//    }

//    return static_cast< int >( linkEndType - receiver1 );
//}

//int getNWayLinkIndexFromEnum( const LinkEndType linkEndType, const int numberOfLinkEnds  )
//{
//    int linkEndIndex;
//    if( linkEndType == transmitter )
//    {
//        linkEndIndex = 0;
//    }
//    else if( linkEndType == receiver )
//    {
//        linkEndIndex = numberOfLinkEnds - 1;
//    }
//    else
//    {
//        if( linkEndType < reflector1 || linkEndType > reflector10 )
//        {
//            std::cerr<<"Error, found link end type "<<linkEndType<<" when getting n-way link end index"<<std::endl;
//        }
//        else
//        {
//            linkEndIndex = 1 + static_cast< int >( linkEndType - reflector1 );
//        }
//    }
//    return linkEndIndex;
//}

//std::vector< int > getNWayLinkEndIndicesFromLinkEndId( const LinkEndId& linkEndid, const LinkEnds& linkEnds )
//{
//    std::vector< LinkEndType >  matchingLinkEndTypes = getNWayLinkIndicesFromLinkEndId(
//                linkEndid, linkEnds );
//    return getNWayLinkEndIndicesFromLinkEndId( matchingLinkEndTypes, linkEnds );
//}

//std::vector< int > getNWayLinkEndIndicesFromLinkEndId( const std::vector< LinkEndType >& linkEndTypes, const LinkEnds& linkEnds )
//{
//    std::vector< int > linkEndIndices;
//    for( unsigned int i = 0; i < linkEndTypes.size( ); i++ )
//    {
//        linkEndIndices.push_back( getNWayLinkIndexFromLinkEndType( linkEndTypes.at( i ), linkEnds.size( ) ) );
//    }
//    return linkEndIndices;
//}


//std::vector< LinkEndType > getNWayLinkIndicesFromLinkEndId( const LinkEndId& linkEndid, const LinkEnds& linkEnds )
//{
//    std::vector< LinkEndType > matchingLinkEndTypes;

//    for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
//    {
//        if( linkEndIterator->second == linkEndid || ( ( linkEndIterator->second.first == linkEndid.first ) &&
//                                                                        linkEndid.second == "" ) )
//        {
//            matchingLinkEndTypes.push_back( linkEndIterator->first );
//        }
//    }
//    return matchingLinkEndTypes;
//}

//int getNWayLinkIndexFromLinkEndType( const LinkEndType linkEndType, const int numberOfLinkEnds )
//{
//    int linkEndIndex;
//    if( linkEndType == transmitter )
//    {
//        linkEndIndex = 0;
//    }
//    else if( linkEndType == receiver )
//    {
//        linkEndIndex = numberOfLinkEnds - 1;
//    }
//    else
//    {
//        linkEndIndex = static_cast< int >( linkEndType ) - static_cast< int >( reflector1 ) + 1;
//        if( linkEndIndex > numberOfLinkEnds - 2 )
//        {
//            std::cerr<<"Error when getting n-way link end index; value too large."<<std::endl;
//        }
//    }
//    return linkEndIndex;
//}


//LinkEndType getNWayLinkEnumFromIndex( const int linkEndIndex, const int numberOfLinkEnds )
//{
//    LinkEndType linkEndType;
//    if( linkEndIndex == 0 )
//    {
//        linkEndType = transmitter;
//    }
//    else if( linkEndIndex == numberOfLinkEnds - 1 )
//    {
//        linkEndType = receiver;
//    }
//    else if( linkEndIndex >= numberOfLinkEnds )
//    {
//        std::cerr<<"Error, found link end index "<<linkEndIndex<<" when getting n-way link end index for "<<
//                   numberOfLinkEnds<<" link end total."<<std::endl;
//    }
//    else
//    {
//        linkEndType = static_cast< LinkEndType >( static_cast< int >( reflector1 ) + ( linkEndIndex - 1 ) );
//    }

//    return linkEndType;
//}

//std::pair< LinkEndType, int > getVaryingLinkEndBaseTypeAndIndex( const LinkEndType varyingLinkEndType )
//{
//    LinkEndType baseLinkEndType;
//    int linkEndIndex;

//    if( varyingLinkEndType < receiver1 )
//    {
//        std::cerr<<"Error when getting varying link end base type, input was "<<varyingLinkEndType<<std::endl;
//    }
//    else if( varyingLinkEndType >= receiver1 && varyingLinkEndType <= receiver10 )
//    {
//        baseLinkEndType = receiver;
//        linkEndIndex = static_cast< int >( varyingLinkEndType - receiver1 );
//    }
//    else if( varyingLinkEndType >= reflector1 && varyingLinkEndType <= reflector10 )
//    {
//        baseLinkEndType = reflector;
//        linkEndIndex = static_cast< int >( varyingLinkEndType - reflector1 );
//    }
//    else if( varyingLinkEndType >= transmitter1 && varyingLinkEndType <= transmitter10 )
//    {
//        baseLinkEndType = transmitter;
//        linkEndIndex = static_cast< int >( varyingLinkEndType - transmitter1 );
//    }
//    else
//    {
//        std::cerr<<"Error when getting varying link end base type, could not recognize input "<<varyingLinkEndType<<std::endl;
//    }
//    return std::make_pair( baseLinkEndType, linkEndIndex );
//}

//LinkEndType getVaryingLinkTypeFromBaseTypeAndIndex( const LinkEndType baseLinkEndType, const int varyingLinkEndIndex )
//{
//    if( varyingLinkEndIndex > 10 )
//    {
//        std::cerr<<"Error when getting varying link end type, only 10 varying link ends supported, "<<varyingLinkEndIndex<<
//                   " are requested."<<std::endl;
//    }
//    LinkEndType varyingLinkEndType;
//    if( baseLinkEndType == receiver )
//    {
//        varyingLinkEndType = static_cast< LinkEndType >( static_cast< int >( receiver1 ) + varyingLinkEndIndex );
//    }
//    else if( baseLinkEndType == reflector )
//    {
//        varyingLinkEndType = static_cast< LinkEndType >( static_cast< int >( reflector1 ) + varyingLinkEndIndex );
//    }
//    else if( baseLinkEndType == transmitter )
//    {
//        varyingLinkEndType = static_cast< LinkEndType >( static_cast< int >( transmitter1 ) + varyingLinkEndIndex );
//    }
//    else
//    {
//        std::cerr<<"Error when getting varying link end type, base link end type is: "<<baseLinkEndType<<std::endl;
//    }
//    return varyingLinkEndType;
//}



//}
