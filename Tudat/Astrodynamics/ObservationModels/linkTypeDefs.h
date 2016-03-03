#ifndef LINKTYPEDEFS_H
#define LINKTYPEDEFS_H

#include <map>
#include <string>
#include <vector>

namespace tudat
{

//! Enum defining different link end types.
/*!
 *  Enum defining different link end types. Enum is used for book-keeping purposes when setting link characteristics
 *  so that one does not need to remember the numerical order of types of link ends in input vector.
 */
enum LinkEndType
{
    transmitter = 0,
    reflector = 1,
    receiver = 2,
    observed_body = 3
};

typedef std::pair< std::string, std::string > LinkEndId;

typedef std::map< LinkEndType, LinkEndId > LinkEnds;

////! Function to split a list of link ends into names of vehicles and names of ground stations per celestial body.
///*!
// *  Function to split a list of link ends into names of vehicles and names of ground stations per celestial body.
// *  \param linkEndList List of all link ends ().
// *  \return Pair with First: Map with as key bodies on which ground stations are located and as values list of ground stations on said body. Second:
// *  Vector of names of bodies without ground stations (typically Vehicles) in link end list.
// */
//std::pair< std::map< std::string, std::vector< std::string > >, std::vector< std::string > > getBodiesWithAndWithoutGroundStations(
//        const std::vector< LinkEndId >& linkEndList );

//std::string getLinkEndNames( const LinkEnds linkEndsToPrint );

//void printLinkEnds( const LinkEnds linkEndsToPrint );

//std::string getLinkEndTypeName( const LinkEndType linkEndType );

//LinkEndType getLinkEndType( const std::string linkEndTypeName );

//void printLinkEndType( const LinkEndType linkEndType );

//int convertReceiverEnumToIndex( const LinkEndType linkEndType );

//std::vector< int > getNWayLinkEndIndicesFromLinkEndId( const LinkEndId& linkEndid, const LinkEnds& linkEnds );

//std::vector< int > getNWayLinkEndIndicesFromLinkEndId( const std::vector< LinkEndType >& linkEndTypes, const LinkEnds& linkEnds );

//std::vector< LinkEndType > getNWayLinkIndicesFromLinkEndId( const LinkEndId& linkEndid, const LinkEnds& linkEnds );

//int getNWayLinkIndexFromLinkEndType( const LinkEndType linkEndType, const int numberOfLinkEnds );

//LinkEndType getNWayLinkEnumFromIndex( const int linkEndIndex, const int numberOfLinkEnds );

//std::pair< LinkEndType, int > getVaryingLinkEndBaseTypeAndIndex( const LinkEndType varyingLinkEndType );

//LinkEndType getVaryingLinkTypeFromBaseTypeAndIndex( const LinkEndType baseLinkEndType, const int varyingLinkEndIndex );


}
#endif // LINKTYPEDEFS_H
