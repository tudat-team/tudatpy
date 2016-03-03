#ifndef OBSERVABLETYPES_H
#define OBSERVABLETYPES_H

#include <string>

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

enum ObservableType
{
    oneWayRange = 0,
    twoWayRange = 1,
    threeWayRange = 2,
    oneWayDiffencedRangeRate = 3,
    angular_position = 4,
    stateObservable = 5,
    multiBaselineRangeRate = 6,
    twoWayDifferencedRangeRate = 7,
    nWayRange = 8,
    nWayRangeRate = 9,
    oneWayTimeTransfer = 10,
    differencedBaselineObservation = 11,
    positionObservable = 12,
    oneWayDoppler = 13

};

std::string getObservableName( const ObservableType observableType );

ObservableType getObservableType( const std::string& observableName );

//int getNumberOfLinkEndTimesForBaseObservable( const ObservableType observableType, const int numberOfLinkEnds = -1 );

//std::vector< int > getLinkEndIndicesForLinkEndTypeAtObservable(
//        const ObservableType observableType, const LinkEndType linkEndType, const int numberOfLinkEnds = -1 );

//std::vector< LinkEndId > getActiveLinkEnds(
//        const std::map< ObservableType, std::vector< LinkEnds > >& linkEndMap );

int getObservableSize( const ObservableType observableType );

}

#endif // OBSERVABLETYPES_H
