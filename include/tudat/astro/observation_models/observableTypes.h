/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVABLETYPES_H
#define TUDAT_OBSERVABLETYPES_H

#include <string>

#include <Eigen/Core>

#include "tudat/astro/observation_models/linkTypeDefs.h"

namespace tudat
{

namespace observation_models
{


//! Enum for types of observations
enum ObservableType
{
    one_way_range = 0,
    angular_position = 1,
    position_observable = 2,
    one_way_doppler = 3,
    one_way_differenced_range = 4,
    n_way_range = 5,
    two_way_doppler = 6,
    euler_angle_313_observable = 7,
    velocity_observable = 8
};


//std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::VectorXd,
//std::pair< std::vector< double >, LinkEndType > > > > getTudatCompatibleObservationsAndTimes(
//        const std::vector< std::tuple< ObservableType, LinkEnds, Eigen::VectorXd,
//        std::vector< double >, LinkEndType > >& tudatpyObservationsAndTimes );

//! Function to get the name (string) associated with a given observable type.
/*!
 * Function to get the name (string) associated with a given observable type.
 * \param observableType Type of observable.
 * \param numberOfLinkEnds Number of link ends in observable
 * \return Name of observable
 */
std::string getObservableName( const ObservableType observableType, const int numberOfLinkEnds = 0 );

//! Function to get the observable type.ssociated with the name (string) of observable.
/*!
 * Function to get the observable type.ssociated with the name (string) of observable.
 * \param observableName of observable
 * \return Type of observable.
 */
ObservableType getObservableType( const std::string& observableName );

//! Function to get the size of an observable of a given type.
/*!
 * Function to get the size of an observable of a given type.
 * \param observableType Type of observable.
 * \return Size of observable.
 */
int getObservableSize( const ObservableType observableType );

bool isObservableOfIntegratedType( const ObservableType observableType );

bool areObservableLinksContinuous( const ObservableType observableType );

LinkEndType getDefaultReferenceLinkEndType(
        const ObservableType observableType );

int getNumberOfLinksInObservable(
        const ObservableType observableType, const int numberOfLinkEnds = -1 );

//! Function to get the indices in link end times/states for a given link end type and observable type
/*!
 * Function to get the indices in link end times/states for a given link end type and observable type
 * \param observableType Type of observable for which link end indices are to be returned
 * \param linkEndType Type of link end for which link end indices are to be returned
 * \param numberOfLinkEnds Number of link ends in observable
 * \return Indices in link end times/states for given link end type and observable type
 */
std::vector< int > getLinkEndIndicesForLinkEndTypeAtObservable(
        const ObservableType observableType, const LinkEndType linkEndType, const int numberOfLinkEnds );


//! Function to retrieve the link end indices in link end states/times that are to be used in viability calculation
/*!
 * Function to retrieve the link end indices in link end states/times that are to be used in viability calculation.
 * Return variable is a vector of pairs, where each the first entry denotes the index of the point at which the link is to be
 * checkd. The second entry denotes the index for the opposite end of the link.
 * \param linkEnds Complete set of link ends for which check is to be performed
 * \param observableType Observable type for which check is to be performed
 * \param linkEndToCheck Link end at which check is to be performed
 * \return Link end indices in link end states/times that are to be used in viability calculation
 */
std::vector< std::pair< int, int > > getLinkStateAndTimeIndicesForLinkEnd(
        const LinkEnds& linkEnds,
        const ObservableType observableType,
        const LinkEndId linkEndToCheck );

std::vector< LinkEndType > getLinkEndTypesForGivenLinkEndId(
        const LinkEnds& linkEnds,
        const LinkEndId linkEndToCheck );

void checkObservationResidualDiscontinuities(
        Eigen::Block< Eigen::VectorXd > observationBlock,
        const ObservableType observableType );


} // namespace observation_models

} // namespace tudat

#endif // TUDAT_OBSERVABLETYPES_H
