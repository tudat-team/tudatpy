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
enum ObservableType {
    undefined_observation_model = -1,
    one_way_range = 0,
    angular_position = 1,
    position_observable = 2,
    one_way_doppler = 3,
    one_way_differenced_range = 4,
    n_way_range = 5,
    two_way_doppler = 6,
    euler_angle_313_observable = 7,
    velocity_observable = 8,
    relative_angular_position = 9,
    n_way_differenced_range = 10,
    relative_position_observable = 11,
    dsn_one_way_averaged_doppler = 12,
    dsn_n_way_averaged_doppler = 13,
    doppler_measured_frequency = 14,
    dsn_n_way_range = 15,
    differenced_time_of_arrival = 16,
    one_way_doppler_measured_frequency = 17,
    differenced_frequency_of_arrival = 18,
    azimuth_elevation_angle = 19,
    pixel_coordinates = 20
};

// std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::VectorXd,
// std::pair< std::vector< double >, LinkEndType > > > > getTudatCompatibleObservationsAndTimes(
//         const std::vector< std::tuple< ObservableType, LinkEnds, Eigen::VectorXd,
//         std::vector< double >, LinkEndType > >& tudatpyObservationsAndTimes );

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

ObservableType getUndifferencedObservableType( const ObservableType differencedObservableType );

ObservableType getUnconcatenatedObservableType( const ObservableType observableType );

ObservableType getBaseObservableType( const ObservableType observableType );

std::pair< std::vector< int >, std::vector< int > > getUndifferencedTimeAndStateIndices( const ObservableType differencedObservableType,
                                                                                         const int numberOfLinkEnds );

std::pair< LinkEnds, LinkEnds > getUndifferencedLinkEnds( const ObservableType differencedObservableType,
                                                          const LinkEnds& differencedLinkEnds );

std::vector< LinkEnds > getUnconcatenatedLinkEnds( const ObservableType concatenatedObservableType, const LinkEnds& concatenatedLinkEnds );

//! Function to get the size of an observable of a given type.
/*!
 * Function to get the size of an observable of a given type.
 * \param observableType Type of observable.
 * \return Size of observable.
 */
int getObservableSize( const ObservableType observableType );

bool doesLinkEndTypeDefineId( const ObservableType observableType );

bool isObservableTypeMultiLink( const ObservableType observableType );

bool isObservableOfIntegratedType( const ObservableType observableType );

/*! Function indicating whether observable type requires transmitting ground station.
 *
 * Function indicating whether observable type requires transmitting ground station in system of bodies for the
 * observations to be simulated.
 * @param observableType Type of observable.
 * @return
 */
bool requiresTransmittingStation( const ObservableType observableType );

/*! Function indicating whether observable type requires receiving ground station.
 *
 * Function indicating whether observable type requires receiving ground station in system of bodies for the
 * observations to be simulated.
 * @param observableType Type of observable.
 * @return
 */
bool requiresFirstReceivingStation( const ObservableType observableType );

/*! Function indicating whether observable type requires a second receiving ground station.
 *
 * Function indicating whether observable type requires a second receiving ground station in system of bodies for the
 * observations to be simulated.
 * @param observableType Type of observable.
 * @return
 */
bool requiresSecondReceivingStation( const ObservableType observableType );

bool isRadiometricObservableType( const ObservableType observableType );

bool isPhaseVelocityBasedObservableType( const ObservableType observableType );

bool isGroupVelocityBasedObservableType( const ObservableType observableType );

bool observableCanHaveRetransmissionDelay( const ObservableType observableType );

bool linkEndIdDefinesSingleLink( const ObservableType observableType );

// bool areObservableLinksContinuous( const ObservableType observableType );

LinkEndType getDefaultReferenceLinkEndType( const ObservableType observableType );

int getNumberOfLinksInObservable( const ObservableType observableType, const int numberOfLinkEnds = -1 );

//! Function to get the indices in link end times/states for a given link end type and observable type
/*!
 * Function to get the indices in link end times/states for a given link end type and observable type
 * \param observableType Type of observable for which link end indices are to be returned
 * \param linkEndType Type of link end for which link end indices are to be returned
 * \param numberOfLinkEnds Number of link ends in observable
 * \return Indices in link end times/states for given link end type and observable type
 */
std::vector< int > getLinkEndIndicesForLinkEndTypeAtObservable( const ObservableType observableType,
                                                                const LinkEndType linkEndType,
                                                                const int numberOfLinkEnds );

//! Function to retrieve the link end indices in link end states/times that are to be used in viability calculation
/*!
 * Function to retrieve the link end indices in link end states/times that are to be used in viability calculation.
 * Return variable is a vector of pairs, where eacgetLinkEndTypesForGivenLinkEndIdh the first entry denotes the index of the point at which
 * the link is to be checkd. The second entry denotes the index for the opposite end of the link.
 * \param linkEnds Complete set of link ends for which check is to be performed
 * \param observableType Observable type for which check is to be performed
 * \param linkEndToCheck Link end at which check is to be performed
 * \return Link end indices in link end states/times that are to be used in viability calculation
 */
std::vector< std::pair< int, int > > getLinkStateAndTimeIndicesForLinkEnd( const LinkEnds& linkEnds,
                                                                           const ObservableType observableType,
                                                                           const LinkEndId linkEndToCheck );

std::vector< LinkEndType > getLinkEndTypesForGivenLinkEndId( const LinkEnds& linkEnds, const LinkEndId linkEndToCheck );

std::vector< int > getLinkEndIndicesForLinkEndIdAtObservable( const ObservableType observableType,
                                                              const LinkEnds& linkEnds,
                                                              const LinkEndId linkEndToCheck );

static const std::map< LinkEndType, int > oneWayLinkStateEntries = { { transmitter, 0 }, { receiver, 1 } };

static const std::map< LinkEndType, int > observedObserverBodiesLinkStateEntries = { { observed_body, 0 }, { observer, 1 } };

static const std::map< LinkEndType, int > observedBodyLinkStateEntries = { { observed_body, 0 } };

std::map< LinkEndType, int > getSingleLinkStateEntryIndices( const ObservableType observableType );

//! Function retrieving link ends information for all interlinks for a given observable type and link ends
std::vector< std::pair< std::pair< LinkEndType, LinkEndId >, std::pair< LinkEndType, LinkEndId > > > getInterlinks(
        const ObservableType observableType,
        const LinkEnds& linkEnds );

//! Struct for defining wrapping range for residual components (e.g. angular residuals)
/*!
 * Struct for defining wrapping range for residual components (e.g. angular residuals).
 * Provides the minimum and maximum of the range to which a residual component should be
 * wrapped (e.g. [-pi, pi] or [0, 2*pi]).
 */
struct ResidualWrappingRange {
    //! Default constructor (no wrapping).
    ResidualWrappingRange( ): minimumRange( 0.0 ), maximumRange( 0.0 ) {}

    //! Constructor with explicit range.
    ResidualWrappingRange( const double minimumRange, const double maximumRange ):
        minimumRange( minimumRange ), maximumRange( maximumRange )
    {}

    //! Minimum of wrapping range.
    double minimumRange;

    //! Maximum of wrapping range.
    double maximumRange;

    //! Period of wrapping (convenience).
    double period( ) const
    {
        return maximumRange - minimumRange;
    }

    //! Center of wrapping range (convenience).
    double center( ) const
    {
        return ( maximumRange + minimumRange ) / 2.0;
    }
};

//! Function to check if an observable type requires residual wrapping.
/*!
 * Function to check if an observable type requires residual wrapping (e.g. for angular
 * observable types where residuals may have 2*pi discontinuities).
 * \param observableType Type of observable.
 * \return True if wrapping is required, false otherwise.
 */
bool isResidualWrappingRequired( const ObservableType observableType );

//! Function to get the wrapping ranges per component for an observable type.
/*!
 * Function to get the wrapping ranges per component for an observable type. Returns
 * an empty vector if no wrapping is needed. Each component has its own range
 * specification (period, center) so that different components of the same observable
 * can have different wrapping behavior.
 * \param observableType Type of observable.
 * \return Vector of wrapping ranges, one per component.
 */
std::vector< ResidualWrappingRange > getResidualWrappingRanges( const ObservableType observableType );

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVABLETYPES_H
