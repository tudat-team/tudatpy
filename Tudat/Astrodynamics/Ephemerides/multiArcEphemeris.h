/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_MULTIARCEPHEMERIS_H
#define TUDAT_MULTIARCEPHEMERIS_H

#include <map>
#include <vector>

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Mathematics/Interpolators/lookupScheme.h"

#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace ephemerides
{

//! Class to define an ephemeris in an arc-wise manner
/*!
 *  Class to define an ephemeris in an arc-wise manner, where each arc is time-delimited and a separate ephemeris object
 *  is provided for each of these arcs.
 */
class MultiArcEphemeris: public Ephemeris
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param singleArcEphemerides Map of single arc ephemerides, with map key the minimum time at which the
     * ephemeris is valid. In case of arc overlaps, the arc with the highest start time is used to determine the state.
     *  \param referenceFrameOrigin Origin of reference frame (string identifier).
     *  \param referenceFrameOrientation Orientation of reference frame (string identifier).
     */
    MultiArcEphemeris(
            const std::map< double, boost::shared_ptr< Ephemeris > >& singleArcEphemerides,
            const std::string& referenceFrameOrigin = "",
            const std::string& referenceFrameOrientation = "" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation ),
        singleArcEphemerides_( utilities::createVectorFromMapValues( singleArcEphemerides ) ),
        arcStartTimes_( utilities::createVectorFromMapKeys( singleArcEphemerides ) )
    {
        // Create times at which the look up changes from one arc to the other.
        arcSplitTimes_ = arcStartTimes_;
        arcSplitTimes_.push_back( std::numeric_limits< double >::max( ) );

        // Create lookup scheme to determine which ephemeris to use.
        lookUpscheme_ = boost::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                    arcSplitTimes_ );
    }

    //! Destructor
    ~MultiArcEphemeris( ){ }

    //! Get state from ephemeris.
    /*!
     * Returns state from ephemeris at given Julian date.
     * \param secondsSinceEpoch Seconds since epoch (J2000) at which ephemeris is to be evaluated.
     * \return State from ephemeris.
     */
    Eigen::Vector6d getCartesianState(
            const double secondsSinceEpoch )
    {
        return singleArcEphemerides_.at( lookUpscheme_->findNearestLowerNeighbour( secondsSinceEpoch ) )->
                getCartesianState( double( secondsSinceEpoch ) );
    }

    //! Get state from ephemeris (long double state output).
    /*!
     * Returns state from ephemeris at given Julian date.
     * \param secondsSinceEpoch Seconds since epoch (J2000) at which ephemeris is to be evaluated.
     * \return State from ephemeris.
     */
    Eigen::Matrix< long double, 6, 1 > getCartesianLongState(
            const double secondsSinceEpoch )
    {
        return singleArcEphemerides_.at( lookUpscheme_->findNearestLowerNeighbour( secondsSinceEpoch ) )->
                getCartesianLongState( secondsSinceEpoch );
    }

    //! Get state from ephemeris (Time time input)
    /*!
     * Returns state from ephemeris at given Julian date.
     * \param currentTime Seconds since epoch (J2000) at which ephemeris is to be evaluated.
     * \return State from ephemeris.
     */
    Eigen::Vector6d getCartesianStateFromExtendedTime(
            const Time& currentTime )
    {
        return singleArcEphemerides_.at( lookUpscheme_->findNearestLowerNeighbour( currentTime ) )->
                getCartesianStateFromExtendedTime( currentTime );
    }

    //! Get state from ephemeris (long double state output and Time time input)
    /*!
     * Returns state from ephemeris at given Julian date.
     * \param currentTime Seconds since epoch (J2000) at which ephemeris is to be evaluated.
     * \return State currentTime ephemeris.
     */
    Eigen::Matrix< long double, 6, 1 > getCartesianLongStateFromExtendedTime(
            const Time& currentTime )
    {
        return singleArcEphemerides_.at( lookUpscheme_->findNearestLowerNeighbour( currentTime ) )->
                getCartesianLongStateFromExtendedTime( currentTime );
    }

    //! Function to reset the constituent arc ephemerides
    /*!
     * Function to reset the constituent arc ephemerides
     * \param singleArcEphemerides New list of arc ephemeris objects
     * \param arcStartTimes New list of ephemeris start times
     */
    void resetSingleArcEphemerides(
            const std::vector< boost::shared_ptr< Ephemeris > >& singleArcEphemerides,
            const std::vector< double >& arcStartTimes )
    {        
        singleArcEphemerides_ = singleArcEphemerides;
        arcStartTimes_ = arcStartTimes;

        // Create times at which the look up changes from one arc to the other.
        arcSplitTimes_ = arcStartTimes_;
        arcSplitTimes_.push_back(  std::numeric_limits< double >::max( ) );
        lookUpscheme_ = boost::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                    arcSplitTimes_ );
    }

    //! Function to reset the constituent arc ephemerides
    /*!
     * Function to reset the constituent arc ephemerides
     * \param singleArcEphemerides Map of of ephemeris objects as a function of their start times
     */
    void resetSingleArcEphemerides(
            const std::map< double, boost::shared_ptr< Ephemeris > >& singleArcEphemerides )
    {
        resetSingleArcEphemerides( utilities::createVectorFromMapValues( singleArcEphemerides ),
                                   utilities::createVectorFromMapKeys( singleArcEphemerides ) );
    }

    //! Function to retrieve times at which the look up changes from one arc to the other.
    /*!
     *  Function to retrieve times at which the look up changes from one arc to the other.
     *  \return Times at which the look up changes from one arc to the other.
     */
    std::vector< double > getArcSplitTimes( )
    {
        return arcSplitTimes_;
    }

    //! Function to retrieve the list of arc ephemeris objects
    /*!
     *  Function to retrieve the list of arc ephemeris objects
     *  \return List of arc ephemeris objects
     */
    std::vector< boost::shared_ptr< Ephemeris > > getSingleArcEphemerides( )
    {
        return singleArcEphemerides_;
    }


private:

    //! List of arc ephemeris objects
    std::vector< boost::shared_ptr< Ephemeris > > singleArcEphemerides_;

    //! List of ephemeris start times
    std::vector< double > arcStartTimes_;

    //! Times at which the look up changes from one arc to the other.
    std::vector< double > arcSplitTimes_;

    //! Lookup scheme to determine which ephemeris to use.
    boost::shared_ptr< interpolators::HuntingAlgorithmLookupScheme< double > > lookUpscheme_;


};

}

}

#endif // TUDAT_MULTIARCEPHEMERIS_H
