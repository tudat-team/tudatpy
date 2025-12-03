/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/basics/tudatExceptions.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace ephemerides
{

//! Class to define an ephemeris in an arc-wise manner
/*!
 *  Class to define an ephemeris in an arc-wise manner, where each arc is time-delimited and a separate ephemeris object
 *  is provided for each of these arcs.
 */
class MultiArcEphemeris : public Ephemeris
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
    MultiArcEphemeris( const std::map< double, std::shared_ptr< Ephemeris > >& singleArcEphemerides,
                       const std::string& referenceFrameOrigin = "",
                       const std::string& referenceFrameOrientation = "",
                       const std::shared_ptr< Ephemeris > defaultEphemeris = nullptr ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation ), defaultEphemeris_( defaultEphemeris )
    {
        resetSingleArcEphemerides( singleArcEphemerides );
    }

    //! Destructor
    ~MultiArcEphemeris( ) {}

    std::pair< bool, int > getCurrentEphemerisArc( const double currentTime );
    //! Get state from ephemeris.
    /*!
     * Returns state from ephemeris at given Julian date.
     * \param currentTime Seconds since epoch (J2000) at which ephemeris is to be evaluated.
     * \return State from ephemeris.
     */
    Eigen::Vector6d getCartesianState( const double currentTime );

    //! Get state from ephemeris (long double state output).
    /*!
     * Returns state from ephemeris at given Julian date.
     * \param currentTime Seconds since epoch (J2000) at which ephemeris is to be evaluated.
     * \return State from ephemeris.
     */
    Eigen::Matrix< long double, 6, 1 > getCartesianLongState( const double currentTime );

    //! Get state from ephemeris (Time time input)
    /*!
     * Returns state from ephemeris at given Julian date.
     * \param currentTime Seconds since epoch (J2000) at which ephemeris is to be evaluated.
     * \return State from ephemeris.
     */
    Eigen::Vector6d getCartesianStateFromExtendedTime( const Time& currentTime );

    //! Get state from ephemeris (long double state output and Time time input)
    /*!
     * Returns state from ephemeris at given Julian date.
     * \param currentTime Seconds since epoch (J2000) at which ephemeris is to be evaluated.
     * \return State currentTime ephemeris.
     */
    Eigen::Matrix< long double, 6, 1 > getCartesianLongStateFromExtendedTime( const Time& currentTime );

    //! Function to reset the constituent arc ephemerides
    /*!
     * Function to reset the constituent arc ephemerides
     * \param singleArcEphemerides New list of arc ephemeris objects
     * \param arcStartTimes New list of ephemeris start times
     */
    void resetSingleArcEphemerides( const std::vector< std::shared_ptr< Ephemeris > >& singleArcEphemerides,
                                    const std::vector< double >& arcStartTimes );

    //! Function to reset the constituent arc ephemerides
    /*!
     * Function to reset the constituent arc ephemerides
     * \param singleArcEphemerides Map of of ephemeris objects as a function of their start times
     */
    void resetSingleArcEphemerides( const std::map< double, std::shared_ptr< Ephemeris > >& singleArcEphemerides )
    {
        for( auto it : singleArcEphemerides )
        {
            std::cout<<"Resetting single arc ephemerides "<<it.first<<std::endl;
        }
        resetSingleArcEphemerides( utilities::createVectorFromMapValues( singleArcEphemerides ),
                                   utilities::createVectorFromMapKeys( singleArcEphemerides ) );
    }

    //! Function to retrieve the list of arc ephemeris objects
    /*!
     *  Function to retrieve the list of arc ephemeris objects
     *  \return List of arc ephemeris objects
     */
    std::vector< std::shared_ptr< Ephemeris > > getSingleArcEphemerides( )
    {
        return singleArcEphemerides_;
    }

private:
    void setArcInitialAndFinalTimes( );

    std::shared_ptr< Ephemeris > defaultEphemeris_;

    //! List of arc ephemeris objects
    std::vector< std::shared_ptr< Ephemeris > > singleArcEphemerides_;

    //! List of ephemeris start times
    std::vector< double > arcStartTimes_;

    std::vector< double > arcEndTimes_;

    //! Lookup scheme to determine which ephemeris to use.
    std::shared_ptr< interpolators::HuntingAlgorithmLookupScheme< double > > lookUpscheme_;
};

// Function that retrieves the time interval at which an ephemeris can be safely interrogated
/*
 * Function that retrieves the time interval at which an ephemeris can be safely interrogated. For most ephemeris types,
 * this function returns the full range of double values ( lowest( ) to max( ) ). For the tabulated ephemeris, the interval
 * on which the interpolator inside this object is valid is checked and returned
 * \param ephemerisModel Ephemeris model for which the interval is to be determined.
 * \return The time interval at which the ephemeris can be safely interrogated
 */
std::pair< double, double > getSafeEphemerisEvaluationInterval( const std::shared_ptr< ephemerides::Ephemeris > ephemerisModel );

}  // namespace ephemerides

}  // namespace tudat

#endif  // TUDAT_MULTIARCEPHEMERIS_H
