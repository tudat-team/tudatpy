/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References: Formulation for Observed and Computed Values of Deep Space Network Data Types for Navigation,
 *      T. Moyer (2000), DEEP SPACE COMMUNICATIONS AND NAVIGATION SERIES
 */

#ifndef TUDAT_TRANSMITTINGFREQUENCIES_H
#define TUDAT_TRANSMITTINGFREQUENCIES_H

#include "tudat/math/quadrature/trapezoidQuadrature.h"
#include "tudat/math/interpolators.h"

namespace tudat
{

namespace ground_stations
{

//! Class to compute the transmitted frequency of a ground station and its integral.
class StationFrequencyInterpolator
{
public:
    //! Constructor
    StationFrequencyInterpolator( ) { }

    //! Destructor
    virtual ~StationFrequencyInterpolator( ) { }

    /*! Templated function to compute the transmitted frequency at the specified time.
     *
     * Templated function to compute the transmitted frequency at the specified time.
     *
     * @param lookupTime Time at which to compute the frequency.
     * @return Frequency value.
     */
    template< typename ObservationScalarType = double, typename TimeType = double >
    ObservationScalarType getTemplatedCurrentFrequency( const TimeType& lookupTime );

    /*! Templated function to compute the integral of the transmitted frequency.
     *
     * Templated function to compute the integral of the transmitted frequency.
     *
     * @param quadratureStartTime Start time of integration interval.
     * @param quadratureEndTime End time of integration interval.
     * @return Frequency integral
     */
    template< typename ObservationScalarType = double, typename TimeType = double >
    ObservationScalarType getTemplatedFrequencyIntegral( const TimeType& quadratureStartTime, const TimeType& quadratureEndTime );

private:

    //! Get frequency (with long double as observation scalar type and double as time type).
    virtual double getCurrentFrequency( const double lookupTime ) = 0;

    //! Get frequency (with double as observation scalar type and Time as time type).
    virtual double getCurrentFrequency( const Time& lookupTime ) = 0;

    //! Get frequency (with long double as observation scalar type and double as time type).
    virtual long double getCurrentLongFrequency( const double lookupTime ) = 0;

    //! Get frequency (with long double as observation scalar type and Time as time type).
    virtual long double getCurrentLongFrequency( const Time& lookupTime ) = 0;

    //! Get frequency integral (with long double as observation scalar type and double as time type).
    virtual double getFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime ) = 0;

    //! Get frequency integral (with double as observation scalar type and Time as time type).
    virtual double getFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime ) = 0;

    //! Get frequency integral (with long double as observation scalar type and double as time type).
    virtual long double getLongFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime ) = 0;

    //! Get frequency integral (with long double as observation scalar type and Time as time type).
    virtual long double getLongFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime ) = 0;

};

//class ConstantFrequencyInterpolator: public StationFrequencyInterpolator< ObservationScalarType, TimeType >
//{
//public:
//    //! Constructor
//    ConstantFrequencyInterpolator( ObservationScalarType frequency ):
//        StationFrequencyInterpolator< ObservationScalarType, TimeType >( ),
//        frequency_( frequency )
//    { }
//
//    //! Destructor
//    ~ConstantFrequencyInterpolator( ) { }
//
//    ObservationScalarType getCurrentFrequency( const TimeType lookupTime )
//    {
//        return frequency_;
//    }
//
//    ObservationScalarType getFrequencyIntegral( const TimeType quadratureStartTime,
//                                                const TimeType quadratureEndTime )
//    {
//        return frequency_ * ( quadratureEndTime - quadratureStartTime );
//    }
//
//private:
//
//    ObservationScalarType frequency_;
//};

//! Class to compute the transmitted frequency of a ground station and its integral, for piecewise frequency (e.g. ramped
//! DSN stations)
class PiecewiseLinearFrequencyInterpolator: public StationFrequencyInterpolator
{
public:

    /*! Constructor
     *
     * Constructor. The end time of each ramp should coincide with the start time of the following one.
     *
     * @param startTimes Start time of each ramp
     * @param endTimes End time of each ramp
     * @param rampRates Rate of each ramp
     * @param startFrequency Start frequency of each ramp
     */
    PiecewiseLinearFrequencyInterpolator(
            const std::vector< double >& startTimes,
            const std::vector< double >& endTimes,
            const std::vector< double >& rampRates,
            const std::vector< double >& startFrequency ):
        StationFrequencyInterpolator( ),
        startTimes_( startTimes ), endTimes_( endTimes ), rampRates_( rampRates ), startFrequencies_( startFrequency )
    {
        // Check if dimensions of all vectors are consistent
        if ( startTimes_.size( ) != endTimes_.size( ) || startTimes_.size( ) != rampRates_.size( ) ||
               startTimes_.size( ) != startFrequencies_.size( ) )
        {
            throw std::runtime_error(
                    "Error when creating piecewise linear frequency interpolator: the dimensions of the specified vectors "
                    "are not consistent: start times (" + std::to_string( startTimes_.size( ) ) + "), end times (" +
                    std::to_string( endTimes_.size( ) ) + "), ramp rates (" + std::to_string( rampRates_.size( ) ) +
                    "), start frequencies (" + std::to_string( startFrequencies_.size( ) ) + ")." );
        }

        // Check if there are no discontinuities between end times and subsequent start times
        for ( unsigned int i = 1; i < startTimes_.size( ); ++i )
        {
            if ( startTimes_.at( i ) != endTimes_.at( i - 1 ) )
            {
                throw std::runtime_error(
                        "Error when creating piecewise linear frequency interpolator: discontinuity between ramp end "
                        "time (" + std::to_string( endTimes_.at( i - 1 ) ) + ") and start time of the following ramp (" +
                        std::to_string( startTimes_.at( i ) ) + ")." );
            }
        }

        startTimeLookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                startTimes_ );
    }

    /*! Templated function to compute the transmitted frequency at the specified time.
     *
     * Templated function to compute the transmitted frequency at the specified time. Frequency is computed according to
     * Eq. 13-60 of Moyer (2000).
     *
     * @param lookupTime Time at which to compute the frequency.
     * @return Frequency value.
     */
    template< typename ObservationScalarType = double, typename TimeType = double >
    ObservationScalarType computeCurrentFrequency( const TimeType lookupTime )
    {
        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour( lookupTime );

        if( lookupTime > endTimes_.at( lowerNearestNeighbour ) || lookupTime < startTimes_.at ( lowerNearestNeighbour ) )
        {
            throw std::runtime_error(
                    "Error when interpolating ramp reference frequency: look up time (" + std::to_string(
                            static_cast< double >( lookupTime ) ) +
                    ") is outside the ramp table interval (" + std::to_string( startTimes_.at( 0 ) ) + " to " +
                    std::to_string( startTimes_.back( ) ) + ")." );
        }

        return startFrequencies_.at( lowerNearestNeighbour ) +
               rampRates_.at( lowerNearestNeighbour ) * ( lookupTime - startTimes_.at( lowerNearestNeighbour ) );
    }

     /*! Templated function to compute the integral of the transmitted frequency.
     *
     * Templated function to compute the integral of the transmitted frequency. Integral is computed according to section
     * 13.3.2.2.2 of Moyer (2000). Generally the integral should only be computed using the time type Time, otherwise
     * issues related to numerical cancelation are likely to occur.
     *
     * @param quadratureStartTime Start time of integration interval.
     * @param quadratureEndTime End time of integration interval.
     * @return Frequency integral
     */
    template< typename ObservationScalarType = double, typename TimeType = Time >
    ObservationScalarType computeFrequencyIntegral( const TimeType quadratureStartTime,
                                                    const TimeType quadratureEndTime )
    {
        ObservationScalarType integral = 0;

        int startTimeLowestNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour( quadratureStartTime );
        int endTimeLowestNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour( quadratureEndTime );

        if ( startTimeLowestNearestNeighbour == endTimeLowestNearestNeighbour )
        {
            integral += static_cast< ObservationScalarType > ( quadratureEndTime - quadratureStartTime ) *
                    ( computeCurrentFrequency< ObservationScalarType, TimeType >( quadratureStartTime ) +
                    computeCurrentFrequency< ObservationScalarType, TimeType >( quadratureEndTime ) ) / 2.0;
        }
        else
        {
            ObservationScalarType timeDelta;

            // First partial ramp
            timeDelta = static_cast< ObservationScalarType >( endTimes_.at( startTimeLowestNearestNeighbour ) - quadratureStartTime );
            integral += timeDelta *
                    ( computeCurrentFrequency< ObservationScalarType, TimeType >( quadratureStartTime ) +
                            static_cast< ObservationScalarType >( rampRates_.at( startTimeLowestNearestNeighbour ) ) *
                    timeDelta / 2.0 );

            // Full ramps
            for( unsigned int i = startTimeLowestNearestNeighbour + 1; i < startTimes_.size( ) &&
                    endTimes_.at( i ) < quadratureEndTime; i++ )
            {
                timeDelta = static_cast< ObservationScalarType >( endTimes_.at( i ) ) - static_cast< ObservationScalarType >( startTimes_.at( i ) );
                integral += timeDelta * ( startFrequencies_.at( i ) + rampRates_.at( i ) * timeDelta / 2.0 );
            }

            // Final partial ramp
            timeDelta = static_cast< ObservationScalarType >( quadratureEndTime - startTimes_.at( endTimeLowestNearestNeighbour ) );
            integral += timeDelta *
                    ( startFrequencies_.at( endTimeLowestNearestNeighbour ) + rampRates_.at( endTimeLowestNearestNeighbour ) *
                    timeDelta / 2.0 );
        }

        return integral;
    }

    //! Function to retrieve ramp start times
    std::vector< double > getStartTimes ( )
    {
        return startTimes_;
    }

    //! Function to retrieve end start times
    std::vector< double > getEndTimes ( )
    {
        return endTimes_;
    }

    //! Function to retrieve the ramp rates
    std::vector< double > getRampRates ( )
    {
        return rampRates_;
    }

    //! Function to retrieve the ramp start frequencies
    std::vector< double > getStartFrequencies ( )
    {
        return startFrequencies_;
    }

private:

    //! Get frequency (with long double as observation scalar type and double as time type).
    virtual double getCurrentFrequency( const double lookupTime )
    {
        return computeCurrentFrequency< double, double >( lookupTime );
    }

    //! Get frequency (with double as observation scalar type and Time as time type).
    virtual double getCurrentFrequency( const Time& lookupTime )
    {
        return computeCurrentFrequency< double, Time >( lookupTime );
    }

    //! Get frequency (with long double as observation scalar type and double as time type).
    virtual long double getCurrentLongFrequency( const double lookupTime )
    {
         return computeCurrentFrequency< long double, double >( lookupTime );
    }

    //! Get frequency (with long double as observation scalar type and Time as time type).
    virtual long double getCurrentLongFrequency( const Time& lookupTime )
    {
         return computeCurrentFrequency< long double, Time >( lookupTime );
    }

    //! Get frequency integral (with long double as observation scalar type and double as time type).
    virtual double getFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        return computeFrequencyIntegral< double, double >( quadratureStartTime, quadratureEndTime );
    }

    //! Get frequency integral (with double as observation scalar type and Time as time type).
    virtual double getFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime )
    {
        return computeFrequencyIntegral< double, Time >( quadratureStartTime, quadratureEndTime );
    }

    //! Get frequency integral (with long double as observation scalar type and double as time type).
    virtual long double getLongFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        return computeFrequencyIntegral< long double, double >( quadratureStartTime, quadratureEndTime );
    }

    //! Get frequency integral (with long double as observation scalar type and Time as time type).
    virtual long double getLongFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime )
    {
        return computeFrequencyIntegral< long double, Time >( quadratureStartTime, quadratureEndTime );
    }

    //! Start time of each ramp
    std::vector< double > startTimes_;
    //! End time of each ramp
    std::vector< double > endTimes_;
    //! Rate of each ramp
    std::vector< double > rampRates_;
    //! Start frequency of each ramp
    std::vector< double > startFrequencies_;

    //! Lookup scheme to find the nearest ramp start time for a given time
    std::shared_ptr< interpolators::LookUpScheme< double > > startTimeLookupScheme_;

};

} // namespace ground_stations

} // namespace tudat

#endif //TUDAT_TRANSMITTINGFREQUENCIES_H
