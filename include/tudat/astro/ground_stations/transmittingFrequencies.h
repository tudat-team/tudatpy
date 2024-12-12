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
#include "tudat/astro/basic_astro/dateTime.h"

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
    template< typename ObservationScalarType = double, typename TimeType = Time >
    ObservationScalarType getTemplatedCurrentFrequency( const TimeType& lookupTime );

    /*! Templated function to compute the integral of the transmitted frequency.
     *
     * Templated function to compute the integral of the transmitted frequency.
     *
     * @param quadratureStartTime Start time of integration interval.
     * @param quadratureEndTime End time of integration interval.
     * @return Frequency integral
     */
    template< typename ObservationScalarType = double, typename TimeType = Time >
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

class ConstantFrequencyInterpolator: public StationFrequencyInterpolator
{
public:
    //! Constructor
    ConstantFrequencyInterpolator( double frequency ):
        StationFrequencyInterpolator( ),
        frequency_( frequency )
    { }

    //! Destructor
    ~ConstantFrequencyInterpolator( ) { }

    template< typename ObservationScalarType = double, typename TimeType = Time >
    ObservationScalarType computeCurrentFrequency( const TimeType lookupTime )
    {
        return frequency_;
    }

    template< typename ObservationScalarType = double, typename TimeType = Time >
    ObservationScalarType computeFrequencyIntegral( const TimeType quadratureStartTime,
                                                    const TimeType quadratureEndTime )
    {
        return frequency_ * ( quadratureEndTime - quadratureStartTime );
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

    double frequency_;
};

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
            const std::vector< Time >& startTimes,
            const std::vector< Time >& endTimes,
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
            // If there are discontinuities
            if ( startTimes_.at( i ) != endTimes_.at( i - 1 ) )
            {
//                // If the start and end times are inconsistent throw error
//                if ( endTimes_.at( i - 1 ) > startTimes_.at( i ) )
//                {
//                    throw std::runtime_error(
//                            "Error when creating piecewise linear frequency interpolator: inconsistency between ramp end "
//                            "time (" + std::to_string( double( endTimes_.at( i - 1 ) ) ) + ";" +
//                            basic_astrodynamics::getCalendarDateFromTime( endTimes_.at( i - 1 ) ).isoString( ) + ") and start time of the following ramp (" +
//                            std::to_string( double( startTimes_.at( i ) ) ) + ";" +
//                            basic_astrodynamics::getCalendarDateFromTime( startTimes_.at( i ) ).isoString( ) +
//                            "); the end is smaller than the start time." );
//                }
//                // If there are gaps in the data save that information
//                else
//                {
//                    invalidTimeBlocksStartTimes_.push_back( endTimes_.at( i - 1 ) );
//                    invalidTimeBlocksEndTimes_.push_back( startTimes_.at( i ) );
//                }
            }
        }

        startTimeLookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< Time > >( startTimes_ );

        if ( !invalidTimeBlocksStartTimes_.empty( ) )
        {
            invalidStartTimeLookupScheme_ = std::make_shared< interpolators::BinarySearchLookupScheme< Time > >( invalidTimeBlocksStartTimes_ );
        }
    }

    /*! Templated function to compute the transmitted frequency at the specified time.
     *
     * Templated function to compute the transmitted frequency at the specified time. Frequency is computed according to
     * Eq. 13-60 of Moyer (2000).
     *
     * @param lookupTime Time at which to compute the frequency.
     * @return Frequency value.
     */
    template< typename ObservationScalarType = double, typename TimeType = Time >
    ObservationScalarType computeCurrentFrequency( const TimeType lookupTimeOriginal )
    {
        TimeType lookupTime = lookupTimeOriginal;
        int lowerNearestNeighbour = -1;
        if( lookupTimeOriginal < startTimes_.at( 0 ) )
        {
            lowerNearestNeighbour = 0;
        }
        else
        {
            try
            {
                lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour( lookupTime );
            }
            catch ( const std::exception &caughtException )
            {
                std::string exceptionText = std::string( caughtException.what( ));
                throw std::runtime_error(
                    "Error when interpolating ramp reference frequency: look up time (" + std::to_string(
                        static_cast< double >( lookupTime )) + ", caught exception: " + exceptionText );
            }
        }

        return startFrequencies_.at( lowerNearestNeighbour ) +
               rampRates_.at( lowerNearestNeighbour ) * ( lookupTime - startTimes_.at( lowerNearestNeighbour ) );
    }

     /*! Templated function to compute the integral of the transmitted frequency.
     *
     * Templated function to compute the integral of the transmitted frequency. Integral is computed according to section
     * 13.3.2.2.2 of Moyer (2000). Generally the integral should only be computed using the time type Time, otherwise
     * issues related to numerical cancellation are likely to occur.
     *
     * @param quadratureStartTime Start time of integration interval.
     * @param quadratureEndTime End time of integration interval.
     * @return Frequency integral
     */
    template< typename ObservationScalarType = double, typename TimeType = Time >
    ObservationScalarType computeFrequencyIntegral( const TimeType quadratureStartTime,
                                                    const TimeType quadratureEndTime )
    {
        if ( invalidStartTimeLookupScheme_ != nullptr )
        {
            int invalidStartLowestNearestNeighbour = invalidStartTimeLookupScheme_->findNearestLowerNeighbour( quadratureStartTime );
            // Integral is valid if both quadrature start/end times are before or after the invalid block. Need to check
            // that they aren't before the block because the lookup scheme might return a lowest nearest neighbour to the
            // right of the lookup time
            if ( !(
                    ( quadratureStartTime > invalidTimeBlocksEndTimes_.at( invalidStartLowestNearestNeighbour ) &&
                    quadratureEndTime > invalidTimeBlocksEndTimes_.at( invalidStartLowestNearestNeighbour ) ) ||
                    ( quadratureStartTime < invalidTimeBlocksStartTimes_.at( invalidStartLowestNearestNeighbour ) &&
                    quadratureEndTime < invalidTimeBlocksStartTimes_.at( invalidStartLowestNearestNeighbour ) ) ) )
            {
                throw std::runtime_error(
                        "Error when integrating ramp reference frequency: look up time (" + std::to_string(
                            static_cast< double >( quadratureStartTime ) ) +
                        ") is in time interval without transmitted frequency (" +
                        std::to_string( double( invalidTimeBlocksStartTimes_.at ( invalidStartLowestNearestNeighbour ) ) ) + " to " +
                        std::to_string( double( invalidTimeBlocksEndTimes_.at ( invalidStartLowestNearestNeighbour ) ) ) + ")." );
            }
        }

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
    std::vector< Time > getStartTimes ( )
    {
        return startTimes_;
    }

    //! Function to retrieve ramp start time
    Time getStartTime ( )
    {
        return startTimes_.front( );
    }

    //! Function to retrieve ramp end times
    std::vector< Time > getEndTimes ( )
    {
        return endTimes_;
    }

    //! Function to retrieve ramp end time
    Time getEndTime( )
    {
        return endTimes_.back( );
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
    std::vector< Time > startTimes_;
    //! End time of each ramp
    std::vector< Time > endTimes_;
    //! Rate of each ramp
    std::vector< double > rampRates_;
    //! Start frequency of each ramp
    std::vector< double > startFrequencies_;

    //! Start and end times of blocks where no frequency was transmitted
    std::vector< Time > invalidTimeBlocksStartTimes_;
    std::vector< Time > invalidTimeBlocksEndTimes_;

    //! Lookup scheme to find the nearest ramp start time for a given time
    std::shared_ptr< interpolators::LookUpScheme< Time > > startTimeLookupScheme_;

    //! Lookup scheme to find the nearest start time of the blocks without frequency transmission
    std::shared_ptr< interpolators::LookUpScheme< Time > > invalidStartTimeLookupScheme_;

};

} // namespace ground_stations

} // namespace tudat

#endif //TUDAT_TRANSMITTINGFREQUENCIES_H
