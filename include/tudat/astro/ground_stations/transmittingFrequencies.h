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

#include <iostream>

#include "tudat/math/quadrature/trapezoidQuadrature.h"
#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/astro/basic_astro/dateTime.h"

namespace tudat
{

namespace ground_stations
{

//! Options for handling queries in gaps between two transmitting frequency ramp intervals.
enum FrequencyGapHandling
{
    extrapolate_at_gaps,
    throw_exception_at_gaps,
    print_error_at_gaps,
    print_error_once_at_gaps
};

//! Class to compute the transmitted frequency of a ground station and its integral.
class StationFrequencyInterpolator
{
public:
    //! Constructor
    StationFrequencyInterpolator( ) {}

    //! Destructor
    virtual ~StationFrequencyInterpolator( ) {}

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

class ConstantFrequencyInterpolator : public StationFrequencyInterpolator
{
public:
    //! Constructor
    ConstantFrequencyInterpolator( double frequency ): StationFrequencyInterpolator( ), frequency_( frequency ) {}

    //! Destructor
    ~ConstantFrequencyInterpolator( ) {}

    template< typename ObservationScalarType = double, typename TimeType = Time >
    ObservationScalarType computeCurrentFrequency( const TimeType lookupTime )
    {
        return frequency_;
    }

    template< typename ObservationScalarType = double, typename TimeType = Time >
    ObservationScalarType computeFrequencyIntegral( const TimeType quadratureStartTime, const TimeType quadratureEndTime )
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
class PiecewiseLinearFrequencyInterpolator : public StationFrequencyInterpolator
{
public:
    /*! Constructor
     *
     * Constructor. By default, gaps between ramps are filled by extrapolating the previous ramp until the next
     * ramp starts.
     *
     * @param startTimes Start time of each ramp
     * @param endTimes End time of each ramp
     * @param rampRates Rate of each ramp
     * @param startFrequency Start frequency of each ramp
     * @param gapHandling Option for handling gaps between consecutive ramp intervals.
     */
    PiecewiseLinearFrequencyInterpolator( const std::vector< Time >& startTimes,
                                          const std::vector< Time >& endTimes,
                                          const std::vector< double >& rampRates,
                                          const std::vector< double >& startFrequency,
                                          const FrequencyGapHandling gapHandling = extrapolate_at_gaps ):
        StationFrequencyInterpolator( ), startTimes_( startTimes ), endTimes_( endTimes ), rampRates_( rampRates ),
        startFrequencies_( startFrequency ), gapHandling_( gapHandling ), hasPrintedGapWarning_( false )
    {
        initialize( );
    }

    void initialize( );

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
            catch( const std::exception& caughtException )
            {
                std::string exceptionText = std::string( caughtException.what( ) );
                throw std::runtime_error( "Error when interpolating ramp reference frequency: look up time (" +
                                          std::to_string( static_cast< double >( lookupTime ) ) + ", caught exception: " + exceptionText );
            }
        }

        if( invalidStartTimeLookupScheme_ != nullptr )
        {
            for( unsigned int i = 0; i < invalidTimeBlocksStartTimes_.size( ); i++ )
            {
                if( lookupTimeOriginal > invalidTimeBlocksStartTimes_.at( i ) &&
                    lookupTimeOriginal < invalidTimeBlocksEndTimes_.at( i ) )
                {
                    const std::string gapMessage = "Error when interpolating ramp reference frequency: look up time (" +
                            std::to_string( static_cast< double >( lookupTimeOriginal ) ) +
                            ") is in time interval without transmitted frequency (" +
                            std::to_string( double( invalidTimeBlocksStartTimes_.at( i ) ) ) + " to " +
                            std::to_string( double( invalidTimeBlocksEndTimes_.at( i ) ) ) + ").";
                    handleGap( gapMessage );
                }
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
    ObservationScalarType computeFrequencyIntegral( const TimeType quadratureStartTime, const TimeType quadratureEndTime )
    {
        if( quadratureEndTime < quadratureStartTime )
        {
            return -computeFrequencyIntegral< ObservationScalarType, TimeType >( quadratureEndTime, quadratureStartTime );
        }

        if( invalidStartTimeLookupScheme_ != nullptr )
        {
            for( unsigned int i = 0; i < invalidTimeBlocksStartTimes_.size( ); i++ )
            {
                if( quadratureStartTime < invalidTimeBlocksEndTimes_.at( i ) &&
                    quadratureEndTime > invalidTimeBlocksStartTimes_.at( i ) )
                {
                    const std::string gapMessage = "Error when integrating ramp reference frequency: time interval (" +
                            std::to_string( static_cast< double >( quadratureStartTime ) ) + " to " +
                            std::to_string( static_cast< double >( quadratureEndTime ) ) +
                            ") overlaps with time interval without transmitted frequency (" +
                            std::to_string( double( invalidTimeBlocksStartTimes_.at( i ) ) ) + " to " +
                            std::to_string( double( invalidTimeBlocksEndTimes_.at( i ) ) ) + ").";
                    handleGap( gapMessage );
                }
            }
        }

        ObservationScalarType integral = 0;

        int currentRamp = -1;
        try
        {
            if( quadratureStartTime < startTimes_.at( 0 ) )
            {
                currentRamp = 0;
            }
            else
            {
                currentRamp = startTimeLookupScheme_->findNearestLowerNeighbour( quadratureStartTime );
            }
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error(
                    "Error when determining current ramp in frequency integral: " + std::string( caughtException.what( ) ) +
                    ", possibly the ground station does not have data"
                    "defined at the requested time." );
        }

        TimeType currentTime = quadratureStartTime;
        while( currentTime < quadratureEndTime )
        {
            TimeType nextRampStartTime = quadratureEndTime;
            if( currentRamp + 1 < static_cast< int >( startTimes_.size( ) ) &&
                startTimes_.at( currentRamp + 1 ) < quadratureEndTime )
            {
                nextRampStartTime = startTimes_.at( currentRamp + 1 );
            }

            ObservationScalarType timeDelta = static_cast< ObservationScalarType >( nextRampStartTime - currentTime );
            ObservationScalarType frequencyAtCurrentTime = startFrequencies_.at( currentRamp ) +
                    static_cast< ObservationScalarType >( rampRates_.at( currentRamp ) ) *
                            static_cast< ObservationScalarType >( currentTime - startTimes_.at( currentRamp ) );
            ObservationScalarType frequencyAtNextRampStartTime = startFrequencies_.at( currentRamp ) +
                    static_cast< ObservationScalarType >( rampRates_.at( currentRamp ) ) *
                            static_cast< ObservationScalarType >( nextRampStartTime - startTimes_.at( currentRamp ) );
            integral += timeDelta * ( frequencyAtCurrentTime + frequencyAtNextRampStartTime ) / 2.0;

            currentTime = nextRampStartTime;
            if( currentRamp + 1 < static_cast< int >( startTimes_.size( ) ) &&
                currentTime == startTimes_.at( currentRamp + 1 ) )
            {
                currentRamp++;
            }
        }

        return integral;
    }

    //! Function to retrieve ramp start times
    std::vector< Time > getStartTimes( )
    {
        return startTimes_;
    }

    //! Function to retrieve ramp start time
    Time getStartTime( )
    {
        return startTimes_.front( );
    }

    //! Function to retrieve ramp end times
    std::vector< Time > getEndTimes( )
    {
        return endTimes_;
    }

    //! Function to retrieve ramp end time
    Time getEndTime( )
    {
        return endTimes_.back( );
    }

    //! Function to retrieve the ramp rates
    std::vector< double > getRampRates( )
    {
        return rampRates_;
    }

    //! Function to retrieve the ramp start frequencies
    std::vector< double > getStartFrequencies( )
    {
        return startFrequencies_;
    }

    void addFrequencyInterpolator( const std::shared_ptr< PiecewiseLinearFrequencyInterpolator > rampsToAdd );

private:
    void handleGap( const std::string& gapMessage )
    {
        switch( gapHandling_ )
        {
            case extrapolate_at_gaps:
                break;
            case throw_exception_at_gaps:
                throw std::runtime_error( gapMessage );
            case print_error_at_gaps:
                std::cerr << gapMessage << std::endl;
                break;
            case print_error_once_at_gaps:
                if( !hasPrintedGapWarning_ )
                {
                    std::cerr << gapMessage << std::endl;
                    hasPrintedGapWarning_ = true;
                }
                break;
            default:
                throw std::runtime_error( "Error when handling transmitting frequency gap: unknown gap handling option." );
        }
    }

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

    //! Option for handling queries in gaps between ramp intervals
    FrequencyGapHandling gapHandling_;

    //! Boolean denoting whether a gap warning has already been printed
    bool hasPrintedGapWarning_;

    //! Lookup scheme to find the nearest ramp start time for a given time
    std::shared_ptr< interpolators::LookUpScheme< Time > > startTimeLookupScheme_;

    //! Lookup scheme to find the nearest start time of the blocks without frequency transmission
    std::shared_ptr< interpolators::LookUpScheme< Time > > invalidStartTimeLookupScheme_;
};

}  // namespace ground_stations

}  // namespace tudat

#endif  // TUDAT_TRANSMITTINGFREQUENCIES_H
