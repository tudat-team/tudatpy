/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRANSMITTINGFREQUENCIES_H
#define TUDAT_TRANSMITTINGFREQUENCIES_H

#include "tudat/math/quadrature/trapezoidQuadrature.h"
#include "tudat/math/interpolators.h"

namespace tudat
{

namespace ground_stations
{

class StationFrequencyInterpolator
{
public:
    //! Constructor
    StationFrequencyInterpolator( ) { }

    //! Destructor
    virtual ~StationFrequencyInterpolator( ) { }

    template< typename ObservationScalarType = double, typename TimeType = double >
    ObservationScalarType getTemplatedCurrentFrequency( const TimeType& lookupTime );

    template< typename ObservationScalarType = double, typename TimeType = double >
    ObservationScalarType getTemplatedFrequencyIntegral( const TimeType& quadratureStartTime, const TimeType& quadratureEndTime );

private:

    virtual double getCurrentFrequency( const double lookupTime ) = 0;

    virtual double getCurrentFrequency( const Time& lookupTime ) = 0;

    virtual long double getCurrentLongFrequency( const double lookupTime ) = 0;

    virtual long double getCurrentLongFrequency( const Time& lookupTime ) = 0;

    virtual double getFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime ) = 0;

    virtual double getFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime ) = 0;

    virtual long double getLongFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime ) = 0;

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

class PiecewiseLinearFrequencyInterpolator: public StationFrequencyInterpolator
{
public:

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

    template< typename ObservationScalarType = double, typename TimeType = double >
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
//            integral += ( quadratureEndTime - quadratureStartTime ) *
//                    ( computeCurrentFrequency< ObservationScalarType, TimeType >( quadratureStartTime ) + rampRates_.at( startTimeLowestNearestNeighbour ) *
//                    ( quadratureEndTime - quadratureStartTime ) / 2.0 );
            std::cout << integral << std::endl << quadratureStartTime << std::endl << quadratureEndTime << std::endl << computeCurrentFrequency< ObservationScalarType, TimeType >( quadratureStartTime )
                << std::endl << computeCurrentFrequency< ObservationScalarType, TimeType >( quadratureEndTime ) << std::endl <<
                static_cast< ObservationScalarType > ( quadratureEndTime - quadratureStartTime ) << std::endl <<
                ( computeCurrentFrequency< ObservationScalarType, TimeType >( quadratureStartTime ) + computeCurrentFrequency< ObservationScalarType, TimeType >( quadratureEndTime ) ) / 2.0 << std::endl << std::endl;
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

    std::vector< double > getStartTimes ( )
    {
        return startTimes_;
    }

    std::vector< double > getEndTimes ( )
    {
        return endTimes_;
    }

    std::vector< double > getRampRates ( )
    {
        return rampRates_;
    }

    std::vector< double > getStartFrequencies ( )
    {
        return startFrequencies_;
    }

private:

    virtual double getCurrentFrequency( const double lookupTime )
    {
        return computeCurrentFrequency< double, double >( lookupTime );
    }

    virtual double getCurrentFrequency( const Time& lookupTime )
    {
        return computeCurrentFrequency< double, Time >( lookupTime );
    }

    virtual long double getCurrentLongFrequency( const double lookupTime )
    {
         return computeCurrentFrequency< long double, double >( lookupTime );
    }

    virtual long double getCurrentLongFrequency( const Time& lookupTime )
    {
         return computeCurrentFrequency< long double, Time >( lookupTime );
    }

        virtual double getFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        return computeFrequencyIntegral< double, double >( quadratureStartTime, quadratureEndTime );
    }

    virtual double getFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime )
    {
        return computeFrequencyIntegral< double, Time >( quadratureStartTime, quadratureEndTime );
    }

    virtual long double getLongFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        return computeFrequencyIntegral< long double, double >( quadratureStartTime, quadratureEndTime );
    }

    virtual long double getLongFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime )
    {
        return computeFrequencyIntegral< long double, Time >( quadratureStartTime, quadratureEndTime );
    }

    std::vector< double > startTimes_;
    std::vector< double > endTimes_;
    std::vector< double > rampRates_;
    std::vector< double > startFrequencies_;

    std::shared_ptr< interpolators::LookUpScheme< double > > startTimeLookupScheme_;

};

} // namespace ground_stations

} // namespace tudat

#endif //TUDAT_TRANSMITTINGFREQUENCIES_H
