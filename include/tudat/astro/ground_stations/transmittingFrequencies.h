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
#include "tudat/io/readOdfFile.h"

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

    virtual double getCurrentFrequency( const double lookupTime ) = 0;

    virtual double getFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime ) = 0;

    virtual double getAveragedFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        return getFrequencyIntegral( quadratureStartTime, quadratureEndTime ) / ( quadratureEndTime - quadratureStartTime );
    }

private:

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

    double getCurrentFrequency( const double lookupTime )
    {
        return frequency_;
    }

    double getFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        return frequency_ * ( quadratureEndTime - quadratureStartTime );
    }

    double getAveragedFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        return getFrequencyIntegral( quadratureStartTime, quadratureEndTime ) / ( quadratureEndTime - quadratureStartTime );
    }

private:

    double frequency_;
};

class PiecewiseLinearFrequencyInterpolator: public StationFrequencyInterpolator
{
public:
    PiecewiseLinearFrequencyInterpolator(
            std::vector< std::shared_ptr< input_output::OdfRampBlock > > rampBlock ):
        StationFrequencyInterpolator( )
    {
        for( unsigned int i = 0; i < rampBlock.size( ); i++ )
        {
            startTimes_.push_back( rampBlock.at( i )->getRampStartTime( ) );
            endTimes_.push_back( rampBlock.at( i )->getRampEndTime( ) );
            rampRates_.push_back( rampBlock.at( i )->getRampRate( ) );
            startFrequencies_.push_back( rampBlock.at( i )->getRampStartFrequency( ) );
        }

        startTimeLookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                startTimes_ );
    }

    PiecewiseLinearFrequencyInterpolator(
            const std::vector< double >& startTimes_,
            const std::vector< double >& endTimes_,
            const std::vector< double >& rampRates_,
            const std::vector< double >& startFrequency_ ):
        StationFrequencyInterpolator( ),
        startTimes_( startTimes_ ), endTimes_( endTimes_ ), rampRates_( rampRates_ ), startFrequencies_( startFrequency_ )
    {
        startTimeLookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                startTimes_ );
    }

    double getFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        std::vector< double > quadratureTimes;
        std::vector < double > quadratureFrequencies;

        // Point corresponding to first partial ramp
        quadratureTimes.push_back ( quadratureStartTime );
        quadratureFrequencies.push_back ( getCurrentFrequency( quadratureStartTime ) );

        // Points corresponding to full ramps
        for( unsigned int i = 1; i < startTimes_.size( ) && startTimes_.at( i ) < quadratureEndTime; i++ )
        {
            quadratureTimes.push_back( startTimes_.at( i ) );
            quadratureFrequencies.push_back( startFrequencies_.at( i ) );
        }

        // Point corresponding to final partial ramp
        quadratureTimes.push_back ( quadratureEndTime );
        quadratureFrequencies.push_back ( getCurrentFrequency( quadratureEndTime ) );

        return numerical_quadrature::performTrapezoidalQuadrature( quadratureTimes, quadratureFrequencies );
    }

    double getCurrentFrequency( const double lookupTime )
    {
        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour( lookupTime );

        if( lookupTime > endTimes_.at( lowerNearestNeighbour ) || lookupTime < startTimes_.at ( lowerNearestNeighbour ) )
        {
            throw std::runtime_error(
                    "Error when interpolating ramp reference frequency: look up time (" + std::to_string( lookupTime ) +
                    ") is outside the ramp table interval (" + std::to_string( startTimes_.at( 0 ) ) + " to " +
                    std::to_string( startTimes_.back( ) ) + ")." );
        }

        return startFrequencies_.at( lowerNearestNeighbour ) +
               rampRates_.at( lowerNearestNeighbour ) * ( lookupTime - startTimes_.at( lowerNearestNeighbour ) );
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

    std::vector< double > startTimes_;
    std::vector< double > endTimes_;
    std::vector< double > rampRates_;
    std::vector< double > startFrequencies_;

    std::shared_ptr< interpolators::LookUpScheme< double > > startTimeLookupScheme_;

};


// All time intervals are assumed to have the same size
class PiecewiseConstantFrequencyInterpolator: public StationFrequencyInterpolator
{
public:
    //! Constructor
    PiecewiseConstantFrequencyInterpolator( std::vector< double > frequencies,
                                            std::vector< double > referenceTimes,
                                            double timeIntervalsSize ):
        StationFrequencyInterpolator( ),
        frequencies_( frequencies ),
        referenceTimes_( referenceTimes ),
        timeIntervalsSize_( timeIntervalsSize )
    {
        if ( frequencies.size( ) != referenceTimes.size( ) )
        {
            throw std::runtime_error("Error when creating piecewise constant frequency interpolator: size of time stamps and "
                                     "frequencies are not consistent.");
        }

        startTimeLookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                referenceTimes_ );
    }

    //! Destructor
    ~PiecewiseConstantFrequencyInterpolator( ) { }

    double getCurrentFrequency( const double lookupTime )
    {
        unsigned int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour( lookupTime );
        unsigned int higherNearestNeighbour = lowerNearestNeighbour + 1;

        // Look-up time closer to lower nearest neighbour
        if ( lookupTime - referenceTimes_.at( lowerNearestNeighbour ) <=  referenceTimes_.at( higherNearestNeighbour ) - lookupTime ||
            lowerNearestNeighbour == referenceTimes_.size( ) - 1 )
        {
            return frequencies_.at( lowerNearestNeighbour );
        }
        // Look-up time closer to higher nearest neighbour
        else
        {
            return frequencies_.at( higherNearestNeighbour );
        }
    }

    double getFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        throw std::runtime_error("Computation of integral not implemented for piecewise constant frequency.");
    }

    double getAveragedFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
    {
        double referenceTime = quadratureStartTime + ( quadratureEndTime - quadratureStartTime ) / 2.0;

        if ( ( referenceTime - quadratureStartTime ) / ( timeIntervalsSize_ / 2.0 ) - 1.0 > 1e-12 ||
            ( quadratureEndTime - referenceTime ) / ( timeIntervalsSize_ / 2.0 ) - 1.0 > 1e-12 )
        {
            throw std::runtime_error("Error when computing the averaged integral of piecewise constant frequency: "
                                     "the specified time interval does not coincide with any piecewise interval.");
        }

        return getCurrentFrequency( referenceTime );
    }

private:

    std::vector< double > frequencies_;
    std::vector< double > referenceTimes_;

    double timeIntervalsSize_;

    std::shared_ptr< interpolators::LookUpScheme< double > > startTimeLookupScheme_;
};

} // namespace ground_stations

} // namespace tudat

#endif //TUDAT_TRANSMITTINGFREQUENCIES_H
