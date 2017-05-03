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

class MultiArcEphemeris: public Ephemeris
{
public:
    MultiArcEphemeris(
            const std::map< std::pair< double, double >, boost::shared_ptr< Ephemeris > >& singleArcEphemerides,
            const std::string& referenceFrameOrigin = "",
            const std::string& referenceFrameOrientation = "" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation ),
        singleArcEphemerides_( utilities::createVectorFromMapValues( singleArcEphemerides ) ),
        arcStartAndEndTimes_( utilities::createVectorFromMapKeys( singleArcEphemerides ) )
    {
        arcSplitTimes_.clear( );
        arcSplitTimes_.push_back( arcStartAndEndTimes_[ 0 ].first );
        for( unsigned int i = 1; i < arcStartAndEndTimes_.size( ); i++ )
        {
            arcSplitTimes_.push_back( arcStartAndEndTimes_.at( i ).first - 0.1 );
        }
        arcSplitTimes_.push_back( arcStartAndEndTimes_.at( arcStartAndEndTimes_.size( ) - 1 ).second );

        lookUpscheme_ = boost::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                    arcSplitTimes_ );
    }

    ~MultiArcEphemeris( ){ }

    Eigen::Vector6d getCartesianState(
            const double secondsSinceEpoch )
    {
        return singleArcEphemerides_.at( lookUpscheme_->findNearestLowerNeighbour( secondsSinceEpoch ) )->
                getCartesianState( double( secondsSinceEpoch ) );
    }

    Eigen::Matrix< long double, 6, 1 > getCartesianLongState(
            const double secondsSinceEpoch )
    {
        return singleArcEphemerides_.at( lookUpscheme_->findNearestLowerNeighbour( secondsSinceEpoch ) )->
                getCartesianLongState( secondsSinceEpoch );
    }

    Eigen::Vector6d getCartesianStateFromExtendedTime(
            const Time& currentTime )
    {
        return singleArcEphemerides_.at( lookUpscheme_->findNearestLowerNeighbour( currentTime ) )->
                getCartesianStateFromExtendedTime( currentTime );
    }

    Eigen::Matrix< long double, 6, 1 > getCartesianLongStateFromExtendedTime(
            const Time& currentTime )
    {
        return singleArcEphemerides_.at( lookUpscheme_->findNearestLowerNeighbour( currentTime ) )->
                getCartesianLongStateFromExtendedTime( currentTime );
    }

    void resetSingleArcEphemerides(
            const std::vector< boost::shared_ptr< Ephemeris > >& singleArcEphemerides,
            const std::vector< std::pair< double, double > >& arcStartAndEndTimes )
    {
        singleArcEphemerides_ = singleArcEphemerides;
        arcStartAndEndTimes_ = arcStartAndEndTimes;

        arcSplitTimes_.clear( );
        arcSplitTimes_.push_back( arcStartAndEndTimes_[ 0 ].first );
        for( unsigned int i = 1; i < arcStartAndEndTimes_.size( ); i++ )
        {
            arcSplitTimes_.push_back( arcStartAndEndTimes_.at( i ).first - 0.1 );
        }
        arcSplitTimes_.push_back( arcStartAndEndTimes_.at( arcStartAndEndTimes_.size( ) - 1 ).second );
        lookUpscheme_ = boost::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
                    arcSplitTimes_ );
    }

    void resetSingleArcEphemerides(
            const std::map< std::pair< double, double >, boost::shared_ptr< Ephemeris > >& singleArcEphemerides )
    {
        resetSingleArcEphemerides( utilities::createVectorFromMapValues( singleArcEphemerides ),
                                   utilities::createVectorFromMapKeys( singleArcEphemerides ) );
    }

    std::vector< double > getArcSplitTimes( )
    {
        return arcSplitTimes_;
    }

    std::vector< boost::shared_ptr< Ephemeris > > getSingleArcEphemerides( )
    {
        return singleArcEphemerides_;
    }


private:
    std::vector< boost::shared_ptr< Ephemeris > > singleArcEphemerides_;

    std::vector< std::pair< double, double > > arcStartAndEndTimes_;

    std::vector< double > arcSplitTimes_;

    boost::shared_ptr< interpolators::HuntingAlgorithmLookupScheme< double > > lookUpscheme_;


};

}

}

#endif // TUDAT_MULTIARCEPHEMERIS_H
