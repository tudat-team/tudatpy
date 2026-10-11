#ifndef TUDAT_BARYCENTER_H
#define TUDAT_BARYCENTER_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/ephemerides/customEphemeris.h"

namespace tudat
{
namespace simulation_setup
{

//! A translational point with prescribed members; deliberately not a Body.
class Barycenter
{
public:
    using State = Eigen::Matrix< long double, 6, 1 >;
    Barycenter( const std::string& name,
                const std::vector< std::string >& members,
                const std::vector< double >& weights,
                const std::vector< std::shared_ptr< ephemerides::Ephemeris > >& nominalEphemerides,
                const std::string& orientation,
                const double offsetAccelerationStep = 10.0 ):
        name_( name ), members_( members ), weights_( weights ), nominalEphemerides_( nominalEphemerides ), orientation_( orientation ),
        offsetAccelerationStep_( offsetAccelerationStep )
    {
        if( members.empty( ) || members.size( ) != weights.size( ) || members.size( ) != nominalEphemerides.size( ) )
        {
            throw std::runtime_error( "Barycenter members, weights and ephemerides must have the same nonzero size." );
        }
        double totalWeight = 0.0;
        for( unsigned int i = 0; i < weights_.size( ); ++i )
        {
            if( !std::isfinite( weights_[ i ] ) || weights_[ i ] <= 0.0 || !nominalEphemerides_[ i ] )
            {
                throw std::runtime_error( "Barycenter weights must be positive and finite and ephemerides must exist." );
            }
            totalWeight += weights_[ i ];
        }
        if( !std::isfinite( totalWeight ) || !std::isfinite( offsetAccelerationStep ) || offsetAccelerationStep <= 0.0 )
        {
            throw std::runtime_error( "Invalid barycenter weight sum or offset acceleration differentiation step." );
        }
        for( auto& weight : weights_ )
        {
            weight /= totalWeight;
        }
        // Capture immutable inputs, not this: replacing a propagated ephemeris must not alter its nominal definition.
        const auto normalizedWeights = weights_;
        nominalEphemeris_ = std::make_shared< ephemerides::CustomEphemeris< Time, long double > >(
                [ nominalEphemerides, normalizedWeights ]( const Time& time ) -> State {
                    const State reference = nominalEphemerides.front( )->getTemplatedStateFromEphemeris< long double, Time >( time );
                    State offset = State::Zero( );
                    for( unsigned int i = 1; i < nominalEphemerides.size( ); ++i )
                    {
                        offset += static_cast< long double >( normalizedWeights[ i ] ) *
                                ( nominalEphemerides[ i ]->getTemplatedStateFromEphemeris< long double, Time >( time ) - reference );
                    }
                    return reference + offset;
                },
                "SSB",
                orientation );
        ephemeris_ = nominalEphemeris_;
    }

    std::string getName( ) const
    {
        return name_;
    }
    const std::vector< std::string >& getMembers( ) const
    {
        return members_;
    }
    const std::vector< double >& getWeights( ) const
    {
        return weights_;
    }
    std::shared_ptr< ephemerides::Ephemeris > getEphemeris( ) const
    {
        return ephemeris_;
    }
    std::shared_ptr< ephemerides::Ephemeris > getNominalEphemeris( ) const
    {
        return nominalEphemeris_;
    }
    void setEphemeris( const std::shared_ptr< ephemerides::Ephemeris >& ephemeris )
    {
        if( !ephemeris || ephemeris->getReferenceFrameOrigin( ) != "SSB" || ephemeris->getReferenceFrameOrientation( ) != orientation_ )
        {
            throw std::runtime_error( "Barycenter ephemeris must use SSB and the environment orientation." );
        }
        ephemeris_ = ephemeris;
        currentStateValid_ = false;
    }

    State getMemberOffset( const std::string& member, const Time& time ) const
    {
        auto iterator = std::find( members_.begin( ), members_.end( ), member );
        if( iterator == members_.end( ) )
        {
            throw std::runtime_error( member + " is not a member of barycenter " + name_ );
        }
        return nominalEphemerides_.at( std::distance( members_.begin( ), iterator ) )
                       ->getTemplatedStateFromEphemeris< long double, Time >( time ) -
                nominalEphemeris_->getTemplatedStateFromEphemeris< long double, Time >( time );
    }

    Eigen::Vector3d getMemberOffsetAcceleration( const std::string& member, const double time ) const
    {
        const long double step = offsetAccelerationStep_;
        const Time epoch( time );
        const State derivative = ( -getMemberOffset( member, epoch + 2.0L * step ) + 8.0L * getMemberOffset( member, epoch + step ) -
                                   8.0L * getMemberOffset( member, epoch - step ) + getMemberOffset( member, epoch - 2.0L * step ) ) /
                ( 12.0L * step );
        return derivative.tail< 3 >( ).cast< double >( );
    }

    void beginStateUpdate( const Time& time )
    {
        currentTime_ = time;
        currentStateValid_ = false;
    }
    void clearCurrentState( )
    {
        currentStateValid_ = false;
    }
    template< typename Scalar >
    void setCurrentState( const Eigen::Matrix< Scalar, 6, 1 >& state )
    {
        currentState_ = state.template cast< long double >( );
        currentStateValid_ = true;
    }
    bool hasCurrentState( ) const
    {
        return currentStateValid_;
    }
    State getCurrentLongState( ) const
    {
        if( !currentStateValid_ )
        {
            throw std::runtime_error( "Current state of barycenter " + name_ + " has not been set." );
        }
        return currentState_;
    }
    Eigen::Vector6d getState( ) const
    {
        return getCurrentLongState( ).cast< double >( );
    }
    State getStateForMember( const Time& time ) const
    {
        if( currentStateValid_ && time == currentTime_ )
        {
            return currentState_;
        }
        return ephemeris_->getTemplatedStateFromEphemeris< long double, Time >( time );
    }

private:
    std::string name_;
    std::vector< std::string > members_;
    std::vector< double > weights_;
    std::vector< std::shared_ptr< ephemerides::Ephemeris > > nominalEphemerides_;
    std::string orientation_;
    double offsetAccelerationStep_;
    std::shared_ptr< ephemerides::Ephemeris > nominalEphemeris_, ephemeris_;
    State currentState_ = State::Zero( );
    Time currentTime_ = Time( TUDAT_NAN );
    bool currentStateValid_ = false;
};

}  // namespace simulation_setup
}  // namespace tudat
#endif
