#ifndef NWAYRANGEOBSERVATIONMODEL_H
#define NWAYRANGEOBSERVATIONMODEL_H

#include <map>

#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/Bodies/celestialBody.h"
#include "Astrodynamics/Ephemerides/createLinkEndEphemeris.h"
#include "Astrodynamics/ObservationModels/observationModel.h"
#include "Astrodynamics/ObservationModels/lightTimeSolution.h"

namespace tudat
{

namespace observation_models
{

//! Class to calculate the n-way range observable, i.e. with n-1 reflections or retransmissions
/*!
 *  Class to calculate the n-way range observable, i.e. with n-1 reflections or retransmissions, with n+1 link ends. The class creates
 *  a 1-way range observable for each constituent range segment. The range calculation may be initiated from any of the n+1 link ends
 */
template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
class NWayRangeObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType, StateScalarType >
{
public:

    //! Typedefs for state type.
    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    //! Constructor.
    /*!
     *  Constructor, creates the constituent 1-way range models, including light time corrections. A global (i.e. for the full n-way range)
     *  observation bias can be set.
     */
    NWayRangeObservationModel( const LinkEnds& linkEndsMap,
                               const std::map< std::string, boost::shared_ptr< bodies::Body > >& bodyMap,
                               const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
                               const boost::function< std::vector< ObservationScalarType >( ) > retransmissionDelaysFunction = NULL,
                               const boost::shared_ptr< ObservationBiasInterface > observationBiasCalculator =
            boost::make_shared< ObservationBiasInterface >( 1 ) ):
        ObservationModel< 1, ObservationScalarType, TimeType, StateScalarType >( nWayRange, observationBiasCalculator ),
        retransmissionDelaysFunction_( retransmissionDelaysFunction )
    {
        numberOfLinkEnds_ = linkEndsMap.size( );

        std::vector< boost::function< StateType( const TimeType ) > > linkEndCompleteEphemerides;
        linkEndCompleteEphemerides.resize( numberOfLinkEnds_ );

        int currentLinkEndIndex;
        for( LinkEnds::const_iterator linkEndIterator = linkEndsMap.begin( ); linkEndIterator != linkEndsMap.end( );
             linkEndIterator++ )
        {
            currentLinkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndIterator->first, numberOfLinkEnds_ );
            linkEndCompleteEphemerides[ currentLinkEndIndex ] =
                    ephemerides::getLinkEndCompleteEphemerisFunction< TimeType, StateScalarType >(
                        linkEndIterator->second, bodyMap );
        }

        std::vector< std::vector< LightTimeCorrectionFunction > > linkCorrectionFunctions;
        linkCorrectionFunctions.resize( numberOfLinkEnds_ - 1 );
        LinkEndType currentTransmitter, currentReceiver;

        for( int i = 0; i < numberOfLinkEnds_ - 1; i++ )
        {
            currentTransmitter = getNWayLinkEnumFromIndex( i, numberOfLinkEnds_ );
            currentReceiver = getNWayLinkEnumFromIndex( i + 1, numberOfLinkEnds_ );

            lightTimeCalculators_.push_back( createLightTimeCalculator< ObservationScalarType, TimeType, StateScalarType >(
                                                 linkEndCompleteEphemerides[ i ], linkEndCompleteEphemerides[ i + 1 ], bodyMap, lightTimeCorrections,
                                                 linkEndsMap.at( currentTransmitter ), linkEndsMap.at( currentReceiver ) ) );
        }

        if( retransmissionDelaysFunction_== NULL )
        {
            retransmissionDelaysFunction_ = boost::lambda::constant(
                        utilities::createConstantValueVector( mathematical_constants::getZero< ObservationScalarType >( ), numberOfLinkEnds_ - 2 ) );
        }
    }

    ~NWayRangeObservationModel( ){ }

    //! Function to compute two-way range observation at given time.
    /*!
             *  This function computes the two-way observation at a given time.
             *  Currently, the time argument can be either the reception or transmission time.
             *  \param time Time at which observation is to be simulated
             *  \param isTimAtRececption True if given time is to be the reception time, false if it is transmission time.
             *  \return Calculated observed two-way range value.
             */
    ObservationScalarType computeObservation( const TimeType time,
                                              const LinkEndType linkEndAssociatedWithTime ) const
    {
        std::vector< double > linkEndTimes;
        std::vector< basic_mathematics::Vector6d > linkEndStates;
        return this->computeObservationAndLinkEndData( time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );
    }

    ObservationScalarType computeObservationAndFullPrecisionLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< TimeType >& linkEndTimes,
            std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) const
    {
        int startLinkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndAssociatedWithTime, numberOfLinkEnds_ );

        ObservationScalarType totalLightTime =
                mathematical_constants::getZero< ObservationScalarType >( );
        ObservationScalarType currentLightTime;

        StateType currentReceiverStateOutput, currentTransmitterStateOutput;

        linkEndTimes.clear( );
        linkEndStates.clear( );

        linkEndTimes.resize( 2 * ( numberOfLinkEnds_ - 1 ) );
        linkEndStates.resize( 2 * ( numberOfLinkEnds_ - 1 ) );

        currentRetransmissionDelays_ = retransmissionDelaysFunction_( );

        if( currentRetransmissionDelays_.size( ) != static_cast< unsigned int >( numberOfLinkEnds_ - 2 ) )
        {
            std::cerr<<"Error when calculating n-way range, retransmission delay vector size is inconsistent"<<std::endl;
        }

        int currentDownIndex = startLinkEndIndex;

        TimeType currentLinkEndStartTime = time;
        while( currentDownIndex > 0 )
        {
            currentLightTime = lightTimeCalculators_.at( currentDownIndex - 1 )->calculateLightTimeWithLinkEndsStates(
                        currentReceiverStateOutput, currentTransmitterStateOutput,
                        currentLinkEndStartTime, 1 );

            linkEndStates[ 2 * ( currentDownIndex - 1 ) + 1 ] = currentReceiverStateOutput;
            linkEndStates[ 2 * ( currentDownIndex - 1 ) ] = currentTransmitterStateOutput;

            linkEndTimes[ 2 * ( currentDownIndex - 1 ) + 1 ] = currentLinkEndStartTime;
            linkEndTimes[ 2 * ( currentDownIndex - 1 )] = currentLinkEndStartTime - currentLightTime;

            if( currentDownIndex > 1 )
            {
                currentLightTime += currentRetransmissionDelays_.at( currentDownIndex - 2 );
            }

            currentLinkEndStartTime -= currentLightTime;

            totalLightTime += currentLightTime;
            currentDownIndex--;
        }

        int currentUpIndex = startLinkEndIndex;
        if( ( startLinkEndIndex != 0 ) && ( startLinkEndIndex != numberOfLinkEnds_ - 1 ) )
        {
            currentLinkEndStartTime = time + currentRetransmissionDelays_.at( startLinkEndIndex - 1 );
            totalLightTime += currentRetransmissionDelays_.at( startLinkEndIndex - 1 );
        }

        while( currentUpIndex < static_cast< int >( lightTimeCalculators_.size( ) ) )
        {
            currentLightTime = lightTimeCalculators_.at( currentUpIndex )->calculateLightTimeWithLinkEndsStates(
                        currentReceiverStateOutput, currentTransmitterStateOutput,
                        currentLinkEndStartTime, 0 );


            linkEndStates[ 2 * currentUpIndex + 1 ] = currentReceiverStateOutput;
            linkEndStates[ 2 * currentUpIndex ] = currentTransmitterStateOutput;

            linkEndTimes[ 2 * currentUpIndex + 1 ] = currentLinkEndStartTime + currentLightTime;
            linkEndTimes[ 2 * currentUpIndex ] = currentLinkEndStartTime;

            if( currentUpIndex < static_cast< int >( lightTimeCalculators_.size( ) ) - 1 )
            {
                currentLightTime += currentRetransmissionDelays_.at( currentUpIndex );
            }

            currentLinkEndStartTime += currentLightTime;

            totalLightTime += currentLightTime;
            currentUpIndex++;
        }

        return totalLightTime * physical_constants::getSpeedOfLight< ObservationScalarType >( );
    }

    void resetRetransmissionDelaysFunction( const boost::function< std::vector< ObservationScalarType >( ) > retransmissionDelaysFunction )
    {
        retransmissionDelaysFunction_ = retransmissionDelaysFunction;
    }

    std::vector< boost::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > > > getLightTimeCalculators( )
    {
        return lightTimeCalculators_;
    }


private:

    std::vector< boost::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > > > lightTimeCalculators_;

    boost::function< std::vector< ObservationScalarType >( ) > retransmissionDelaysFunction_;

    mutable std::vector< ObservationScalarType > currentRetransmissionDelays_;

    int numberOfLinkEnds_;
};

}

}

#endif // NWAYRANGEOBSERVATIONMODEL_H
