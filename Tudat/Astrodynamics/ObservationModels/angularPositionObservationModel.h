#ifndef ANGULARPOSITIONOBSERVATIONMODEL_H
#define ANGULARPOSITIONOBSERVATIONMODEL_H

#include <map>
#include <string>

#include <Eigen/Core>

#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/Ephemerides/createLinkEndEphemeris.h"
#include "Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Astrodynamics/ObservationModels/observationModel.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
class AngularPositionObservationModel: public ObservationModel< 2, ObservationScalarType, TimeType, StateScalarType >
{
public:
    AngularPositionObservationModel( const boost::shared_ptr< observation_models::LightTimeCalculator< > > lightTimeCalculator ):
        ObservationModel< 2, ObservationScalarType, TimeType, StateScalarType >( ),
        lightTimeCalculator_( lightTimeCalculator ) { }

    AngularPositionObservationModel( const std::pair< std::string, std::string >& observationReferencePoint,
                                     const std::pair< std::string, std::string >& transmittingSystem,
                                     const std::map< std::string, boost::shared_ptr< bodies::Body > >& bodyMap,
                                     const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections =
                                     std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( )  ):
        ObservationModel< 2, ObservationScalarType, TimeType, StateScalarType >( angular_position )
    {
        using namespace tudat::observation_models;
        using namespace tudat::ephemerides;
        using namespace tudat::bodies;

        typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;


        // Declare functions to retrieve states of transmitter and receiever.
        boost::function< StateType( const TimeType ) > transmitterCompleteEphemeris =
                ephemerides::getLinkEndCompleteEphemerisFunction< TimeType, StateScalarType >( transmittingSystem, bodyMap );
        boost::function< StateType( const TimeType ) > observationReferencePointCompleteEphemeris =
                ephemerides::getLinkEndCompleteEphemerisFunction< TimeType, StateScalarType >( observationReferencePoint, bodyMap );

        lightTimeCalculator_ = createLightTimeCalculator< ObservationScalarType, TimeType, StateScalarType >(
                            transmitterCompleteEphemeris, observationReferencePointCompleteEphemeris, bodyMap, lightTimeCorrections,
                    transmittingSystem, observationReferencePoint );
    }

    Eigen::Matrix< ObservationScalarType, 2, 1 > computeObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime ) const
    {
        std::vector< double > linkEndTimes;
        std::vector< basic_mathematics::Vector6d > linkEndStates;

        return this->computeObservationsAndLinkEndData( time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates );
    }



    Eigen::Matrix< ObservationScalarType, 2, 1 > computeObservationsAndFullPrecisionLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< TimeType >& linkEndTimes,
            std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) const

    {
        bool isTimeAtReception;
        if( linkEndAssociatedWithTime == receiver )
        {
            isTimeAtReception = 1;
        }
        else if( linkEndAssociatedWithTime == transmitter )
        {
            isTimeAtReception = 0;
        }
        else
        {
            std::cerr<<"Error when calculating angular position observation, link end is not transmitter or receiver"<<std::endl;
        }

        Eigen::Matrix< StateScalarType, 6, 1 > receiverState;
        Eigen::Matrix< StateScalarType, 6, 1 > transmitterState;

        ObservationScalarType lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                    receiverState, transmitterState, time, isTimeAtReception );

        Eigen::Matrix< StateScalarType, 3, 1 > vectorToTransmitter = transmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 );
        Eigen::Matrix< ObservationScalarType, 3, 1 > sphericalRelativeCoordinates =
                coordinate_conversions::convertCartesianToSphericalTemplate< StateScalarType >( vectorToTransmitter ).
                template cast< ObservationScalarType >( );

        linkEndTimes.clear( );
        linkEndStates.clear( );
        linkEndStates.push_back( transmitterState );
        linkEndStates.push_back( receiverState );

        if( isTimeAtReception )
        {
            linkEndTimes.push_back( time - lightTime );
            linkEndTimes.push_back( time );
        }
        else
        {
            linkEndTimes.push_back( time );
            linkEndTimes.push_back( time + lightTime );
        }

        return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << sphericalRelativeCoordinates.z( ),
                 mathematical_constants::PI / 2.0 - sphericalRelativeCoordinates.y( ) ).finished( );
    }

    boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > > getLightTimeCalculator( )
    {
        return lightTimeCalculator_;
    }

private:

    //! Object to calculate light time.
    /*!
     *  Object to calculate light time, including possible corrections from troposphere, relativistic corrections, etc.
     */
    boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > > lightTimeCalculator_;
};

}

}

#endif // ANGULARPOSITIONOBSERVATIONMODEL_H
