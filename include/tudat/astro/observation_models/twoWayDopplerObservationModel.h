/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TWOWAYDOPPLEROBSERVATIONMODEL_H
#define TUDAT_TWOWAYDOPPLEROBSERVATIONMODEL_H

#include "tudat/astro/observation_models/oneWayDopplerObservationModel.h"

namespace tudat
{

namespace observation_models
{


//! Computes observable the (simplified) one-way Doppler observation between two link ends, omitting proper time rates and
//! light time corrections.
/*!
 *  Computes observable the (simplified) one-way Doppler observation between two link ends, omitting proper time rates and
 *  light time corrections. The observable is defined as d f_{B}/d_f_{A} - 1, with A the transmitter and B the receiver, and
 *  f the frequency of the signal.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class TwoWayDopplerObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;


    //! Constructor.
    /*!
     *  Constructor,
     *  \param uplinkDopplerCalculator Object that computes the one-way Doppler observable for the uplink
     *  \param downlinkDopplerCalculator Object that computes the one-way Doppler observable for the downlink
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    TwoWayDopplerObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< observation_models::MultiLegLightTimeCalculator< ObservationScalarType, TimeType > >
            multiLegLightTimeCalculator,
            const std::shared_ptr< observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > >
            uplinkDopplerCalculator,
            const std::shared_ptr< observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > >
            downlinkDopplerCalculator,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const bool normalizeWithSpeedOfLight = false ):
        ObservationModel< 1, ObservationScalarType, TimeType >( two_way_doppler, linkEnds, observationBiasCalculator ),
        multiLegLightTimeCalculator_( multiLegLightTimeCalculator ),
        uplinkDopplerCalculator_( uplinkDopplerCalculator ),
        downlinkDopplerCalculator_( downlinkDopplerCalculator )
    {
        setNormalizeWithSpeedOfLight( normalizeWithSpeedOfLight );
        uplinkDopplerCalculator_->setNormalizeWithSpeedOfLight( true );
        downlinkDopplerCalculator_->setNormalizeWithSpeedOfLight( true );
    }


    //! Destructor
    ~TwoWayDopplerObservationModel( ){ }

    //! Function to compute one-way Doppler observable without any corrections.
    /*!
     *  Function to compute one-way Doppler  observable without any corrections, i.e. the true physical Doppler as computed
     *  from the defined link ends. It does not include system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in double precision. These states and times are returned by
     *  reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal one-way Doppler ObservationAncilliarySimulationSettings<observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySetings = nullptr )
    {
        linkEndTimes.clear( );
        linkEndStates.clear( );
        multiLegLightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
            time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancilliarySetings );

        std::vector< double > uplinkLinkEndTimes = { linkEndTimes.at( 0 ), linkEndTimes.at( 1 ) };
        std::vector< Eigen::Matrix< double, 6, 1 > > uplinkLinkEndStates = { linkEndStates.at( 0 ), linkEndStates.at( 1 ) };

        std::vector< double > downlinkLinkEndTimes = { linkEndTimes.at( 2 ), linkEndTimes.at( 3 ) };
        std::vector< Eigen::Matrix< double, 6, 1 > > downlinkLinkEndStates = { linkEndStates.at( 2 ), linkEndStates.at( 3 ) };

        Eigen::Matrix< ObservationScalarType, 1, 1 > uplinkDoppler, downlinkDoppler;

        LinkEndType uplinkReferenceLinkEnd, downlinkReferenceLinkEnd;
        switch( linkEndAssociatedWithTime )
        {
        case receiver:
            uplinkReferenceLinkEnd = receiver;
            downlinkReferenceLinkEnd = receiver;
            break;
        case reflector1:
            uplinkReferenceLinkEnd = receiver;
            downlinkReferenceLinkEnd = transmitter;
            break;
        case transmitter:
            uplinkReferenceLinkEnd = transmitter;
            downlinkReferenceLinkEnd = transmitter;
            break;
        default:
            throw std::runtime_error(
                        "Error when calculating two way Doppler observation, link end is not transmitter or receiver" );
        }

        uplinkDoppler = uplinkDopplerCalculator_->computeIdealDopplerWithLinkEndData(
            uplinkReferenceLinkEnd, uplinkLinkEndTimes, uplinkLinkEndStates );
        downlinkDoppler = downlinkDopplerCalculator_->computeIdealDopplerWithLinkEndData(
            downlinkReferenceLinkEnd, downlinkLinkEndTimes, downlinkLinkEndStates );


        return multiplicationTerm_ * ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << downlinkDoppler( 0 ) * uplinkDoppler( 0 ) +
                 downlinkDoppler( 0 ) + uplinkDoppler( 0 ) ).finished( );
    }

    //! Function to retrieve the object that computes the one-way Doppler observable for the uplink
    /*!
     * Function to retrieve the object that computes the one-way Doppler observable for the uplink
     * \return Object that computes the one-way Doppler observable for the uplink
     */
    std::shared_ptr< observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > >
    getUplinkDopplerCalculator( )
    {
        return uplinkDopplerCalculator_;
    }

    //! Function to retrieve the object that computes the one-way Doppler observable for the downlink
    /*!
     * Function to retrieve the object that computes the one-way Doppler observable for the downlink
     * \return Object that computes the one-way Doppler observable for the downlink
     */
    std::shared_ptr< observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > >
    getDownlinkDopplerCalculator( )
    {
        return downlinkDopplerCalculator_;
    }

    void setNormalizeWithSpeedOfLight( const bool normalizeWithSpeedOfLight )
    {
        normalizeWithSpeedOfLight_ = normalizeWithSpeedOfLight;
        multiplicationTerm_ = normalizeWithSpeedOfLight ? mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 ) :
                                                          physical_constants::getSpeedOfLight< ObservationScalarType >( );
    }

    bool getNormalizeWithSpeedOfLight( )
    {
        return normalizeWithSpeedOfLight_;
    }

    ObservationScalarType getMultiplicationTerm( )
    {
        return multiplicationTerm_;
    }

private:

    std::shared_ptr< observation_models::MultiLegLightTimeCalculator< ObservationScalarType, TimeType > > multiLegLightTimeCalculator_;

    //! Object that computes the one-way Doppler observable for the uplink
    std::shared_ptr< observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > > uplinkDopplerCalculator_;

    //! Object that computes the one-way Doppler observable for the downlink
    std::shared_ptr< observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > > downlinkDopplerCalculator_;

    ObservationScalarType multiplicationTerm_;

    bool normalizeWithSpeedOfLight_;
};

}

}

#endif // TUDAT_TWOWAYDOPPLEROBSERVATIONMODEL_H
