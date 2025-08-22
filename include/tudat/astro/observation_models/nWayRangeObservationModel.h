/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NWAYRANGEOBSERVATIONMODEL_H
#define TUDAT_NWAYRANGEOBSERVATIONMODEL_H

#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
namespace tudat
{

namespace observation_models
{

//! Class for simulating n-way range observables.
/*!
 *  Class for simulating n-way range observations. It connects an arbitrary number of links by rane observables to obtain the
 *  n-way range. In most practical cases (e.g. DSN radio ranging, SLR), n will be equal to 2. Note that here the number of 'ways'
 *  represents the number of legs that the signal travels along. As a result, for both DSN 2-way and 3-way observables, n equals
 *  2 in this model. The difference is that the first and last link end will be the same for the former case, and different for
 *  the latter case. The retransmission of a signal at the intermediate link ends can be done with a (negative or positive) delay.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class NWayRangeObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;

    //! Constructor.
    /*!
     *  Constructor,
     *  \param lightTimeCalculators List of objects to compute the light-times (including any corrections w.r.t. Euclidean case)
     *  for each leg of the n-way range. First entry starts at transmitter; last entry is to receiver.
     *  \param observationBiasCalculator Object for calculating (system-dependent) errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    NWayRangeObservationModel(
            const LinkEnds& linkEnds,
            const std::vector< std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >
                    lightTimeCalculators,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                    std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModel< 1, ObservationScalarType, TimeType >( n_way_range, linkEnds, observationBiasCalculator ),
        linkEnds_( linkEnds )
    {
        multiLegLightTimeCalculator_ =
                std::make_shared< observation_models::MultiLegLightTimeCalculator< ObservationScalarType, TimeType > >(
                        lightTimeCalculators, lightTimeConvergenceCriteria );
    }

    NWayRangeObservationModel( const LinkEnds& linkEnds,
                               const std::shared_ptr< observation_models::MultiLegLightTimeCalculator< ObservationScalarType, TimeType > >
                                       multiLegLightTimeCalculator,
                               const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( n_way_range, linkEnds, observationBiasCalculator ),
        linkEnds_( linkEnds ), multiLegLightTimeCalculator_( multiLegLightTimeCalculator )
    { }

    //! Destructor
    ~NWayRangeObservationModel( ) { }

    //! Function to compute n-way range observable without any corrections.
    /*!
     *  Function to compute n-way range  observable without any corrections, i.e. the true physical range as computed
     *  from the defined link ends. The time argument can be at any of the link ends
     *  involved in the onbservation (including the intermediate link ends) by the linkEndAssociatedWithTime input).
     *  In the case where the reference link end is an intermediate link end, the inpit time denotes the reception time
     *  of the signal at this station (which need not be the same as its retransmission time).
     *  Note that this observable does include light-time corrections, which represent physically true corrections. It does not
     *  include e.g. system-dependent measurement errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal n-way range observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySetings = nullptr )
    {
        std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySetingsToUse;
        setFrequencyProperties( time, linkEndAssociatedWithTime, ancilliarySetings, ancilliarySetingsToUse );

        ObservationScalarType totalLightTime = multiLegLightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancilliarySetingsToUse );

        // Return total range observation.
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( )
                 << totalLightTime * physical_constants::getSpeedOfLight< ObservationScalarType >( ) )
                .finished( );
    }

    std::vector< std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > > getLightTimeCalculators( )
    {
        return multiLegLightTimeCalculator_->getLightTimeCalculators( );
    }

    std::shared_ptr< MultiLegLightTimeCalculator< ObservationScalarType, TimeType > > getMultiLegLightTimeCalculator( )
    {
        return multiLegLightTimeCalculator_;
    }

    void setFrequencyInterpolatorAndTurnaroundRatio(
        std::shared_ptr< ground_stations::StationFrequencyInterpolator > frequencyInterpolator,
        std::function< double( FrequencyBands uplinkBand, FrequencyBands downlinkBand ) > turnaroundRatio )
    {
        frequencyInterpolator_ = frequencyInterpolator;
        turnaroundRatio_ = turnaroundRatio;
        timeScaleConverter_ = earth_orientation::createDefaultTimeConverter( );
    }
private:

    bool setFrequencyProperties( const TimeType time,
                                 const LinkEndType linkEndAssociatedWithTime,
                                 const std::shared_ptr< ObservationAncilliarySimulationSettings > inputAncilliarySetings,
                                 std::shared_ptr< ObservationAncilliarySimulationSettings >& ancilliarySetingsToUse )
    {
        if( frequencyInterpolator_ != nullptr )
        {
            if( linkEndAssociatedWithTime != receiver )
            {
                throw std::runtime_error(
                    "Error when computing n-way range, frequency interpolator use is only compatible with transmitter reference "
                    "frequency at present" );
            }
            else
            {
                if( inputAncilliarySetings == nullptr )
                {
                    ancilliarySetingsToUse = std::make_shared< ObservationAncilliarySimulationSettings >( );
                }
                else
                {
                    ancilliarySetingsToUse = inputAncilliarySetings;
                }


                setTransmissionReceptionFrequencies( multiLegLightTimeCalculator_,
                                          timeScaleConverter_,
                                          frequencyInterpolator_,
                                          time,
                                          linkEndAssociatedWithTime,
                                          ancilliarySetingsToUse,
                                          getTurnaroundRatio( ancilliarySetingsToUse ) );
            }
            return true;
        }
        else
        {
            ancilliarySetingsToUse = inputAncilliarySetings;
            return false;
        }
    }

    ObservationScalarType getTurnaroundRatio( const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings )
    {
        std::vector< FrequencyBands > frequencyBands;
        FrequencyBands referenceUplinkBand;
        try
        {
            frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error( "Error when retrieving ancillary settings for N-way range observable: " +
                                      std::string( caughtException.what( ) ) );
        }

        if( frequencyBands.size( ) != linkEnds_.size( ) - 1 )
        {
            throw std::runtime_error(
                "Error when retrieving frequency bands ancillary settings for N-way range observable: "
                "size (" +
                std::to_string( frequencyBands.size( ) ) + ") is inconsistent with number of links (" +
                std::to_string( linkEnds_.size( ) - 1 ) + ")." );
        }
        FrequencyBands uplinkBand = frequencyBands.at( 0 );
        FrequencyBands downlinkBand = frequencyBands.at( 1 );

        // Set approximate up- and down-link frequencies.
        return static_cast< ObservationScalarType >( turnaroundRatio_( uplinkBand, downlinkBand ) );
    }

    LinkEnds linkEnds_;

    // Object that iteratively computes the light time of multiple legs
    std::shared_ptr< MultiLegLightTimeCalculator< ObservationScalarType, TimeType > > multiLegLightTimeCalculator_;

    std::shared_ptr< ground_stations::StationFrequencyInterpolator > frequencyInterpolator_;

    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter_;

    std::function< double( FrequencyBands uplinkBand, FrequencyBands downlinkBand ) > turnaroundRatio_;

};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_NWAYRANGEOBSERVATIONMODEL_H
