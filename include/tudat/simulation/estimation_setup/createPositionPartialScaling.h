/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_POSITIONPARTIALSCALING_H
#define TUDAT_POSITIONPARTIALSCALING_H

#include <vector>
#include <map>

#include "tudat/astro/orbit_determination/observation_partials/directObservationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/angularPositionPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayRangePartial.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayDopplerPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/nWayRangePartial.h"
#include "tudat/astro/orbit_determination/observation_partials/differencedObservationPartial.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace observation_partials
{
//! Function to create an object that computes the scaling of the state partials to obtain proper time rate partials
/*!
 * Function to create an object that computes the scaling of the state partials to obtain proper time rate partials. A single
 * scaling object is used for a single link end of the one-way Doppler partials
 * \param dopplerProperTimeInterface Object that is used to computed proper-time rate in one-way Doppler modelkkl
 * \param oneWayDopplerLinkEnds Link ends of observable
 * \param linkEndAtWhichPartialIsComputed Link end for which proper-time partials are to be created
 * \return Scaling object for proper-time rate partials
 */
inline std::shared_ptr< OneWayDopplerProperTimeComponentScaling > createDopplerProperTimePartials(
        const std::shared_ptr< observation_models::DopplerProperTimeRateInterface > dopplerProperTimeInterface,
        const observation_models::LinkEnds oneWayDopplerLinkEnds,
        const observation_models::LinkEndType linkEndAtWhichPartialIsComputed )
{
    std::shared_ptr< OneWayDopplerProperTimeComponentScaling >  properTimeRateDopplerPartial = nullptr;
    if( dopplerProperTimeInterface == nullptr )
    {
        properTimeRateDopplerPartial = nullptr;
    }
    else if( std::dynamic_pointer_cast< observation_models::DirectFirstOrderDopplerProperTimeRateInterface >(
                 dopplerProperTimeInterface ) != nullptr )
    {
        properTimeRateDopplerPartial = std::make_shared< OneWayDopplerDirectFirstOrderProperTimeComponentScaling >(
                    std::dynamic_pointer_cast< observation_models::DirectFirstOrderDopplerProperTimeRateInterface >(
                        dopplerProperTimeInterface ), linkEndAtWhichPartialIsComputed, oneWayDopplerLinkEnds, true );
    }
    else
    {
        std::cerr << "Warning, proper time contribution to Doppler observable not incorporated into Doppler partial " << std::endl;
        properTimeRateDopplerPartial = nullptr;
    }
    return properTimeRateDopplerPartial;

}

template< int ObservationSize >
class ObservationPartialScalingCreator
{
public:

    template< typename ParameterType = double, typename TimeType = double >
    static std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > createPositionScalingObject(
            const observation_models::LinkEnds& linkEnds,
            const observation_models::ObservableType observableType,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > observationModel = nullptr );

    template< typename ParameterType = double, typename TimeType = double >
    static std::shared_ptr< PositionPartialScaling > createDifferencedPositionPartialScalingObject(
            const observation_models::ObservableType differencedObservableType,
            const std::shared_ptr< PositionPartialScaling > firstPositionPartialScaling,
            const std::shared_ptr< PositionPartialScaling > secondPositionPartialScaling,
            const simulation_setup::SystemOfBodies& bodies );
};

template< >
class ObservationPartialScalingCreator< 1 >
{
public:

    template< typename ParameterType = double, typename TimeType = double >
    static std::shared_ptr< DirectPositionPartialScaling< 1 > > createPositionScalingObject(
            const observation_models::LinkEnds& linkEnds,
            const observation_models::ObservableType observableType,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< observation_models::ObservationModel< 1, ParameterType, TimeType > > observationModel = nullptr  )
    {
        std::shared_ptr< DirectPositionPartialScaling< 1 > > positionPartialScaler;

        switch( observableType )
        {
        case observation_models::one_way_range:
            positionPartialScaler = std::make_shared< OneWayRangeScaling >( );
            break;
        case observation_models::one_way_doppler:
        {
            std::shared_ptr< observation_models::OneWayDopplerObservationModel< ParameterType, TimeType > > dopplerObservationModel =
                    std::dynamic_pointer_cast< observation_models::OneWayDopplerObservationModel< ParameterType, TimeType > >(
                        observationModel );
            if( dopplerObservationModel == nullptr )
            {
                throw std::runtime_error( "Error when creating one-way Doppler partial scaling object, input observation model is incompatible" );
            }
            // Create scaling object, to be used for all one-way doppler partials in current link end.
            std::shared_ptr< OneWayDopplerProperTimeComponentScaling > transmitterProperTimePartials =
                    createDopplerProperTimePartials( dopplerObservationModel->getTransmitterProperTimeRateCalculator( ), linkEnds,
                                                     observation_models::transmitter );
            std::shared_ptr< OneWayDopplerProperTimeComponentScaling > receiverProperTimePartials =
                    createDopplerProperTimePartials( dopplerObservationModel->getReceiverProperTimeRateCalculator( ), linkEnds,
                                                     observation_models::receiver  );

            std::function< Eigen::Vector6d( const double )> transmitterNumericalStateDerivativeFunction =
                    std::bind( &numerical_derivatives::computeCentralDifferenceFromFunction< Eigen::Vector6d, double >,
                                 simulation_setup::getLinkEndCompleteEphemerisFunction< double, double >(
                                     linkEnds.at( observation_models::transmitter ), bodies ), std::placeholders::_1, 100.0,
                                 numerical_derivatives::order8 );
            std::function< Eigen::Vector6d( const double )> receiverNumericalStateDerivativeFunction =
                    std::bind( numerical_derivatives::computeCentralDifferenceFromFunction< Eigen::Vector6d, double >,
                                 simulation_setup::getLinkEndCompleteEphemerisFunction< double, double >(
                                     linkEnds.at( observation_models::receiver ), bodies ), std::placeholders::_1, 100.0,
                                 numerical_derivatives::order8 );

            positionPartialScaler =
                    std::make_shared< OneWayDopplerScaling >(
                        std::bind( &linear_algebra::evaluateSecondBlockInStateVector, transmitterNumericalStateDerivativeFunction, std::placeholders::_1 ),
                        std::bind( &linear_algebra::evaluateSecondBlockInStateVector, receiverNumericalStateDerivativeFunction, std::placeholders::_1 ),
                        dopplerObservationModel->getNormalizeWithSpeedOfLight( ) ? physical_constants::SPEED_OF_LIGHT : 1.0,
                        transmitterProperTimePartials,
                        receiverProperTimePartials );

            break;
        }
        default:
            throw std::runtime_error( "Error when creating partial scaler for " +
                                      observation_models::getObservableName( observableType, linkEnds.size( ) ) +
                                      ", type not yet rezognized. " );
        }

        return positionPartialScaler;
    }    

    template< typename ParameterType = double, typename TimeType = double >
    static std::shared_ptr< PositionPartialScaling > createDifferencedPositionPartialScalingObject(
            const observation_models::ObservableType differencedObservableType,
            const std::shared_ptr< PositionPartialScaling > firstPositionPartialScaling,
            const std::shared_ptr< PositionPartialScaling > secondPositionPartialScaling,
            const simulation_setup::SystemOfBodies& bodies )

    {
        std::shared_ptr< PositionPartialScaling > positionPartialScaler;

        switch( differencedObservableType )
        {
        case observation_models::one_way_differenced_range:
        {
            if( std::dynamic_pointer_cast< OneWayRangeScaling >( firstPositionPartialScaling ) == nullptr )
            {
                throw std::runtime_error( "Error when creating one-way differenced range partial scaling object, first range partial is of incompatible type" );
            }
            if( std::dynamic_pointer_cast< OneWayRangeScaling >( secondPositionPartialScaling ) == nullptr )
            {
                throw std::runtime_error( "Error when creating one-way differenced range partial scaling object, second range partial is of incompatible type" );
            }
            positionPartialScaler = std::make_shared< DifferencedObservablePartialScaling >(
                        firstPositionPartialScaling, secondPositionPartialScaling,
                        observation_models::getUndifferencedTimeAndStateIndices( observation_models::one_way_differenced_range, 2 ) );
            break;
        }
        case observation_models::n_way_differenced_range:
        case observation_models::dsn_n_way_averaged_doppler:
        {
            if( std::dynamic_pointer_cast< NWayRangeScaling >( firstPositionPartialScaling ) == nullptr )
            {
                throw std::runtime_error( "Error when creating " + observation_models::getObservableName( differencedObservableType ) +
                    " partial scaling object, first range partial is of incompatible type" );
            }
            if( std::dynamic_pointer_cast< NWayRangeScaling >( secondPositionPartialScaling ) == nullptr )
            {
                throw std::runtime_error( "Error when creating " + observation_models::getObservableName( differencedObservableType ) +
                    " partial scaling object, second range partial is of incompatible type" );
            }

            int numberOfLinkEnds = std::dynamic_pointer_cast< NWayRangeScaling >( firstPositionPartialScaling )->getNumberOfLinkEnds( );
            if( std::dynamic_pointer_cast< NWayRangeScaling >( firstPositionPartialScaling )->getNumberOfLinkEnds( ) !=
                    std::dynamic_pointer_cast< NWayRangeScaling >( secondPositionPartialScaling )->getNumberOfLinkEnds( ) )
            {
                throw std::runtime_error( "Error when creating " + observation_models::getObservableName( differencedObservableType ) +
                    " partial scaling object, first and second range partials are incompatible" );
            }
            positionPartialScaler = std::make_shared< DifferencedObservablePartialScaling >(
                        firstPositionPartialScaling, secondPositionPartialScaling,
                        observation_models::getUndifferencedTimeAndStateIndices( differencedObservableType, numberOfLinkEnds ) );
            break;
        }
        default:
            throw std::runtime_error( "Error when creating differenced observable partial scaler for " +
                                      observation_models::getObservableName( differencedObservableType ) +
                                      ", type not recognized. " );
        }
        return positionPartialScaler;
    }
};


template< >
class ObservationPartialScalingCreator< 2 >
{
public:

    template< typename ParameterType = double, typename TimeType = double >
    static std::shared_ptr< DirectPositionPartialScaling< 2 > > createPositionScalingObject(
            const observation_models::LinkEnds& linkEnds,
            const observation_models::ObservableType observableType,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< observation_models::ObservationModel< 2, ParameterType, TimeType > > observationModel = nullptr  )
    {
        std::shared_ptr< DirectPositionPartialScaling< 2 > > positionPartialScaler;

        switch( observableType )
        {
        case observation_models::angular_position:
            positionPartialScaler = std::make_shared< AngularPositionScaling >( );
            break;
        default:
            throw std::runtime_error( "Error when creating partial scaler for " +
                                      observation_models::getObservableName( observableType, linkEnds.size( ) ) +
                                      ", type not recognized. " );
        }

        return positionPartialScaler;
    }

    template< typename ParameterType = double, typename TimeType = double >
    static std::shared_ptr< PositionPartialScaling > createDifferencedPositionPartialScalingObject(
            const observation_models::ObservableType differencedObservableType,
            const std::shared_ptr< PositionPartialScaling > firstPositionPartialScaling,
            const std::shared_ptr< PositionPartialScaling > secondPositionPartialScaling,
            const simulation_setup::SystemOfBodies& bodies )

    {
        std::shared_ptr< PositionPartialScaling > positionPartialScaler;

        switch( differencedObservableType )
        {
        case observation_models::relative_angular_position:
        {
            if( std::dynamic_pointer_cast< AngularPositionScaling >( firstPositionPartialScaling ) == nullptr )
            {
                throw std::runtime_error( "Error when creating relative angular position partial scaling object, first range partial is of incompatible type" );
            }
            if( std::dynamic_pointer_cast< AngularPositionScaling >( secondPositionPartialScaling ) == nullptr )
            {
                throw std::runtime_error( "Error when creating relative angular position partial scaling object, second range partial is of incompatible type" );
            }
            std::function< void( const observation_models::LinkEndType ) > customCheckFunction =
                    []( const observation_models::LinkEndType fixedLinkEnd )
            {
                if ( fixedLinkEnd != observation_models::receiver )
                {
                    throw std::runtime_error( "Error when updating a relative angular position scaling object, fixed link end time different from receiver." );
                }
            };
            positionPartialScaler = std::make_shared< DifferencedObservablePartialScaling >(
                        firstPositionPartialScaling, secondPositionPartialScaling,
                        observation_models::getUndifferencedTimeAndStateIndices( observation_models::relative_angular_position, 3 ), customCheckFunction );
            break;
        }
        default:
            throw std::runtime_error( "Error when creating differenced observable partial scaler for " +
                                      observation_models::getObservableName( differencedObservableType ) +
                                      ", type not recognized. " );
        }
        return positionPartialScaler;
    }
};

template< >
class ObservationPartialScalingCreator< 3 >
{
public:

    template< typename ParameterType = double, typename TimeType = double >
    static std::shared_ptr< DirectPositionPartialScaling< 3 > > createPositionScalingObject(
            const observation_models::LinkEnds& linkEnds,
            const observation_models::ObservableType observableType,
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< observation_models::ObservationModel< 3, ParameterType, TimeType > > observationModel = nullptr  )
    {
        std::shared_ptr< DirectPositionPartialScaling< 3 > > positionPartialScaler;

        switch( observableType )
        {
        case observation_models::position_observable:
            positionPartialScaler = std::make_shared< PositionObservationScaling >( );
            break;
        case observation_models::velocity_observable:
            positionPartialScaler = std::make_shared< VelocityObservationScaling >( );
            break;
        case observation_models::relative_position_observable:
            positionPartialScaler = std::make_shared< RelativePositionObservationScaling >( );
            break;
        default:
            throw std::runtime_error( "Error when creating partial scaler for " +
                                      observation_models::getObservableName( observableType, linkEnds.size( ) ) +
                                      ", type not yet rezognized. " );
        }

        return positionPartialScaler;
    }
};


}

}


#endif // TUDAT_POSITIONPARTIALSCALING_H
