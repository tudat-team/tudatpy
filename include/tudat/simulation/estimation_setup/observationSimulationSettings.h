/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONSIMULATIONSETTINGS_H
#define TUDAT_OBSERVATIONSIMULATIONSETTINGS_H

#include <memory>
#include <boost/bind.hpp>
#include <functional>

#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/simulation/estimation_setup/observations.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/statistics/randomVariableGenerator.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"

namespace tudat
{

namespace simulation_setup
{

extern int noiseSeed;

std::function< Eigen::VectorXd( const double ) > getNoiseFunctionForObservable(
        const std::function< double( const double ) > singleNoiseFunction,
        const observation_models::ObservableType observableType );

//! Base struct for defining times at which observations are to be simulated.
/*!
 *  Base struct for defining times at which observations are to be simulated. Here, only the link end from which the
 *  observation is to be calculated is defined. Derived classes are used for defining the times themselves
 *  (either directly or through some algorithm).
 */
template< typename TimeType = double >
struct ObservationSimulationSettings
{
    //! Constructor, defines link end type.
    /*!
     *  Constructor, defines link end type from which observations are to be simulated.
     *  \param linkEndType Link end type from which observations are to be simulated.
     */
    ObservationSimulationSettings(
            const observation_models::ObservableType observableType,
            const observation_models::LinkEnds linkEnds,
            const observation_models::LinkEndType linkEndType = observation_models::unidentified_link_end,
            const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
            std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ),
            const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction = nullptr ):
        observableType_( observableType ), linkEnds_( linkEnds ),
        linkEndType_( linkEndType == observation_models::unidentified_link_end ? observation_models::getDefaultReferenceLinkEndType( observableType ) : linkEndType ),
        viabilitySettingsList_( viabilitySettingsList ), observationNoiseFunction_( observationNoiseFunction )
    {
        dependentVariableCalculator_ = std::make_shared< ObservationDependentVariableCalculator >(
                    observableType_, linkEnds_ );
    }

    //! Destructor.
    virtual ~ObservationSimulationSettings( ){ }

    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

    observation_models::LinkEnds getLinkEnds( )
    {
        return linkEnds_;
    }

    observation_models::LinkEndType getReferenceLinkEndType( )
    {
        return linkEndType_;
    }

    std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > > getViabilitySettingsList( )
    {
        return viabilitySettingsList_;
    }

    void setViabilitySettingsList(
            const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList )
    {
        viabilitySettingsList_ = viabilitySettingsList;
    }

    std::function< Eigen::VectorXd( const double ) > getObservationNoiseFunction( )
    {
        return observationNoiseFunction_;
    }

    void setObservationNoiseFunction(
            const std::function< Eigen::VectorXd( const double ) >& observationNoiseFunction )
    {
        observationNoiseFunction_ = observationNoiseFunction;
    }

    void setObservationNoiseFunction(
            const std::function< double( const double ) >& observationNoiseFunction )
    {
        observationNoiseFunction_ = getNoiseFunctionForObservable(
                    observationNoiseFunction, observableType_ );
    }

    std::shared_ptr< ObservationDependentVariableCalculator > getDependentVariableCalculator( )
    {
        return dependentVariableCalculator_;
    }

protected:

    // Type of observable to be simulated
    observation_models::ObservableType observableType_;

    // List of link ends for the observations to be simulated
    observation_models::LinkEnds linkEnds_;

    // Reference link end type from which observations are to be simulated.
    observation_models::LinkEndType linkEndType_;

    // Settings used to check whether observtion is possible (non-viable observations are not simulated)
    std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > > viabilitySettingsList_;

    // Settings for variables that are to be saved along with the observables.
    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > observationDependentVariableSettings_;

    // Function to generate noise to add to observations that are to be simulated
    std::function< Eigen::VectorXd( const double ) > observationNoiseFunction_;

    std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator_;
};


//! Struct to define a list of observation times, fully defined before simulating the observations
/*!
 *  Struct to define a list of observation times, fully defined before simulating the observations. Simulations are simulated
 *  at the times stored in this struct. Some may be discarded due to the use of vaibility settins
 */
template< typename TimeType = double >
struct TabulatedObservationSimulationSettings: public ObservationSimulationSettings< TimeType >
{
    //! Constructor
    /*!
     * Constructor
     * \param linkEndType Link end type from which observations are to be simulated.
     * \param simulationTimes List of times at which to perform the observation simulation
     */
    TabulatedObservationSimulationSettings(
            const observation_models::ObservableType observableType,
            const observation_models::LinkEnds linkEnds,
            const std::vector< TimeType >& simulationTimes,
            const observation_models::LinkEndType linkEndType = observation_models::unidentified_link_end,
            const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
            std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ),
            const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction = nullptr  ):
        ObservationSimulationSettings< TimeType >(
            observableType, linkEnds, linkEndType, viabilitySettingsList, observationNoiseFunction ),
        simulationTimes_( simulationTimes ){ }

    //! Destructor
    ~TabulatedObservationSimulationSettings( ){ }

    //! List of times at which to perform the observation simulation
    std::vector< TimeType > simulationTimes_;
};

template< typename TimeType = double >
inline std::shared_ptr< ObservationSimulationSettings< TimeType > > tabulatedObservationSimulationSettings(
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds,
        const std::vector< TimeType >& simulationTimes,
        const observation_models::LinkEndType linkEndType = observation_models::receiver,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
        std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ),
        const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction = nullptr  )
{
    return std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                observableType, linkEnds, simulationTimes, linkEndType, viabilitySettingsList,
                observationNoiseFunction );
}

template< typename TimeType = double >
std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > createTabulatedObservationSimulationSettingsList(
        const std::map< observation_models::ObservableType, std::vector< observation_models::LinkEnds > > linkEndsPerObservable,
        const std::vector< TimeType >& simulationTimes,
        const observation_models::LinkEndType linkEndType = observation_models::receiver,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
        std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( )
        )
{
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > observationSimulationSettingsList;
    for( auto observableIterator : linkEndsPerObservable )
    {
        for( unsigned int i = 0; i < observableIterator.second.size( ); i++ )
        {
            observationSimulationSettingsList.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                            observableIterator.first, observableIterator.second.at( i ), simulationTimes, linkEndType,
                            viabilitySettingsList) );
        }
    }
    return observationSimulationSettingsList;
}

template< typename TimeType = double >
void clearNoiseFunctionFromObservationSimulationSettings(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationSimulationSettings )
{
    for( unsigned int i = 0; i < observationSimulationSettings.size( ); i++ )
    {
        observationSimulationSettings.at( i )->setObservationNoiseFunction( std::function< Eigen::VectorXd( const double ) >( ) );
    }
}

void addDependentVariables(
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList,
        const SystemOfBodies& bodies );

template< typename TimeType = double >
void addDependentVariableToSingleObservationSimulationSettings(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >& dependentVariableList,
        const SystemOfBodies& bodies )
{
    observationSimulationSettings->getDependentVariableCalculator( )->addDependentVariables(
                dependentVariableList, bodies );
}

template< typename TimeType = double >
void addViabilityToSingleObservationSimulationSettings(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList )
{
    std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > > viabilitySettingsToAdd;
    for( unsigned int j = 0; j < viabilitySettingsList.size( ); j++ )
    {
        viabilitySettingsToAdd.push_back( viabilitySettingsList.at( j ) );
    }

    if( viabilitySettingsToAdd.size( ) > 0 )
    {
        std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > > currentViabilitySettingsList =
                observationSimulationSettings->getViabilitySettingsList( );
        currentViabilitySettingsList.insert(
                    currentViabilitySettingsList.end( ), viabilitySettingsToAdd.begin( ), viabilitySettingsToAdd.end( ) );
        observationSimulationSettings->setViabilitySettingsList( currentViabilitySettingsList );
    }
}

template< typename TimeType = double, typename DataType >
void addNoiseToSingleObservationSimulationSettings(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > > observationSimulationSettings,
        const std::function< DataType( const double ) > observationNoiseFunction )
{
    observationSimulationSettings->setObservationNoiseFunction( observationNoiseFunction );
}

template< typename TimeType = double >
void addGaussianNoiseToSingleObservationSimulationSettings(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > >& observationSimulationSettings,
        const double observationNoiseAmplitude )
{
    std::function< Eigen::VectorXd( const double ) > noiseFunction =
            statistics::getIndependentGaussianNoiseFunction(
                observationNoiseAmplitude, 0.0, noiseSeed, observation_models::getObservableSize(
                    observationSimulationSettings->getObservableType( ) ) );
    noiseSeed++;

    observationSimulationSettings->setObservationNoiseFunction( noiseFunction );
}

template< typename TimeType = double >
void modifyObservationSimulationSettings(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationSimulationSettings,
        const std::function< void( const std::shared_ptr< ObservationSimulationSettings< TimeType > > ) > modificationFunction )
{
    for( unsigned int i = 0; i < observationSimulationSettings.size( ); i++ )
    {
        modificationFunction( observationSimulationSettings.at( i ) );
    }
}

template< typename TimeType = double >
void modifyObservationSimulationSettings(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationSimulationSettings,
        const std::function< void( const std::shared_ptr< ObservationSimulationSettings< TimeType > > ) > modificationFunction,
        const observation_models::ObservableType observableType )
{
    for( unsigned int i = 0; i < observationSimulationSettings.size( ); i++ )
    {
        if( observationSimulationSettings.at( i )->getObservableType( ) == observableType )
        {
            modificationFunction( observationSimulationSettings.at( i ) );
        }
    }
}

template< typename TimeType = double >
void modifyObservationSimulationSettings(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationSimulationSettings,
        const std::function< void( const std::shared_ptr< ObservationSimulationSettings< TimeType > > ) > modificationFunction,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds& linkEnds )
{
    for( unsigned int i = 0; i < observationSimulationSettings.size( ); i++ )
    {
        if( observationSimulationSettings.at( i )->getObservableType( ) == observableType &&
                observationSimulationSettings.at( i )->getLinkEnds( ) == linkEnds )
        {
            modificationFunction( observationSimulationSettings.at( i ) );
        }
    }
}

template< typename TimeType = double, typename... ArgTypes >
void addViabilityToObservationSimulationSettings(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList,
        ArgTypes... args )
{
    std::function< void( const std::shared_ptr< ObservationSimulationSettings< TimeType > > ) > modificationFunction =
            std::bind( &addViabilityToSingleObservationSimulationSettings< TimeType >,
                       std::placeholders::_1, viabilitySettingsList );
    modifyObservationSimulationSettings(
                observationSimulationSettings,
                modificationFunction, args ... );
}


template< typename TimeType = double, typename DataType = double, typename... ArgTypes >
void addNoiseFunctionToObservationSimulationSettings(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationSimulationSettings,
        const std::function< DataType( const double ) > observationNoiseFunction,
        ArgTypes... args )
{
    std::function< void( const std::shared_ptr< ObservationSimulationSettings< TimeType > > ) > modificationFunction =
            std::bind( &addNoiseToSingleObservationSimulationSettings< TimeType, DataType >,
                       std::placeholders::_1, observationNoiseFunction );
    modifyObservationSimulationSettings(
                observationSimulationSettings, modificationFunction, args ... );
}

template< typename TimeType = double, typename... ArgTypes  >
void addGaussianNoiseFunctionToObservationSimulationSettings(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationSimulationSettings,
        const double observationNoiseAmplitude,
        ArgTypes... args )
{
    std::function< void( const std::shared_ptr< ObservationSimulationSettings< TimeType > > ) > modificationFunction =
            std::bind( &addGaussianNoiseToSingleObservationSimulationSettings< TimeType >,
                       std::placeholders::_1, observationNoiseAmplitude );
    modifyObservationSimulationSettings(
                observationSimulationSettings, modificationFunction, args ... );
}

template< typename TimeType = double, typename... ArgTypes  >
void addDependentVariablesToObservationSimulationSettings(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >& dependentVariableList,
        const SystemOfBodies& bodies,
        ArgTypes... args )
{
    std::function< void( const std::shared_ptr< ObservationSimulationSettings< TimeType > > ) > modificationFunction =
            std::bind( &addDependentVariableToSingleObservationSimulationSettings< TimeType >,
                       std::placeholders::_1, dependentVariableList, bodies );
    modifyObservationSimulationSettings(
                observationSimulationSettings, modificationFunction, args ... );
}


////! Function to simulate a fixed number of simulations, in an arcwise manner, taking into account viability settings
///*!
// *  Function to simulate a fixed number of simulations, in an arcwise manner, taking into account viability settings. This class
// *  defines an observation arc starting every X seconds, with observations simulated every M seconds from the start of this
// *  interval until N observations have been simulated, taking into account that some are discarded due to vaibiloty settings.
// */
//template< typename TimeType = double >
//struct ArcLimitedObservationSimulationSettings: public ObservationSimulationSettings< TimeType >
//{
//    //! Constructor
//    /*!
//     * Constructor
//     * \param linkEndType Reference link end type
//     * \param startTime Time at which to start simulating observations
//     * \param endTime Time at which to end simulating observations
//     * \param observationInterval Time between two subsequent observations (e.g. integration time)
//     * \param arcDuration Duration of an arc over which observationLimitPerArc observations are to be simulated
//     * \param observationLimitPerArc Number of observations that are to be simulated per arc
//     */
//    ArcLimitedObservationSimulationSettings(
//            const observation_models::ObservableType observableType,
//            const observation_models::LinkEnds linkEnds,
//            const TimeType startTime, const TimeType endTime, const TimeType observationInterval,
//            const TimeType arcDuration, const int observationLimitPerArc,
//            const observation_models::LinkEndType linkEndType = observation_models::receiver,
//            const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
//            std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ) ):
//        ObservationSimulationSettings< TimeType >(
//            observableType, linkEnds, linkEndType, viabilitySettingsList ),
//        startTime_( startTime ), endTime_( endTime ), observationInterval_( observationInterval ),
//        arcDuration_( arcDuration ), observationLimitPerArc_( observationLimitPerArc ){ }

//    ~ArcLimitedObservationSimulationSettings( ){ }

//    //! Time at which to start simulating observations
//    TimeType startTime_;

//    //! Time at which to end simulating observations
//    TimeType endTime_;

//    //! Time between two subsequent observations (e.g. integration time)
//    TimeType observationInterval_;

//    //! Duration of an arc over which observationLimitPerArc observations are to be simulated
//    TimeType arcDuration_;

//    //! Number of observations that are to be simulated per arc
//    int observationLimitPerArc_;
//};

template< typename TimeType >
std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >  getObservationSimulationSettings(
        const std::map< observation_models::ObservableType, std::vector< observation_models::LinkEnds > >& linkEndsPerObservable,
        const std::vector< TimeType >& observationTimes,
        const observation_models::LinkEndType referenceLinkEnd = observation_models::receiver )
{
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;
    for( auto it : linkEndsPerObservable )
    {
        observation_models::ObservableType currentObservable = it.first;
        std::vector< observation_models::LinkEnds > currentLinkEndsList = it.second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                            currentObservable, currentLinkEndsList.at( i ), observationTimes, referenceLinkEnd ) );
        }
    }
    return measurementSimulationInput;
}

}

}
#endif // TUDAT_OBSERVATIONSIMULATIONSETTINGS_H
