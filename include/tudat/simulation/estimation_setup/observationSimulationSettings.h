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

#include <functional>

#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/statistics/randomVariableGenerator.h"
#include "tudat/math/statistics/multiVariateGaussianProbabilityDistributions.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"

namespace tudat
{

namespace simulation_setup
{

extern int noiseSeed;

int getDefaultNoiseSeed( );

std::function< Eigen::VectorXd( const double ) > getNoiseFunctionForObservable(
        const std::function< double( const double ) > singleNoiseFunction,
        const observation_models::ObservableType observableType );


struct ObservationNoiseModel
{

public:
    ObservationNoiseModel( ): observationSize_( -1 ){ }

    virtual ~ObservationNoiseModel( ){ }

    virtual Eigen::VectorXd getObservationNoise(
        const double& observationTime,
        const Eigen::VectorXd& calculatedObservation,
        const std::vector< Eigen::Vector6d >& vectorOfStates,
        const std::vector< double >& vectorOfTimes ) = 0;

    virtual void setObservationSize( const int observationSize )
    {
        observationSize_ = observationSize;
    }

    int getObservationSize( )
    {
        return observationSize_;
    }

protected:

    int observationSize_;
};

struct UnivariateGaussianObservationNoiseModel: public ObservationNoiseModel
{

public:
    UnivariateGaussianObservationNoiseModel(
        const double noiseAmplitude,
        const double noiseMean = 0.0,
        const int gaussianNoiseSeed = getDefaultNoiseSeed( ) ):
            ObservationNoiseModel( ),
            noiseAmplitude_( noiseAmplitude ),
            noiseMean_( noiseMean )
        {
            gaussianNoiseFunction_ = statistics::createBoostContinuousRandomVariableGeneratorFunction(
                statistics::normal_boost_distribution, { noiseMean, noiseAmplitude }, gaussianNoiseSeed );
        }

    virtual ~UnivariateGaussianObservationNoiseModel( ){ }

    Eigen::VectorXd getObservationNoise(
        const double& observationTime,
        const Eigen::VectorXd& calculatedObservation,
        const std::vector< Eigen::Vector6d >& vectorOfStates,
        const std::vector< double >& vectorOfTimes )
    {
        double currentNoise = gaussianNoiseFunction_( );
        currentNoiseValue_.setConstant( currentNoise );
        return currentNoiseValue_;
    }

    virtual void setObservationSize( const int observationSize )
    {
        observationSize_ = observationSize;
        currentNoiseValue_ = Eigen::VectorXd::Zero( observationSize_ );
    }

protected:

    Eigen::VectorXd currentNoiseValue_;

    std::function< double( ) > gaussianNoiseFunction_;

    double noiseAmplitude_;

    double noiseMean_;
};

//! Base struct for defining times at which observations are to be simulated.
/*!
 *  Base struct for defining times at which observations are to be simulated. Here, only the link end from which the
 *  observation is to be calculated is defined. Derived classes are used for defining the times themselves
 *  (either directly or through some algorithm).
 */
template< typename TimeType = double >
class ObservationSimulationSettings
{
public:
    //! Constructor, defines link end type.
    /*!
     *  Constructor, defines link end type from which observations are to be simulated.
     *  \param linkEndType Link end type from which observations are to be simulated.
     */
    ObservationSimulationSettings(
            const observation_models::ObservableType observableType,
            const observation_models::LinkDefinition& linkEnds,
            const observation_models::LinkEndType linkEndType = observation_models::unidentified_link_end,
            const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
            std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ),
            const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction = nullptr,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr ):
        observableType_( observableType ), linkEnds_( linkEnds ),
        linkEndType_( linkEndType == observation_models::unidentified_link_end ? observation_models::getDefaultReferenceLinkEndType( observableType ) : linkEndType ),
        viabilitySettingsList_( viabilitySettingsList ), observationNoiseFunction_( observationNoiseFunction ),
        ancilliarySettings_( ancilliarySettings )
    {
        if( ancilliarySettings_ == nullptr )
        {
            ancilliarySettings_ = observation_models::getDefaultAncilliaryObservationSettings( observableType );
        }
        dependentVariableCalculator_ = std::make_shared< ObservationDependentVariableCalculator >(
                    observableType_, linkEnds_ );

    }

    //! Destructor.
    virtual ~ObservationSimulationSettings( ){ }

    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

    void setObservableType( observation_models::ObservableType observableType )
    {
        observableType_ = observableType;
    }

    observation_models::LinkDefinition getLinkEnds( )
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

    void setAncilliarySettings(
            std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >& ancilliarySettings )
    {
        ancilliarySettings_ = ancilliarySettings;
    }


    std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > getAncilliarySettings( )
    {
        return ancilliarySettings_;
    }

protected:

    // Type of observable to be simulated
    observation_models::ObservableType observableType_;

    // List of link ends for the observations to be simulated
    observation_models::LinkDefinition linkEnds_;

    // Reference link end type from which observations are to be simulated.
    observation_models::LinkEndType linkEndType_;

    // Settings used to check whether observtion is possible (non-viable observations are not simulated)
    std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > > viabilitySettingsList_;

    // Settings for variables that are to be saved along with the observables.
    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > observationDependentVariableSettings_;

    // Function to generate noise to add to observations that are to be simulated
    std::function< Eigen::VectorXd( const double ) > observationNoiseFunction_;

    std::shared_ptr< ObservationDependentVariableCalculator > dependentVariableCalculator_;

    std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings_;
};


//! Struct to define a list of observation times, fully defined before simulating the observations
/*!
 *  Struct to define a list of observation times, fully defined before simulating the observations. Simulations are simulated
 *  at the times stored in this struct. Some may be discarded due to the use of vaibility settins
 */
template< typename TimeType = double >
class TabulatedObservationSimulationSettings: public ObservationSimulationSettings< TimeType >
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param linkEndType Link end type from which observations are to be simulated.
     * \param simulationTimes List of times at which to perform the observation simulation
     */
    TabulatedObservationSimulationSettings(
            const observation_models::ObservableType observableType,
            const observation_models::LinkDefinition& linkEnds,
            const std::vector< TimeType >& simulationTimes,
            const observation_models::LinkEndType linkEndType = observation_models::unidentified_link_end,
            const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
            std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ),
            const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction = nullptr,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr  ):
        ObservationSimulationSettings< TimeType >(
            observableType, linkEnds, linkEndType, viabilitySettingsList, observationNoiseFunction, ancilliarySettings ),
        simulationTimes_( simulationTimes ){ }

    //! Destructor
    ~TabulatedObservationSimulationSettings( ){ }

    //! List of times at which to perform the observation simulation
    std::vector< TimeType > simulationTimes_;
};

template< typename TimeType = double >
class PerArcObservationSimulationSettings: public ObservationSimulationSettings< TimeType >
{
public:
    PerArcObservationSimulationSettings(
            const observation_models::ObservableType observableType,
            const observation_models::LinkDefinition& linkEnds,
            const TimeType startTime, const TimeType endTime, const TimeType intervalBetweenObservations,
            const std::shared_ptr< observation_models::ObservationViabilitySettings > arcDefiningConstraint,
            const TimeType minimumArcDuration = TUDAT_NAN, TimeType maximumArcDuration = TUDAT_NAN,
            const TimeType minimumTimeBetweenArcs = TUDAT_NAN,
            const observation_models::LinkEndType linkEndType = observation_models::unidentified_link_end,
            const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >&
            additionalViabilitySettingsList =
            std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ),
            const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction = nullptr,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr  ):
        ObservationSimulationSettings< TimeType >(
            observableType, linkEnds, linkEndType, additionalViabilitySettingsList, observationNoiseFunction, ancilliarySettings ),
        startTime_( startTime ), endTime_( endTime ), intervalBetweenObservations_( intervalBetweenObservations ),
        arcDefiningConstraint_( arcDefiningConstraint ),
        minimumArcDuration_( minimumArcDuration ), maximumArcDuration_( maximumArcDuration ),
        minimumTimeBetweenArcs_( minimumTimeBetweenArcs ),
    additionalViabilitySettingsList_( additionalViabilitySettingsList ){ }

    ~PerArcObservationSimulationSettings( ){ }

    TimeType startTime_;

    TimeType endTime_;

    TimeType intervalBetweenObservations_;

    std::shared_ptr< observation_models::ObservationViabilitySettings > arcDefiningConstraint_;

    TimeType minimumArcDuration_;

    TimeType maximumArcDuration_;

    TimeType minimumTimeBetweenArcs_;

    std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > > additionalViabilitySettingsList_;

};

template< typename TimeType = double >
std::shared_ptr< ObservationSimulationSettings< TimeType > > perturbObservationTime(
    const std::shared_ptr< ObservationSimulationSettings< TimeType > > originalSettings, const double timePerturbation )
{
    std::shared_ptr< ObservationSimulationSettings< TimeType > > newSettings;

    std::shared_ptr< TabulatedObservationSimulationSettings< TimeType > > originalTabulatedSettings =
        std::dynamic_pointer_cast< TabulatedObservationSimulationSettings< TimeType > >( originalSettings );
    if( std::dynamic_pointer_cast< TabulatedObservationSimulationSettings< TimeType > >( originalSettings ) != nullptr )
    {
        std::vector< TimeType > perturbedObservationTimes = originalTabulatedSettings->simulationTimes_;
        std::transform(perturbedObservationTimes.begin(), perturbedObservationTimes.end(), perturbedObservationTimes.begin(),
                  std::bind(std::plus<double>(), std::placeholders::_1, timePerturbation) );
        newSettings = std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
            originalTabulatedSettings->getObservableType( ),
            originalTabulatedSettings->getLinkEnds( ),
            perturbedObservationTimes,
            originalTabulatedSettings->getReferenceLinkEndType( ),
            originalTabulatedSettings->getViabilitySettingsList( ),
            originalTabulatedSettings->getObservationNoiseFunction( ),
            originalTabulatedSettings->getAncilliarySettings( ) );
    }
    else
    {
        throw std::runtime_error( "Error, could not perturb observation time; settings are not tabulated." );
    }
    return newSettings;
}

template< typename TimeType = double >
inline std::shared_ptr< ObservationSimulationSettings< TimeType > > tabulatedObservationSimulationSettings(
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition& linkEnds,
        const std::vector< TimeType >& simulationTimes,
        const observation_models::LinkEndType linkEndType = observation_models::receiver,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList =
        std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ),
        const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction = nullptr,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr)
{
    return std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                observableType, linkEnds, simulationTimes, linkEndType, viabilitySettingsList,
                observationNoiseFunction, ancilliarySettings );
}

template< typename TimeType = double >
inline std::shared_ptr< ObservationSimulationSettings< TimeType > > perArcObservationSimulationSettings(
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition& linkEnds,
        const TimeType startTime, const TimeType endTime, const TimeType intervalBetweenObservations,
        const std::shared_ptr< observation_models::ObservationViabilitySettings > arcDefiningConstraint,
        const TimeType minimumArcDuration = TUDAT_NAN, TimeType maximumArcDuration = TUDAT_NAN,
        const TimeType minimumTimeBetweenArcs = TUDAT_NAN,
        const observation_models::LinkEndType linkEndType = observation_models::unidentified_link_end,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >&
        additionalViabilitySettingsList =
        std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ),
        const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction = nullptr  )
{
    return std::make_shared< PerArcObservationSimulationSettings< TimeType > >(
                observableType, linkEnds,  startTime, endTime, intervalBetweenObservations,
                arcDefiningConstraint, minimumArcDuration, maximumArcDuration, minimumTimeBetweenArcs, linkEndType,
                additionalViabilitySettingsList, observationNoiseFunction );
}

template< typename TimeType = double >
std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > createTabulatedObservationSimulationSettingsList(
        const std::map< observation_models::ObservableType, std::vector< observation_models::LinkDefinition > > linkEndsPerObservable,
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
std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >
    perArcObservationSimulationSettingsList(
        const std::map< observation_models::ObservableType, std::vector< observation_models::LinkDefinition > > linkEndsPerObservable,
        const TimeType startTime, const TimeType endTime, const TimeType intervalBetweenObservations,
        const std::shared_ptr< observation_models::ObservationViabilitySettings > arcDefiningConstraint,
        const TimeType minimumArcDuration = TUDAT_NAN, TimeType maximumArcDuration = TUDAT_NAN,
        const TimeType minimumTimeBetweenArcs = TUDAT_NAN,
        const observation_models::LinkEndType linkEndType = observation_models::unidentified_link_end,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >&
        additionalViabilitySettingsList =
        std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ) )
{
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > observationSimulationSettingsList;
    for( auto observableIterator : linkEndsPerObservable )
    {
        for( unsigned int i = 0; i < observableIterator.second.size( ); i++ )
        {
            observationSimulationSettingsList.push_back(
                        std::make_shared< PerArcObservationSimulationSettings< TimeType > >(
                            observableIterator.first, observableIterator.second.at( i ),  startTime, endTime, intervalBetweenObservations,
                            arcDefiningConstraint, minimumArcDuration, maximumArcDuration, minimumTimeBetweenArcs, linkEndType,
                            additionalViabilitySettingsList ) );
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

template< typename TimeType = double >
void addDependentVariableToSingleObservationSimulationSettings(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >& dependentVariableList,
        const SystemOfBodies& bodies )
{
    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > extendedDependentVariablesList;

    // Retrieve interlinks information for current observable type and link ends
    ObservableType observableType = observationSimulationSettings->getObservableType( );
    LinkEnds linkEnds = observationSimulationSettings->getLinkEnds( ).linkEnds_;
    std::vector< std::pair< std::pair< LinkEndType, LinkEndId >, std::pair< LinkEndType, LinkEndId > > > interlinksInSet = getInterlinks( observableType, linkEnds );

    // Parse all dependent variable settings
    for ( auto settings : dependentVariableList )
    {
        // Create complete list of all dependent variable settings compatible with the original settings (possibly not fully defined, i.e. with
        // missing information on link ends, etc.) for the given observable type and link ends
        std::vector< std::shared_ptr< ObservationDependentVariableSettings > > allSettingsToCreate =
                createAllCompatibleDependentVariableSettings( observableType, linkEnds, settings );
        for ( auto it : allSettingsToCreate )
        {
            extendedDependentVariablesList.push_back( it );
        }
    }
    // Add all relevant dependent variables
    observationSimulationSettings->getDependentVariableCalculator( )->addDependentVariables( extendedDependentVariablesList, bodies );
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
void addAncilliarySettingsToSingleObservationSimulationSettings(
        const std::shared_ptr< ObservationSimulationSettings< TimeType > >& observationSimulationSettings,
        std::shared_ptr<observation_models::ObservationAncilliarySimulationSettings> &ancilliarySettings)
{
    observationSimulationSettings->setAncilliarySettings(
            const_cast<std::shared_ptr<observation_models::ObservationAncilliarySimulationSettings> &>(ancilliarySettings));
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
        const observation_models::LinkDefinition& linkEnds )
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

template< typename TimeType = double, typename... ArgTypes  >
void addAncilliarySettingsToObservationSimulationSettings(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >& observationSimulationSettings,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >& ancilliarySettings,
        ArgTypes... args)
{

    std::function< void(  std::shared_ptr< ObservationSimulationSettings< TimeType > > ) > modificationFunction =
            std::bind( &addAncilliarySettingsToSingleObservationSimulationSettings< TimeType >,
                       std::placeholders::_1, ancilliarySettings );
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
//            const observation_models::LinkDefinition& linkEnds,
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
        const std::map< observation_models::ObservableType, std::vector< observation_models::LinkDefinition > >& linkEndsPerObservable,
        const std::vector< TimeType >& observationTimes,
        const observation_models::LinkEndType referenceLinkEnd = observation_models::receiver )
{
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;
    for( auto it : linkEndsPerObservable )
    {
        observation_models::ObservableType currentObservable = it.first;
        std::vector< observation_models::LinkDefinition > currentLinkEndsList = it.second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                            currentObservable, currentLinkEndsList.at( i ), observationTimes, referenceLinkEnd ) );
        }
    }
    return measurementSimulationInput;
}

template< typename TimeType >
std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > >  getObservationSimulationSettings(
        const std::map< observation_models::ObservableType, std::vector< observation_models::LinkEnds > >& linkEndsPerObservable,
        const std::vector< TimeType >& observationTimes,
        const observation_models::LinkEndType referenceLinkEnd = observation_models::receiver )
{
    std::map< observation_models::ObservableType, std::vector< observation_models::LinkDefinition > > linkDefsPerObservable;
    for( auto it : linkEndsPerObservable )
    {
        for( unsigned int i = 0; i < it.second.size( ); i++ )
        {
            linkDefsPerObservable[ it.first ].push_back( it.second.at( i ) );
        }
    }
    return getObservationSimulationSettings( linkDefsPerObservable, observationTimes, referenceLinkEnd );
}

//extern template class ObservationSimulationSettings< double >;
//extern template class TabulatedObservationSimulationSettings< double >;
//extern template class PerArcObservationSimulationSettings< double >;


}

}
#endif // TUDAT_OBSERVATIONSIMULATIONSETTINGS_H
