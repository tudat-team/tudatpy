/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONOUTPUT
#define TUDAT_OBSERVATIONOUTPUT

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/observationInterfacesForwardDeclarations.h"
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"

namespace tudat
{

namespace simulation_setup
{

typedef std::function< Eigen::VectorXd( const std::vector< double >&,
                                        const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                                        const Eigen::VectorXd&,
                                        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ) >
        ObservationDependentVariableFunction;

typedef std::function< void( Eigen::VectorXd&,
                             const std::vector< double >&,
                             const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                             const Eigen::VectorXd&,
                             const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ) >
        ObservationDependentVariableAddFunction;

void checkObservationDependentVariableEnvironment( const SystemOfBodies& bodies,
                                                   const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

std::pair< int, int > getLinkEndStateTimeIndices(
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition linkEnds,
        const observation_models::LinkEndId linkEndId,
        const observation_models::LinkEndType linkEndRole = observation_models::unidentified_link_end,
        const observation_models::LinkEndType originatingLinkEndRole = observation_models::unidentified_link_end,
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined );

ObservationDependentVariableFunction getStationObservationAngleFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition linkEnds );

ObservationDependentVariableFunction getObservationDoubleDependentVariableFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition linkEnds );

ObservationDependentVariableFunction getObservationVectorDependentVariableFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition linkEnds );

class ObservationDependentVariableBookkeeping
{
public:
    ObservationDependentVariableBookkeeping( const observation_models::ObservableType observableType,
                                             const observation_models::LinkDefinition& linkEnds ):
        observableType_( observableType ), linkEnds_( linkEnds )
    {
        totalDependentVariableSize_ = 0;
    }

    //! Register a new dependent variable entry. Returns (startIndex, size).
    //!
    //! When `sizeOverride` is non-negative it is used verbatim (needed for dependent variables
    //! whose size is known by the Calculator but not by `getObservationDependentVariableSize`,
    //! e.g. `light_time_correction_components`). Otherwise the size is resolved via
    //! `getObservationDependentVariableSize( settings, linkEnds )`.
    std::pair< int, int > addDependentVariable( const std::shared_ptr< ObservationDependentVariableSettings > settings,
                                                const int sizeOverride = -1 );

    void addDependentVariables( const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList );

    std::pair< int, int > getDependentVariableIndices( const std::shared_ptr< ObservationDependentVariableSettings > dependentVariables );

    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

    observation_models::LinkDefinition getLinkEnds( )
    {
        return linkEnds_;
    }

    std::map< std::pair< int, int >, std::shared_ptr< ObservationDependentVariableSettings > > getSettingsIndicesAndSizes( ) const;

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > getDependentVariableSettings( ) const
    {
        return settingsList_;
    }

    int getTotalDependentVariableSize( ) const
    {
        return totalDependentVariableSize_;
    }

    void clearSettings( )
    {
        settingsList_.clear( );
        dependentVariableStartIndices_.clear( );
        dependentVariableSizes_.clear( );
        deferredSettings_.clear( );
        totalDependentVariableSize_ = 0;
    }

    //! Store a setting whose size depends on state not yet available (e.g.
    //! `light_time_correction_components` with an empty type-filter — its size equals the
    //! number of registered light-time corrections, which is only known once the observation
    //! model has been built). `ObservationDependentVariableCalculator` flushes this list and
    //! registers the settings properly once constructed with a populated leg map.
    void addDeferredSetting( const std::shared_ptr< ObservationDependentVariableSettings > setting )
    {
        for( const auto& existingSetting : settingsList_ )
        {
            if( existingSetting->areSettingsCompatible( setting ) )
            {
                return;
            }
        }
        for( const auto& existingSetting : deferredSettings_ )
        {
            if( existingSetting->areSettingsCompatible( setting ) )
            {
                return;
            }
        }
        deferredSettings_.push_back( setting );
    }

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > takeDeferredSettings( )
    {
        std::vector< std::shared_ptr< ObservationDependentVariableSettings > > out;
        out.swap( deferredSettings_ );
        return out;
    }

    const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >& getDeferredSettings( ) const
    {
        return deferredSettings_;
    }

private:
    observation_models::ObservableType observableType_;

    observation_models::LinkDefinition linkEnds_;

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList_;

    std::vector< int > dependentVariableStartIndices_;

    std::vector< int > dependentVariableSizes_;

    int totalDependentVariableSize_;

    //! Settings whose layout cannot be resolved yet (see `addDeferredSetting`). Picked up and
    //! turned into real entries by `ObservationDependentVariableCalculator` once the leg map is known.
    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > deferredSettings_;
};

class ObservationDependentVariableCalculator
{
public:
    ObservationDependentVariableCalculator(
            const std::shared_ptr< ObservationDependentVariableBookkeeping > dependentVariableBookkeeping,
            const SystemOfBodies& bodies,
            const std::map< std::pair< observation_models::LinkEndType, observation_models::LinkEndType >,
                            std::vector< std::shared_ptr< observation_models::LightTimeCalculatorBase > > >& legLightTimeCalculators ):
        dependentVariableBookkeeping_( dependentVariableBookkeeping ), legLightTimeCalculators_( legLightTimeCalculators )
    {
        // First, build add-functions for whatever settings already live on the bookkeeping. This
        // covers everything previously added through `addDependentVariable` plus any
        // `light_time_correction_components` settings whose layout was resolved on a previous
        // Calculator instance.
        for( unsigned int i = 0; i < dependentVariableBookkeeping_->getDependentVariableSettings( ).size( ); i++ )
        {
            std::pair< int, int > indices = dependentVariableBookkeeping_->getDependentVariableIndices(
                    dependentVariableBookkeeping_->getDependentVariableSettings( ).at( i ) );
            addDependentVariableFunction(
                    dependentVariableBookkeeping_->getDependentVariableSettings( ).at( i ), bodies, indices.first, indices.second );
        }

        // Then, drain any deferred `light_time_correction_components` settings. We always drain,
        // even if the leg map is empty (i.e. the observation model is one that the simulator does
        // not extract leg calculators from): `registerLightTimeCorrectionComponents` will throw
        // a clear "no LightTimeCalculator for leg X" error per entry, surfacing the unsupported
        // observable rather than silently dropping the user's request.
        for( const auto& settings : dependentVariableBookkeeping_->takeDeferredSettings( ) )
        {
            registerLightTimeCorrectionComponents( settings );
        }
    }

    Eigen::VectorXd calculateDependentVariables( const std::vector< double >& linkEndTimes,
                                                 const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                                                 const Eigen::VectorXd& observation,
                                                 const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > );

    void addDependentVariable( const std::shared_ptr< ObservationDependentVariableSettings > settings, const SystemOfBodies& bodies );

    void addDependentVariables( const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList,
                                const SystemOfBodies& bodies );

    std::pair< int, int > getDependentVariableIndices( const std::shared_ptr< ObservationDependentVariableSettings > dependentVariables );

    std::shared_ptr< ObservationDependentVariableBookkeeping > getDependentVariableBookkeeping( )
    {
        return dependentVariableBookkeeping_;
    }

    const std::map< std::pair< observation_models::LinkEndType, observation_models::LinkEndType >,
                    std::vector< std::shared_ptr< observation_models::LightTimeCalculatorBase > > >&
    getLegLightTimeCalculators( ) const
    {
        return legLightTimeCalculators_;
    }

private:
    void addDependentVariableFunction( const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
                                       const SystemOfBodies& bodies,
                                       const int currentIndex,
                                       const int parameterSize );

    std::shared_ptr< ObservationDependentVariableBookkeeping > dependentVariableBookkeeping_;

    std::vector< std::function< void( Eigen::VectorXd&,
                                      const std::vector< double >&,
                                      const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                                      const Eigen::VectorXd&,
                                      const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ) > >
            dependentVariableAddFunctions_;

    //! Per-leg light-time calculators used only for `light_time_correction_components`. Populated
    //! at simulate-time by the simulator (or left empty if no leg-specific variables are requested).
    std::map< std::pair< observation_models::LinkEndType, observation_models::LinkEndType >,
              std::vector< std::shared_ptr< observation_models::LightTimeCalculatorBase > > >
            legLightTimeCalculators_;

    //! Worker function that registers a single `light_time_correction_components` setting. Assumes
    //! `legLightTimeCalculators_` is already populated.
    void registerLightTimeCorrectionComponents( const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );
};

}  // namespace simulation_setup

}  // namespace tudat
#endif  // TUDAT_OBSERVATIONOUTPUT
