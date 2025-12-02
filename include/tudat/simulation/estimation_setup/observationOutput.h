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

#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"

namespace tudat
{

namespace simulation_setup
{

typedef std::function< Eigen::VectorXd( const std::vector< double >&,
                                        const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                                        const Eigen::VectorXd&,
                                        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ) >
        ObservationDependentVariableFunction;

typedef std::function< void( Eigen::VectorXd&,
                             const std::vector< double >&,
                             const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                             const Eigen::VectorXd&,
                             const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ) >
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

    std::pair< int, int > addDependentVariable( const std::shared_ptr< ObservationDependentVariableSettings > settings );

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
        totalDependentVariableSize_ = 0;
    }

private:
    observation_models::ObservableType observableType_;

    observation_models::LinkDefinition linkEnds_;

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList_;

    std::vector< int > dependentVariableStartIndices_;

    std::vector< int > dependentVariableSizes_;

    int totalDependentVariableSize_;
};

class ObservationDependentVariableCalculator
{
public:
    ObservationDependentVariableCalculator( const observation_models::ObservableType observableType,
                                            const observation_models::LinkDefinition& linkEnds ):
        dependentVariableBookkeeping_( std::make_shared< ObservationDependentVariableBookkeeping >( observableType, linkEnds ) )
    {}

    ObservationDependentVariableCalculator( const std::shared_ptr< ObservationDependentVariableBookkeeping > dependentVariableBookkeeping,
                                            const SystemOfBodies& bodies ): dependentVariableBookkeeping_( dependentVariableBookkeeping )
    {
        for( unsigned int i = 0; i < dependentVariableBookkeeping_->getDependentVariableSettings( ).size( ); i++ )
        {
            std::pair< int, int > indices = dependentVariableBookkeeping_->getDependentVariableIndices(
                    dependentVariableBookkeeping_->getDependentVariableSettings( ).at( i ) );
            addDependentVariableFunction(
                    dependentVariableBookkeeping_->getDependentVariableSettings( ).at( i ), bodies, indices.first, indices.second );
        }
    }

    Eigen::VectorXd calculateDependentVariables( const std::vector< double >& linkEndTimes,
                                                 const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                                                 const Eigen::VectorXd& observation,
                                                 const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > );

    void addDependentVariable( const std::shared_ptr< ObservationDependentVariableSettings > settings, const SystemOfBodies& bodies );

    void addDependentVariables( const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList,
                                const SystemOfBodies& bodies );

    std::pair< int, int > getDependentVariableIndices( const std::shared_ptr< ObservationDependentVariableSettings > dependentVariables );

    std::shared_ptr< ObservationDependentVariableBookkeeping > getDependentVariableBookkeeping( )
    {
        return dependentVariableBookkeeping_;
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
                                      const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ) > >
            dependentVariableAddFunctions_;
};

}  // namespace simulation_setup

}  // namespace tudat
#endif  // TUDAT_OBSERVATIONOUTPUT
