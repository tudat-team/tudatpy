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

typedef std::function< Eigen::VectorXd( const std::vector< double >& ,
                                        const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                                        const Eigen::VectorXd&,
                                        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ) >
                                        ObservationDependentVariableFunction;

typedef std::function< void( Eigen::VectorXd&,
                             const std::vector< double >&,
                             const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                             const Eigen::VectorXd&,
                             const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ) > ObservationDependentVariableAddFunction;

void checkObservationDependentVariableEnvironment(
        const SystemOfBodies& bodies,
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

class ObservationDependentVariableCalculator
{
public:

    ObservationDependentVariableCalculator( const observation_models::ObservableType observableType,
                                            const observation_models::LinkDefinition& linkEnds ):
    observableType_( observableType ), linkEnds_( linkEnds )
    {
        totalDependentVariableSize_ = 0;
    }

    Eigen::VectorXd calculateDependentVariables(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::VectorXd& observation,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > );

    void addDependentVariable(
            const std::shared_ptr< ObservationDependentVariableSettings > settings,
            const SystemOfBodies& bodies );

    void addDependentVariables(
            const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList,
            const SystemOfBodies& bodies );

    std::pair< int, int > getDependentVariableIndices(
            const std::shared_ptr< ObservationDependentVariableSettings > dependentVariables );

    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

    observation_models::LinkDefinition getLinkEnds( )
    {
        return linkEnds_;
    }

    std::map< std::pair< int, int >, std::shared_ptr< ObservationDependentVariableSettings > > getSettingsIndicesAndSizes( ) const
    {
        std::map< std::pair< int, int >, std::shared_ptr< ObservationDependentVariableSettings > > settingsStartIndices;
        for ( unsigned int i = 0 ; i < dependentVariableStartIndices_.size( ) ; i++ )
        {
            settingsStartIndices[ std::make_pair( dependentVariableStartIndices_[ i ], dependentVariableSizes_[ i ] ) ] = settingsList_[ i ];
        }
        return settingsStartIndices;
    }

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > getDependentVariableSettings( ) const
    {
        return settingsList_;
    }

    int getTotalDependentVariableSize( ) const
    {
        return totalDependentVariableSize_;
    }

private:

    observation_models::ObservableType observableType_;

    observation_models::LinkDefinition linkEnds_;

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList_;

    std::vector< std::function< void(
            Eigen::VectorXd&,
            const std::vector< double >&,
            const std::vector< Eigen::Matrix< double, 6, 1 > >&,
            const Eigen::VectorXd&,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ) > > dependentVariableAddFunctions_;

    std::vector< int > dependentVariableStartIndices_;

    std::vector< int > dependentVariableSizes_;

    int totalDependentVariableSize_;


};



}

}
#endif // TUDAT_OBSERVATIONOUTPUT
