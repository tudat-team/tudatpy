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

typedef std::function< double( const std::vector< double >& ,
                               const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                               const Eigen::VectorXd& ) > DoubleObservationDependentVariableFunction;

typedef std::function< Eigen::VectorXd( const std::vector< double >& ,
                                        const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                                        const Eigen::VectorXd& ) > VectorObservationDependentVariableFunction;

typedef std::function< void( Eigen::VectorXd&,
                             const std::vector< double >&,
                             const std::vector< Eigen::Matrix< double, 6, 1 > >&,
                             const Eigen::VectorXd& ) > ObservationDependentVariableAddFunction;

void checkObservationDependentVariableEnvironment(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );


DoubleObservationDependentVariableFunction getBodyAvoidanceFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< BodyAvoidanceObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds );

DoubleObservationDependentVariableFunction getTargetRangeFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< InterlinkObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds );

DoubleObservationDependentVariableFunction getStationObservationAngleFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds );

DoubleObservationDependentVariableFunction getObservationDoubleDependentVariableFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds );

VectorObservationDependentVariableFunction getObservationVectorDependentVariableFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds );

class ObservationDependentVariableCalculator
{
public:

    ObservationDependentVariableCalculator( const observation_models::ObservableType observableType,
                                            const observation_models::LinkEnds& linkEnds ):
    observableType_( observableType ), linkEnds_( linkEnds )
    {
        totalDependentVariableSize_ = 0.0;
    }

    Eigen::VectorXd calculateDependentVariables(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::VectorXd& observation );

    void addDependentVariable(
            const std::shared_ptr< ObservationDependentVariableSettings > settings,
            const SystemOfBodies& bodies );

    void addDependentVariables(
            const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList,
            const SystemOfBodies& bodies );

    std::pair< int, int > getDependentVariableIndices(
            const std::shared_ptr< ObservationDependentVariableSettings > dependentVariables );

private:

    observation_models::ObservableType observableType_;

    observation_models::LinkEnds linkEnds_;

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList_;

    std::vector< std::function< void(
            Eigen::VectorXd&,
            const std::vector< double >&,
            const std::vector< Eigen::Matrix< double, 6, 1 > >&,
            const Eigen::VectorXd& ) > > dependentVariableAddFunctions_;

    std::vector< int > dependentVariableStartIndices_;

    std::vector< int > dependentVariableSizes_;

    int totalDependentVariableSize_;


};



}

}
#endif // TUDAT_OBSERVATIONOUTPUT
