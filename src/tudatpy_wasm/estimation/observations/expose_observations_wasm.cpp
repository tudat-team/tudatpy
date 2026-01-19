/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "../../wasm_module.h"
#include "../../eigen_wasm.h"
#include "../../stl_wasm.h"
#include "../../shared_ptr_wasm.h"

#include <tudat/simulation/estimation_setup/singleObservationSet.h>
#include <tudat/simulation/estimation_setup/observationCollection.h>
#include <tudat/simulation/estimation_setup/simulateObservations.h>

namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

// Helper function to create SingleObservationSet without dependent variables (simpler signature)
namespace tudat {
namespace simulation_setup {

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< tom::SingleObservationSet< ObservationScalarType, TimeType > >
singleObservationSetWithoutDependentVariables(
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                observations,
        const std::vector< TimeType > observationTimes,
        const tom::LinkEndType referenceLinkEnd,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >
                ancilliarySettings = nullptr )
{
    return std::make_shared< tom::SingleObservationSet< ObservationScalarType, TimeType > >(
            observableType,
            linkEnds,
            observations,
            observationTimes,
            referenceLinkEnd,
            std::vector< Eigen::VectorXd >( ),
            nullptr,
            ancilliarySettings );
}

} // namespace simulation_setup
} // namespace tudat

WASM_MODULE_PATH("estimation_observations")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observations) {
    using namespace emscripten;

    // ========================================================================
    // SingleObservationSet class - Full implementation
    // ========================================================================
    class_<tom::SingleObservationSet<double, double>>("estimation_observations_SingleObservationSet")
        .smart_ptr<std::shared_ptr<tom::SingleObservationSet<double, double>>>(
            "shared_ptr_SingleObservationSet")
        // Constructor
        .constructor<const tom::ObservableType,
                     const tom::LinkDefinition,
                     const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>>,
                     const std::vector<double>,
                     const tom::LinkEndType,
                     const std::vector<Eigen::VectorXd>,
                     const std::shared_ptr<tss::ObservationDependentVariableCalculator>,
                     const std::shared_ptr<tom::ObservationAncilliarySimulationSettings>>()
        // Getter methods
        .function("getObservableType", &tom::SingleObservationSet<double, double>::getObservableType)
        .function("getLinkEnds", &tom::SingleObservationSet<double, double>::getLinkEnds)
        .function("setLinkEnds", &tom::SingleObservationSet<double, double>::setLinkEnds)
        .function("getReferenceLinkEnd", &tom::SingleObservationSet<double, double>::getReferenceLinkEnd)
        .function("getNumberOfObservables", &tom::SingleObservationSet<double, double>::getNumberOfObservables)
        .function("getSingleObservableSize", &tom::SingleObservationSet<double, double>::getSingleObservableSize)
        .function("getTotalObservationSetSize", &tom::SingleObservationSet<double, double>::getTotalObservationSetSize)
        .function("getTimeBounds", &tom::SingleObservationSet<double, double>::getTimeBounds)
        .function("getObservations", &tom::SingleObservationSet<double, double>::getObservations)
        .function("getObservationTimes", &tom::SingleObservationSet<double, double>::getObservationTimes)
        .function("getObservationsVector", &tom::SingleObservationSet<double, double>::getObservationsVector)
        .function("getComputedObservations", &tom::SingleObservationSet<double, double>::getComputedObservations)
        .function("getComputedObservationsVector", &tom::SingleObservationSet<double, double>::getComputedObservationsVector)
        .function("getResiduals", &tom::SingleObservationSet<double, double>::getResiduals)
        .function("getResidualsVector", &tom::SingleObservationSet<double, double>::getResidualsVector)
        .function("getRmsResiduals", &tom::SingleObservationSet<double, double>::getRmsResiduals)
        .function("getMeanResiduals", &tom::SingleObservationSet<double, double>::getMeanResiduals)
        .function("getWeights", &tom::SingleObservationSet<double, double>::getWeights)
        .function("getWeightsVector", &tom::SingleObservationSet<double, double>::getWeightsVector)
        .function("getObservationsDependentVariables", &tom::SingleObservationSet<double, double>::getObservationsDependentVariables)
        .function("getDependentVariableHistory", &tom::SingleObservationSet<double, double>::getDependentVariableHistory)
        .function("getObservationsHistory", &tom::SingleObservationSet<double, double>::getObservationsHistory)
        .function("getAncilliarySettings", &tom::SingleObservationSet<double, double>::getAncilliarySettings)
        .function("getFilteredObservationSet", &tom::SingleObservationSet<double, double>::getFilteredObservationSet)
        .function("getNumberOfFilteredObservations", &tom::SingleObservationSet<double, double>::getNumberOfFilteredObservations)
        .function("getObservationsDependentVariablesMatrix", &tom::SingleObservationSet<double, double>::getObservationsDependentVariablesMatrix)
        // Setter methods - using overload resolution with select_overload
        .function("setObservationsFromVector",
            select_overload<void(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>>&)>(
                &tom::SingleObservationSet<double, double>::setObservations))
        .function("setObservationsFromEigen",
            select_overload<void(const Eigen::Matrix<double, Eigen::Dynamic, 1>&)>(
                &tom::SingleObservationSet<double, double>::setObservations))
        .function("setResidualsFromVector",
            select_overload<void(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>>&)>(
                &tom::SingleObservationSet<double, double>::setResiduals))
        .function("setResidualsFromEigen",
            select_overload<void(const Eigen::Matrix<double, Eigen::Dynamic, 1>&)>(
                &tom::SingleObservationSet<double, double>::setResiduals))
        .function("setConstantWeightScalar",
            select_overload<void(const double)>(
                &tom::SingleObservationSet<double, double>::setConstantWeight))
        .function("setConstantWeightVector",
            select_overload<void(const Eigen::Matrix<double, Eigen::Dynamic, 1>&)>(
                &tom::SingleObservationSet<double, double>::setConstantWeight))
        .function("setTabulatedWeights",
            select_overload<void(const Eigen::VectorXd&)>(
                &tom::SingleObservationSet<double, double>::setTabulatedWeights))
        .function("setObservationsDependentVariables", &tom::SingleObservationSet<double, double>::setObservationsDependentVariables)
        .function("filterObservations", &tom::SingleObservationSet<double, double>::filterObservations)
        .function("getCompatibleDependentVariablesSettingsList", &tom::SingleObservationSet<double, double>::getCompatibleDependentVariablesSettingsList)
        .function("getAllCompatibleDependentVariables", &tom::SingleObservationSet<double, double>::getAllCompatibleDependentVariables)
        .function("getSingleDependentVariableHistory", &tom::SingleObservationSet<double, double>::getSingleDependentVariableHistory);

    // Factory function for creating SingleObservationSet without dependent variables
    function("estimation_observations_single_observation_set",
        &tss::singleObservationSetWithoutDependentVariables<double, double>);

    // ========================================================================
    // ObservationCollection class - Full implementation
    // ========================================================================
    class_<tom::ObservationCollection<double, double>>("estimation_observations_ObservationCollection")
        .smart_ptr<std::shared_ptr<tom::ObservationCollection<double, double>>>(
            "shared_ptr_ObservationCollection")
        // Constructor
        .constructor<std::vector<std::shared_ptr<tom::SingleObservationSet<double, double>>>>()
        // Basic properties
        .function("getTotalObservableSize", &tom::ObservationCollection<double, double>::getTotalObservableSize)
        .function("getConcatenatedTimeVector", &tom::ObservationCollection<double, double>::getConcatenatedTimeVector)
        .function("getConcatenatedDoubleTimeVector", &tom::ObservationCollection<double, double>::getConcatenatedDoubleTimeVector)
        .function("getUnparsedConcatenatedWeights", &tom::ObservationCollection<double, double>::getUnparsedConcatenatedWeights)
        .function("getObservationVector", &tom::ObservationCollection<double, double>::getObservationVector)
        .function("getObservationTypeStartAndSize", &tom::ObservationCollection<double, double>::getObservationTypeStartAndSize)
        .function("getObservationSetStartAndSizePerLinkEndIndex", &tom::ObservationCollection<double, double>::getObservationSetStartAndSizePerLinkEndIndex)
        .function("getSortedObservationSets", &tom::ObservationCollection<double, double>::getSortedObservationSets)
        .function("getInverseLinkEndIdentifierMap", &tom::ObservationCollection<double, double>::getInverseLinkEndIdentifierMap)
        .function("getLinkEndsPerObservableType", &tom::ObservationCollection<double, double>::getLinkEndsPerObservableType)
        .function("getLinkDefinitionsPerObservable", &tom::ObservationCollection<double, double>::getLinkDefinitionsPerObservable)
        .function("getTimeBoundsDouble", &tom::ObservationCollection<double, double>::getTimeBoundsDouble)
        .function("getTimeBounds", &tom::ObservationCollection<double, double>::getTimeBounds)
        .function("getSortedObservationSetsTimeBounds", &tom::ObservationCollection<double, double>::getSortedObservationSetsTimeBounds)
        // Set methods
        .function("setObservationsFromEigen",
            select_overload<void(const Eigen::Matrix<double, Eigen::Dynamic, 1>&)>(
                &tom::ObservationCollection<double, double>::setObservations))
        .function("setObservationsWithParser",
            select_overload<void(const Eigen::Matrix<double, Eigen::Dynamic, 1>&, const std::shared_ptr<tom::ObservationCollectionParser>)>(
                &tom::ObservationCollection<double, double>::setObservations))
        .function("setResidualsFromEigen",
            select_overload<void(const Eigen::Matrix<double, Eigen::Dynamic, 1>&)>(
                &tom::ObservationCollection<double, double>::setResiduals))
        .function("setResidualsWithParser",
            select_overload<void(const Eigen::Matrix<double, Eigen::Dynamic, 1>&, const std::shared_ptr<tom::ObservationCollectionParser>)>(
                &tom::ObservationCollection<double, double>::setResiduals))
        // Get methods with parser
        .function("getLinkDefinitionsForSingleObservable", &tom::ObservationCollection<double, double>::getLinkDefinitionsForSingleObservable)
        .function("getSingleLinkAndTypeObservationSets", &tom::ObservationCollection<double, double>::getSingleLinkAndTypeObservationSets)
        .function("getObservableTypes", &tom::ObservationCollection<double, double>::getObservableTypes)
        .function("getBodiesInLinkEnds", &tom::ObservationCollection<double, double>::getBodiesInLinkEnds)
        .function("getReferencePointsInLinkEnds", &tom::ObservationCollection<double, double>::getReferencePointsInLinkEnds)
        .function("getTimeBoundsListDouble", &tom::ObservationCollection<double, double>::getTimeBoundsListDouble)
        .function("getTimeBoundsList", &tom::ObservationCollection<double, double>::getTimeBoundsList)
        .function("getTimeBoundsPerSetDouble", &tom::ObservationCollection<double, double>::getTimeBoundsPerSetDouble)
        .function("getTimeBoundsPerSet", &tom::ObservationCollection<double, double>::getTimeBoundsPerSet)
        .function("getObservations", &tom::ObservationCollection<double, double>::getObservations)
        .function("getConcatenatedObservations", &tom::ObservationCollection<double, double>::getConcatenatedObservations)
        .function("getObservationTimesDouble", &tom::ObservationCollection<double, double>::getObservationTimesDouble)
        .function("getObservationTimes", &tom::ObservationCollection<double, double>::getObservationTimes)
        .function("getConcatenatedObservationTimes", &tom::ObservationCollection<double, double>::getConcatenatedObservationTimes)
        .function("getConcatenatedDoubleObservationTimes", &tom::ObservationCollection<double, double>::getConcatenatedDoubleObservationTimes)
        .function("getObservationsAndTimesDouble", &tom::ObservationCollection<double, double>::getObservationsAndTimesDouble)
        .function("getObservationsAndTimes", &tom::ObservationCollection<double, double>::getObservationsAndTimes)
        .function("getConcatenatedObservationsAndTimesDouble", &tom::ObservationCollection<double, double>::getConcatenatedObservationsAndTimesDouble)
        .function("getConcatenatedObservationsAndTimes", &tom::ObservationCollection<double, double>::getConcatenatedObservationsAndTimes)
        .function("getWeights", &tom::ObservationCollection<double, double>::getWeights)
        .function("getConcatenatedWeights", &tom::ObservationCollection<double, double>::getConcatenatedWeights)
        .function("getResiduals", &tom::ObservationCollection<double, double>::getResiduals)
        .function("getConcatenatedResiduals", &tom::ObservationCollection<double, double>::getConcatenatedResiduals)
        .function("getRmsResiduals", &tom::ObservationCollection<double, double>::getRmsResiduals)
        .function("getMeanResiduals", &tom::ObservationCollection<double, double>::getMeanResiduals)
        .function("getComputedObservations", &tom::ObservationCollection<double, double>::getComputedObservations)
        .function("getConcatenatedComputedObservations", &tom::ObservationCollection<double, double>::getConcatenatedComputedObservations)
        // Weight setting methods
        .function("setConstantWeightScalar",
            select_overload<void(const double, const std::shared_ptr<tom::ObservationCollectionParser>)>(
                &tom::ObservationCollection<double, double>::setConstantWeight))
        .function("setConstantWeightVector",
            select_overload<void(const Eigen::VectorXd, const std::shared_ptr<tom::ObservationCollectionParser>)>(
                &tom::ObservationCollection<double, double>::setConstantWeight))
        .function("setTabulatedWeightsWithParser",
            select_overload<void(const Eigen::VectorXd, const std::shared_ptr<tom::ObservationCollectionParser>)>(
                &tom::ObservationCollection<double, double>::setTabulatedWeights))
        // Collection manipulation
        .function("appendObservationCollection", &tom::ObservationCollection<double, double>::appendObservationCollection)
        .function("getSingleObservationSets", &tom::ObservationCollection<double, double>::getSingleObservationSets)
        .function("printObservationSetsStartAndSize", &tom::ObservationCollection<double, double>::printObservationSetsStartAndSize)
        .function("removeSingleObservationSets",
            select_overload<void(std::shared_ptr<tom::ObservationCollectionParser>)>(
                &tom::ObservationCollection<double, double>::removeSingleObservationSets))
        .function("removeEmptySingleObservationSets", &tom::ObservationCollection<double, double>::removeEmptySingleObservationSets)
        // Dependent variables
        .function("addDependentVariable", &tom::ObservationCollection<double, double>::addDependentVariable)
        .function("getDependentVariables", &tom::ObservationCollection<double, double>::getDependentVariables)
        .function("getConcatenatedDependentVariables", &tom::ObservationCollection<double, double>::getConcatenatedDependentVariables)
        .function("getCompatibleDependentVariablesSettingsList", &tom::ObservationCollection<double, double>::getCompatibleDependentVariablesSettingsList)
        .function("getAllCompatibleDependentVariables", &tom::ObservationCollection<double, double>::getAllCompatibleDependentVariables)
        .function("getDependentVariableHistoryPerObservationSetDouble", &tom::ObservationCollection<double, double>::getDependentVariableHistoryPerObservationSetDouble)
        .function("getDependentVariableHistoryPerObservationSet", &tom::ObservationCollection<double, double>::getDependentVariableHistoryPerObservationSet)
        .function("getDependentVariableHistoryDouble", &tom::ObservationCollection<double, double>::getDependentVariableHistoryDouble)
        .function("getDependentVariableHistory", &tom::ObservationCollection<double, double>::getDependentVariableHistory);

    // ========================================================================
    // ObservationCollectionParser class
    // ========================================================================
    class_<tom::ObservationCollectionParser>("estimation_observations_ObservationCollectionParser")
        .smart_ptr<std::shared_ptr<tom::ObservationCollectionParser>>(
            "shared_ptr_ObservationCollectionParser")
        .constructor<>();

    // ========================================================================
    // ObservationFilterBase class
    // ========================================================================
    class_<tom::ObservationFilterBase>("estimation_observations_ObservationFilterBase")
        .smart_ptr<std::shared_ptr<tom::ObservationFilterBase>>(
            "shared_ptr_ObservationFilterBase");

    // ========================================================================
    // ObservationSetSplitterBase class
    // ========================================================================
    class_<tom::ObservationSetSplitterBase>("estimation_observations_ObservationSetSplitterBase")
        .smart_ptr<std::shared_ptr<tom::ObservationSetSplitterBase>>(
            "shared_ptr_ObservationSetSplitterBase");

    // ========================================================================
    // Free functions
    // ========================================================================

    // compute_residuals_and_dependent_variables
    function("estimation_observations_compute_residuals_and_dependent_variables",
        &tss::computeResidualsAndDependentVariables<double, double>);

    // merge_observation_collections
    function("estimation_observations_merge_observation_collections",
        &tss::mergeObservationCollections<double, double>);

    // filter_observations for SingleObservationSet
    function("estimation_observations_filter_single_observation_set",
        select_overload<std::shared_ptr<tom::SingleObservationSet<double, double>>(
            const std::shared_ptr<tom::SingleObservationSet<double, double>>,
            const std::shared_ptr<tom::ObservationFilterBase>,
            const bool)>(&tom::filterObservations<double, double>));

    // split_observation_set
    function("estimation_observations_split_observation_set",
        select_overload<std::vector<std::shared_ptr<tom::SingleObservationSet<double, double>>>(
            const std::shared_ptr<tom::SingleObservationSet<double, double>>,
            const std::shared_ptr<tom::ObservationSetSplitterBase>,
            const bool)>(&tom::splitObservationSet<double, double>));

    // create_filtered_observation_collection (with parser)
    function("estimation_observations_create_filtered_observation_collection",
        select_overload<std::shared_ptr<tom::ObservationCollection<double, double>>(
            const std::shared_ptr<tom::ObservationCollection<double, double>>,
            const std::shared_ptr<tom::ObservationFilterBase>,
            const std::shared_ptr<tom::ObservationCollectionParser>)>(
                &tom::filterObservations<double, double>));

    // split_observation_collection
    function("estimation_observations_split_observation_collection",
        &tom::splitObservationSets<double, double>);

    // create_new_observation_collection
    function("estimation_observations_create_new_observation_collection",
        &tom::createNewObservationCollection<double, double>);
}

#endif
