/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_observations.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <map>
#include <sstream>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"
#include "tudat/simulation/estimation_setup/simulateObservationsLegacy.h"
#include "observations_processing/expose_observations_processing.h"
#include "observations_geometry/expose_observations_geometry.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
namespace te = tudat::ephemerides;

namespace
{

const char* legacyObservationDeprecationGuide =
        "https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-dataset-deprecation.html";

std::string getObservationApiReferenceLink( const std::string& replacementApi )
{
    const std::string apiAnchor = replacementApi.substr( 0, replacementApi.find( ' ' ) );
    return "https://py.api.tudat.space/en/latest/estimation/observations.html#tudatpy.estimation.observations." + apiAnchor;
}

std::string getSingleObservationSetReplacement( const std::string& memberName )
{
    const static std::map< std::string, std::string > replacements = {
        { "set_observations", "ObservationDataset.set_observations_for_set" },
        { "set_residuals", "ObservationDataset.set_residuals_for_set" },
        { "set_constant_weight", "ObservationDataset.set_constant_single_observation_scalar_weight_for_set" },
        { "set_tabulated_weights", "ObservationDataset.set_weight_vector_for_set" },
        { "filter_observations",
          "ObservationDataset.reject_observations, ObservationDataset.remove_observations or ObservationDataset.create_new_and_drop" },
        { "observable_type", "ObservationDataset.get_observation_set_metadata" },
        { "link_definition", "ObservationDataset.get_observation_set_metadata and ObservationDataset.get_link_definition" },
        { "reference_link_end", "ObservationDataset.get_observation_set_metadata" },
        { "number_of_observables", "ObservationDataset.number_of_observations_for_set" },
        { "single_observable_size", "ObservationDataset.get_observation_set_metadata" },
        { "total_observation_set_size", "ObservationDataset.total_scalar_size_for_set" },
        { "time_bounds", "ObservationDataset.time_bounds_for_set" },
        { "list_of_observations", "ObservationDataset.observations_for_set" },
        { "observation_times", "ObservationDataset.observation_times_for_set" },
        { "concatenated_observations", "ObservationDataset.observation_vector_for_set" },
        { "computed_observations", "ObservationDataset.computed_observations_for_set" },
        { "concatenated_computed_observations", "ObservationDataset.computed_observation_vector_for_set" },
        { "residuals", "ObservationDataset.residuals_for_set" },
        { "concatenated_residuals", "ObservationDataset.residual_vector_for_set" },
        { "rms_residuals", "ObservationDataset.rms_residuals_for_set" },
        { "mean_residuals", "ObservationDataset.mean_residuals_for_set" },
        { "weights", "ObservationDataset.weights_for_set" },
        { "concatenad_weights", "ObservationDataset.weight_vector_for_set" },
        { "weights_vector", "ObservationDataset.weight_vector_for_set" },
        { "dependent_variables", "ObservationDataset.dependent_variables_for_set" },
        { "dependent_variables_history", "ObservationDataset.dependent_variables_for_set" },
        { "observations_history", "ObservationDataset.observations_for_set" },
        { "ancillary_settings", "ObservationDataset.ancillary_settings" },
        { "filtered_observation_set", "ObservationDataset.create_new_and_drop" },
        { "number_filtered_observations", "ObservationDataset.observation_ids_matching_condition" },
        { "single_dependent_variable", "ObservationDataset.single_dependent_variable_for_set" },
        { "compatible_dependent_variable_settings", "ObservationDataset.compatible_dependent_variable_settings_for_set" },
        { "compatible_dependent_variables_list", "ObservationDataset.all_compatible_dependent_variables_for_set" },
        { "single_dependent_variable_history", "ObservationDataset.single_dependent_variable_for_set" },
        { "dependent_variables_matrix", "ObservationDataset.dependent_variables_for_set" }
    };

    const auto replacementIterator = replacements.find( memberName );
    return replacementIterator == replacements.end( ) ? "ObservationDataset" : replacementIterator->second;
}

std::string getObservationCollectionReplacement( const std::string& memberName )
{
    const static std::map< std::string, std::string > replacements = {
        { "concatenated_times", "ObservationDataset.ordered_flattened_observation_data.times" },
        { "concatenated_times_objects", "ObservationDataset.ordered_flattened_observation_data.times" },
        { "concatenated_weights", "ObservationDataset.ordered_flattened_observation_data.weight_vector" },
        { "concatenated_observations", "ObservationDataset.ordered_flattened_observation_data.observation_vector" },
        { "concatenated_link_definition_ids", "ObservationDataset.ordered_flattened_observation_data.set_ids" },
        { "link_definition_ids", "ObservationDataset.link_definition" },
        { "observable_type_start_index_and_size", "ObservationDataset.ordered_flattened_observation_data" },
        { "observation_set_start_index_and_size", "ObservationDataset.ordered_flattened_observation_data" },
        { "observation_vector_size", "ObservationDataset.total_scalar_size" },
        { "sorted_observation_sets", "ObservationDataset.observation_set_metadata" },
        { "filter_observations",
          "ObservationDataset.reject_observations, ObservationDataset.remove_observations or ObservationDataset.create_new_and_drop" },
        { "split_observation_sets", "ObservationDataset.create_new_and_keep plus add_observation_set_from_dataset" },
        { "create_new_observation_collection", "ObservationDataset.create_new_and_keep" },
        { "set_constant_weight", "ObservationDataset.set_constant_single_observation_scalar_weight" },
        { "set_tabulated_weights", "ObservationDataset.set_weight_vector_for_set" },
        { "get_observation_dataset", "ObservationDataset" }
    };

    const auto replacementIterator = replacements.find( memberName );
    return replacementIterator == replacements.end( ) ? "ObservationDataset" : replacementIterator->second;
}

void warnLegacyObservationInterface( const std::string& interfaceName, const std::string& replacementApi = "ObservationDataset" )
{
    const std::string message = interfaceName + " is deprecated and kept only for backwards compatibility. Use " + replacementApi +
            " instead. API reference: " + getObservationApiReferenceLink( replacementApi ) +
            ". Migration guide: " + legacyObservationDeprecationGuide;
    if( PyErr_WarnEx( PyExc_DeprecationWarning, message.c_str( ), 1 ) < 0 )
    {
        throw py::error_already_set( );
    }
}

py::object getLegacyAttributeWithWarning( const py::object& self,
                                          const std::string& attributeName,
                                          const std::string& interfaceName,
                                          const std::string& replacementApi )
{
    if( attributeName.rfind( "__", 0 ) != 0 )
    {
        warnLegacyObservationInterface( interfaceName + "." + attributeName, replacementApi );
    }
    PyObject* attribute = PyObject_GenericGetAttr( self.ptr( ), py::str( attributeName ).ptr( ) );
    if( attribute == nullptr )
    {
        throw py::error_already_set( );
    }
    return py::reinterpret_steal< py::object >( attribute );
}

Eigen::VectorXd castObservationSelectionConditionVectorLimit( const py::object& value )
{
    if( py::isinstance< py::float_ >( value ) || py::isinstance< py::int_ >( value ) )
    {
        return ( Eigen::VectorXd( 1 ) << value.cast< double >( ) ).finished( );
    }
    return value.cast< Eigen::VectorXd >( );
}

bool isProvided( const py::object& value )
{
    return !value.is_none( );
}

template< typename ObservationScalarType, typename TimeType >
struct ObservationWeightSelector {
    bool useSetId = false;
    unsigned int setId = 0;
    tom::ObservationSelectionCondition< ObservationScalarType, TimeType > condition =
            tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::all( );
};

template< typename ObservationScalarType, typename TimeType >
void addObservationWeightCondition( tom::ObservationSelectionCondition< ObservationScalarType, TimeType >& combinedCondition,
                                    bool& hasCondition,
                                    const tom::ObservationSelectionCondition< ObservationScalarType, TimeType >& condition )
{
    combinedCondition = hasCondition ? ( combinedCondition && condition ) : condition;
    hasCondition = true;
}

template< typename ObservationScalarType, typename TimeType >
ObservationWeightSelector< ObservationScalarType, TimeType > createObservationWeightSelector(
        const py::object& condition,
        const py::object& setId,
        const py::object& observableType,
        const py::object& linkDefinition,
        const py::object& linkEndType,
        const py::object& linkEnd,
        const py::object& timeBounds,
        const py::object& active,
        const py::object& rejected,
        const py::object& residualAbsoluteValueGreaterThan,
        const py::object& observationAbsoluteValueGreaterThan,
        const py::object& dependentVariableGreaterThan )
{
    const bool hasCondition = isProvided( condition );
    const bool hasSetId = isProvided( setId );
    const bool hasObservableType = isProvided( observableType );
    const bool hasLinkDefinition = isProvided( linkDefinition );
    const bool hasLinkEndType = isProvided( linkEndType );
    const bool hasLinkEnd = isProvided( linkEnd );
    const bool hasTimeBounds = isProvided( timeBounds );
    const bool hasActive = isProvided( active );
    const bool hasRejected = isProvided( rejected );
    const bool hasResidualLimit = isProvided( residualAbsoluteValueGreaterThan );
    const bool hasObservationLimit = isProvided( observationAbsoluteValueGreaterThan );
    const bool hasDependentVariableLimit = isProvided( dependentVariableGreaterThan );

    const bool hasNonSetSelector = hasObservableType || hasLinkDefinition || hasLinkEndType || hasLinkEnd || hasTimeBounds || hasActive ||
            hasRejected || hasResidualLimit || hasObservationLimit || hasDependentVariableLimit;

    if( hasCondition && ( hasSetId || hasNonSetSelector ) )
    {
        throw py::value_error( "condition cannot be combined with selector keywords." );
    }
    if( hasSetId && hasNonSetSelector )
    {
        throw py::value_error( "set_id cannot be combined with other selector keywords." );
    }
    if( hasActive && !active.cast< bool >( ) )
    {
        throw py::value_error( "active only accepts True." );
    }
    if( hasRejected && !rejected.cast< bool >( ) )
    {
        throw py::value_error( "rejected only accepts True." );
    }
    if( hasActive && hasRejected )
    {
        throw py::value_error( "active and rejected cannot both be selected." );
    }
    if( !( hasCondition || hasSetId || hasNonSetSelector ) )
    {
        throw py::value_error( "At least one observation selector must be provided." );
    }

    ObservationWeightSelector< ObservationScalarType, TimeType > selector;
    if( hasSetId )
    {
        selector.useSetId = true;
        selector.setId = setId.cast< unsigned int >( );
        return selector;
    }
    if( hasCondition )
    {
        selector.condition = condition.cast< tom::ObservationSelectionCondition< ObservationScalarType, TimeType > >( );
        return selector;
    }

    bool hasCombinedCondition = false;
    if( hasObservableType )
    {
        addObservationWeightCondition( selector.condition,
                                       hasCombinedCondition,
                                       tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::observableType(
                                               observableType.cast< tom::ObservableType >( ) ) );
    }
    if( hasLinkDefinition )
    {
        addObservationWeightCondition( selector.condition,
                                       hasCombinedCondition,
                                       tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::linkDefinition(
                                               linkDefinition.cast< tom::LinkDefinition >( ) ) );
    }
    if( hasLinkEndType )
    {
        addObservationWeightCondition( selector.condition,
                                       hasCombinedCondition,
                                       tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::linkEndType(
                                               linkEndType.cast< tom::LinkEndType >( ) ) );
    }
    if( hasLinkEnd )
    {
        const py::tuple linkEndTuple = linkEnd.cast< py::tuple >( );
        if( linkEndTuple.size( ) != 2 )
        {
            throw py::value_error( "link_end must be a tuple of (link_end_type, link_end_id)." );
        }
        addObservationWeightCondition(
                selector.condition,
                hasCombinedCondition,
                tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::linkEnd(
                        linkEndTuple[ 0 ].cast< tom::LinkEndType >( ), linkEndTuple[ 1 ].cast< tom::LinkEndId >( ) ) );
    }
    if( hasTimeBounds )
    {
        const py::tuple timeBoundsTuple = timeBounds.cast< py::tuple >( );
        if( timeBoundsTuple.size( ) != 2 )
        {
            throw py::value_error( "time_bounds must be a tuple of (start_time, end_time)." );
        }
        addObservationWeightCondition( selector.condition,
                                       hasCombinedCondition,
                                       tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::timeBounds(
                                               timeBoundsTuple[ 0 ].cast< TimeType >( ), timeBoundsTuple[ 1 ].cast< TimeType >( ) ) );
    }
    if( hasActive )
    {
        addObservationWeightCondition( selector.condition,
                                       hasCombinedCondition,
                                       tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::active( ) );
    }
    if( hasRejected )
    {
        addObservationWeightCondition( selector.condition,
                                       hasCombinedCondition,
                                       tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::rejected( ) );
    }
    if( hasResidualLimit )
    {
        addObservationWeightCondition(
                selector.condition,
                hasCombinedCondition,
                tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::residualAbsoluteValueGreaterThan(
                        castObservationSelectionConditionVectorLimit( residualAbsoluteValueGreaterThan ) ) );
    }
    if( hasObservationLimit )
    {
        addObservationWeightCondition(
                selector.condition,
                hasCombinedCondition,
                tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::observationAbsoluteValueGreaterThan(
                        castObservationSelectionConditionVectorLimit( observationAbsoluteValueGreaterThan ) ) );
    }
    if( hasDependentVariableLimit )
    {
        const py::tuple dependentVariableTuple = dependentVariableGreaterThan.cast< py::tuple >( );
        if( dependentVariableTuple.size( ) != 2 && dependentVariableTuple.size( ) != 3 )
        {
            throw py::value_error(
                    "dependent_variable_greater_than must be a tuple of (dependent_variable_settings, limit) or "
                    "(dependent_variable_settings, limit, return_first_compatible_settings)." );
        }
        const bool returnFirstCompatibleSettings =
                dependentVariableTuple.size( ) == 3 ? dependentVariableTuple[ 2 ].cast< bool >( ) : false;
        addObservationWeightCondition(
                selector.condition,
                hasCombinedCondition,
                tom::ObservationSelectionCondition< ObservationScalarType, TimeType >::dependentVariableGreaterThan(
                        dependentVariableTuple[ 0 ].cast< std::shared_ptr< tss::ObservationDependentVariableSettings > >( ),
                        castObservationSelectionConditionVectorLimit( dependentVariableTuple[ 1 ].cast< py::object >( ) ),
                        returnFirstCompatibleSettings ) );
    }
    return selector;
}

struct ObservationDocSpec {
    std::string signature;
    std::string returnType;
    std::string summary;
};

struct ObservationDocArgument {
    std::string name;
    std::string defaultValue;
};

std::string observationArgumentType( const std::string& argumentName )
{
    const static std::map< std::string, std::string > argumentTypes = {
        { "observable_type", "tudatpy.estimation.observable_models_setup.model_settings.ObservableType" },
        { "link_definition", "tudatpy.estimation.observable_models_setup.links.LinkDefinition" },
        { "observations", "list[numpy.ndarray]" },
        { "times", "list[tudatpy.astro.time_representation.Time] or list[float]" },
        { "reference_link_end", "tudatpy.estimation.observable_models_setup.links.LinkEndType" },
        { "dependent_variables", "list[numpy.ndarray], optional" },
        { "dependent_variable_bookkeeping", "object, optional" },
        { "ancillary_settings", "object, optional" },
        { "weights", "list[numpy.ndarray], optional" },
        { "residuals", "list[numpy.ndarray], optional" },
        { "sort_observations", "bool, optional" },
        { "erase_duplicate_observations", "bool, optional" },
        { "weight_settings", "tudatpy.estimation.observations.ObservationWeightSettings" },
        { "source_dataset", "tudatpy.estimation.observations.ObservationDataset" },
        { "source_set_id", "int" },
        { "set_id", "int" },
        { "condition", "tudatpy.estimation.observations.ObservationSelectionCondition" },
        { "reason", "str, optional" },
        { "include_rejected", "bool, optional" },
        { "include_inactive", "bool, optional" },
        { "weight", "float or numpy.ndarray" },
        { "weight_vector", "numpy.ndarray" },
        { "weight_matrix", "numpy.ndarray" },
        { "observation_id", "int" },
        { "row_observation_ids", "list[int]" },
        { "column_observation_ids", "list[int]" },
        { "weight_block", "numpy.ndarray" },
        { "row_components", "list[int], optional" },
        { "column_components", "list[int], optional" },
        { "make_symmetric", "bool, optional" },
        { "indices_to_remove", "list[int]" },
        { "print_warning", "bool, optional" },
        { "observation_vector", "numpy.ndarray" },
        { "residual_vector", "numpy.ndarray" },
        { "scalar_component_id", "int" },
        { "dependent_variable_settings",
          "tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings" },
        { "return_first_compatible_settings", "bool, optional" },
        { "dependent_variable_index_and_size", "tuple[int, int]" },
        { "link_definition_id", "int" },
        { "ancillary_settings_id", "int" },
        { "dependent_variable_layout_id", "int" },
        { "observation_dataset", "tudatpy.estimation.observations.ObservationDataset" },
        { "observation_simulators", "list[tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator]" },
        { "bodies", "tudatpy.dynamics.environment.SystemOfBodies" },
        { "viewer_index", "int" },
        { "component_index", "int" }
    };

    const auto typeIterator = argumentTypes.find( argumentName );
    return typeIterator == argumentTypes.end( ) ? "object" : typeIterator->second;
}

std::string observationArgumentDescription( const std::string& argumentName )
{
    const static std::map< std::string, std::string > argumentDescriptions = {
        { "observable_type", "Observable type stored in the observation set." },
        { "link_definition", "Link definition shared by the observations." },
        { "observations", "Vector-valued observations, one entry per observation event." },
        { "times", "Observation times, one entry per observation event." },
        { "reference_link_end", "Link end at which the observation time is defined." },
        { "dependent_variables", "Dependent-variable vectors, one entry per observation event." },
        { "dependent_variable_bookkeeping", "Bookkeeping that describes the dependent-variable vector layout." },
        { "ancillary_settings", "Ancillary settings associated with the observations." },
        { "weights", "Per-observation scalar-component weight vectors." },
        { "residuals", "Residual vectors, one entry per observation event." },
        { "sort_observations", "Whether observations should be sorted by time after insertion." },
        { "erase_duplicate_observations", "Whether duplicate times in the new set should be removed." },
        { "weight_settings", "Compact scalar, per-observation block or set-level block weight policy." },
        { "source_dataset", "Dataset containing the observation set to copy." },
        { "source_set_id", "Identifier of the set in ``source_dataset``." },
        { "set_id", "Observation set identifier." },
        { "condition", "Row-level condition selecting observations." },
        { "reason", "Text stored with each rejected row." },
        { "include_rejected", "Whether rejected observations should be included." },
        { "include_inactive", "Whether inactive or rejected observations should be included." },
        { "weight", "Weight value, vector or matrix to assign." },
        { "weight_vector", "Concatenated scalar-component weights for one set." },
        { "weight_matrix", "Dense weight matrix." },
        { "observation_id", "Observation row identifier." },
        { "row_observation_ids", "Observation row identifiers selecting block rows." },
        { "column_observation_ids", "Observation row identifiers selecting block columns." },
        { "weight_block", "Dense block to insert into the sparse weight matrix." },
        { "row_components", "Component indices selected from each row observation. Empty selects all components." },
        { "column_components", "Component indices selected from each column observation. Empty selects all components." },
        { "make_symmetric", "Whether to also store the transposed block." },
        { "indices_to_remove", "Indices within the selected observation set to remove." },
        { "print_warning", "Whether to print a warning when duplicate observations are removed." },
        { "observation_vector", "Concatenated observation vector for one set." },
        { "residual_vector", "Concatenated residual vector for one set." },
        { "scalar_component_id", "Scalar-component row identifier." },
        { "dependent_variable_settings", "Settings identifying the dependent variable to retrieve." },
        { "return_first_compatible_settings", "Whether to use the first compatible stored dependent variable." },
        { "dependent_variable_index_and_size", "Pair containing the stored dependent-variable start index and size." },
        { "link_definition_id", "Index of the link definition in the dataset registry." },
        { "ancillary_settings_id", "Index of the ancillary settings in the dataset registry." },
        { "dependent_variable_layout_id", "Index of the dependent-variable bookkeeping in the dataset registry." },
        { "observation_dataset", "Dataset containing observed values and metadata." },
        { "observation_simulators", "Observation simulators used to recompute observables represented by the dataset." },
        { "bodies", "Environment used by the observation simulators." },
        { "viewer_index", "Index of the selected observation inside the viewer." },
        { "component_index", "Scalar component index inside the vector-valued observation." }
    };

    const auto descriptionIterator = argumentDescriptions.find( argumentName );
    return descriptionIterator == argumentDescriptions.end( ) ? "Input argument." : descriptionIterator->second;
}

std::string replaceAll( std::string input, const std::string& oldText, const std::string& newText )
{
    std::size_t position = 0;
    while( ( position = input.find( oldText, position ) ) != std::string::npos )
    {
        input.replace( position, oldText.length( ), newText );
        position += newText.length( );
    }
    return input;
}

std::string observationSignatureType( const std::string& argumentName )
{
    if( argumentName == "dependent_variable_index_and_size" )
    {
        return "tuple";
    }
    std::string type = observationArgumentType( argumentName );
    type = replaceAll( type, ", optional", "" );
    type = replaceAll( type, " or ", " | " );
    return type;
}

std::string observationReturnSignatureType( const std::string& returnType )
{
    return replaceAll( returnType, " or ", " | " );
}

std::vector< ObservationDocArgument > observationDocArguments( const std::string& signature )
{
    std::vector< ObservationDocArgument > arguments;
    const std::size_t openParenthesis = signature.find( '(' );
    const std::size_t closeParenthesis = signature.rfind( ')' );
    if( openParenthesis == std::string::npos || closeParenthesis == std::string::npos || closeParenthesis <= openParenthesis + 1 )
    {
        return arguments;
    }

    std::string argumentList = signature.substr( openParenthesis + 1, closeParenthesis - openParenthesis - 1 );
    std::stringstream stream( argumentList );
    std::string argument;
    while( std::getline( stream, argument, ',' ) )
    {
        argument.erase( 0, argument.find_first_not_of( " \t" ) );
        argument.erase( argument.find_last_not_of( " \t" ) + 1 );
        const std::size_t defaultPosition = argument.find( '=' );
        std::string defaultValue;
        if( defaultPosition != std::string::npos )
        {
            defaultValue = argument.substr( defaultPosition + 1 );
            defaultValue.erase( 0, defaultValue.find_first_not_of( " \t" ) );
            defaultValue.erase( defaultValue.find_last_not_of( " \t" ) + 1 );
            argument = argument.substr( 0, defaultPosition );
            argument.erase( argument.find_last_not_of( " \t" ) + 1 );
        }
        if( !argument.empty( ) )
        {
            arguments.push_back( { argument, defaultValue } );
        }
    }
    return arguments;
}

std::string makeObservationSignature( const std::string& functionName, const ObservationDocSpec& docSpec )
{
    std::ostringstream signature;
    signature << functionName << "(";
    const std::vector< ObservationDocArgument > arguments = observationDocArguments( docSpec.signature );
    for( unsigned int i = 0; i < arguments.size( ); ++i )
    {
        if( i > 0 )
        {
            signature << ", ";
        }
        signature << arguments.at( i ).name << ": " << observationSignatureType( arguments.at( i ).name );
        if( !arguments.at( i ).defaultValue.empty( ) )
        {
            signature << " = " << arguments.at( i ).defaultValue;
        }
    }
    signature << ") -> " << observationReturnSignatureType( docSpec.returnType );
    return signature.str( );
}

std::string makeObservationDocString( const std::string& functionName, const ObservationDocSpec& docSpec )
{
    std::ostringstream doc;
    doc << makeObservationSignature( functionName, docSpec ) << "\n\n";
    doc << docSpec.summary << "\n\n";
    doc << "Parameters\n----------\n";
    const std::vector< ObservationDocArgument > arguments = observationDocArguments( docSpec.signature );
    if( arguments.empty( ) )
    {
        doc << "None\n    This function has no input arguments.\n";
    }
    else
    {
        for( const ObservationDocArgument& argument : arguments )
        {
            doc << argument.name << " : " << observationArgumentType( argument.name ) << "\n";
            doc << "    " << observationArgumentDescription( argument.name ) << "\n";
        }
    }
    doc << "\nReturns\n-------\n" << docSpec.returnType << "\n";
    if( docSpec.returnType == "None" )
    {
        doc << "    This function modifies the dataset in place or performs the requested action without returning a value.";
    }
    else
    {
        doc << "    " << docSpec.summary;
    }
    return doc.str( );
}

const char* observationDatasetDoc( const std::string& methodName )
{
    const static std::map< std::string, ObservationDocSpec > docSpecs = {
        { "add_observation_set",
          { "(observable_type, link_definition, observations, times, reference_link_end, dependent_variables=[], "
            "dependent_variable_bookkeeping=None, ancillary_settings=None, weights=[], residuals=[], sort_observations=False, "
            "erase_duplicate_observations=False)",
            "int",
            "Add a logical observation set and return its dataset set identifier." } },
        { "add_observation_set_with_weights",
          { "(observable_type, link_definition, observations, times, reference_link_end, weight_settings, dependent_variables=[], "
            "dependent_variable_bookkeeping=None, ancillary_settings=None, residuals=[])",
            "int",
            "Add a logical observation set and initialize its weights from settings." } },
        { "add_observation_set_from_dataset",
          { "(source_dataset, source_set_id)", "int", "Copy one observation set from another dataset." } },
        { "get_observation_set_metadata",
          { "(set_id)", "tudatpy.estimation.observations.ObservationSetMetadata", "Return metadata for one observation set." } },
        { "create_viewer",
          { "(condition)",
            "tudatpy.estimation.observations.ObservationDatasetViewer",
            "Create a read-only viewer over selected observations." } },
        { "create_new_and_keep",
          { "(condition)",
            "tudatpy.estimation.observations.ObservationDataset",
            "Create a new dataset containing selected observations." } },
        { "create_new_and_drop",
          { "(condition)",
            "tudatpy.estimation.observations.ObservationDataset",
            "Create a new dataset excluding selected observations." } },
        { "reject_observations", { "(condition, reason='')", "None", "Mark selected observations as rejected." } },
        { "restore_observations", { "(condition)", "None", "Restore selected observations to active status." } },
        { "remove_observations", { "(condition)", "None", "Physically remove selected observations." } },
        { "delete_rejected_observations", { "()", "None", "Physically remove all currently rejected observations." } },
        { "estimation_flattened_observation_data",
          { "(include_rejected=False)",
            "tudatpy.estimation.observations.FlattenedObservationData",
            "Return flattened data for estimation." } },
        { "computation_flattened_observation_data",
          { "(include_rejected=True)",
            "tudatpy.estimation.observations.FlattenedObservationData",
            "Return flattened data for recomputation." } },
        { "ordered_flattened_observation_data",
          { "(include_inactive=True)",
            "tudatpy.estimation.observations.FlattenedObservationData",
            "Return flattened data in ordered output order." } },
        { "rms_residuals_for_set", { "(set_id)", "numpy.ndarray", "Return root-mean-square residuals for one set." } },
        { "mean_residuals_for_set", { "(set_id)", "numpy.ndarray", "Return mean residuals for one set." } },
        { "set_constant_weight", { "(condition, weight)", "None", "Set one scalar weight for all observations matching a condition." } },
        { "set_constant_diagonal_weight",
          { "(condition, weight)", "None", "Set one diagonal weight vector for all matching observations." } },
        { "set_constant_matrix_weight",
          { "(condition, weight)", "None", "Set one dense observable-size weight matrix for all matching observations." } },
        { "set_weight_vector_for_set",
          { "(set_id, weight_vector)", "None", "Replace the concatenated scalar-component weight vector for one set." } },
        { "set_weight_matrix_for_set", { "(set_id, weight_matrix)", "None", "Store one full dense set-level weight matrix." } },
        { "has_weight_matrix_for_set", { "(set_id)", "bool", "Return whether one set has an explicitly stored set-level weight matrix." } },
        { "set_weight_matrix_for_observation",
          { "(observation_id, weight_matrix)", "None", "Store one dense observable-size weight matrix for an observation row." } },
        { "has_weight_matrix_for_observation",
          { "(observation_id)", "bool", "Return whether one observation row has an explicitly stored dense weight matrix." } },
        { "set_weight_block",
          { "(row_observation_ids, column_observation_ids, weight_block, row_components=[], column_components=[], make_symmetric=False)",
            "None",
            "Store an advanced dense weight block selected by observation ids." } },
        { "observation_row",
          { "(observation_id)", "tudatpy.estimation.observations.ObservationDatasetRow", "Return row metadata for one observation." } },
        { "scalar_component_row",
          { "(scalar_component_id)",
            "tudatpy.estimation.observations.ObservationScalarComponentRow",
            "Return metadata for one scalar component." } },
        { "observation_ids_for_set", { "(set_id)", "list[int]", "Return observation row identifiers belonging to one set." } },
        { "observation_ids_matching_condition",
          { "(condition)", "list[int]", "Return observation row identifiers selected by a condition." } },
        { "observations_for_set", { "(set_id)", "list[numpy.ndarray]", "Return all observation vectors in one set." } },
        { "observation_vector_for_set", { "(set_id)", "numpy.ndarray", "Return all observations in one set as a concatenated vector." } },
        { "computed_observations_for_set",
          { "(set_id)", "list[numpy.ndarray]", "Return computed observations for one set, defined as observed values minus residuals." } },
        { "computed_observation_vector_for_set",
          { "(set_id)", "numpy.ndarray", "Return computed observations for one set as a concatenated vector." } },
        { "observation_value", { "(observation_id)", "numpy.ndarray", "Return the vector-valued observation for one observation row." } },
        { "observation_times_for_set", { "(set_id)", "list", "Return all observation times in one set." } },
        { "observation_time",
          { "(observation_id)", "tudatpy.astro.time_representation.Time or float", "Return the time of one observation row." } },
        { "weights_for_set", { "(set_id)", "list[numpy.ndarray]", "Return all per-observation weight vectors in one set." } },
        { "weight_vector_for_set",
          { "(set_id)", "numpy.ndarray", "Return all scalar-component weights in one set as a concatenated vector." } },
        { "weight_value", { "(observation_id)", "numpy.ndarray", "Return the scalar-component weight vector for one observation row." } },
        { "weight_matrix_for_observation",
          { "(observation_id)", "numpy.ndarray", "Return the dense observable-size weight matrix for one observation row." } },
        { "weight_matrix_for_set", { "(set_id)", "numpy.ndarray", "Return the full dense weight matrix for one set." } },
        { "residuals_for_set", { "(set_id)", "list[numpy.ndarray]", "Return all residual vectors in one set." } },
        { "residual_vector_for_set", { "(set_id)", "numpy.ndarray", "Return all residuals in one set as a concatenated vector." } },
        { "residual_value", { "(observation_id)", "numpy.ndarray", "Return the residual vector for one observation row." } },
        { "dependent_variables_for_set",
          { "(set_id)", "list[numpy.ndarray]", "Return all dependent-variable vectors stored for one set." } },
        { "dependent_variables", { "(observation_id)", "numpy.ndarray", "Return the dependent-variable vector for one observation row." } },
        { "single_dependent_variable_for_set",
          { "(set_id, dependent_variable_settings, return_first_compatible_settings=False)",
            "numpy.ndarray",
            "Return values of one dependent variable stored for one set." } },
        { "single_dependent_variable_for_set_by_index",
          { "(set_id, dependent_variable_index_and_size)",
            "numpy.ndarray",
            "Return dependent-variable values by stored index and size." } },
        { "compatible_dependent_variable_settings_for_set",
          { "(set_id, dependent_variable_settings)", "list", "Return compatible dependent-variable settings in one set." } },
        { "all_compatible_dependent_variables_for_set",
          { "(set_id, dependent_variable_settings)",
            "list[numpy.ndarray]",
            "Return all compatible dependent-variable values in one set." } },
        { "set_dependent_variables_for_set",
          { "(set_id, dependent_variables)", "None", "Replace all dependent-variable vectors in one set." } },
        { "clear_dependent_variables_for_set", { "(set_id)", "None", "Clear all dependent-variable vectors in one set." } },
        { "set_observations_for_set", { "(set_id, observations)", "None", "Replace all observation vectors in one set." } },
        { "set_observation_vector_for_set",
          { "(set_id, observation_vector)", "None", "Replace all observations in one set from a vector." } },
        { "set_residuals_for_set", { "(set_id, residuals)", "None", "Replace all residual vectors in one set." } },
        { "set_residual_vector_for_set", { "(set_id, residual_vector)", "None", "Replace all residuals in one set from a vector." } },
        { "add_observations_to_set",
          { "(set_id, observations, times, dependent_variables=[], weights=[], residuals=[], sort_observations=True)",
            "None",
            "Append observations to an existing set." } },
        { "remove_observations_from_set", { "(set_id, indices_to_remove)", "None", "Remove observations from one set by index." } },
        { "erase_duplicate_observations_from_set",
          { "(set_id, print_warning=True)", "None", "Remove duplicate observations from one set." } },
        { "time_bounds_for_set", { "(set_id)", "tuple", "Return the minimum and maximum observation time in one set." } },
        { "number_of_observations_for_set", { "(set_id)", "int", "Return the number of observation rows in one set." } },
        { "total_scalar_size_for_set", { "(set_id)", "int", "Return the number of scalar components in one set." } },
        { "link_definition",
          { "(link_definition_id)",
            "tudatpy.estimation.observable_models_setup.links.LinkDefinition",
            "Return a link definition from the dataset registry." } },
        { "ancillary_settings", { "(ancillary_settings_id)", "object", "Return ancillary settings from the dataset registry." } },
        { "dependent_variable_bookkeeping",
          { "(dependent_variable_layout_id)", "object", "Return dependent-variable bookkeeping from the dataset registry." } },
        { "set_constant_single_observation_scalar_weight",
          { "(condition, weight)", "None", "Set one scalar weight for all observations matching a condition." } },
        { "set_constant_single_observation_diagonal_weight",
          { "(condition, weight)", "None", "Set one diagonal weight vector for all observations matching a condition." } },
        { "set_constant_single_observation_matrix_weight",
          { "(condition, weight)", "None", "Set one dense observable-size weight matrix for all observations matching a condition." } },
        { "set_constant_single_observation_scalar_weight_for_set",
          { "(set_id, weight)", "None", "Set one scalar weight for every scalar component in an observation set." } },
        { "set_constant_single_observation_diagonal_weight_for_set",
          { "(set_id, weight)", "None", "Set one diagonal weight vector for every observation in a set." } },
        { "set_constant_single_observation_matrix_weight_for_set",
          { "(set_id, weight)", "None", "Set one dense observable-size weight matrix for every observation in a set." } }
    };

    const static std::map< std::string, std::string > docStrings = [ & ]( ) {
        std::map< std::string, std::string > docs;
        for( const auto& docSpec : docSpecs )
        {
            docs[ docSpec.first ] = makeObservationDocString( docSpec.first, docSpec.second );
        }
        return docs;
    }( );

    return docStrings.at( methodName ).c_str( );
}

const char* observationModuleFunctionDoc( const std::string& functionName )
{
    const static std::map< std::string, ObservationDocSpec > docSpecs = {
        { "compute_residuals_and_dependent_variables_for_dataset",
          { "(observation_dataset, observation_simulators, bodies)",
            "None",
            "Compute simulated observations, residuals and dependent variables for a dataset." } },
        { "observation_simulation_settings_from_dataset",
          { "(observation_dataset, bodies)", "list", "Create observation simulation settings from a dataset." } }
    };

    const static std::map< std::string, std::string > docStrings = [ & ]( ) {
        std::map< std::string, std::string > docs;
        for( const auto& docSpec : docSpecs )
        {
            docs[ docSpec.first ] = makeObservationDocString( docSpec.first, docSpec.second );
        }
        return docs;
    }( );

    return docStrings.at( functionName ).c_str( );
}

const char* supportingObservationObjectDoc( const std::string& methodName )
{
    const static std::map< std::string, ObservationDocSpec > docSpecs = {
        { "viewer_observation_row",
          { "(viewer_index)",
            "tudatpy.estimation.observations.ObservationDatasetRow",
            "Return row metadata for one selected observation." } },
        { "viewer_observation_value",
          { "(viewer_index)", "numpy.ndarray", "Return the vector-valued observation at one selected viewer index." } },
        { "viewer_observation_time",
          { "(viewer_index)",
            "tudatpy.astro.time_representation.Time or float",
            "Return the observation time at one selected viewer index." } },
        { "viewer_create_viewer",
          { "(condition)",
            "tudatpy.estimation.observations.ObservationDatasetViewer",
            "Create a narrower read-only viewer selected from this viewer." } },
        { "viewer_estimation_flattened_observation_data",
          { "(include_rejected=False)",
            "tudatpy.estimation.observations.FlattenedObservationData",
            "Return flattened data for estimation." } },
        { "viewer_ordered_flattened_observation_data",
          { "(include_inactive=True)",
            "tudatpy.estimation.observations.FlattenedObservationData",
            "Return flattened data in ordered output order." } },
        { "unique_observation_ids_for_set",
          { "(set_id)", "list[int]", "Return unique observation row identifiers for one set in flattened-data row order." } },
        { "flattened_row",
          { "(observation_id, component_index)", "int", "Return the flattened scalar row for one observation row and component index." } }
    };

    const static std::map< std::string, std::string > docStrings = [ & ]( ) {
        std::map< std::string, std::string > docs;
        for( const auto& docSpec : docSpecs )
        {
            std::string publicName = docSpec.first;
            publicName = replaceAll( publicName, "viewer_", "" );
            docs[ docSpec.first ] = makeObservationDocString( publicName, docSpec.second );
        }
        return docs;
    }( );

    return docStrings.at( methodName ).c_str( );
}

}  // namespace

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< tom::SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSetWithoutDependentVariables(
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType > observationTimes,
        const tom::LinkEndType referenceLinkEnd,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
{
    warnLegacyObservationInterface( "single_observation_set", "ObservationDataset.add_observation_set" );
    return std::make_shared< tom::SingleObservationSet< ObservationScalarType, TimeType > >( observableType,
                                                                                             linkEnds,
                                                                                             observations,
                                                                                             observationTimes,
                                                                                             referenceLinkEnd,
                                                                                             std::vector< Eigen::VectorXd >( ),
                                                                                             nullptr,
                                                                                             ancillarySettings );
}

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace estimation
{
namespace observations
{

void expose_observations( py::module& m )
{
    auto observations_processing = m.def_submodule( "observations_processing" );
    observations_processing::expose_observations_processing( observations_processing );

    auto observations_geometry = m.def_submodule( "observations_geometry" );
    observations_geometry::expose_observations_geometry( observations_geometry );

    // OBSERVATION DATASET

    py::class_< tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                               "ObservationSetMetadata",
                                                                               R"doc(
Metadata describing one logical observation set stored in an :class:`ObservationDataset`.

The metadata identifies the observable type, link definition, reference link end,
observable size and the registered ancillary/dependent-variable layouts used by
all observations in the set.
)doc" )
            .def_readonly(
                    "observable_type",
                    &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::observableType_,
                    R"doc(tudatpy.estimation.observable_models_setup.model_settings.ObservableType: Observable type stored in this set.)doc" )
            .def_readonly( "link_definition_id",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::linkDefinitionId_,
                           R"doc(int: Index of the link definition in the dataset link-definition registry.)doc" )
            .def_readonly(
                    "reference_link_end",
                    &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::referenceLinkEnd_,
                    R"doc(tudatpy.estimation.observable_models_setup.links.LinkEndType: Reference link end used for all observations in this set.)doc" )
            .def_readonly( "observable_size",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::observableSize_,
                           R"doc(int: Number of scalar components in one observation of this set.)doc" )
            .def_readonly( "ancillary_settings_id",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::ancillarySettingsId_,
                           R"doc(int: Index of the ancillary settings in the dataset registry.)doc" )
            .def_readonly( "dependent_variable_layout_id",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::dependentVariableLayoutId_,
                           R"doc(int: Index of the dependent-variable bookkeeping layout in the dataset registry.)doc" );

    py::class_< tom::ObservationDatasetRow< TIME_TYPE > >( m,
                                                           "ObservationDatasetRow",
                                                           R"doc(
Row-level storage metadata for one observation inside an :class:`ObservationDataset`.

Each row points to the first scalar value of the observation in the dataset-wide
scalar-value storage and records the observation time, owning set and index
within that set.
)doc" )
            .def_readonly( "time",
                           &tom::ObservationDatasetRow< TIME_TYPE >::time_,
                           R"doc(tudatpy.astro.time_representation.Time or float: Observation time.)doc" )
            .def_readonly(
                    "set_id", &tom::ObservationDatasetRow< TIME_TYPE >::setId_, R"doc(int: Identifier of the owning observation set.)doc" )
            .def_readonly( "first_scalar_component",
                           &tom::ObservationDatasetRow< TIME_TYPE >::firstScalarComponent_,
                           R"doc(int: Index of this observation's first scalar value in the dataset-wide scalar-value storage.)doc" )
            .def_readonly( "scalar_size",
                           &tom::ObservationDatasetRow< TIME_TYPE >::scalarSize_,
                           R"doc(int: Observable size of this observation: the number of scalar values it contributes.)doc" )
            .def_readonly( "index_in_set",
                           &tom::ObservationDatasetRow< TIME_TYPE >::indexInSet_,
                           R"doc(int: Index of this observation within its observation set.)doc" )
            .def_readonly( "is_active",
                           &tom::ObservationDatasetRow< TIME_TYPE >::isActive_,
                           R"doc(bool: Whether this row is active in estimation/covariance flattened data.)doc" )
            .def_readonly( "rejection_reason",
                           &tom::ObservationDatasetRow< TIME_TYPE >::rejectionReason_,
                           R"doc(str: Optional text describing why this observation was rejected.)doc" );

    py::class_< tom::ObservationScalarComponentRow >( m,
                                                      "ObservationScalarComponentRow",
                                                      R"doc(
Storage metadata for one scalar component of an observation.

The row records the owning observation and component index inside that
observation.
)doc" )
            .def_readonly( "observation_id",
                           &tom::ObservationScalarComponentRow::observationId_,
                           R"doc(int: Identifier of the owning observation.)doc" )
            .def_readonly( "component_index",
                           &tom::ObservationScalarComponentRow::componentIndex_,
                           R"doc(int: Component index within the owning vector-valued observation.)doc" );

    py::class_< tom::ObservationWeightBlock >( m,
                                               "ObservationWeightBlock",
                                               R"doc(
Advanced dense weight block over selected scalar components.

This type is intended for rare off-diagonal correlations that are not naturally
represented as per-observation weights or as a full set-level weight block.
)doc" )
            .def( py::init<>( ), R"doc(Create an empty observation weight block.)doc" )
            .def_readwrite( "row_scalar_component_ids",
                            &tom::ObservationWeightBlock::rowScalarComponentIds_,
                            R"doc(Scalar component ids corresponding to the block rows.)doc" )
            .def_readwrite( "column_scalar_component_ids",
                            &tom::ObservationWeightBlock::columnScalarComponentIds_,
                            R"doc(Scalar component ids corresponding to the block columns.)doc" )
            .def_readwrite( "weight_block",
                            &tom::ObservationWeightBlock::weightBlock_,
                            R"doc(Dense weight block for the selected scalar components.)doc" );

    py::class_< tom::ObservationWeightSettings >( m,
                                                  "ObservationWeightSettings",
                                                  R"doc(
Weight policy used when adding a new observation set.

Use the static constructors to request compact scalar weights, per-observation
scalar weights, per-observation matrix blocks or a full set-level block.
)doc" )
            .def( py::init<>( ), R"doc(Create settings for default unit weights.)doc" )
            .def_static( "default_weights",
                         &tom::ObservationWeightSettings::defaultWeights,
                         R"doc(Return settings for default unit weights.)doc" )
            .def_static( "constant_scalar",
                         &tom::ObservationWeightSettings::constantScalar,
                         py::arg( "weight" ),
                         R"doc(Return settings using one scalar weight for every observation.)doc" )
            .def_static( "scalar_per_observation",
                         &tom::ObservationWeightSettings::scalarPerObservation,
                         py::arg( "weights" ),
                         R"doc(Return settings using one scalar weight per observation.)doc" )
            .def_static( "constant_block",
                         &tom::ObservationWeightSettings::constantBlock,
                         py::arg( "weight_block" ),
                         R"doc(Return settings using one observable-size matrix block for every observation.)doc" )
            .def_static( "block_per_observation",
                         &tom::ObservationWeightSettings::blockPerObservation,
                         py::arg( "weight_blocks" ),
                         R"doc(Return settings using one observable-size matrix block per observation.)doc" )
            .def_static( "set_block",
                         &tom::ObservationWeightSettings::setBlock,
                         py::arg( "weight_block" ),
                         R"doc(Return settings using one full set-level matrix block.)doc" );

    {
        py::options flattenedObservationDataDocOptions;
        flattenedObservationDataDocOptions.disable_function_signatures( );

        py::class_< tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                     "FlattenedObservationData",
                                                                                     R"doc(
Flattened vector data created from an :class:`ObservationDataset`.

This object contains the concatenated observation, residual and weight vectors,
together with the scalar-component provenance needed to map each entry back to
a dataset observation row and observation set. The diagonal weights are always
available through :attr:`weight_vector`. The full matrix is returned as a sparse
matrix and is only needed when off-diagonal terms are present.
)doc" )
                .def_property_readonly( "observation_vector",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationVector,
                                        R"doc(numpy.ndarray: Concatenated vector of observed values.)doc" )
                .def_property_readonly( "residual_vector",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualVector,
                                        R"doc(numpy.ndarray: Concatenated vector of residual values.)doc" )
                .def_property_readonly( "weight_vector",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightVector,
                                        R"doc(numpy.ndarray: Concatenated vector of scalar observation weights.)doc" )
                .def_property_readonly( "weight_matrix",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightMatrix,
                                        R"doc(
scipy.sparse.spmatrix: Sparse weight matrix in the same order as :attr:`observation_vector`.

For diagonal-only weights this matrix is generated from :attr:`weight_vector`.
For off-diagonal weights it contains the materialized sparse matrix assembled
from per-observation blocks, set-level blocks and advanced scalar-component
blocks.
)doc" )
                .def_property_readonly( "sparse_weight_matrix",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getSparseWeightMatrix,
                                        R"doc(
scipy.sparse.spmatrix: Sparse weight matrix in the same row order as :attr:`observation_vector`.

For diagonal-only flattened data, prefer :attr:`weight_vector`; requesting this
property materializes the sparse diagonal matrix.
)doc" )
                .def_property_readonly( "is_diagonal_weight_only",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::isDiagonalWeightOnly,
                                        R"doc(bool: True when the weight matrix contains no off-diagonal entries.)doc" )
                .def_property_readonly( "has_off_diagonal_weights",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::hasOffDiagonalWeights,
                                        R"doc(bool: True when the weight matrix contains off-diagonal entries.)doc" )
                .def_property_readonly(
                        "times",
                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getTimes,
                        R"doc(list[tudatpy.astro.time_representation.Time] or list[float]: Observation time associated with each scalar component.)doc" )
                .def_property_readonly( "observation_ids",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationIds,
                                        R"doc(list[int]: Observation row identifier associated with each scalar component.)doc" )
                .def_property_readonly( "set_ids",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getSetIds,
                                        R"doc(list[int]: Observation set identifier associated with each scalar component.)doc" )
                .def_property_readonly( "scalar_component_ids",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getScalarComponentIds,
                                        R"doc(list[int]: Scalar-component row identifier for each flattened scalar entry.)doc" )
                .def_property_readonly( "set_ids_in_row_order",
                                        &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getSetIdsInRowOrder,
                                        R"doc(list[int]: Unique observation set identifiers in the order in which they first appear.)doc" )
                .def( "unique_observation_ids_for_set",
                      &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getUniqueObservationIdsForSetInRowOrder,
                      py::arg( "set_id" ),
                      supportingObservationObjectDoc( "unique_observation_ids_for_set" ) )
                .def( "flattened_row",
                      &tom::FlattenedObservationData< STATE_SCALAR_TYPE, TIME_TYPE >::getFlattenedRow,
                      py::arg( "observation_id" ),
                      py::arg( "component_index" ),
                      supportingObservationObjectDoc( "flattened_row" ) );
    }

    py::enum_< tom::ObservationSelectionConditionType >(
            m, "ObservationSelectionConditionType", R"doc(Type of an observation-selection condition node.)doc" )
            .value( "all", tom::ObservationSelectionConditionType::all )
            .value( "observable_type", tom::ObservationSelectionConditionType::observable_type )
            .value( "link_definition", tom::ObservationSelectionConditionType::link_definition )
            .value( "link_end_type", tom::ObservationSelectionConditionType::link_end_type )
            .value( "link_end", tom::ObservationSelectionConditionType::link_end )
            .value( "set_id", tom::ObservationSelectionConditionType::set_id )
            .value( "time_bounds", tom::ObservationSelectionConditionType::time_bounds )
            .value( "time_greater_equal", tom::ObservationSelectionConditionType::time_greater_equal )
            .value( "time_greater_than", tom::ObservationSelectionConditionType::time_greater_than )
            .value( "time_less_equal", tom::ObservationSelectionConditionType::time_less_equal )
            .value( "time_less_than", tom::ObservationSelectionConditionType::time_less_than )
            .value( "active", tom::ObservationSelectionConditionType::active )
            .value( "residual_absolute_value_greater_than", tom::ObservationSelectionConditionType::residual_absolute_value_greater_than )
            .value( "observation_absolute_value_greater_than",
                    tom::ObservationSelectionConditionType::observation_absolute_value_greater_than )
            .value( "dependent_variable_greater_than", tom::ObservationSelectionConditionType::dependent_variable_greater_than )
            .value( "and_condition", tom::ObservationSelectionConditionType::and_condition )
            .value( "or_condition", tom::ObservationSelectionConditionType::or_condition )
            .value( "not_condition", tom::ObservationSelectionConditionType::not_condition )
            .value( "custom", tom::ObservationSelectionConditionType::custom );

    py::class_< tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                      "ObservationSelectionCondition",
                                                                                      R"doc(
Composable row-level condition used to select observations in an ObservationDataset.

Conditions operate on individual observation rows, not complete observation
sets. Combine conditions with ``&`` and ``|`` and negate them with ``~``. The
condition stores an inspectable query tree for conditions created through the
public builders.
)doc" )
            .def( py::init<>( ), R"doc(Create a condition that selects all observations.)doc" )
            .def_property_readonly( "condition_type",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getConditionType,
                                    R"doc(Type of this condition node.)doc" )
            .def_property_readonly( "condition_type_string",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getConditionTypeString,
                                    R"doc(String name of this condition node type.)doc" )
            .def_property_readonly( "child_conditions",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getChildConditions,
                                    py::return_value_policy::reference_internal,
                                    R"doc(Child conditions for logical AND, OR and NOT nodes.)doc" )
            .def_property_readonly( "observable_type_value",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableType,
                                    R"doc(Observable type stored by observable-type condition nodes.)doc" )
            .def_property_readonly( "link_definition_value",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinition,
                                    R"doc(Link definition stored by link-definition condition nodes.)doc" )
            .def_property_readonly( "link_end_type_value",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEndType,
                                    R"doc(Link-end type stored by link-end condition nodes.)doc" )
            .def_property_readonly( "link_end_id_value",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEndId,
                                    R"doc(Link-end identifier stored by link-end condition nodes.)doc" )
            .def_property_readonly( "set_id_value",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getSetId,
                                    R"doc(Set identifier stored by set-id condition nodes.)doc" )
            .def_property_readonly( "time_bounds_value",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBounds,
                                    R"doc(Inclusive time bounds stored by time-bound condition nodes.)doc" )
            .def_property_readonly( "time_value",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeValue,
                                    R"doc(Time value stored by one-sided time-comparison condition nodes.)doc" )
            .def_property_readonly( "vector_limit",
                                    &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::getVectorLimit,
                                    R"doc(Component limit stored by value-threshold condition nodes.)doc" )
            .def_static( "all",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::all,
                         R"doc(Return a condition that selects all observations.)doc" )
            .def_static( "observable_type",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::observableType,
                         py::arg( "observable_type" ),
                         R"doc(Return a condition selecting observations of one observable type.)doc" )
            .def_static( "link_definition",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::linkDefinition,
                         py::arg( "link_definition" ),
                         R"doc(Return a condition selecting observations with a matching link definition.)doc" )
            .def_static( "link_end_type",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::linkEndType,
                         py::arg( "link_end_type" ),
                         R"doc(Return a condition selecting observations whose link definition contains a link-end type.)doc" )
            .def_static( "link_end",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::linkEnd,
                         py::arg( "link_end_type" ),
                         py::arg( "link_end_id" ),
                         R"doc(Return a condition selecting observations with a specific link end.)doc" )
            .def_static( "set_id",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::setId,
                         py::arg( "set_id" ),
                         R"doc(Return a condition selecting observations from one observation set.)doc" )
            .def_static( "time_bounds",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::timeBounds,
                         py::arg( "start_time" ),
                         py::arg( "end_time" ),
                         R"doc(Return a condition selecting rows with start_time <= observation_time <= end_time.)doc" )
            .def_static( "time_greater_equal",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::timeGreaterEqual,
                         py::arg( "time" ),
                         R"doc(Return a condition selecting rows with observation_time >= time.)doc" )
            .def_static( "time_greater_than",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::timeGreaterThan,
                         py::arg( "time" ),
                         R"doc(Return a condition selecting rows with observation_time > time.)doc" )
            .def_static( "time_less_equal",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::timeLessEqual,
                         py::arg( "time" ),
                         R"doc(Return a condition selecting rows with observation_time <= time.)doc" )
            .def_static( "time_less_than",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::timeLessThan,
                         py::arg( "time" ),
                         R"doc(Return a condition selecting rows with observation_time < time.)doc" )
            .def_static( "active",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::active,
                         R"doc(Return a condition selecting active observations.)doc" )
            .def_static( "rejected",
                         &tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::rejected,
                         R"doc(Return a condition selecting rejected observations.)doc" )
            .def_static(
                    "residual_absolute_value_greater_than",
                    []( const py::object& limit ) {
                        return tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::residualAbsoluteValueGreaterThan(
                                castObservationSelectionConditionVectorLimit( limit ) );
                    },
                    py::arg( "limit" ),
                    R"doc(Return a condition selecting rows where any residual component exceeds the supplied absolute limit.)doc" )
            .def_static(
                    "observation_absolute_value_greater_than",
                    []( const py::object& limit ) {
                        return tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::observationAbsoluteValueGreaterThan(
                                castObservationSelectionConditionVectorLimit( limit ) );
                    },
                    py::arg( "limit" ),
                    R"doc(Return a condition selecting rows where any observation component exceeds the supplied absolute limit.)doc" )
            .def_static(
                    "dependent_variable_greater_than",
                    []( const std::shared_ptr< tss::ObservationDependentVariableSettings >& dependentVariableSettings,
                        const py::object& limit,
                        const bool returnFirstCompatibleSettings ) {
                        return tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::dependentVariableGreaterThan(
                                dependentVariableSettings,
                                castObservationSelectionConditionVectorLimit( limit ),
                                returnFirstCompatibleSettings );
                    },
                    py::arg( "dependent_variable_settings" ),
                    py::arg( "limit" ),
                    py::arg( "return_first_compatible_settings" ) = false,
                    R"doc(Return a condition selecting rows where a compatible dependent-variable component exceeds the limit.)doc" )
            .def(
                    "__and__",
                    []( const tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >& lhs,
                        const tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >& rhs ) { return lhs && rhs; },
                    py::is_operator( ),
                    R"doc(Return the logical AND of two conditions.)doc" )
            .def(
                    "__or__",
                    []( const tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >& lhs,
                        const tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >& rhs ) { return lhs || rhs; },
                    py::is_operator( ),
                    R"doc(Return the logical OR of two conditions.)doc" )
            .def(
                    "__invert__",
                    []( const tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >& condition ) { return !condition; },
                    py::is_operator( ),
                    R"doc(Return the logical negation of a condition.)doc" )
            .def(
                    "__bool__",
                    []( const tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >& ) {
                        throw py::type_error(
                                "ObservationSelectionCondition objects cannot be converted to bool. Use &, |, and ~ to combine conditions, "
                                "not and/or/not." );
                    },
                    R"doc(Always raise; use &, | and ~ instead of and/or/not.)doc" );

    {
        py::options observationDatasetViewerDocOptions;
        observationDatasetViewerDocOptions.disable_function_signatures( );

        py::class_< tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                     "ObservationDatasetViewer",
                                                                                     R"doc(
Read-only view on a selected subset of an ObservationDataset.

The viewer stores observation row identifiers selected from a parent dataset and
exposes only inspection and flattened-data methods. It is invalidated if the parent
dataset is structurally modified.
)doc" )
                .def_property_readonly( "number_of_observations",
                                        &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservations,
                                        R"doc(int: Number of selected observation rows.)doc" )
                .def_property_readonly( "observation_ids",
                                        &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationIds,
                                        R"doc(list[int]: Selected observation row identifiers.)doc" )
                .def( "observation_row",
                      &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationRow,
                      py::arg( "viewer_index" ),
                      supportingObservationObjectDoc( "viewer_observation_row" ) )
                .def( "observation_value",
                      &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationValue,
                      py::arg( "viewer_index" ),
                      supportingObservationObjectDoc( "viewer_observation_value" ) )
                .def( "observation_time",
                      &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTime,
                      py::arg( "viewer_index" ),
                      supportingObservationObjectDoc( "viewer_observation_time" ) )
                .def( "create_viewer",
                      &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::createViewer,
                      py::keep_alive< 0, 1 >( ),
                      py::arg( "condition" ),
                      supportingObservationObjectDoc( "viewer_create_viewer" ) )
                .def( "estimation_flattened_observation_data",
                      &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::createEstimationFlattenedObservationData,
                      py::arg( "include_rejected" ) = false,
                      supportingObservationObjectDoc( "viewer_estimation_flattened_observation_data" ) )
                .def( "ordered_flattened_observation_data",
                      &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::createOrderedFlattenedObservationData,
                      py::arg( "include_inactive" ) = true,
                      supportingObservationObjectDoc( "viewer_ordered_flattened_observation_data" ) );
    }

    {
        py::options observationDatasetDocOptions;
        observationDatasetDocOptions.disable_function_signatures( );

        py::class_< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >,
                    std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                  "ObservationDataset",
                                                                                                  R"doc(
Backend storage object for observation data.

An ``ObservationDataset`` stores observation values, residuals, weights,
dependent variables, set metadata and link/ancillary registries in a single
dataset-centric representation.
)doc" )
                .def( py::init<>( ), R"doc(Create an empty observation dataset.)doc" )
                .def( "add_observation_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSet,
                      py::arg( "observable_type" ),
                      py::arg( "link_definition" ),
                      py::arg( "observations" ),
                      py::arg( "times" ),
                      py::arg( "reference_link_end" ),
                      py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                      py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                      py::arg( "ancillary_settings" ) = nullptr,
                      py::arg( "weights" ) = std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                      py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                      py::arg( "sort_observations" ) = false,
                      py::arg( "erase_duplicate_observations" ) = false,
                      observationDatasetDoc( "add_observation_set" ) )
                .def( "add_observation_set_from_dataset",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSetFromDataset,
                      py::arg( "source_dataset" ),
                      py::arg( "source_set_id" ),
                      observationDatasetDoc( "add_observation_set_from_dataset" ) )
                .def( "add_observation_set_with_weights",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSetWithWeights,
                      py::arg( "observable_type" ),
                      py::arg( "link_definition" ),
                      py::arg( "observations" ),
                      py::arg( "times" ),
                      py::arg( "reference_link_end" ),
                      py::arg( "weight_settings" ),
                      py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                      py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                      py::arg( "ancillary_settings" ) = nullptr,
                      py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                      observationDatasetDoc( "add_observation_set_with_weights" ) )
                .def( "set_observations_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationsForSet,
                      py::arg( "set_id" ),
                      py::arg( "observations" ),
                      observationDatasetDoc( "set_observations_for_set" ) )
                .def( "set_observation_vector_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationVectorForSet,
                      py::arg( "set_id" ),
                      py::arg( "observation_vector" ),
                      observationDatasetDoc( "set_observation_vector_for_set" ) )
                .def( "set_residuals_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setResidualsForSet,
                      py::arg( "set_id" ),
                      py::arg( "residuals" ),
                      observationDatasetDoc( "set_residuals_for_set" ) )
                .def( "set_residual_vector_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setResidualVectorForSet,
                      py::arg( "set_id" ),
                      py::arg( "residual_vector" ),
                      observationDatasetDoc( "set_residual_vector_for_set" ) )
                .def( "set_constant_single_observation_scalar_weight_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationScalarWeightForSet,
                      py::arg( "set_id" ),
                      py::arg( "weight" ),
                      observationDatasetDoc( "set_constant_single_observation_scalar_weight_for_set" ) )
                .def( "set_constant_single_observation_diagonal_weight_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationDiagonalWeightForSet,
                      py::arg( "set_id" ),
                      py::arg( "weight" ),
                      observationDatasetDoc( "set_constant_single_observation_diagonal_weight_for_set" ) )
                .def( "set_constant_single_observation_matrix_weight_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationMatrixWeightForSet,
                      py::arg( "set_id" ),
                      py::arg( "weight" ),
                      observationDatasetDoc( "set_constant_single_observation_matrix_weight_for_set" ) )
                .def( "set_constant_single_observation_scalar_weight",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationScalarWeight,
                      py::arg( "condition" ),
                      py::arg( "weight" ),
                      observationDatasetDoc( "set_constant_single_observation_scalar_weight" ) )
                .def(
                        "set_constant_single_observation_scalar_weight",
                        []( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& dataset,
                            const double weight,
                            const py::object& condition,
                            const py::object& setId,
                            const py::object& observableType,
                            const py::object& linkDefinition,
                            const py::object& linkEndType,
                            const py::object& linkEnd,
                            const py::object& timeBounds,
                            const py::object& active,
                            const py::object& rejected,
                            const py::object& residualAbsoluteValueGreaterThan,
                            const py::object& observationAbsoluteValueGreaterThan,
                            const py::object& dependentVariableGreaterThan ) {
                            const auto selector =
                                    createObservationWeightSelector< STATE_SCALAR_TYPE, TIME_TYPE >( condition,
                                                                                                     setId,
                                                                                                     observableType,
                                                                                                     linkDefinition,
                                                                                                     linkEndType,
                                                                                                     linkEnd,
                                                                                                     timeBounds,
                                                                                                     active,
                                                                                                     rejected,
                                                                                                     residualAbsoluteValueGreaterThan,
                                                                                                     observationAbsoluteValueGreaterThan,
                                                                                                     dependentVariableGreaterThan );
                            if( selector.useSetId )
                            {
                                dataset.setConstantSingleObservationScalarWeightForSet( selector.setId, weight );
                            }
                            else
                            {
                                dataset.setConstantSingleObservationScalarWeight( selector.condition, weight );
                            }
                        },
                        py::arg( "weight" ),
                        py::kw_only( ),
                        py::arg( "condition" ) = py::none( ),
                        py::arg( "set_id" ) = py::none( ),
                        py::arg( "observable_type" ) = py::none( ),
                        py::arg( "link_definition" ) = py::none( ),
                        py::arg( "link_end_type" ) = py::none( ),
                        py::arg( "link_end" ) = py::none( ),
                        py::arg( "time_bounds" ) = py::none( ),
                        py::arg( "active" ) = py::none( ),
                        py::arg( "rejected" ) = py::none( ),
                        py::arg( "residual_absolute_value_greater_than" ) = py::none( ),
                        py::arg( "observation_absolute_value_greater_than" ) = py::none( ),
                        py::arg( "dependent_variable_greater_than" ) = py::none( ) )
                .def( "set_constant_weight",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationScalarWeight,
                      py::arg( "condition" ),
                      py::arg( "weight" ),
                      observationDatasetDoc( "set_constant_weight" ) )
                .def( "set_constant_single_observation_diagonal_weight",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationDiagonalWeight,
                      py::arg( "condition" ),
                      py::arg( "weight" ),
                      observationDatasetDoc( "set_constant_single_observation_diagonal_weight" ) )
                .def(
                        "set_constant_single_observation_diagonal_weight",
                        []( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& dataset,
                            const Eigen::VectorXd& weight,
                            const py::object& condition,
                            const py::object& setId,
                            const py::object& observableType,
                            const py::object& linkDefinition,
                            const py::object& linkEndType,
                            const py::object& linkEnd,
                            const py::object& timeBounds,
                            const py::object& active,
                            const py::object& rejected,
                            const py::object& residualAbsoluteValueGreaterThan,
                            const py::object& observationAbsoluteValueGreaterThan,
                            const py::object& dependentVariableGreaterThan ) {
                            const auto selector =
                                    createObservationWeightSelector< STATE_SCALAR_TYPE, TIME_TYPE >( condition,
                                                                                                     setId,
                                                                                                     observableType,
                                                                                                     linkDefinition,
                                                                                                     linkEndType,
                                                                                                     linkEnd,
                                                                                                     timeBounds,
                                                                                                     active,
                                                                                                     rejected,
                                                                                                     residualAbsoluteValueGreaterThan,
                                                                                                     observationAbsoluteValueGreaterThan,
                                                                                                     dependentVariableGreaterThan );
                            if( selector.useSetId )
                            {
                                dataset.setConstantSingleObservationDiagonalWeightForSet( selector.setId, weight );
                            }
                            else
                            {
                                dataset.setConstantSingleObservationDiagonalWeight( selector.condition, weight );
                            }
                        },
                        py::arg( "weight" ),
                        py::kw_only( ),
                        py::arg( "condition" ) = py::none( ),
                        py::arg( "set_id" ) = py::none( ),
                        py::arg( "observable_type" ) = py::none( ),
                        py::arg( "link_definition" ) = py::none( ),
                        py::arg( "link_end_type" ) = py::none( ),
                        py::arg( "link_end" ) = py::none( ),
                        py::arg( "time_bounds" ) = py::none( ),
                        py::arg( "active" ) = py::none( ),
                        py::arg( "rejected" ) = py::none( ),
                        py::arg( "residual_absolute_value_greater_than" ) = py::none( ),
                        py::arg( "observation_absolute_value_greater_than" ) = py::none( ),
                        py::arg( "dependent_variable_greater_than" ) = py::none( ) )
                .def( "set_constant_diagonal_weight",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationDiagonalWeight,
                      py::arg( "condition" ),
                      py::arg( "weight" ),
                      observationDatasetDoc( "set_constant_diagonal_weight" ) )
                .def( "set_constant_single_observation_matrix_weight",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationMatrixWeight,
                      py::arg( "condition" ),
                      py::arg( "weight" ),
                      observationDatasetDoc( "set_constant_single_observation_matrix_weight" ) )
                .def(
                        "set_constant_single_observation_matrix_weight",
                        []( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& dataset,
                            const Eigen::MatrixXd& weight,
                            const py::object& condition,
                            const py::object& setId,
                            const py::object& observableType,
                            const py::object& linkDefinition,
                            const py::object& linkEndType,
                            const py::object& linkEnd,
                            const py::object& timeBounds,
                            const py::object& active,
                            const py::object& rejected,
                            const py::object& residualAbsoluteValueGreaterThan,
                            const py::object& observationAbsoluteValueGreaterThan,
                            const py::object& dependentVariableGreaterThan ) {
                            const auto selector =
                                    createObservationWeightSelector< STATE_SCALAR_TYPE, TIME_TYPE >( condition,
                                                                                                     setId,
                                                                                                     observableType,
                                                                                                     linkDefinition,
                                                                                                     linkEndType,
                                                                                                     linkEnd,
                                                                                                     timeBounds,
                                                                                                     active,
                                                                                                     rejected,
                                                                                                     residualAbsoluteValueGreaterThan,
                                                                                                     observationAbsoluteValueGreaterThan,
                                                                                                     dependentVariableGreaterThan );
                            if( selector.useSetId )
                            {
                                dataset.setConstantSingleObservationMatrixWeightForSet( selector.setId, weight );
                            }
                            else
                            {
                                dataset.setConstantSingleObservationMatrixWeight( selector.condition, weight );
                            }
                        },
                        py::arg( "weight" ),
                        py::kw_only( ),
                        py::arg( "condition" ) = py::none( ),
                        py::arg( "set_id" ) = py::none( ),
                        py::arg( "observable_type" ) = py::none( ),
                        py::arg( "link_definition" ) = py::none( ),
                        py::arg( "link_end_type" ) = py::none( ),
                        py::arg( "link_end" ) = py::none( ),
                        py::arg( "time_bounds" ) = py::none( ),
                        py::arg( "active" ) = py::none( ),
                        py::arg( "rejected" ) = py::none( ),
                        py::arg( "residual_absolute_value_greater_than" ) = py::none( ),
                        py::arg( "observation_absolute_value_greater_than" ) = py::none( ),
                        py::arg( "dependent_variable_greater_than" ) = py::none( ) )
                .def( "set_constant_matrix_weight",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationMatrixWeight,
                      py::arg( "condition" ),
                      py::arg( "weight" ),
                      observationDatasetDoc( "set_constant_matrix_weight" ) )
                .def( "set_weight_vector_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightVectorForSet,
                      py::arg( "set_id" ),
                      py::arg( "weight_vector" ),
                      observationDatasetDoc( "set_weight_vector_for_set" ) )
                .def( "set_weight_matrix_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightMatrixForSet,
                      py::arg( "set_id" ),
                      py::arg( "weight_matrix" ),
                      observationDatasetDoc( "set_weight_matrix_for_set" ) )
                .def( "has_weight_matrix_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::hasWeightMatrixForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "has_weight_matrix_for_set" ) )
                .def( "set_weight_matrix_for_observation",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightMatrixForObservation,
                      py::arg( "observation_id" ),
                      py::arg( "weight_matrix" ),
                      observationDatasetDoc( "set_weight_matrix_for_observation" ) )
                .def( "has_weight_matrix_for_observation",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::hasWeightMatrixForObservation,
                      py::arg( "observation_id" ),
                      observationDatasetDoc( "has_weight_matrix_for_observation" ) )
                .def( "set_weight_block",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightBlock,
                      py::arg( "row_observation_ids" ),
                      py::arg( "column_observation_ids" ),
                      py::arg( "weight_block" ),
                      py::arg( "row_components" ) = std::vector< unsigned int >( ),
                      py::arg( "column_components" ) = std::vector< unsigned int >( ),
                      py::arg( "make_symmetric" ) = false,
                      observationDatasetDoc( "set_weight_block" ) )
                .def( "extra_weight_blocks",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getExtraWeightBlocks,
                      R"doc(Return the advanced scalar-component weight blocks stored on this dataset.)doc" )
                .def_property_readonly( "has_extra_weight_blocks",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::hasExtraWeightBlocks,
                                        R"doc(True when the dataset stores advanced scalar-component weight blocks.)doc" )
                .def( "add_observations_to_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationsToSet,
                      py::arg( "set_id" ),
                      py::arg( "observations" ),
                      py::arg( "times" ),
                      py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                      py::arg( "weights" ) = std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                      py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                      py::arg( "sort_observations" ) = true,
                      observationDatasetDoc( "add_observations_to_set" ) )
                .def( "remove_observations_from_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::removeObservationsFromSet,
                      py::arg( "set_id" ),
                      py::arg( "indices_to_remove" ),
                      observationDatasetDoc( "remove_observations_from_set" ) )
                .def( "remove_observations",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::removeObservations,
                      py::arg( "condition" ),
                      observationDatasetDoc( "remove_observations" ) )
                .def( "delete_rejected_observations",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::removeRejectedObservations,
                      observationDatasetDoc( "delete_rejected_observations" ) )
                .def( "filtered_observation_indices",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getFilteredObservationIndices,
                      py::arg( "set_id" ),
                      py::arg( "observation_filter" ) )
                .def( "time_bounds_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "time_bounds_for_set" ) )
                .def( "erase_duplicate_observations_from_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::eraseDuplicateObservationsFromSet,
                      py::arg( "set_id" ),
                      py::arg( "print_warning" ) = true,
                      observationDatasetDoc( "erase_duplicate_observations_from_set" ) )
                .def( "computed_observations_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservationsForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "computed_observations_for_set" ) )
                .def( "computed_observation_vector_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservationVectorForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "computed_observation_vector_for_set" ) )
                .def( "rms_residuals_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResidualsForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "rms_residuals_for_set" ) )
                .def( "mean_residuals_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResidualsForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "mean_residuals_for_set" ) )
                .def_property_readonly( "number_of_observation_sets",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservationSets,
                                        R"doc(Number of logical observation sets stored in the dataset.)doc" )
                .def_property_readonly( "number_of_observations",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservations,
                                        R"doc(Number of observation rows stored in the dataset.)doc" )
                .def_property_readonly( "total_scalar_size",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalScalarSize,
                                        R"doc(Number of scalar observation components stored in the dataset.)doc" )
                .def_property_readonly(
                        "observation_set_metadata",
                        py::overload_cast<>( &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSetMetadata,
                                             py::const_ ),
                        R"doc(List of metadata entries for all observation sets.)doc" )
                .def( "get_observation_set_metadata",
                      py::overload_cast< const unsigned int >(
                              &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSetMetadata, py::const_ ),
                      py::arg( "set_id" ),
                      observationDatasetDoc( "get_observation_set_metadata" ) )
                .def_property_readonly( "observation_rows",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationRows,
                                        R"doc(Row metadata for all observations.)doc" )
                .def( "observation_row",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationRow,
                      py::arg( "observation_id" ),
                      observationDatasetDoc( "observation_row" ) )
                .def_property_readonly( "scalar_component_rows",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getScalarComponentRows,
                                        R"doc(Row metadata for all scalar components.)doc" )
                .def( "scalar_component_row",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getScalarComponentRow,
                      py::arg( "scalar_component_id" ),
                      observationDatasetDoc( "scalar_component_row" ) )
                .def( "observation_ids_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationIdsForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "observation_ids_for_set" ) )
                .def( "observations_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "observations_for_set" ) )
                .def( "observation_vector_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationVectorForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "observation_vector_for_set" ) )
                .def( "observation_value",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationValue,
                      py::arg( "observation_id" ),
                      observationDatasetDoc( "observation_value" ) )
                .def( "observation_times_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimesForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "observation_times_for_set" ) )
                .def( "observation_time",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTime,
                      py::arg( "observation_id" ),
                      observationDatasetDoc( "observation_time" ) )
                .def( "weights_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "weights_for_set" ) )
                .def( "weight_vector_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightVectorForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "weight_vector_for_set" ) )
                .def( "weight_matrix_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightMatrixForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "weight_matrix_for_set" ) )
                .def( "weight_value",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightValue,
                      py::arg( "observation_id" ),
                      observationDatasetDoc( "weight_value" ) )
                .def( "weight_matrix_for_observation",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightMatrixForObservation,
                      py::arg( "observation_id" ),
                      observationDatasetDoc( "weight_matrix_for_observation" ) )
                .def( "residuals_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualsForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "residuals_for_set" ) )
                .def( "residual_vector_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualVectorForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "residual_vector_for_set" ) )
                .def( "residual_value",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualValue,
                      py::arg( "observation_id" ),
                      observationDatasetDoc( "residual_value" ) )
                .def( "dependent_variables_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariablesForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "dependent_variables_for_set" ) )
                .def( "dependent_variables",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariables,
                      py::arg( "observation_id" ),
                      observationDatasetDoc( "dependent_variables" ) )
                .def( "single_dependent_variable_for_set",
                      py::overload_cast< const unsigned int,
                                         const std::shared_ptr< tss::ObservationDependentVariableSettings >&,
                                         const bool >(
                              &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariableForSet, py::const_ ),
                      py::arg( "set_id" ),
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "return_first_compatible_settings" ) = false,
                      observationDatasetDoc( "single_dependent_variable_for_set" ) )
                .def( "single_dependent_variable_for_set_by_index",
                      py::overload_cast< const unsigned int, const std::pair< int, int >& >(
                              &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariableForSet, py::const_ ),
                      py::arg( "set_id" ),
                      py::arg( "dependent_variable_index_and_size" ),
                      observationDatasetDoc( "single_dependent_variable_for_set_by_index" ) )
                .def( "compatible_dependent_variable_settings_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariableSettingsForSet,
                      py::arg( "set_id" ),
                      py::arg( "dependent_variable_settings" ),
                      observationDatasetDoc( "compatible_dependent_variable_settings_for_set" ) )
                .def( "all_compatible_dependent_variables_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariablesForSet,
                      py::arg( "set_id" ),
                      py::arg( "dependent_variable_settings" ),
                      observationDatasetDoc( "all_compatible_dependent_variables_for_set" ) )
                .def( "set_dependent_variables_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setDependentVariablesForSet,
                      py::arg( "set_id" ),
                      py::arg( "dependent_variables" ),
                      observationDatasetDoc( "set_dependent_variables_for_set" ) )
                .def( "clear_dependent_variables_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::clearDependentVariablesForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "clear_dependent_variables_for_set" ) )
                .def( "number_of_observations_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservationsForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "number_of_observations_for_set" ) )
                .def( "total_scalar_size_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalScalarSizeForSet,
                      py::arg( "set_id" ),
                      observationDatasetDoc( "total_scalar_size_for_set" ) )
                .def( "link_definition",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinition,
                      py::arg( "link_definition_id" ),
                      observationDatasetDoc( "link_definition" ) )
                .def_property_readonly( "number_of_link_definitions",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfLinkDefinitions,
                                        R"doc(Number of unique link definitions registered in the dataset.)doc" )
                .def( "ancillary_settings",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettings,
                      py::arg( "ancillary_settings_id" ),
                      observationDatasetDoc( "ancillary_settings" ) )
                .def( "dependent_variable_bookkeeping",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableBookkeeping,
                      py::arg( "dependent_variable_layout_id" ),
                      observationDatasetDoc( "dependent_variable_bookkeeping" ) )
                .def(
                        "observation_set_ids",
                        []( const tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& dataset,
                            const std::shared_ptr< tom::ObservationCollectionParser >& observationParser ) {
                            warnLegacyObservationInterface( "ObservationDataset.observation_set_ids(observation_parser)",
                                                            "ObservationDataset.observation_ids_matching_condition" );
                            return dataset.getObservationSetIds( observationParser );
                        },
                        py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ) )
                .def( "observation_ids_matching_condition",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationIdsMatchingCondition,
                      py::arg( "condition" ),
                      observationDatasetDoc( "observation_ids_matching_condition" ) )
                .def( "create_viewer",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createViewer,
                      py::keep_alive< 0, 1 >( ),
                      py::arg( "condition" ),
                      observationDatasetDoc( "create_viewer" ) )
                .def( "create_new_and_keep",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createNewAndKeep,
                      py::arg( "condition" ),
                      observationDatasetDoc( "create_new_and_keep" ) )
                .def( "create_new_and_drop",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createNewAndDrop,
                      py::arg( "condition" ),
                      observationDatasetDoc( "create_new_and_drop" ) )
                .def( "reject_observations",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::rejectObservations,
                      py::arg( "condition" ),
                      py::arg( "reason" ) = "",
                      observationDatasetDoc( "reject_observations" ) )
                .def( "restore_observations",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::restoreObservations,
                      py::arg( "condition" ),
                      observationDatasetDoc( "restore_observations" ) )
                .def( "ordered_flattened_observation_data",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createOrderedFlattenedObservationData,
                      py::arg( "include_inactive" ) = true,
                      observationDatasetDoc( "ordered_flattened_observation_data" ) )
                .def( "estimation_flattened_observation_data",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createEstimationFlattenedObservationData,
                      py::arg( "include_rejected" ) = false,
                      observationDatasetDoc( "estimation_flattened_observation_data" ) )
                .def( "computation_flattened_observation_data",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createComputationFlattenedObservationData,
                      py::arg( "include_rejected" ) = true,
                      observationDatasetDoc( "computation_flattened_observation_data" ) );
    }

    {
        py::options legacyDocOptions;
        legacyDocOptions.disable_user_defined_docstrings( );
        legacyDocOptions.disable_function_signatures( );

        m.def(
                "create_observation_dataset_from_single_observation_set",
                []( const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >& observationSet ) {
                    warnLegacyObservationInterface( "create_observation_dataset_from_single_observation_set",
                                                    "ObservationDataset.add_observation_set" );
                    return tom::createObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >( observationSet );
                },
                py::arg( "observation_set" ) );

        m.def(
                "create_observation_dataset_from_collection",
                []( const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >& observationCollection ) {
                    warnLegacyObservationInterface( "create_observation_dataset_from_collection", "ObservationDataset" );
                    return tom::createObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >( observationCollection );
                },
                py::arg( "observation_collection" ) );

        m.def(
                "create_single_observation_set_from_dataset",
                []( const std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >& observationDataset,
                    const unsigned int setId ) {
                    warnLegacyObservationInterface( "create_single_observation_set_from_dataset",
                                                    "ObservationDataset.observation_ids_for_set" );
                    return tom::createSingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >( observationDataset, setId );
                },
                py::arg( "observation_dataset" ),
                py::arg( "set_id" ) = 0 );

        m.def(
                "create_observation_collection_from_dataset",
                []( const std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >& observationDataset ) {
                    warnLegacyObservationInterface( "create_observation_collection_from_dataset", "ObservationDataset" );
                    return tom::createObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >( observationDataset );
                },
                py::arg( "observation_dataset" ) );
    }

    // SINGLE OBSERVATION SET

    {
        py::options legacyDocOptions;
        legacyDocOptions.disable_user_defined_docstrings( );
        legacyDocOptions.disable_function_signatures( );

        py::class_< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >,
                    std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                    "SingleObservationSet",
                                                                                                    R"doc(
        Class for storing and handling a single set of observations.

        This class stores a set of observations of a single observable type and a single link definition.
        The observations are stored as a vector of Eigen matrices, where each entry represents a single
        observation at a given time. The class also stores the observation times, reference link end,
        and other metadata.

        The pybind ``SingleObservationSet`` object is a list-like object, which can be iterated over, and from which specific
        observation data can be retrieved by index.

        When iterating, a tuple is returned with the following entries:
        - ``int``: index of the observation.
        - ``float``: time of the observation.
        - ``numpy.ndarray``: value of the observation.

        When using the ``[]`` operator, a tuple is returned with the following entries:
        - ``float``: time of the observation.
        - ``numpy.ndarray``: value of the observation.

        Parameters
        ----------
        observable_type : tudatpy.estimation.observable_models_setup.model_settings.ObservableType
            Type of observable.
        link_ends : tudatpy.estimation.observable_models_setup.links.LinkDefinition
            Definition of the link ends for the observation.
        observations : list[numpy.ndarray]
            List of observations. Each entry is a vector representing a single observation.
        observation_epochs : list[float]
            List of observation times.
        reference_link_end : tudatpy.estimation.observable_models_setup.links.LinkEndType
            Reference link end for the observation.
        observation_dependent_variables : list[numpy.ndarray], optional
            List of dependent variables for each observation.
        dependent_variable_calculator : tudatpy.estimation.observations.ObservationDependentVariableCalculator, optional
            Calculator for dependent variables.
        ancillary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings, optional
            Ancillary settings for the observation.
      )doc" )
                .def( py::init( []( const tom::ObservableType observableType,
                                    const tom::LinkDefinition linkEnds,
                                    const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > > observations,
                                    const std::vector< TIME_TYPE > observationEpochs,
                                    const tom::LinkEndType referenceLinkEnd,
                                    const std::vector< Eigen::VectorXd > observationDependentVariables,
                                    const std::shared_ptr< tss::ObservationDependentVariableBookkeeping > dependentVariableBookkeeping,
                                    const std::shared_ptr< tom::ObservationAncillarySimulationSettings > ancillarySettings ) {
                          warnLegacyObservationInterface( "SingleObservationSet", "ObservationDataset.add_observation_set" );
                          return std::make_shared< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >(
                                  observableType,
                                  linkEnds,
                                  observations,
                                  observationEpochs,
                                  referenceLinkEnd,
                                  observationDependentVariables,
                                  dependentVariableBookkeeping,
                                  ancillarySettings );
                      } ),
                      py::arg( "observable_type" ),
                      py::arg( "link_ends" ),
                      py::arg( "observations" ),
                      py::arg( "observation_epochs" ),
                      py::arg( "reference_link_end" ),
                      py::arg( "observation_dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                      py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                      py::arg( "ancillary_settings" ) = nullptr )
                .def(
                        "__getattribute__",
                        []( const py::object& self, const std::string& attributeName ) {
                            return getLegacyAttributeWithWarning(
                                    self, attributeName, "SingleObservationSet", getSingleObservationSetReplacement( attributeName ) );
                        },
                        py::arg( "name" ) )
                .def( "set_observations",
                      py::overload_cast< const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                              &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                      py::arg( "observations" ),
                      R"doc(
Sets the observation values from a list of vectors.

Parameters
----------
observations : list[numpy.ndarray]
    The new list of observations.
)doc" )
                .def( "set_observations",
                      py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                              &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                      py::arg( "observations" ),
                      R"doc(
Sets the observation values from a single concatenated vector.

Parameters
----------
observations : numpy.ndarray
    A single vector containing all observations concatenated.
)doc" )
                .def( "set_residuals",
                      py::overload_cast< const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                              &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                      py::arg( "residuals" ),
                      R"doc(
Sets the residuals for all observations from a list of vectors.

Parameters
----------
residuals : list[numpy.ndarray]
    The new list of residuals.
)doc" )
                .def( "set_residuals",
                      py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                              &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                      py::arg( "residuals" ),
                      R"doc(
Sets the residuals for all observations from a single concatenated vector.

Parameters
----------
residuals : numpy.ndarray
    A single vector containing all residuals concatenated.
)doc" )
                .def( "set_constant_weight",
                      py::overload_cast< const double >( &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                      py::arg( "weight" ),
                      R"doc(
Sets a constant scalar weight for all observations.

Parameters
----------
weight : float
    The constant weight to apply.
)doc" )
                .def( "set_constant_weight",
                      py::overload_cast< const Eigen::Matrix< double, Eigen::Dynamic, 1 >& >(
                              &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                      py::arg( "weight" ),
                      R"doc(
Sets a constant weight vector for all observations.

Parameters
----------
weight : numpy.ndarray
    The constant weight vector to apply to each observation.
)doc" )
                .def( "set_tabulated_weights",
                      py::overload_cast< const Eigen::VectorXd& >(
                              &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                      py::arg( "weights" ),
                      R"doc(
Sets weights for all observations from a single concatenated vector.

Parameters
----------
weights : numpy.ndarray
    A single vector containing all weights concatenated.
)doc" )
                .def(
                        "filter_observations",
                        []( tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >& observationSet,
                            const std::shared_ptr< tom::ObservationFilterBase >& observationFilter,
                            const bool saveFilteredObservations ) {
                            warnLegacyObservationInterface( "SingleObservationSet.filter_observations",
                                                            "ObservationDataset.reject_observations" );
                            observationSet.filterObservations( observationFilter, saveFilteredObservations );
                        },
                        py::arg( "filter" ),
                        py::arg( "save_filtered_obs" ) = true,
                        R"doc(
Filters observations based on a given filter criterion.

Parameters
----------
filter : tudatpy.estimation.observations.observations_processing.ObservationFilterBase
    The filter to apply.
save_filtered_obs : bool, optional
    If true, the filtered observations are stored in a separate set. Defaults to true.
)doc" )
                .def_property_readonly( "observable_type",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableType,
                                        R"doc(

         **read-only**

         Type of observable for which the object stores observations

         :type: ObservableType
      )doc" )
                .def_property( "link_definition",
                               &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEnds,
                               &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setLinkEnds,
                               R"doc(

         Definition of the link ends for which the object stores observations

         :type: LinkDefinition
      )doc" )
                .def_property_readonly( "reference_link_end",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getReferenceLinkEnd,
                                        R"doc(

         **read-only**

         Reference link end for all stored observations

         :type: LinkEndType
      )doc" )
                .def_property_readonly( "number_of_observables",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservables,
                                        R"doc(
Returns the number of observations in the set.

Returns
-------
int
    The number of observations.
)doc" )
                .def_property_readonly( "single_observable_size",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleObservableSize,
                                        R"doc(
Returns the size of a single observation value (e.g., 1 for range, 2 for angular position).

Returns
-------
int
    The size of a single observation.
)doc" )
                .def_property_readonly( "total_observation_set_size",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalObservationSetSize,
                                        R"doc(
Returns the total size of all observation values in the set.

Returns
-------
int
    The total size of the observation set.
)doc" )
                .def_property_readonly( "time_bounds",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBounds,
                                        R"doc(
Returns the start and end time of the observation set.

Returns
-------
tuple[tudatpy.astro.time_representation.Time, tudatpy.astro.time_representation.Time]
    The start and end time of the observations.
)doc" )
                .def_property_readonly( "list_of_observations",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                                        R"doc(

         **read-only**

         List of separate stored observations. Each entry of this list is a vector containing a single observation. In cases where the observation is single-valued (range, Doppler), the vector is size 1, but for multi-valued observations such as angular position, each vector in the list will have size >1

         :type: list[ numpy.ndarray[numpy.float64[m, 1]] ]
      )doc" )
                .def_property_readonly( "observation_times",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimes,
                                        R"doc(

         **read-only**

         List of reference times for each of the observations in ``list_of_observations``

         :type: list[ tudatpy.astro.time_representation.Time ]
      )doc" )
                .def_property_readonly( "concatenated_observations",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsVector,
                                        R"doc(

         **read-only**

         Concatenated vector of all stored observations

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
                .def_property_readonly( "computed_observations",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservations,
                                        R"doc(
Returns the vector of computed observation values (observed - residual).

Returns
-------
list[numpy.ndarray]
    The list of computed observations.
)doc" )
                .def_property_readonly( "concatenated_computed_observations",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservationsVector,
                                        R"doc(
Returns all computed observations concatenated into a single vector.

Returns
-------
numpy.ndarray
    A single vector containing all computed observations.
)doc" )
                .def_property_readonly( "residuals",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getResiduals,
                                        R"doc(
Returns the vector of residuals for all observations.

Returns
-------
list[numpy.ndarray]
    The list of residuals.
)doc" )
                .def_property_readonly( "concatenated_residuals",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualsVector,
                                        R"doc(
Returns all residuals concatenated into a single vector.

Returns
-------
numpy.ndarray
    A single vector containing all residuals.
)doc" )
                .def_property_readonly( "rms_residuals",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResiduals,
                                        R"doc(
Returns the Root Mean Square (RMS) of the residuals.

Returns
-------
numpy.ndarray
    A vector containing the RMS of residuals for each component of the observable.
)doc" )
                .def_property_readonly( "mean_residuals",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResiduals,
                                        R"doc(
Returns the mean of the residuals.

Returns
-------
numpy.ndarray
    A vector containing the mean of residuals for each component of the observable.
)doc" )
                .def_property_readonly( "weights", &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeights, R"doc(
Returns the vector of weights for all observations.

Returns
-------
list[numpy.ndarray]
    The list of weights.
)doc" )
                .def_property_readonly( "concatenad_weights",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsVector,
                                        R"doc(
Returns all weights concatenated into a single vector.

Returns
-------
numpy.ndarray
    A single vector containing all weights.
)doc" )
                .def_property( "dependent_variables",
                               &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsDependentVariables,
                               &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationsDependentVariables,
                               R"doc(List of dependent variables for all observations.)doc" )
                .def_property_readonly( "dependent_variables_history",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistory,
                                        R"doc(
Time history of all observation dependent variables.

Returns
-------
dict[tudatpy.astro.time_representation.Time, numpy.ndarray]
    A map from observation time to the vector of dependent variables.
)doc" )
                .def_property_readonly( "observations_history",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsHistory,
                                        R"doc(

         **read-only**

         Dictionary of observations sorted by time. Created by making a dictionary with ``observation_times`` as keys and ``list_of_observations`` as values.

         :type: dict[ tudatpy.astro.time_representation.Time, numpy.ndarray[numpy.float64[m, 1]] ]
      )doc" )
                .def_property_readonly( "ancillary_settings",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettings,
                                        R"doc(

         **read-only**

         Ancillary settings for all stored observations

         :type: ObservationAncillarySimulationSettings
      )doc" )
                .def_property( "weights_vector",
                               &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsVector,
                               &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights,
                               R"doc(Concatenated vector of weights for all observations.)doc" )
                .def_property_readonly( "filtered_observation_set",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getFilteredObservationSet,
                                        R"doc(
Returns the set of filtered observations.

Returns
-------
SingleObservationSet
    The observation set containing filtered observations.
)doc" )
                .def_property_readonly( "number_filtered_observations",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfFilteredObservations,
                                        R"doc(
Returns the number of observations that have been filtered out.

Returns
-------
int
    The number of filtered observations.
)doc" )
                .def( "single_dependent_variable",
                      py::overload_cast< std::shared_ptr< tss::ObservationDependentVariableSettings >, const bool >(
                              &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariable ),
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "return_first_compatible_settings" ) = false,
                      R"doc(
Returns the values of a single dependent variable (specified by dependent variable settings).

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    Settings for the dependent variable to retrieve.
return_first_compatible_settings : bool, optional
    If true, returns the first compatible variable found. Defaults to false.

Returns
-------
numpy.ndarray
    A matrix of the dependent variable values over all observation times.
)doc" )
                .def( "compatible_dependent_variable_settings",
                      &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                      py::arg( "dependent_variable_settings" ),
                      R"doc(
Returns the list of all dependent variable settings compatible with the settings provided as inputs.

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    The settings to check for compatibility.

Returns
-------
list[tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings]
    A list of compatible dependent variable settings.
)doc" )
                .def( "compatible_dependent_variables_list",
                      &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariables,
                      py::arg( "dependent_variable_settings" ),
                      R"doc(
Returns a vector containing the values of all dependent variables compatible with the settings provided as input.

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    The settings to check for compatibility.

Returns
-------
list[numpy.ndarray]
    A list of matrices, each containing values of a compatible dependent variable.
)doc" )
                .def( "single_dependent_variable_history",
                      &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariableHistory,
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "return_first_compatible_settings" ) = false,
                      R"doc(
Returns the time history of a single dependent variable (specified by settings).

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    Settings for the dependent variable to retrieve.
return_first_compatible_settings : bool, optional
    If true, returns the first compatible variable found. Defaults to false.

Returns
-------
dict[float, numpy.ndarray]
    A map from observation time to the value of the specified dependent variable.
)doc" )
                .def_property_readonly( "dependent_variables_matrix",
                                        &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsDependentVariablesMatrix,
                                        R"doc(
Returns the dependent variables for all observations as a matrix.

Returns
-------
numpy.ndarray
    A matrix where each row corresponds to an observation and columns to dependent variables.
)doc" );
    }

    {
        py::options legacyDocOptions;
        legacyDocOptions.disable_user_defined_docstrings( );
        legacyDocOptions.disable_function_signatures( );

        m.def( "single_observation_set",
               &tss::singleObservationSetWithoutDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >,
               py::arg( "observable_type" ),
               py::arg( "link_definition" ),
               py::arg( "observations" ),
               py::arg( "observation_times" ),
               py::arg( "reference_link_end" ),
               py::arg( "ancillary_settings" ) = nullptr,
               R"doc(

        Deprecated. Use :func:`~tudatpy.estimation.observations.create_single_observation_set` instead.

        )doc" );

        m.def(
                "create_single_observation_set",
                []( const tom::ObservableType observableType,
                    const tom::LinkEnds& linkEnds,
                    const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& observations,
                    const std::vector< TIME_TYPE > observationTimes,
                    const tom::LinkEndType referenceLinkEnd,
                    const std::shared_ptr< tom::ObservationAncillarySimulationSettings > ancillarySettings ) {
                    warnLegacyObservationInterface( "create_single_observation_set", "ObservationDataset.add_observation_set" );
                    return tom::createSingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >(
                            observableType, linkEnds, observations, observationTimes, referenceLinkEnd, ancillarySettings );
                },
                py::arg( "observable_type" ),
                py::arg( "link_ends" ),
                py::arg( "observations" ),
                py::arg( "observation_times" ),
                py::arg( "reference_link_end" ),
                py::arg( "ancillary_settings" ) = nullptr,
                R"doc(
        Factory function to create a `SingleObservationSet` object.

        This function creates a `SingleObservationSet` object from a list of observations and their corresponding times.
        This is a simplified factory function that does not include dependent variables.

        Parameters
        ----------
        observable_type : :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
            Type of observable.
        link_ends : dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`]
            Definition of the link ends for the observation.
        observations : list[numpy.ndarray]
            List of observations. Each entry is a vector representing a single observation.
        observation_times : list[float]
            List of observation times.
        reference_link_end : :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`
            Reference link end for the observation.
        ancillary_settings : :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings`, optional
            Ancillary settings for the observation.

        Returns
        -------
        tudatpy.estimation.observations.SingleObservationSet
            A `SingleObservationSet` object.
        )doc" );
    }

    // OBSERVATION COLLECTION

    {
        py::options legacyDocOptions;
        legacyDocOptions.disable_user_defined_docstrings( );
        legacyDocOptions.disable_function_signatures( );

        py::class_< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
                    std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                     "ObservationCollection",
                                                                                                     R"doc(

         Class collecting all observations and associated data for use in an estimation.

         Class containing the full set of observations and associated data, typically for input into the estimation. When using simulated data,
         this class is instantiated via a call to the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. More information is provided
         on the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__





      )doc" )
                .def( py::init( []( const std::vector< std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                                            observationSets ) {
                          warnLegacyObservationInterface( "ObservationCollection", "ObservationDataset" );
                          return std::make_shared< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >( observationSets );
                      } ),
                      py::arg( "observation_sets" ),
                      R"doc(
Constructor for the ObservationCollection class.

Parameters
----------
observation_sets : list[tudatpy.estimation.observations.SingleObservationSet]
    List of single observation sets to be included in the collection.
)doc" )
                .def(
                        "__getattribute__",
                        []( const py::object& self, const std::string& attributeName ) {
                            return getLegacyAttributeWithWarning(
                                    self, attributeName, "ObservationCollection", getObservationCollectionReplacement( attributeName ) );
                        },
                        py::arg( "name" ) )
                .def_property_readonly( "concatenated_times",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDoubleTimeVector,
                                        R"doc(

         **read-only**

         Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__ for details on storage order

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
                .def_property_readonly( "concatenated_times_objects",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedTimeVector,
                                        R"doc(

         **read-only**

         Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__ for details on storage order

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
                .def_property_readonly( "concatenated_weights",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getUnparsedConcatenatedWeights,
                                        R"doc(
         **read-only**

         Vector containing concatenated observation weights.

         :type: numpy.ndarray[numpy.float64[m, 1]]
)doc" )
                .def_property_readonly( "concatenated_observations",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationVector,
                                        R"doc(

         **read-only**

         Vector containing concatenated observable values. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__ for details on storage order

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
                .def_property_readonly(
                        "concatenated_link_definition_ids",
                        py::overload_cast<>( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedLinkEndIds ),
                        R"doc(

         **read-only**

         Vector containing concatenated indices identifying the link ends. Each set of link ends is assigned a unique integer identifier (for a given instance of this class). The definition of a given integer identifier with the link ends is given by this class' :func:`link_definition_ids` function. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__ for details on storage order of the present vector.

         :type: numpy.ndarray[ int ]
      )doc" )
                .def_property_readonly( "link_definition_ids",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getInverseLinkEndIdentifierMap,
                                        R"doc(

         **read-only**

         Dictionary mapping a link end integer identifier to the specific link ends

         :type: dict[ int, dict[ LinkEndType, LinkEndId ] ]
      )doc" )
                .def_property_readonly( "observable_type_start_index_and_size",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTypeStartAndSize,
                                        R"doc(

         **read-only**

         Dictionary defining per observable type (dict key), the index in the full observation vector (:func:`concatenated_observations`) where the given observable type starts, and the number of subsequent entries in this vector containing a value of an observable of this type

         :type: dict[ ObservableType, [ int, int ] ]
      )doc" )
                .def_property_readonly(
                        "observation_set_start_index_and_size",
                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSetStartAndSizePerLinkEndIndex,
                        R"doc(

         **read-only**

         The nested dictionary/list returned by this property mirrors the structure of the :func:`sorted_observation_sets` property of this class. The present function provides the start index and size of the observables in the full observation vector that come from the corresponding `SingleObservationSet` in the :func:`sorted_observation_sets` Consequently, the present property returns a nested dictionary defining per observable type, link end identifier, and `SingleObservationSet` index (for the given observable type and link end identifier), where the observables in the given `SingleObservationSet` starts, and the number of subsequent entries in this vector containing data from it.

         :type: dict[ ObservableType, dict[ int, list[ int, int ] ] ]
      )doc" )
                .def_property_readonly( "observation_vector_size",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalObservableSize,
                                        R"doc(

         **read-only**

         Length of the total vector of observations

         :type: int
      )doc" )
                .def_property_readonly( "sorted_observation_sets",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSortedObservationSets,
                                        R"doc(

         **read-only**

         The nested dictionary/list contains the list of `SingleObservationSet` objects, in the same method as they are stored internally in the present class. Specifics on the storage order are given in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__

         :type: dict[ ObservableType, dict[ int, list[ SingleObservationSet ] ] ]
      )doc" )
                .def_property_readonly( "link_ends_per_observable_type",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEndsPerObservableType,
                                        R"doc(
         **read-only**

         Dictionary mapping each observable type to a list of link ends for which observations are available.

         :type: dict[ ObservableType, list[LinkEnds] ]
)doc" )
                .def_property_readonly( "link_definitions_per_observable",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinitionsPerObservable,
                                        R"doc(
         **read-only**

         Dictionary mapping each observable type to a list of link definitions.

         :type: dict[ ObservableType, list[LinkDefinition] ]
)doc" )
                .def_property_readonly( "time_bounds",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsDouble,
                                        R"doc(
         **read-only**

         Pair of minimum and maximum observation time in the collection.

         :type: tuple[float, float]
)doc" )
                .def_property_readonly( "time_bounds_time_object",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBounds,
                                        R"doc(
         **read-only**

         Pair of minimum and maximum observation time in the collection, as Time objects.

         :type: tuple[Time, Time]
)doc" )
                .def_property_readonly( "sorted_per_set_time_bounds",
                                        &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSortedObservationSetsTimeBounds,
                                        R"doc(
         **read-only**

         Nested dictionary providing the time bounds for each single observation set, sorted by observable type and link end ID.

         :type: dict[ ObservableType, dict[ int, list[tuple[float, float]] ] ]
)doc" )
                .def( "set_observations",
                      py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                      py::arg( "observations" ),
                      R"doc(
Function to reset the full list of observable values. The order of the observations must be the same as for :attr:`~ObservationCollection.concatenated_observations`

Parameters
----------
observations : np.ndarray
    New list of observable values
     )doc" )
                .def( "set_observations",
                      py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&,
                                         const std::shared_ptr< tom::ObservationCollectionParser > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                      py::arg( "observations" ),
                      py::arg( "observation_parser" ),
                      R"doc(
Function to reset a subset of all observable values, with this subset defined by the ``observation_parser`` input.
The order and size of the new observation vector must be the same as when calling :attr:`~ObservationCollection.concatenated_observations` on
an ``ObservationCollection`` containing only the parsed observation.

Parameters
----------
observations : np.ndarray
    New list of observable values
observation_parser : ObservationCollectionParser
    Observation parser with which to select the subset of observations that is to be reset
     )doc" )
                .def( "set_observations",
                      py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                         Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                      py::arg( "observations_per_parser" ),
                      R"doc(
Function to reset a subset of all observable values, with this subset defined by a list of observation parsers input.
Each observation parser is associated with a new set of observable values.
The order and size of the new observation vector for each parser must be the same as when calling :attr:`~ObservationCollection.concatenated_observations` on
an ``ObservationCollection`` containing only the parsed observation (from a single parser). NOTE: since the multiple parsers
are handled in order (iterating over the keys of ``observations_per_parser``) some observations may be reset several times,
in case.

Parameters
----------
observations : np.ndarray
    New list of observable values
observation_parser : ObservationCollectionParser
    Observation parser with which to select the subset of observations that is to be reset
     )doc" )
                .def( "set_residuals",
                      py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                      py::arg( "residuals" ),
                      R"doc(
Function to reset the full list of observation residuals. The order of the residuals must be the same as for :attr:`~ObservationCollection.concatenated_observations`

Parameters
----------
residuals : np.ndarray
    New list of observation residuals
     )doc" )
                .def( "set_residuals",
                      py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&,
                                         const std::shared_ptr< tom::ObservationCollectionParser > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                      py::arg( "residuals" ),
                      py::arg( "observation_parser" ),
                      R"doc(
Function to reset a subset of all observation residuals, with this subset defined by the ``observation_parser`` input.

Parameters
----------
residuals : np.ndarray
    New list of observation residuals
observation_parser : ObservationCollectionParser
    Observation parser with which to select the subset of residuals that is to be reset
     )doc" )
                .def( "set_residuals",
                      py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                         Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                      py::arg( "residuals_per_parser" ),
                      R"doc(
Function to reset a subset of all observation residuals, with this subset defined by a list of observation parsers input.

Parameters
----------
residuals_per_parser : dict[ObservationCollectionParser, np.ndarray]
    Dictionary mapping observation parsers to new residual values.
     )doc" )
                .def( "get_link_definitions_for_observables",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinitionsForSingleObservable,
                      py::arg( "observable_type" ),
                      R"doc(
         Function to get all link definitions for a given observable type.

         Parameters
         ----------
         observable_type : :class:`tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
             Observable type for which link definitions are to be retrieved.
         Returns
         -------
         list[ LinkDefinition ]
             List of link definitions for the given observable type.
     )doc" )
                .def( "get_single_link_and_type_observations",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleLinkAndTypeObservationSets,
                      py::arg( "observable_type" ),
                      py::arg( "link_definition" ),
                      R"doc(
         Function to get all observation sets for a given observable type and link definition.


         Parameters
         ----------
         observable_type : :class:`tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
             Observable type of which observations are to be simulated.
         link_ends : LinkDefinition
             Link ends for which observations are to be simulated.
         Returns
         -------
         list[ SingleObservationSet ]
             List of observation sets for given observable type and link definition.
     )doc" )
                .def( "get_observable_types",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableTypes,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the observable types for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[ObservableType]
             A list of observable types present in the selected subset.
     )doc" )
                .def( "get_bodies_in_link_ends",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getBodiesInLinkEnds,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the names of bodies present in the link ends of a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[str]
             A list of body names.
     )doc" )
                .def( "get_reference_points_in_link_ends",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getReferencePointsInLinkEnds,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the names of reference points (e.g., ground stations) in the link ends of a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[str]
             A list of reference point names.
     )doc" )
                .def( "get_time_bounds_list",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsListDouble,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the time bounds for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[float, float]]
             A list of time bound pairs (start_time, end_time).
     )doc" )
                .def( "get_time_bounds_list_time_object",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsList,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the time bounds for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[Time, Time]]
             A list of time bound pairs (start_time, end_time) as Time objects.
     )doc" )
                .def( "get_time_bounds_per_set",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsPerSetDouble,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the time bounds for each set in a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[float, float]]
             A list of time bound pairs for each observation set.
     )doc" )
                .def( "get_time_bounds_per_set_time_object",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsPerSet,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the time bounds for each set in a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[Time, Time]]
             A list of time bound pairs for each observation set as Time objects.
     )doc" )
                .def( "get_observations",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of observation vectors, one for each selected observation set.
     )doc" )
                .def( "get_concatenated_observations",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservations,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the concatenated observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of observation values.
     )doc" )
                .def( "get_observation_times",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimesDouble,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the observation times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[list[float]]
             A list of lists of observation times, one for each selected observation set.
     )doc" )
                .def( "get_observation_times_objects",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimes,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the observation times for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[list[Time]]
             A list of lists of observation times as Time objects.
     )doc" )
                .def( "get_concatenated_observation_times_objects",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationTimes,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the concatenated observation times for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[Time]
             A concatenated list of observation times as Time objects.
     )doc" )
                .def( "get_concatenated_observation_times",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDoubleObservationTimes,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the concatenated observation times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[float]
             A concatenated list of observation times.
     )doc" )
                .def( "get_observations_and_times",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsAndTimesDouble,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the observations and times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], list[list[float]]]
             A pair containing a list of observation vectors and a list of lists of observation times.
     )doc" )
                .def( "get_observations_and_times_objects",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsAndTimes,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the observations and times for a subset of observation sets, with times as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], list[list[Time]]]
             A pair containing a list of observation vectors and a list of lists of observation times as Time objects.
     )doc" )
                .def( "get_concatenated_observations_and_times",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationsAndTimesDouble,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the concatenated observations and times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[numpy.ndarray, list[float]]
             A pair containing the concatenated observation vector and list of times.
     )doc" )
                .def( "get_concatenated_observations_and_times_objects",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationsAndTimes,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the concatenated observations and times for a subset of observation sets, with times as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[numpy.ndarray, list[Time]]
             A pair containing the concatenated observation vector and list of times as Time objects.
     )doc" )
                .def( "get_concatenated_link_definition_ids",
                      py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedLinkEndIds ),
                      py::arg( "observation_parser" ),
                      R"doc(
         Get the concatenated link definition IDs for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser
             Object that is used to select a subset of the observation sets.

         Returns
         -------
         list[int]
             A list of concatenated link end IDs.
     )doc" )
                .def( "get_weights",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getWeights,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the weights for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of weight vectors, one for each selected observation set.
     )doc" )
                .def( "get_concatenated_weights",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedWeights,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the concatenated weights for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of weights.
     )doc" )
                .def( "get_residuals",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getResiduals,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of residual vectors, one for each selected observation set.
     )doc" )
                .def( "get_concatenated_residuals",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedResiduals,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the concatenated residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of residuals.
     )doc" )
                .def( "get_rms_residuals",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResiduals,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the RMS of residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of RMS residual vectors, one for each selected observation set.
     )doc" )
                .def( "get_mean_residuals",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResiduals,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the mean of residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of mean residual vectors, one for each selected observation set.
     )doc" )
                .def( "get_computed_observations",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservations,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the computed observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of computed observation vectors, one for each selected observation set.
     )doc" )
                .def( "get_concatenated_computed_observations",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedComputedObservations,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get the concatenated computed observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of computed observation values.
     )doc" )
                .def( "set_constant_weight",
                      py::overload_cast< const double, const std::shared_ptr< tom::ObservationCollectionParser > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                      py::arg( "weight" ),
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Set a constant weight for a subset of observation sets.

         Parameters
         ----------
         weight : float
             The constant weight to set.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
                .def( "set_constant_weight",
                      py::overload_cast< const Eigen::VectorXd, const std::shared_ptr< tom::ObservationCollectionParser > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                      py::arg( "weight" ),
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Set a constant weight vector for a subset of observation sets.

         Parameters
         ----------
         weight : numpy.ndarray
             The constant weight vector to set.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
                .def( "set_constant_weight_per_observation_parser",
                      py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, double > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightPerObservable ),
                      py::arg( "weights_per_observation_parser" ),
                      R"doc(
         Set a constant weight for multiple subsets of observation sets.

         Parameters
         ----------
         weights_per_observation_parser : dict[tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, float]
             A dictionary mapping observation parsers to constant weights.
     )doc" )
                .def( "set_constant_weight_per_observation_parser",
                      py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightPerObservable ),
                      py::arg( "weights_per_observation_parser" ),
                      R"doc(
         Set a constant weight vector for multiple subsets of observation sets.

         Parameters
         ----------
         weights_per_observation_parser : dict[tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, numpy.ndarray]
             A dictionary mapping observation parsers to constant weight vectors.
     )doc" )
                .def( "set_tabulated_weights",
                      py::overload_cast< const Eigen::VectorXd, const std::shared_ptr< tom::ObservationCollectionParser > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                      py::arg( "tabulated_weights" ),
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Set tabulated weights for a subset of observation sets.

         Parameters
         ----------
         tabulated_weights : numpy.ndarray
             The vector of tabulated weights to set.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
                .def( "set_tabulated_weights",
                      py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                      py::arg( "tabulated_weights" ),
                      R"doc(
         Set tabulated weights for multiple subsets of observation sets.

         Parameters
         ----------
         tabulated_weights : dict[tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, numpy.ndarray]
             A dictionary mapping observation parsers to tabulated weight vectors.
     )doc" )
                .def( "append",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::appendObservationCollection,
                      py::arg( "observation_collection_to_append" ) )
                .def(
                        "filter_observations",
                        []( tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >& observationCollection,
                            const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                            std::shared_ptr< tom::ObservationFilterBase > >& observationFilters,
                            const bool saveFilteredObservations ) {
                            warnLegacyObservationInterface( "ObservationCollection.filter_observations",
                                                            "ObservationDataset.reject_observations" );
                            observationCollection.filterObservations( observationFilters, saveFilteredObservations );
                        },
                        py::arg( "observation_filters" ),
                        py::arg( "save_filtered_observations" ) = true,
                        R"doc(
         Filter observations using a set of filters.

         This function filters the observations in the collection based on a map of observation filters, each associated with an observation parser.

         Parameters
         ----------
         observation_filters : dict[tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, tudatpy.estimation.observations.observations_processing.ObservationFilterBase]
             A dictionary mapping observation parsers to observation filters.
         save_filtered_observations : bool, optional
             If true, the filtered-out observations are saved within each observation set, by default True.
     )doc" )
                .def(
                        "filter_observations",
                        []( tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >& observationCollection,
                            const std::shared_ptr< tom::ObservationFilterBase >& observationFilter,
                            const std::shared_ptr< tom::ObservationCollectionParser >& observationParser,
                            const bool saveFilteredObservations ) {
                            warnLegacyObservationInterface( "ObservationCollection.filter_observations",
                                                            "ObservationDataset.reject_observations" );
                            observationCollection.filterObservations( observationFilter, observationParser, saveFilteredObservations );
                        },
                        py::arg( "observation_filters" ),
                        py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                        py::arg( "save_filtered_observations" ) = true,
                        R"doc(
         Filter observations using a single filter.

         This function filters a subset of observations (or all) using a single observation filter.

         Parameters
         ----------
         observation_filter : tudatpy.estimation.observations.observations_processing.ObservationFilterBase
             The observation filter to apply.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
         save_filtered_observations : bool, optional
             If true, the filtered-out observations are saved within each observation set, by default True.
     )doc" )
                .def(
                        "split_observation_sets",
                        []( tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >& observationCollection,
                            const std::shared_ptr< tom::ObservationSetSplitterBase >& observationSetSplitter,
                            const std::shared_ptr< tom::ObservationCollectionParser >& observationParser ) {
                            warnLegacyObservationInterface( "ObservationCollection.split_observation_sets",
                                                            "ObservationDataset.create_new_and_keep" );
                            observationCollection.splitObservationSets( observationSetSplitter, observationParser );
                        },
                        py::arg( "observation_set_splitter" ),
                        py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                        R"doc(
         Split observation sets based on a splitter.

         This function splits a subset of observation sets (or all) into smaller sets based on the criteria defined by the splitter.

         Parameters
         ----------
         observation_set_splitter : tudatpy.estimation.observations.observations_processing.ObservationSetSplitterBase
             The splitter to use for splitting the observation sets.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
                .def( "get_single_observation_sets",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleObservationSets,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get a subset of single observation sets.

         This function returns a list of the single observation sets that are selected by the provided observation parser.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tudatpy.estimation.observations.SingleObservationSet]
             A list of the selected single observation sets.
     )doc" )
                .def( "print_observation_sets_start_and_size",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::printObservationSetsStartAndSize,
                      R"doc(Prints the structure of the observation collection, showing the start index and size of each individual observation set.)doc" )
                .def( "remove_single_observation_sets",
                      py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeSingleObservationSets ),
                      py::arg( "observation_parser" ),
                      R"doc(
         Remove a subset of single observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser
             Object that is used to select the subset of observation sets to remove.
     )doc" )
                .def( "set_reference_point",
                      py::overload_cast< tss::SystemOfBodies&,
                                         const Eigen::Vector3d&,
                                         const std::string&,
                                         const std::string&,
                                         const tom::LinkEndType,
                                         const std::shared_ptr< tom::ObservationCollectionParser > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoint ),
                      py::arg( "bodies" ),
                      py::arg( "antenna_position" ),
                      py::arg( "antenna_name" ),
                      py::arg( "spacecraft_name" ),
                      py::arg( "link_end_type" ),
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Set a fixed reference point for a subset of observations.

         Parameters
         ----------
         bodies : tudatpy.dynamics.environment.SystemOfBodies
             System of bodies.
         antenna_position : numpy.ndarray
             Position of the antenna in the spacecraft body-fixed frame.
         antenna_name : str
             Name of the antenna/reference point.
         spacecraft_name : str
             Name of the spacecraft body.
         link_end_type : LinkEndType
             Link end type to which the reference point should be applied.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select the observation sets to which the reference point should be applied.
     )doc" )
                .def( "set_reference_points",
                      py::overload_cast< tss::SystemOfBodies&,
                                         const std::map< double, Eigen::Vector3d >&,
                                         const std::string&,
                                         const tom::LinkEndType,
                                         const std::shared_ptr< tom::ObservationCollectionParser > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoints ),
                      py::arg( "bodies" ),
                      py::arg( "antenna_switch_history" ),
                      py::arg( "spacecraft_name" ),
                      py::arg( "link_end_type" ),
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Set multiple reference points based on an antenna switch history.

         Parameters
         ----------
         bodies : tudatpy.dynamics.environment.SystemOfBodies
             System of bodies.
         antenna_switch_history : dict[float, numpy.ndarray]
             Dictionary mapping time to antenna position.
         spacecraft_name : str
             Name of the spacecraft body.
         link_end_type : LinkEndType
             Link end type to which the reference points should be applied.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select the observation sets.
     )doc" )
                .def( "set_reference_point",
                      py::overload_cast< tss::SystemOfBodies&,
                                         const std::shared_ptr< te::Ephemeris >,
                                         const std::string&,
                                         const std::string&,
                                         const tom::LinkEndType,
                                         const std::shared_ptr< tom::ObservationCollectionParser > >(
                              &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoint ),
                      py::arg( "bodies" ),
                      py::arg( "antenna_body_fixed_ephemeris" ),
                      py::arg( "antenna_name" ),
                      py::arg( "spacecraft_name" ),
                      py::arg( "link_end_type" ),
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Set a time-varying reference point for a subset of observations.

         Parameters
         ----------
         bodies : tudatpy.dynamics.environment.SystemOfBodies
             System of bodies.
         antenna_body_fixed_ephemeris : tudatpy.dynamics.environment.Ephemeris
             Ephemeris of the antenna in the spacecraft body-fixed frame.
         antenna_name : str
             Name of the antenna/reference point.
         spacecraft_name : str
             Name of the spacecraft body.
         link_end_type : LinkEndType
             Link end type to which the reference point should be applied.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select the observation sets.
     )doc" )
                .def( "set_transponder_delay",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTransponderDelay,
                      py::arg( "spacecraft_name" ),
                      py::arg( "transponder_delay" ),
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Deprecated: set the transponder delay for a subset of observations by modifying the
         retransmission delay in their ancillary settings.

         For new simulations, set the default transponder delay on the spacecraft vehicle systems
         before creating the observation model:
         ``bodies.get_body(spacecraft_name).vehicle_systems.transponder_delay = transponder_delay``.

         Parameters
         ----------
         spacecraft_name : str
             Name of the spacecraft with the transponder.
         transponder_delay : float
             The transponder delay in seconds.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select the observation sets.
     )doc" )
                .def( "remove_empty_observation_sets",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeEmptySingleObservationSets,
                      R"doc(Remove all single observation sets that contain no observations.)doc" )
                .def( "add_dependent_variable",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::addDependentVariable,
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Add an observation dependent variable to a subset of the single observation sets.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to add.
         bodies : tudatpy.dynamics.environment.SystemOfBodies
             System of bodies containing the environment.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select the observation sets to which the variable should be added.
         
         Returns
         -------
         tudatpy.estimation.observations.observations_processing.ObservationCollectionParser
             A parser that can be used to retrieve the added dependent variable.
     )doc" )
                .def( "dependent_variable",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariables,
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "first_compatible_settings" ) = false,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Retrieve the values of a given dependent variable, sorted per single observation set.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], tudatpy.estimation.observations.observations_processing.ObservationCollectionParser]
             A pair containing a list of matrices with the dependent variable values and the parser used.
     )doc" )
                .def( "concatenated_dependent_variable",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDependentVariables,
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "first_compatible_settings" ) = false,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Retrieve the concatenated values of a given dependent variable.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[numpy.ndarray, tudatpy.estimation.observations.observations_processing.ObservationCollectionParser]
             A pair containing a matrix with the concatenated dependent variable values and the parser used.
     )doc" )
                .def( "compatible_dependent_variable_settings",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get a list of all dependent variable settings compatible with the input settings.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[list[list[tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings]], tudatpy.estimation.observations.observations_processing.ObservationCollectionParser]
             A pair containing a list of lists of compatible settings and the parser used.
     )doc" )
                .def( "compatible_dependent_variables_list",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariables,
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Get all dependent variables compatible with the input settings.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[list[list[numpy.ndarray]], tudatpy.estimation.observations.observations_processing.ObservationCollectionParser]
             A pair containing a list of lists of dependent variable values and the parser used.
     )doc" )
                .def( "dependent_variable_history_per_set",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryPerObservationSetDouble,
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "first_compatible_settings" ) = false,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Retrieve the time history of a given dependent variable, sorted per observation set.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         list[dict[float, numpy.ndarray]]
             A list of maps from time to dependent variable value, one for each set.
     )doc" )
                .def( "dependent_variable_history_per_set_time_object",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryPerObservationSet,
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "first_compatible_settings" ) = false,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Retrieve the time history of a given dependent variable, sorted per observation set, with times as Time objects.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         list[dict[Time, numpy.ndarray]]
             A list of maps from time to dependent variable value, one for each set, with times as Time objects.
     )doc" )
                .def( "dependent_variable_history",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryDouble,
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "first_compatible_settings" ) = false,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Retrieve the concatenated time history of a given dependent variable.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         dict[float, numpy.ndarray]
             A map from time to dependent variable value.
     )doc" )
                .def( "dependent_variable_history_time_object",
                      &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistory,
                      py::arg( "dependent_variable_settings" ),
                      py::arg( "first_compatible_settings" ) = false,
                      py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                      R"doc(
         Retrieve the concatenated time history of a given dependent variable, with times as Time objects.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         dict[Time, numpy.ndarray]
             A map from time to dependent variable value, with times as Time objects.
     )doc" );
    }

    {
        py::options legacyDocOptions;
        legacyDocOptions.disable_user_defined_docstrings( );
        legacyDocOptions.disable_function_signatures( );

        m.def(
                "compute_residuals_and_dependent_variables",
                []( std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > observationCollection,
                    const std::vector< std::shared_ptr< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                            observationSimulators,
                    const tss::SystemOfBodies& bodies ) {
                    warnLegacyObservationInterface( "compute_residuals_and_dependent_variables",
                                                    "compute_residuals_and_dependent_variables_for_dataset" );
                    tss::computeResidualsAndDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >(
                            observationCollection, observationSimulators, bodies );
                },
                py::arg( "observation_collection" ),
                py::arg( "observation_simulators" ),
                py::arg( "bodies" ),
                R"doc(
        Computes residuals and dependent variables for a given observation collection.

        This function simulates observations based on the settings of the input `observation_collection`
        (which typically contains real data). It then computes the residuals by subtracting the
        simulated observations from the original observations. The computed residuals and any
        observation-dependent variables are then stored in the input `observation_collection`.

        Parameters
        ----------
        observation_collection : tudatpy.estimation.observations.ObservationCollection
            The collection of observations for which to compute residuals and dependent variables.
            This object is modified in-place.
        observation_simulators : list[tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator]
            List of observation simulators to be used for computing the observations.
        bodies : tudatpy.dynamics.environment.SystemOfBodies
            The system of bodies required for the observation simulation.
        )doc" );
    }

    {
        py::options observationFunctionDocOptions;
        observationFunctionDocOptions.disable_function_signatures( );

        m.def( "observation_simulation_settings_from_dataset",
               &tss::getObservationSimulationSettingsFromObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >,
               py::arg( "observation_dataset" ),
               py::arg( "bodies" ),
               observationModuleFunctionDoc( "observation_simulation_settings_from_dataset" ) );

        m.def(
                "compute_residuals_and_dependent_variables_for_dataset",
                []( const std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >& observationDataset,
                    const std::vector< std::shared_ptr< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                            observationSimulators,
                    const tss::SystemOfBodies& bodies ) {
                    if( observationDataset == nullptr )
                    {
                        throw std::runtime_error(
                                "Error when computing residuals and dependent variables for dataset, input dataset is None." );
                    }
                    tss::computeResidualsAndDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >(
                            observationDataset, observationSimulators, bodies );
                },
                py::arg( "observation_dataset" ),
                py::arg( "observation_simulators" ),
                py::arg( "bodies" ),
                observationModuleFunctionDoc( "compute_residuals_and_dependent_variables_for_dataset" ) );
    }

    {
        py::options legacyDocOptions;
        legacyDocOptions.disable_user_defined_docstrings( );
        legacyDocOptions.disable_function_signatures( );

        m.def(
                "filter_observations",
                []( const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > originalObservationSet,
                    const std::shared_ptr< tom::ObservationFilterBase > observationFilter,
                    const bool saveFilteredObservations ) {
                    warnLegacyObservationInterface( "filter_observations", "ObservationDataset.reject_observations" );
                    return tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE >(
                            originalObservationSet, observationFilter, saveFilteredObservations );
                },
                py::arg( "original_observation_set" ),
                py::arg( "observation_filter" ),
                py::arg( "save_filtered_observations" ) = false,
                R"doc(

Deprecated. Use :func:`~tudatpy.estimation.observations.create_filtered_observation_set` instead.


        )doc" );

        m.def(
                "create_filtered_observation_set",
                []( const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > originalObservationSet,
                    const std::shared_ptr< tom::ObservationFilterBase > observationFilter,
                    const bool saveFilteredObservations ) {
                    warnLegacyObservationInterface( "create_filtered_observation_set", "ObservationDataset.create_new_and_drop" );
                    return tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE >(
                            originalObservationSet, observationFilter, saveFilteredObservations );
                },
                py::arg( "original_observation_set" ),
                py::arg( "observation_filter" ),
                py::arg( "save_filtered_observations" ) = false,
                R"doc(

        Filters a single observation set and returns a new set containing the filtered observations.

        This function creates a copy of the input observation set and applies the given filter.
        The original observation set is not modified.

        Parameters
        ----------
        original_observation_set : tudatpy.estimation.observations.SingleObservationSet
            The observation set to filter.
        observation_filter : tudatpy.estimation.observations.observations_processing.ObservationFilterBase
            The filter to apply.
        save_filtered_observations : bool, optional
            If true, the observations that are filtered out are stored within the new observation set. Defaults to false.

        Returns
        -------
        tudatpy.estimation.observations.SingleObservationSet
            A new observation set containing only the observations that passed the filter.
        )doc" );

        m.def(
                "split_observation_set",
                []( const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > originalObservationSet,
                    const std::shared_ptr< tom::ObservationSetSplitterBase > observationSplitter,
                    const bool printWarning ) {
                    warnLegacyObservationInterface( "split_observation_set", "ObservationDataset.create_new_and_keep" );
                    return tom::splitObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >(
                            originalObservationSet, observationSplitter, printWarning );
                },
                py::arg( "original_observation_set" ),
                py::arg( "observation_splitter" ),
                py::arg( "print_warning" ) = true,
                R"doc(
        Splits a single observation set into multiple sets based on a splitter.

        This function takes an observation set and divides it into a list of smaller observation sets
        according to the criteria defined in the observation splitter.

        Parameters
        ----------
        original_observation_set : tudatpy.estimation.observations.SingleObservationSet
            The observation set to split.
        observation_splitter : tudatpy.estimation.observations.observations_processing.ObservationSetSplitterBase
            The splitter defining the splitting criteria.
        print_warning : bool, optional
            If true, a warning is printed if the original set contains filtered observations that will be lost. Defaults to true.

        Returns
        -------
        list[tudatpy.estimation.observations.SingleObservationSet]
            A list of new observation sets resulting from the split.
        )doc" );

        m.def(
                "merge_observation_collections",
                []( const std::vector< std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                            observationCollectionList ) {
                    warnLegacyObservationInterface( "merge_observation_collections",
                                                    "ObservationDataset.add_observation_set_from_dataset" );
                    return tss::mergeObservationCollections< STATE_SCALAR_TYPE, TIME_TYPE >( observationCollectionList );
                },
                py::arg( "observation_collection_list" ) );

        // The following functions create a new ObservationCollection object from an existing one

        m.def(
                "create_filtered_observation_collection",
                []( const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > originalObservationCollection,
                    const std::map< std::shared_ptr< tom::ObservationCollectionParser >, std::shared_ptr< tom::ObservationFilterBase > >&
                            observationFiltersMap ) {
                    warnLegacyObservationInterface( "create_filtered_observation_collection", "ObservationDataset.create_new_and_drop" );
                    return tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE >( originalObservationCollection, observationFiltersMap );
                },
                py::arg( "original_observation_collection" ),
                py::arg( "observation_filters_map" ),
                R"doc(

        Creates a new, filtered observation collection from an existing one using multiple filters.

        This function applies a set of filters, each associated with a parser, to an observation collection
        and returns a new collection with the filtered observations. The original collection is not modified.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The observation collection to filter.
        observation_filters_map : dict[tudatpy.estimation.observations.ObservationCollectionParser, tudatpy.estimation.observations.observations_processing.ObservationFilterBase]
            A dictionary mapping parsers to filters.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the filtered observations.

        )doc" );

        m.def(
                "create_filtered_observation_collection",
                []( const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > originalObservationCollection,
                    const std::shared_ptr< tom::ObservationFilterBase > observationFilter,
                    const std::shared_ptr< tom::ObservationCollectionParser > observationParser ) {
                    warnLegacyObservationInterface( "create_filtered_observation_collection", "ObservationDataset.create_new_and_drop" );
                    return tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE >(
                            originalObservationCollection, observationFilter, observationParser );
                },
                py::arg( "original_observation_collection" ),
                py::arg( "observation_filter" ),
                py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                R"doc(

        Creates a new, filtered observation collection from an existing one using a single filter.

        This function applies a single filter to a subset of an observation collection (selected by a parser)
        and returns a new collection with the filtered observations. The original collection is not modified.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The observation collection to filter.
        observation_filter : tudatpy.estimation.observations.observations_processing.ObservationFilterBase
            The filter to apply.
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, optional
            Parser to select the subset of observations to filter. Defaults to an empty parser (all observations).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the filtered observations.

        )doc" );

        m.def(
                "split_observation_collection",
                []( const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > originalObservationCollection,
                    const std::shared_ptr< tom::ObservationSetSplitterBase > observationSetSplitter,
                    const std::shared_ptr< tom::ObservationCollectionParser > observationParser ) {
                    warnLegacyObservationInterface( "split_observation_collection", "ObservationDataset.create_new_and_keep" );
                    return tom::splitObservationSets< STATE_SCALAR_TYPE, TIME_TYPE >(
                            originalObservationCollection, observationSetSplitter, observationParser );
                },
                py::arg( "original_observation_collection" ),
                py::arg( "observation_set_splitter" ),
                py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                R"doc(
        Creates a new observation collection by splitting sets from an existing collection.

        This function splits observation sets from the original collection based on a splitter and
        returns a new collection containing the resulting sets. The original collection is not modified.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The observation collection from which to split sets.
        observation_set_splitter : tudatpy.estimation.observations.observations_processing.ObservationSetSplitterBase
            The splitter defining how to split the sets.
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, optional
            Parser to select which observation sets to split. Defaults to an empty parser (all sets).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection with the split observation sets.
        )doc" );

        m.def(
                "create_new_observation_collection",
                []( const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > originalObservationCollection,
                    const std::shared_ptr< tom::ObservationCollectionParser > observationParser ) {
                    warnLegacyObservationInterface( "create_new_observation_collection", "ObservationDataset.create_new_and_keep" );
                    return tom::createNewObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >( originalObservationCollection,
                                                                                                observationParser );
                },
                py::arg( "original_observation_collection" ),
                py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                R"doc(

        Creates a new observation collection containing a subset of an existing collection.

        This function selects a subset of observation sets from an original collection using a parser
        and creates a new collection containing only those sets.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The collection from which to extract a subset.
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, optional
            Parser to select the observation sets to include in the new collection. Defaults to an empty parser (all sets).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the selected subset of observation sets.
        )doc" );
    }
}

}  // namespace observations
}  // namespace estimation
}  // namespace tudatpy
