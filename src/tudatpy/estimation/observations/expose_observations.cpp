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

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"
#include "tudat/basics/deprecationWarnings.h"
#include "tudat/io/serialization/pybind_helpers.h"
#include "tudat/io/serialization/registrations_estimation.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/createObservationCollection.h"
#include "tudat/simulation/estimation_setup/compressDopplerObservationCollection.h"
#include "tudat/simulation/estimation_setup/simulateObservationsLegacy.h"
#include "observations_processing/expose_observations_processing.h"
#include "observations_geometry/expose_observations_geometry.h"
#include "expose_observations_bindings.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
namespace te = tudat::ephemerides;
namespace tdat = tudat::data;

namespace
{

using InspectionDataset = tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >;
using InspectionCondition = tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >;

std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > > evaluateTrackingDataAngularCorrections(
        const std::shared_ptr< tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE > >& trackingData,
        const tss::SystemOfBodies& bodies,
        const std::vector< TIME_TYPE >& epochsTdb )
{
    if( bodies.getFrameOrientation( ) != "J2000" )
    {
        throw std::runtime_error( "Angular observation corrections require a J2000 global frame orientation." );
    }
    if( bodies.getFrameOrigin( ) != "SSB" )
    {
        throw std::runtime_error( "Angular observation corrections require an SSB global frame origin." );
    }
    const auto& settings = trackingData->getObservationCorrectionSettings( ).value( );
    const auto linkEnds = tom::getLinkEndsFromTrackingData( trackingData->getLinkEnds( ) );
    const auto& target = linkEnds.at( tom::transmitter );
    const auto& observer = linkEnds.at( tom::receiver );
    if( !target.stationName_.empty( ) )
    {
        throw std::runtime_error( "Angular observation corrections require a body-centre transmitter." );
    }
    py::object observerReferenceName = py::none( );
    if( !observer.stationName_.empty( ) )
    {
        observerReferenceName = py::str( observer.stationName_ );
    }

    const auto& sourceObservations = trackingData->getObservations( );
    py::array_t< double > observations( { static_cast< py::ssize_t >( sourceObservations.size( ) ), py::ssize_t( 3 ) } );
    auto values = observations.mutable_unchecked< 2 >( );
    for( unsigned int i = 0; i < sourceObservations.size( ); ++i )
    {
        values( i, 0 ) = static_cast< double >( epochsTdb.at( i ) );
        values( i, 1 ) = sourceObservations.at( i )( 0 );
        values( i, 2 ) = sourceObservations.at( i )( 1 );
    }

    const auto module = py::module_::import( "tudatpy.estimation.observations.observation_corrections" );
    Eigen::MatrixXd corrections = Eigen::MatrixXd::Zero( sourceObservations.size( ), 2 );
    if( settings.photocenterBodyDimensions_.size( ) != 0 )
    {
        py::object dimensions = settings.photocenterBodyDimensions_.size( ) == 1
                ? py::object( py::float_( settings.photocenterBodyDimensions_( 0 ) ) )
                : py::cast( settings.photocenterBodyDimensions_ );
        corrections += module.attr( "photocenter_correction_angular_observations" )(
                                     py::arg( "observations" ) = observations,
                                     py::arg( "body_dimensions" ) = dimensions,
                                     py::arg( "bodies" ) = py::cast( &bodies, py::return_value_policy::reference ),
                                     py::arg( "body_name" ) = target.bodyName_,
                                     py::arg( "observer_body_name" ) = observer.bodyName_,
                                     py::arg( "observer_reference_name" ) = observerReferenceName )
                               .cast< Eigen::MatrixXd >( );
    }
    if( !settings.lightDeflectionBodies_.empty( ) )
    {
        corrections += module.attr( "light_deflection_correction_angular_observations" )(
                                     py::arg( "observations" ) = observations,
                                     py::arg( "bodies" ) = py::cast( &bodies, py::return_value_policy::reference ),
                                     py::arg( "body_name" ) = target.bodyName_,
                                     py::arg( "observer_body_name" ) = observer.bodyName_,
                                     py::arg( "observer_reference_name" ) = observerReferenceName,
                                     py::arg( "perturbing_bodies_list" ) = settings.lightDeflectionBodies_ )
                               .cast< Eigen::MatrixXd >( );
    }
    std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > > result;
    result.reserve( sourceObservations.size( ) );
    for( Eigen::Index i = 0; i < corrections.rows( ); ++i )
    {
        result.push_back( corrections.row( i ).transpose( ).cast< STATE_SCALAR_TYPE >( ) );
    }
    return result;
}

std::shared_ptr< InspectionDataset > createObservationDatasetFromTrackingData(
        const std::vector< std::shared_ptr< tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE > > >& trackingData,
        tss::SystemOfBodies& bodies,
        const bool applyCorrections )
{
    return tom::createObservationDatasetFromTrackingData< STATE_SCALAR_TYPE, TIME_TYPE >(
            trackingData, bodies, applyCorrections, evaluateTrackingDataAngularCorrections );
}

tom::ObservationOrdering inspectionOrdering( const std::string& ordering )
{
    if( ordering == "internal" )
    {
        return tom::ObservationOrdering::internal;
    }
    if( ordering == "estimation" )
    {
        return tom::ObservationOrdering::estimation;
    }
    throw py::value_error( "ordering must be 'internal' or 'estimation'" );
}

py::list inspectionRows( const std::vector< tom::ObservationDatasetRow< TIME_TYPE > >& rows )
{
    py::list result;
    for( const auto& row : rows )
    {
        py::dict value;
        value[ "observation_id" ] = row.observationId_;
        value[ "time" ] = py::cast( row.time_, py::return_value_policy::copy );
        value[ "set_id" ] = row.setId_;
        value[ "first_scalar_component" ] = row.firstScalarComponent_;
        value[ "scalar_size" ] = row.scalarSize_;
        value[ "index_in_set" ] = row.indexInSet_;
        value[ "dependent_variables" ] = py::cast( row.dependentVariableValues_, py::return_value_policy::copy );
        value[ "is_active" ] = row.isActive_;
        value[ "rejection_reason" ] = row.rejectionReason_;
        result.append( std::move( value ) );
    }
    return result;
}

py::dict inspectionMetadata( const InspectionDataset::InspectionMetadata& metadata )
{
    py::dict result;
    for( const auto& entry : metadata )
    {
        const auto& set = std::get< 0 >( entry.second );
        py::dict value;
        value[ "observable_type" ] = set.observableType_;
        value[ "reference_link_end" ] = set.referenceLinkEnd_;
        value[ "observable_size" ] = set.observableSize_;
        value[ "link_definition_id" ] = set.linkDefinitionId_;
        value[ "ancillary_settings_id" ] = set.ancillarySettingsId_;
        value[ "dependent_variable_layout_id" ] = set.dependentVariableLayoutId_;
        value[ "link_definition" ] = py::cast( std::get< 1 >( entry.second ), py::return_value_policy::copy );
        // These settings are copies, so that modifying them does not change the dataset.
        value[ "ancillary_settings" ] = py::cast( std::get< 2 >( entry.second ) );
        value[ "dependent_variable_layout" ] = py::cast( std::get< 3 >( entry.second ) );
        result[ py::int_( entry.first ) ] = std::move( value );
    }
    return result;
}

py::dict inspectionData( const InspectionDataset& dataset,
                         const InspectionCondition& condition,
                         const std::vector< std::string >& fields,
                         const std::string& ordering )
{
    const auto order = inspectionOrdering( ordering );
    // Check the requested fields before retrieving any data.
    const std::set< std::string > supported = { "times",   "observations",        "residuals",         "observation_ids", "set_ids",
                                                "rows",    "dependent_variables", "scalar_components", "weight_diagonal", "weight_matrix",
                                                "metadata" };
    std::set< std::string > requested;
    for( const auto& field : fields )
    {
        if( supported.count( field ) == 0 )
        {
            throw py::value_error( "Unknown observation field: " + field );
        }
        if( !requested.insert( field ).second )
        {
            throw py::value_error( "Duplicate observation field: " + field );
        }
    }
    py::dict result;
    for( const auto& field : fields )
    {
        if( field == "times" )
        {
            result[ py::str( field ) ] = py::cast( dataset.getTimes( condition, order ) );
        }
        else if( field == "observations" )
        {
            result[ py::str( field ) ] = py::cast( dataset.getObservations( condition, order ) );
        }
        else if( field == "residuals" )
        {
            result[ py::str( field ) ] = py::cast( dataset.getResiduals( condition, order ) );
        }
        else if( field == "observation_ids" )
        {
            result[ py::str( field ) ] = py::cast( dataset.getObservationIds( condition, order ) );
        }
        else if( field == "set_ids" )
        {
            result[ py::str( field ) ] = py::cast( dataset.getSetIds( condition, order ) );
        }
        else if( field == "rows" )
        {
            result[ py::str( field ) ] = inspectionRows( dataset.getRows( condition, order ) );
        }
        else if( field == "dependent_variables" )
        {
            result[ py::str( field ) ] = py::cast( dataset.getDependentVariableValues( condition, order ) );
        }
        else if( field == "scalar_components" )
        {
            result[ py::str( field ) ] = py::cast( dataset.getScalarComponents( condition, order ) );
        }
        else if( field == "weight_diagonal" )
        {
            result[ py::str( field ) ] = py::cast( dataset.getWeightDiagonal( condition, order ) );
        }
        else if( field == "weight_matrix" )
        {
            result[ py::str( field ) ] = py::cast( dataset.getWeightMatrix( condition, order ) );
        }
        else if( field == "metadata" )
        {
            result[ py::str( field ) ] = inspectionMetadata( dataset.getMetadata( condition, order ) );
        }
    }
    return result;
}

void computeDatasetResiduals(
        const std::shared_ptr< InspectionDataset >& observationDataset,
        const std::vector< std::shared_ptr< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >& observationSimulators,
        const tss::SystemOfBodies& bodies )
{
    if( observationDataset == nullptr )
    {
        throw std::runtime_error( "Error when computing residuals and dependent variables for dataset, input dataset is None." );
    }
    tss::computeResidualsAndDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >( observationDataset, observationSimulators, bodies );
}

Eigen::VectorXd castObservationSelectionConditionVectorLimit( const py::object& value )
{
    Eigen::VectorXd limit;
    if( !py::isinstance< py::str >( value ) )
    {
        const auto array = py::array_t< double, py::array::forcecast >::ensure( value );
        if( array && array.ndim( ) <= 1 )
        {
            limit = Eigen::Map< const Eigen::VectorXd >( array.data( ), array.size( ) );
        }
    }
    if( limit.size( ) == 0 )
    {
        throw py::type_error( "limit must be a number or a one-dimensional array of numbers (one per component)." );
    }
    if( limit.hasNaN( ) )
    {
        throw py::value_error( "limit must not be NaN." );
    }
    return limit;
}

template< typename ScalarType >
std::vector< Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > > castObservationVectorList( const py::object& value,
                                                                                         const std::string& argumentName,
                                                                                         const int expectedObservationCount = -1 )
{
    if( value.is_none( ) )
    {
        return std::vector< Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > >( );
    }

    py::array_t< ScalarType, py::array::c_style | py::array::forcecast > array =
            py::array_t< ScalarType, py::array::c_style | py::array::forcecast >::ensure( value );
    if( array )
    {
        const py::buffer_info buffer = array.request( );
        std::vector< Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > > result;

        if( buffer.ndim == 1 )
        {
            const auto scalarCount = static_cast< int >( buffer.shape.at( 0 ) );
            const ScalarType* data = static_cast< const ScalarType* >( buffer.ptr );
            if( scalarCount == 0 )
            {
                return result;
            }
            if( expectedObservationCount >= 0 && scalarCount == expectedObservationCount )
            {
                result.reserve( scalarCount );
                for( int i = 0; i < scalarCount; ++i )
                {
                    Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > scalarObservation( 1 );
                    scalarObservation( 0 ) = data[ i ];
                    result.push_back( scalarObservation );
                }
                return result;
            }

            Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > observation( scalarCount );
            for( int i = 0; i < scalarCount; ++i )
            {
                observation( i ) = data[ i ];
            }
            result.push_back( observation );
            return result;
        }
        else if( buffer.ndim == 2 || ( buffer.ndim == 3 && buffer.shape.at( 2 ) == 1 ) )
        {
            const auto rowCount = static_cast< int >( buffer.shape.at( 0 ) );
            const auto columnCount = static_cast< int >( buffer.shape.at( 1 ) );
            const ScalarType* data = static_cast< const ScalarType* >( buffer.ptr );
            result.reserve( rowCount );
            for( int i = 0; i < rowCount; ++i )
            {
                Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > observation( columnCount );
                for( int j = 0; j < columnCount; ++j )
                {
                    observation( j ) = data[ i * columnCount + j ];
                }
                result.push_back( observation );
            }
            return result;
        }

        throw py::type_error( argumentName +
                              " must be a 1D array, a 2D array with one observation per row, "
                              "a 3D array with trailing singleton dimension, or a sequence of vectors." );
    }

    try
    {
        return value.cast< std::vector< Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > > >( );
    }
    catch( const py::cast_error& )
    {
        throw py::type_error( argumentName + " must be convertible to a sequence of Eigen-compatible vectors." );
    }
}

template< typename TimeType >
TimeType castObservationTime( const py::object& value )
{
    if( py::isinstance< py::float_ >( value ) || py::isinstance< py::int_ >( value ) )
    {
        return TimeType( value.cast< double >( ) );
    }
    return value.cast< TimeType >( );
}

template< typename TimeType >
std::vector< TimeType > castObservationTimeList( const py::object& value, const std::string& argumentName )
{
    if( value.is_none( ) )
    {
        return std::vector< TimeType >( );
    }

    // Preserve Time objects, including their split-second precision. Only numeric
    // NumPy arrays use the contiguous double fast path; object arrays and Python
    // sequences are converted one element at a time below.
    if( py::isinstance< py::array >( value ) && py::reinterpret_borrow< py::array >( value ).dtype( ).kind( ) != 'O' )
    {
        const py::dtype dtype = py::reinterpret_borrow< py::array >( value ).dtype( );
        if( dtype.kind( ) != 'f' && dtype.kind( ) != 'i' && dtype.kind( ) != 'u' )
        {
            throw py::type_error( argumentName +
                                  " must be an array of numbers (seconds since J2000 in TDB); numpy datetime64, timedelta64 and "
                                  "boolean arrays are not supported." );
        }
        if( dtype.kind( ) == 'f' && dtype.itemsize( ) < static_cast< py::ssize_t >( sizeof( double ) ) )
        {
            throw py::type_error( argumentName + " must not be a single-precision array, which cannot represent epochs accurately." );
        }
        auto array = py::array_t< double, py::array::c_style | py::array::forcecast >::ensure( value );
        if( !array || array.ndim( ) != 1 )
        {
            throw py::type_error( argumentName + " must be a one-dimensional sequence of epochs." );
        }
        const auto buffer = array.request( );
        const double* data = static_cast< const double* >( buffer.ptr );
        std::vector< TimeType > result;
        result.reserve( buffer.size );
        for( py::ssize_t i = 0; i < buffer.size; ++i )
        {
            result.emplace_back( data[ i ] );
        }
        return result;
    }

    try
    {
        // A list owns its items (e.g. those of a pandas Series, which are created on access), and is iterated by value
        const py::list items( value );
        std::vector< TimeType > result;
        result.reserve( items.size( ) );
        for( py::handle item : items )
        {
            result.push_back( castObservationTime< TimeType >( py::reinterpret_borrow< py::object >( item ) ) );
        }
        return result;
    }
    catch( const py::cast_error& )
    {
        throw py::type_error( argumentName + " must be convertible to a sequence of Time objects or numeric epochs." );
    }
    catch( const py::error_already_set& )
    {
        throw py::type_error( argumentName + " must be convertible to a sequence of Time objects or numeric epochs." );
    }
}

}  // namespace

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
Class holding the properties that are shared by all observations in an observation set.

Objects of this class are retrieved from an :class:`~tudatpy.estimation.observations.ObservationDataset`, using its
:attr:`~tudatpy.estimation.observations.ObservationDataset.observation_set_metadata` attribute or its
:meth:`~tudatpy.estimation.observations.ObservationDataset.get_observation_set_metadata` function. The link definition,
ancillary settings and dependent variable settings of the set are stored in the dataset, and are identified here by an
integer id.
)doc" )
            .def_readonly( "observable_type",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::observableType_,
                           R"doc(
**read-only**

Observable type of the observations in the set.

:type: tudatpy.estimation.observable_models_setup.model_settings.ObservableType
)doc" )
            .def_readonly( "link_definition_id",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::linkDefinitionId_,
                           R"doc(
**read-only**

Id of the link definition of the set, which can be retrieved with
:meth:`~tudatpy.estimation.observations.ObservationDataset.link_definition`.

:type: int
)doc" )
            .def_readonly( "reference_link_end",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::referenceLinkEnd_,
                           R"doc(
**read-only**

Link end at which the observation times of the set are defined.

:type: tudatpy.estimation.observable_models_setup.links.LinkEndType
)doc" )
            .def_readonly( "observable_size",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::observableSize_,
                           R"doc(
**read-only**

Number of components of a single observation in the set (e.g. 1 for range, 2 for angular position).

:type: int
)doc" )
            .def_readonly( "ancillary_settings_id",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::ancillarySettingsId_,
                           R"doc(
**read-only**

Id of the ancillary settings of the set, which can be retrieved with
:meth:`~tudatpy.estimation.observations.ObservationDataset.ancillary_settings`.

:type: int
)doc" )
            .def_readonly( "dependent_variable_layout_id",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::dependentVariableLayoutId_,
                           R"doc(
**read-only**

Id of the object defining the dependent variables stored for the set, which can be retrieved with
:meth:`~tudatpy.estimation.observations.ObservationDataset.dependent_variable_bookkeeping`.

:type: int
)doc" );

    py::class_< tom::ObservationDatasetRow< TIME_TYPE > >( m,
                                                           "ObservationDatasetRow",
                                                           R"doc(
Class holding the properties of a single observation in an :class:`~tudatpy.estimation.observations.ObservationDataset`.

Objects of this class are retrieved with :meth:`~tudatpy.estimation.observations.ObservationDataset.observation_row`
and :attr:`~tudatpy.estimation.observations.ObservationDataset.observation_rows`.
)doc" )
            .def_readonly( "observation_id", &tom::ObservationDatasetRow< TIME_TYPE >::observationId_ )
            .def_readonly( "time",
                           &tom::ObservationDatasetRow< TIME_TYPE >::time_,
                           R"doc(
**read-only**

Time of the observation, in seconds since J2000 (TDB), defined at the reference link end of its observation set.

:type: tudatpy.astro.time_representation.Time
)doc" )
            .def_readonly( "set_id", &tom::ObservationDatasetRow< TIME_TYPE >::setId_, R"doc(
**read-only**

Id of the observation set to which the observation belongs.

:type: int
)doc" )
            .def_readonly( "first_scalar_component",
                           &tom::ObservationDatasetRow< TIME_TYPE >::firstScalarComponent_,
                           R"doc(
**read-only**

Index of the first component of the observation in the storage of all observation components in the dataset
(see :meth:`~tudatpy.estimation.observations.ObservationDataset.scalar_component_row`).

:type: int
)doc" )
            .def_readonly( "scalar_size",
                           &tom::ObservationDatasetRow< TIME_TYPE >::scalarSize_,
                           R"doc(
**read-only**

Number of components of the observation.

:type: int
)doc" )
            .def_readonly( "index_in_set",
                           &tom::ObservationDatasetRow< TIME_TYPE >::indexInSet_,
                           R"doc(
**read-only**

Index of the observation within its observation set.

:type: int
)doc" )
            .def_readonly( "is_active",
                           &tom::ObservationDatasetRow< TIME_TYPE >::isActive_,
                           R"doc(
**read-only**

Whether the observation is used in the estimation (``False`` if the observation is rejected).

:type: bool
)doc" )
            .def_readonly( "rejection_reason",
                           &tom::ObservationDatasetRow< TIME_TYPE >::rejectionReason_,
                           R"doc(
**read-only**

Reason provided when the observation was rejected (empty if no reason was provided, or if the observation is not
rejected).

:type: str
)doc" );

    py::class_< tom::ObservationScalarComponentRow >( m,
                                                      "ObservationScalarComponentRow",
                                                      R"doc(
Class identifying a single component of an observation in an :class:`~tudatpy.estimation.observations.ObservationDataset`.

Objects of this class are retrieved with
:meth:`~tudatpy.estimation.observations.ObservationDataset.scalar_component_row`.
)doc" )
            .def_readonly( "observation_id",
                           &tom::ObservationScalarComponentRow::observationId_,
                           R"doc(
**read-only**

Id of the observation to which the component belongs.

:type: int
)doc" )
            .def_readonly( "component_index",
                           &tom::ObservationScalarComponentRow::componentIndex_,
                           R"doc(
**read-only**

Index of the component within the observation (e.g. 0 for the right ascension and 1 for the declination of an
angular position observation).

:type: int
)doc" );

    {
        py::class_< tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                  "ObservationVectorData",
                                                                                  R"doc(
Class holding the concatenated observations, residuals and weights of an observation dataset, as used in the estimation.

Objects of this class are created with
:meth:`~tudatpy.estimation.observations.ObservationDataset.observation_vector_data`. Each entry of the vectors in this
object is a single component of an observation, with the observations in the ``"estimation"`` order (see
:class:`~tudatpy.estimation.observations.ObservationDataset`). For each entry, the attributes :attr:`observation_ids`, :attr:`set_ids` and :attr:`times` give
the observation, observation set and time to which it belongs. The data in this object is a copy, which is not updated
when the dataset is modified.
)doc" )
                .def_property_readonly( "observation_vector",
                                        &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationVector,
                                        R"doc(
**read-only**

Concatenated vector of all observation components.

:type: numpy.ndarray[numpy.float64[m, 1]]
)doc" )
                .def_property_readonly( "residual_vector",
                                        &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualVector,
                                        R"doc(
**read-only**

Concatenated vector of the residuals, in the same order as :attr:`observation_vector`.

:type: numpy.ndarray[numpy.float64[m, 1]]
)doc" )
                .def_property_readonly( "weight_vector",
                                        &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightVector,
                                        R"doc(
**read-only**

Diagonal of the weight matrix, in the same order as :attr:`observation_vector`.

:type: numpy.ndarray[numpy.float64[m, 1]]
)doc" )
                .def_property_readonly( "sparse_weight_matrix",
                                        &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getSparseWeightMatrix,
                                        R"doc(
**read-only**

Full weight matrix, with rows and columns in the same order as :attr:`observation_vector`. If
:attr:`has_off_diagonal_weights` is ``False``, this matrix is diagonal, and contains the same information as
:attr:`weight_vector`.

:type: scipy.sparse.csc_matrix
)doc" )
                .def( "inverse_weight_matrix_for_observation",
                      &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getInverseWeightMatrixForObservation,
                      py::arg( "observation_id" ),
                      py::return_value_policy::copy,
                      R"doc(
Function to retrieve the covariance matrix (inverse of the weight matrix) of a single observation.

If the observation is correlated with other observations in its set, this is the block of the inverse of the weight
matrix of the complete set (including rejected observations) that corresponds to the observation. All weights must
be positive. Correlations between observations in different sets are not supported.

Parameters
----------
observation_id : int
    Id of the observation, which may be any observation in the dataset at the time this object was created
    (including rejected observations).

Returns
-------
numpy.ndarray[numpy.float64[m, m]]
    Covariance matrix of the observation, with one row and column per observation component.
)doc" )
                .def_property_readonly( "has_off_diagonal_weights",
                                        &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::hasOffDiagonalWeights,
                                        R"doc(
**read-only**

Whether the weight matrix contains off-diagonal entries (correlations between observation components).

:type: bool
)doc" )
                .def_property_readonly(
                        "times",
                        []( const tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >& observationVectorData ) {
                            const std::vector< TIME_TYPE > rawTimes = observationVectorData.getTimes( );
                            std::vector< double > convertedTimes;
                            convertedTimes.reserve( rawTimes.size( ) );
                            for( const TIME_TYPE& time : rawTimes )
                            {
                                convertedTimes.push_back( static_cast< double >( time ) );
                            }
                            return convertedTimes;
                        },
                        R"doc(
**read-only**

Observation time of each entry of :attr:`observation_vector`, in seconds since J2000 (TDB).

:type: list[float]
)doc" )
                .def_property_readonly( "observation_ids",
                                        &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationIds,
                                        R"doc(
**read-only**

Id of the observation to which each entry of :attr:`observation_vector` belongs.

:type: list[int]
)doc" )
                .def_property_readonly( "set_ids",
                                        &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getSetIds,
                                        R"doc(
**read-only**

Id of the observation set to which each entry of :attr:`observation_vector` belongs.

:type: list[int]
)doc" )
                .def_property_readonly( "link_definition_ids",
                                        &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinitionIds,
                                        R"doc(
**read-only**

Id of the link definition (see :meth:`~tudatpy.estimation.observations.ObservationDataset.link_definition`) of each
entry of :attr:`observation_vector`.

:type: list[int]
)doc" )
                .def_property_readonly( "scalar_component_ids",
                                        &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getScalarComponentIds,
                                        R"doc(
**read-only**

Index of each entry of :attr:`observation_vector` in the storage of all observation components in the dataset. The
corresponding observation id and component index are given by
:meth:`~tudatpy.estimation.observations.ObservationDataset.scalar_component_row`.

:type: list[int]
)doc" )
                .def_property_readonly( "set_ids_in_row_order",
                                        &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getSetIdsInRowOrder,
                                        R"doc(
**read-only**

Ids of the observation sets in this object, in the order in which they first appear in :attr:`observation_vector`.

:type: list[int]
)doc" )
                .def( "unique_observation_ids_for_set",
                      &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getUniqueObservationIdsForSetInRowOrder,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the ids of the observations of one observation set, in the order in which they appear in
:attr:`observation_vector`.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
list[int]
    Ids of the observations of the set.
)doc" )
                .def( "vector_row",
                      &tom::ObservationVectorData< STATE_SCALAR_TYPE, TIME_TYPE >::getVectorRow,
                      py::arg( "observation_id" ),
                      py::arg( "component_index" ),
                      R"doc(
Function to retrieve the index in :attr:`observation_vector` of a single component of an observation.

Parameters
----------
observation_id : int
    Id of the observation.
component_index : int
    Index of the component within the observation.

Returns
-------
int
    Index of the entry in :attr:`observation_vector`.
)doc" );
    }

    py::enum_< InspectionCondition::Type >( m,
                                            "ObservationSelectionConditionType",
                                            R"doc(
Enumeration of the kinds of conditions that an
:class:`~tudatpy.estimation.observations.ObservationSelectionCondition` can represent, see
:attr:`~tudatpy.estimation.observations.ObservationSelectionCondition.condition_type`. The entries correspond to the
functions of :class:`~tudatpy.estimation.observations.ObservationSelectionCondition` with the same name, with
``and_condition``, ``or_condition`` and ``not_condition`` denoting conditions that are created with ``&``, ``|`` and
``~``, respectively.
)doc" )
            .value( "all", InspectionCondition::Type::all )
            .value( "observable_type", InspectionCondition::Type::observable_type )
            .value( "link_definition", InspectionCondition::Type::link_definition )
            .value( "link_end_type", InspectionCondition::Type::link_end_type )
            .value( "link_end", InspectionCondition::Type::link_end )
            .value( "set_id", InspectionCondition::Type::set_id )
            .value( "time_greater_equal", InspectionCondition::Type::time_greater_equal )
            .value( "time_greater_than", InspectionCondition::Type::time_greater_than )
            .value( "time_less_equal", InspectionCondition::Type::time_less_equal )
            .value( "time_less_than", InspectionCondition::Type::time_less_than )
            .value( "active", InspectionCondition::Type::active )
            .value( "residual_absolute_value_greater_than", InspectionCondition::Type::residual_absolute_value_greater_than )
            .value( "observation_absolute_value_greater_than", InspectionCondition::Type::observation_absolute_value_greater_than )
            .value( "dependent_variable_greater_than", InspectionCondition::Type::dependent_variable_greater_than )
            .value( "dependent_variable_less_than", InspectionCondition::Type::dependent_variable_less_than )
            .value( "and_condition", InspectionCondition::Type::and_condition )
            .value( "or_condition", InspectionCondition::Type::or_condition )
            .value( "not_condition", InspectionCondition::Type::not_condition )
            .value( "custom", InspectionCondition::Type::custom );

    py::class_< InspectionCondition >( m,
                                       "ObservationSelectionCondition",
                                       R"doc(
Class defining a condition that selects observations in an :class:`~tudatpy.estimation.observations.ObservationDataset`.

A condition is evaluated for each individual observation (not for complete observation sets). Conditions are normally
created with :data:`~tudatpy.estimation.observations.observation_query`, but can also be created with the static
functions of this class. Conditions can be combined with ``&`` (and) and ``|`` (or), and negated with ``~``.

The contents of a condition can be retrieved from its ``condition_type``, from its ``child_conditions`` (for conditions
that are combined with ``&``, ``|`` or ``~``), and from the attributes holding the value that it compares to.
)doc" )
            .def( py::init<>( ), R"doc(Class constructor, creating a condition that selects all observations.)doc" )
            .def_static( "all", &InspectionCondition::all, R"doc(
Function to create a condition that selects all observations.

Returns
-------
ObservationSelectionCondition
    Condition selecting all observations.
)doc" )
            .def_static( "observable_type",
                         &InspectionCondition::observableType,
                         py::arg( "observable_type" ),
                         R"doc(
Function to create a condition that selects the observations of a given observable type.

Parameters
----------
observable_type : tudatpy.estimation.observable_models_setup.model_settings.ObservableType
    Observable type of the observations to select.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations of the given observable type.
)doc" )
            .def_static( "link_definition",
                         &InspectionCondition::linkDefinition,
                         py::arg( "link_definition" ),
                         R"doc(
Function to create a condition that selects the observations with a given link definition.

Parameters
----------
link_definition : tudatpy.estimation.observable_models_setup.links.LinkDefinition
    Link definition (all link ends) of the observations to select.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with the given link definition.
)doc" )
            .def_static( "link_end_type",
                         &InspectionCondition::linkEndType,
                         py::arg( "link_end_type" ),
                         R"doc(
Function to create a condition that selects the observations that have a link end with a given role, for any body and
reference point of this link end.

Parameters
----------
link_end_type : tudatpy.estimation.observable_models_setup.links.LinkEndType
    Role of the link end in the observation (e.g. receiver or transmitter).

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations that have a link end with the given role.
)doc" )
            .def_static( "link_end",
                         &InspectionCondition::linkEnd,
                         py::arg( "link_end_type" ),
                         py::arg( "link_end_id" ),
                         R"doc(
Function to create a condition that selects the observations with a given link end.

Parameters
----------
link_end_type : tudatpy.estimation.observable_models_setup.links.LinkEndType
    Role of the link end in the observation (e.g. receiver or transmitter).
link_end_id : tudatpy.estimation.observable_models_setup.links.LinkEndId
    Body and reference point of the link end.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with the given link end.
)doc" )
            .def_static( "set_id",
                         &InspectionCondition::setId,
                         py::arg( "set_id" ),
                         R"doc(
Function to create a condition that selects the observations of a given observation set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations of the given observation set.
)doc" )
            .def_static( "time_bounds",
                         &InspectionCondition::timeBounds,
                         py::arg( "start_time" ),
                         py::arg( "end_time" ),
                         R"doc(
Function to create a condition that selects the observations with a time between (and including) two bounds.

Parameters
----------
start_time : tudatpy.astro.time_representation.Time
    Lower bound of the observation time, in seconds since J2000 (TDB).
end_time : tudatpy.astro.time_representation.Time
    Upper bound of the observation time, in seconds since J2000 (TDB).

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations within the time bounds.
)doc" )
            .def_static( "time_greater_equal",
                         &InspectionCondition::timeGreaterEqual,
                         py::arg( "time" ),
                         R"doc(
Function to create a condition that selects the observations with a time greater than or equal to a given value.

Parameters
----------
time : tudatpy.astro.time_representation.Time
    Limit of the observation time, in seconds since J2000 (TDB).

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with a time greater than or equal to ``time``.
)doc" )
            .def_static( "time_greater_than",
                         &InspectionCondition::timeGreaterThan,
                         py::arg( "time" ),
                         R"doc(
Function to create a condition that selects the observations with a time greater than a given value.

Parameters
----------
time : tudatpy.astro.time_representation.Time
    Limit of the observation time, in seconds since J2000 (TDB).

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with a time greater than ``time``.
)doc" )
            .def_static( "time_less_equal",
                         &InspectionCondition::timeLessEqual,
                         py::arg( "time" ),
                         R"doc(
Function to create a condition that selects the observations with a time less than or equal to a given value.

Parameters
----------
time : tudatpy.astro.time_representation.Time
    Limit of the observation time, in seconds since J2000 (TDB).

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with a time less than or equal to ``time``.
)doc" )
            .def_static( "time_less_than",
                         &InspectionCondition::timeLessThan,
                         py::arg( "time" ),
                         R"doc(
Function to create a condition that selects the observations with a time less than a given value.

Parameters
----------
time : tudatpy.astro.time_representation.Time
    Limit of the observation time, in seconds since J2000 (TDB).

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with a time less than ``time``.
)doc" )
            .def_static( "active", &InspectionCondition::active, R"doc(
Function to create a condition that selects the observations that are not rejected.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations that are not rejected.
)doc" )
            .def_static( "rejected", &InspectionCondition::rejected, R"doc(
Function to create a condition that selects the rejected observations.

Returns
-------
ObservationSelectionCondition
    Condition selecting the rejected observations.
)doc" )
            .def_static(
                    "residual_absolute_value_greater_than",
                    []( const py::object& limit ) {
                        return InspectionCondition::residualAbsoluteValueGreaterThan(
                                castObservationSelectionConditionVectorLimit( limit ) );
                    },
                    py::arg( "limit" ),
                    R"doc(
Function to create a condition that selects the observations for which the absolute value of any component of the
residual is larger than a limit.

Parameters
----------
limit : float or numpy.ndarray[numpy.float64[m, 1]]
    Limit on the absolute value of the residual, either as a single value for all components, or with one value per
    observation component.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with a residual larger than the limit.
)doc" )
            .def_static(
                    "observation_absolute_value_greater_than",
                    []( const py::object& limit ) {
                        return InspectionCondition::observationAbsoluteValueGreaterThan(
                                castObservationSelectionConditionVectorLimit( limit ) );
                    },
                    py::arg( "limit" ),
                    R"doc(
Function to create a condition that selects the observations for which the absolute value of any component of the
observation is larger than a limit.

Parameters
----------
limit : float or numpy.ndarray[numpy.float64[m, 1]]
    Limit on the absolute value of the observation, either as a single value for all components, or with one value
    per observation component.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with a value larger than the limit.
)doc" )
            .def_static(
                    "dependent_variable_greater_than",
                    []( const std::shared_ptr< tss::ObservationDependentVariableSettings >& dependentVariableSettings,
                        const py::object& limit,
                        const bool returnFirstCompatibleSettings ) {
                        return InspectionCondition::dependentVariableGreaterThan( dependentVariableSettings,
                                                                                  castObservationSelectionConditionVectorLimit( limit ),
                                                                                  returnFirstCompatibleSettings );
                    },
                    py::arg( "dependent_variable_settings" ),
                    py::arg( "limit" ),
                    py::arg( "return_first_compatible_settings" ) = false,
                    R"doc(
Function to create a condition that selects the observations for which any component of a dependent variable is
larger than a limit.

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    Settings of the dependent variable. Observations in sets for which no compatible dependent variable is stored
    are not selected.
limit : float or numpy.ndarray[numpy.float64[m, 1]]
    Limit on the (signed) value of the dependent variable, either as a single value for all components, or with one
    value per dependent variable component.
return_first_compatible_settings : bool, default = False
    Whether to use the first compatible dependent variable if a set contains several dependent variables that are
    compatible with ``dependent_variable_settings``. If ``False``, an error is raised in this case.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with a dependent variable larger than the limit.
)doc" )
            .def_static( "rejection_reason",
                         &InspectionCondition::rejectionReason,
                         py::arg( "reason" ),
                         R"doc(
Function to create a condition that selects the observations with a given rejection reason (see
:attr:`ObservationDatasetRow.rejection_reason`). The outlier rejection during an estimation uses the reasons
``"simple_outlier_rejection"`` and ``"carpino_outlier_rejection"``.

Parameters
----------
reason : str
    Rejection reason. An empty string also selects the observations that are not rejected.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with the given rejection reason.
)doc" )
            .def_static( "observation_ids",
                         &InspectionCondition::observationIds,
                         py::arg( "observation_ids" ),
                         R"doc(
Function to create a condition that selects the observations with the given ids (see
:meth:`ObservationDataset.get_observation_ids`).

Parameters
----------
observation_ids : list[int]
    Ids of the observations to select. Ids that do not exist in the dataset are ignored.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with the given ids.
)doc" )
            .def_static(
                    "dependent_variable_less_than",
                    []( const std::shared_ptr< tss::ObservationDependentVariableSettings >& dependentVariableSettings,
                        const py::object& limit,
                        const bool returnFirstCompatibleSettings ) {
                        return InspectionCondition::dependentVariableLessThan( dependentVariableSettings,
                                                                               castObservationSelectionConditionVectorLimit( limit ),
                                                                               returnFirstCompatibleSettings );
                    },
                    py::arg( "dependent_variable_settings" ),
                    py::arg( "limit" ),
                    py::arg( "return_first_compatible_settings" ) = false,
                    R"doc(
Function to create a condition that selects the observations for which any component of a dependent variable is
smaller than a limit.

Unlike the negation of :meth:`dependent_variable_greater_than`, this condition does not select observations in
sets for which no compatible dependent variable is stored.

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    Settings of the dependent variable. Observations in sets for which no compatible dependent variable is stored
    are not selected.
limit : float or numpy.ndarray[numpy.float64[m, 1]]
    Limit on the (signed) value of the dependent variable, either as a single value for all components, or with one
    value per dependent variable component.
return_first_compatible_settings : bool, default = False
    Whether to use the first compatible dependent variable if a set contains several dependent variables that are
    compatible with ``dependent_variable_settings``. If ``False``, an error is raised in this case.

Returns
-------
ObservationSelectionCondition
    Condition selecting the observations with a dependent variable smaller than the limit.
)doc" )
            .def_property_readonly( "condition_type", &InspectionCondition::getType, R"doc(
**read-only**

Kind of the condition, for instance a condition on the observable type, or the combination of two conditions with
``&``.

:type: ObservationSelectionConditionType
)doc" )
            .def_property_readonly( "child_conditions", &InspectionCondition::getChildConditions, R"doc(
**read-only**

Conditions that are combined by this condition: two for an ``and_condition`` or ``or_condition`` (created with ``&``
or ``|``), one for a ``not_condition`` (created with ``~``), and none for other kinds of conditions.

:type: list[ObservationSelectionCondition]
)doc" )
            .def_property_readonly( "observable_type_value", &InspectionCondition::getObservableType, R"doc(
**read-only**

Observable type that is selected by an ``observable_type`` condition (not used by other kinds of conditions).

:type: tudatpy.estimation.observable_models_setup.model_settings.ObservableType
)doc" )
            .def_property_readonly( "link_definition_value", &InspectionCondition::getLinkDefinition, R"doc(
**read-only**

Link definition that is selected by a ``link_definition`` condition (not used by other kinds of conditions).

:type: tudatpy.estimation.observable_models_setup.links.LinkDefinition
)doc" )
            .def_property_readonly( "link_end_type_value", &InspectionCondition::getLinkEndType, R"doc(
**read-only**

Role of the link end that is selected by a ``link_end_type`` or ``link_end`` condition (not used by other kinds of
conditions).

:type: tudatpy.estimation.observable_models_setup.links.LinkEndType
)doc" )
            .def_property_readonly( "link_end_id_value", &InspectionCondition::getLinkEndId, R"doc(
**read-only**

Body and reference point of the link end that is selected by a ``link_end`` condition (not used by other kinds of
conditions).

:type: tudatpy.estimation.observable_models_setup.links.LinkEndId
)doc" )
            .def_property_readonly( "set_id_value", &InspectionCondition::getSetId, R"doc(
**read-only**

Id of the observation set that is selected by a ``set_id`` condition (not used by other kinds of conditions).

:type: int
)doc" )
            .def_property_readonly( "time_value", &InspectionCondition::getTimeValue, R"doc(
**read-only**

Time (in seconds since J2000, TDB) to which the observation time is compared by a ``time_greater_equal``,
``time_greater_than``, ``time_less_equal`` or ``time_less_than`` condition (not used by other kinds of conditions).

:type: tudatpy.astro.time_representation.Time
)doc" )
            .def_property_readonly( "vector_limit", &InspectionCondition::getVectorLimit, R"doc(
**read-only**

Limit used by a ``residual_absolute_value_greater_than``, ``observation_absolute_value_greater_than`` or
dependent-variable condition (not used by other kinds of conditions).

:type: numpy.ndarray[numpy.float64[m, 1]]
)doc" )
            .def_property_readonly( "dependent_variable_settings", &InspectionCondition::getDependentVariableSettings, R"doc(
**read-only**

Settings of the dependent variable used by a dependent-variable condition (``None`` for other kinds
of conditions).

:type: tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
)doc" )
            .def_property_readonly( "return_first_compatible_settings", &InspectionCondition::getReturnFirstCompatibleSettings, R"doc(
**read-only**

Whether a dependent-variable condition uses the first compatible dependent variable if a set
contains several (not used by other kinds of conditions).

:type: bool
)doc" )
            .def(
                    "__and__",
                    []( const InspectionCondition& lhs, const InspectionCondition& rhs ) { return lhs && rhs; },
                    py::is_operator( ),
                    R"doc(Returns a condition that selects the observations selected by both conditions (used through ``&``).)doc" )
            .def(
                    "__or__",
                    []( const InspectionCondition& lhs, const InspectionCondition& rhs ) { return lhs || rhs; },
                    py::is_operator( ),
                    R"doc(Returns a condition that selects the observations selected by either condition (used through ``|``).)doc" )
            .def(
                    "__invert__",
                    []( const InspectionCondition& condition ) { return !condition; },
                    py::is_operator( ),
                    R"doc(Returns a condition that selects the observations not selected by this condition (used through ``~``).)doc" )
            .def(
                    "__bool__",
                    []( const InspectionCondition& ) {
                        throw py::type_error(
                                "ObservationSelectionCondition objects cannot be converted to bool. Use &, |, and ~ to combine conditions, "
                                "not and/or/not; use observation_query.time.between( t0, t1 ) instead of a chained comparison, and "
                                ".isin( [...] ) instead of 'in'." );
                    },
                    R"doc(Raises an error, since conditions cannot be converted to a boolean: combine them with ``&``, ``|`` and ``~`` instead of ``and``, ``or`` and ``not``.)doc" );

    {
        py::class_< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >,
                    std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                  "ObservationDataset",
                                                                                                  R"doc(
Class holding a set of observations and associated data, for use in an estimation.

An observation dataset stores observations together with their times, residuals, weights and dependent variables. Each
observation is a measurement at a single time, which may have several components (e.g. right ascension and declination).
Observations with the same observable type, link ends, reference link end and ancillary settings are grouped into
observation sets. Observations and observation sets are identified by integer ids, which do not change when
observations are sorted, rejected or removed. When using simulated data, objects of this class are created with
:func:`~tudatpy.estimation.observations.simulate_observation_dataset`; real data is loaded with
:func:`~tudatpy.estimation.observations.create_observation_dataset_from_tracking_data`.

Most functions of this class act on the observations selected by a ``condition``, created with
:data:`~tudatpy.estimation.observations.observation_query`. Observations can be rejected (see
:meth:`reject_observations`): rejected observations remain in the dataset, but are not used in the estimation. The
weights may be diagonal, or may contain correlations between the components of an observation, between observations
in a set, or between observations in different sets. If the weights correlate an observation that is rejected with one
that is not, the estimation uses the rows and columns of the weight matrix that belong to the observations that are not
rejected; the weights of these observations are not recomputed from the covariance of the remaining observations.

Functions that retrieve data take an ``ordering`` input, which is one of:

* ``"internal"``: the order in which the observations were added to the dataset;
* ``"estimation"``: the order used in the estimation (and in :meth:`observation_vector_data`), in which the observations
  are sorted by observable type, then by link ends, then by observation set, and then by their order within the set.

The data returned by these functions is a copy, so that modifying it does not change the dataset.

This class replaces the (deprecated) :class:`~tudatpy.estimation.observations.ObservationCollection` and
:class:`~tudatpy.estimation.observations.SingleObservationSet` classes.
)doc" )
                .def( py::init<>( ), R"doc(Class constructor, creating an empty observation dataset.)doc" )
                .def(
                        "get_times",
                        []( const InspectionDataset& self, const InspectionCondition& condition, const std::string& ordering ) {
                            return self.getTimes( condition, inspectionOrdering( ordering ) );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        R"doc(
Function to retrieve the times of the selected observations.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).

Returns
-------
list[tudatpy.astro.time_representation.Time]
    Time of each selected observation, in seconds since J2000 (TDB).
)doc" )
                .def(
                        "get_observations",
                        []( const InspectionDataset& self, const InspectionCondition& condition, const std::string& ordering ) {
                            return self.getObservations( condition, inspectionOrdering( ordering ) );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        R"doc(
Function to retrieve the values of the selected observations.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).

Returns
-------
list[numpy.ndarray[numpy.float64[m, 1]]]
    Value of each selected observation, with one entry per observation component.
)doc" )
                .def(
                        "get_residuals",
                        []( const InspectionDataset& self, const InspectionCondition& condition, const std::string& ordering ) {
                            return self.getResiduals( condition, inspectionOrdering( ordering ) );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        R"doc(
Function to retrieve the residuals of the selected observations.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).

Returns
-------
list[numpy.ndarray[numpy.float64[m, 1]]]
    Residual of each selected observation, with one entry per observation component (zero if not computed).
)doc" )
                .def(
                        "get_observation_ids",
                        []( const InspectionDataset& self, const InspectionCondition& condition, const std::string& ordering ) {
                            return self.getObservationIds( condition, inspectionOrdering( ordering ) );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        R"doc(
Function to retrieve the ids of the selected observations.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).

Returns
-------
list[int]
    Id of each selected observation.
)doc" )
                .def(
                        "get_set_ids",
                        []( const InspectionDataset& self, const InspectionCondition& condition, const std::string& ordering ) {
                            return self.getSetIds( condition, inspectionOrdering( ordering ) );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        R"doc(
Function to retrieve the observation set ids of the selected observations.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).

Returns
-------
list[int]
    Id of the observation set of each selected observation (so a set id is repeated for every selected observation of
    the set; use, e.g., ``numpy.unique`` for the distinct set ids).
)doc" )
                .def(
                        "get_rows",
                        []( const InspectionDataset& self, const InspectionCondition& condition, const std::string& ordering ) {
                            return inspectionRows( self.getRows( condition, inspectionOrdering( ordering ) ) );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        R"doc(
Function to retrieve the properties of the selected observations, as dictionaries.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).

Returns
-------
list[dict]
    Dictionary for each selected observation, with the keys ``observation_id``, ``time``, ``set_id``,
    ``first_scalar_component``, ``scalar_size``, ``index_in_set``, ``dependent_variables``, ``is_active`` and
    ``rejection_reason`` (see :class:`~tudatpy.estimation.observations.ObservationDatasetRow`).
)doc" )
                .def(
                        "get_dependent_variables",
                        []( const InspectionDataset& self,
                            const InspectionCondition& condition,
                            const std::string& ordering,
                            const std::shared_ptr< tss::ObservationDependentVariableSettings >& dependentVariableSettings,
                            const bool returnFirstCompatibleSettings ) {
                            return self.getDependentVariableValues(
                                    condition, inspectionOrdering( ordering ), dependentVariableSettings, returnFirstCompatibleSettings );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        py::arg( "dependent_variable_settings" ) = nullptr,
                        py::arg( "return_first_compatible_settings" ) = false,
                        R"doc(
Function to retrieve the dependent variables of the selected observations.

By default, all dependent variables of each observation are retrieved. To retrieve a single dependent variable, its
settings are provided as ``dependent_variable_settings``.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings, default = None
    Settings of the single dependent variable that is to be retrieved. If ``None``, all dependent variables are
    retrieved.
return_first_compatible_settings : bool, default = False
    Whether to use the first compatible dependent variable if a set contains several dependent variables that are
    compatible with ``dependent_variable_settings``. If ``False``, an error is raised in this case.

Returns
-------
list[numpy.ndarray[numpy.float64[m, 1]]]
    Dependent variables of each selected observation, in the same order as the output of :meth:`get_times` for the
    same ``condition`` and ``ordering``. An entry is empty if no dependent variables are stored for the observation
    or, when ``dependent_variable_settings`` is provided, if its observation set has no compatible dependent variable.
    An error is raised if ``dependent_variable_settings`` is provided, and the dependent variable is defined for an
    observation set but its values have not been computed (see
    :func:`~tudatpy.estimation.observations.compute_residuals_and_dependent_variables`).
    If ``dependent_variable_settings`` is ``None``, the content of these vectors is given by the
    ``dependent_variable_layout`` entry of :meth:`get_metadata`.
)doc" )
                .def(
                        "get_scalar_components",
                        []( const InspectionDataset& self, const InspectionCondition& condition, const std::string& ordering ) {
                            return self.getScalarComponents( condition, inspectionOrdering( ordering ) );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        R"doc(
Function to retrieve the observation id and component index of each component of the selected observations.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).

Returns
-------
list[tuple[int, int]]
    Observation id and component index of each component. This is the order of the entries of
    :meth:`get_weight_diagonal`, and of the rows and columns of :meth:`get_weight_matrix`.
)doc" )
                .def(
                        "get_weight_diagonal",
                        []( const InspectionDataset& self, const InspectionCondition& condition, const std::string& ordering ) {
                            return self.getWeightDiagonal( condition, inspectionOrdering( ordering ) );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        R"doc(
Function to retrieve the diagonal of the weight matrix of the selected observations.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).

Returns
-------
numpy.ndarray[numpy.float64[m, 1]]
    Weight of each component of the selected observations, in the order given by :meth:`get_scalar_components`.
)doc" )
                .def(
                        "get_weight_matrix",
                        []( const InspectionDataset& self, const InspectionCondition& condition, const std::string& ordering ) {
                            return self.getWeightMatrix( condition, inspectionOrdering( ordering ) );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        R"doc(
Function to retrieve the weight matrix of the selected observations, including the correlations between them.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).

Returns
-------
scipy.sparse.csc_matrix
    Weight matrix, with rows and columns in the order given by :meth:`get_scalar_components`.
)doc" )
                .def(
                        "get_metadata",
                        []( const InspectionDataset& self, const InspectionCondition& condition, const std::string& ordering ) {
                            return inspectionMetadata( self.getMetadata( condition, inspectionOrdering( ordering ) ) );
                        },
                        py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                        py::arg( "ordering" ) = "internal",
                        R"doc(
Function to retrieve the properties of the observation sets that contain the selected observations.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
ordering : str, default = "internal"
    Not used by this function, since the result is sorted by set id.

Returns
-------
dict[int, dict]
    Dictionary with the set id as key and, as value, a dictionary with the keys ``observable_type``,
    ``reference_link_end``, ``observable_size``, ``link_definition_id``, ``ancillary_settings_id``,
    ``dependent_variable_layout_id`` (see :class:`~tudatpy.estimation.observations.ObservationSetMetadata`),
    ``link_definition``, ``ancillary_settings`` (``None`` if the set has none) and ``dependent_variable_layout``. The
    latter is a dictionary with, for each dependent variable, its start index and size in the dependent variable vector
    as key and its settings as value.
)doc" )
                .def( "get_data",
                      &inspectionData,
                      py::arg_v( "condition", InspectionCondition::all( ), "observation_query.all" ),
                      py::arg( "fields" ) = std::vector< std::string >{ "times", "observations" },
                      py::arg( "ordering" ) = "internal",
                      R"doc(
Function to retrieve several types of data of the selected observations at once.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations (see :data:`~tudatpy.estimation.observations.observation_query`). By default,
    all observations are selected, including rejected observations.
fields : list[str], default = ["times", "observations"]
    Data to retrieve, each given by the name of a ``get_...`` function of this class without the ``get_`` prefix:
    ``"times"``, ``"observations"``, ``"residuals"``, ``"observation_ids"``, ``"set_ids"``, ``"rows"``,
    ``"dependent_variables"``, ``"scalar_components"``, ``"weight_diagonal"``, ``"weight_matrix"`` or ``"metadata"``.
    Each field may be given only once.
ordering : str, default = "internal"
    Order of the returned data, either ``"internal"`` or ``"estimation"`` (see :class:`ObservationDataset`).

Returns
-------
dict[str, object]
    Dictionary with the requested fields as keys, and the output of the corresponding ``get_...`` function as values.
)doc" )
                .def(
                        "add_observation_set",
                        []( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& self,
                            const tom::ObservableType observableType,
                            const tom::LinkDefinition& linkDefinition,
                            py::object observations,
                            py::object times,
                            const tom::LinkEndType referenceLinkEnd,
                            py::object dependentVariables,
                            std::shared_ptr< tss::ObservationDependentVariableBookkeeping > dependentVariableBookkeeping,
                            std::shared_ptr< tom::ObservationAncillarySimulationSettings > ancillarySettings,
                            py::object weights,
                            py::object residuals,
                            const bool sortObservations,
                            const bool eraseDuplicateObservations ) {
                            const std::vector< TIME_TYPE > convertedTimes = castObservationTimeList< TIME_TYPE >( times, "times" );
                            const int observationCount = static_cast< int >( convertedTimes.size( ) );
                            return self.addObservationSet(
                                    observableType,
                                    linkDefinition,
                                    castObservationVectorList< STATE_SCALAR_TYPE >( observations, "observations", observationCount ),
                                    convertedTimes,
                                    referenceLinkEnd,
                                    castObservationVectorList< double >( dependentVariables, "dependent_variables", observationCount ),
                                    dependentVariableBookkeeping,
                                    ancillarySettings,
                                    castObservationVectorList< double >( weights, "weights", observationCount ),
                                    castObservationVectorList< STATE_SCALAR_TYPE >( residuals, "residuals", observationCount ),
                                    sortObservations,
                                    eraseDuplicateObservations );
                        },
                        py::arg( "observable_type" ),
                        py::arg( "link_definition" ),
                        py::arg( "observations" ),
                        py::arg( "times" ),
                        py::arg( "reference_link_end" ),
                        py::arg( "dependent_variables" ) = py::list( ),
                        py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                        py::arg( "ancillary_settings" ) = nullptr,
                        py::arg( "weights" ) = py::list( ),
                        py::arg( "residuals" ) = py::list( ),
                        py::arg( "sort_observations" ) = false,
                        py::arg( "erase_duplicate_observations" ) = false,
                        R"doc(
Function to add an observation set to the dataset.

Parameters
----------
observable_type : tudatpy.estimation.observable_models_setup.model_settings.ObservableType
    Observable type of the observations.
link_definition : tudatpy.estimation.observable_models_setup.links.LinkDefinition
    Link ends of the observations.
observations : list[numpy.ndarray[numpy.float64[m, 1]]] or numpy.ndarray[numpy.float64[n, m]]
    Values of the observations, as a list with one vector per observation, or as an array with one row per
    observation. For an observable with a single component, a one-dimensional array with one value per observation
    may also be used.
times : list[tudatpy.astro.time_representation.Time] or numpy.ndarray[numpy.float64[n]]
    Times of the observations, in seconds since J2000 (TDB).
reference_link_end : tudatpy.estimation.observable_models_setup.links.LinkEndType
    Link end at which the observation times are defined.
dependent_variables : list[numpy.ndarray[numpy.float64[k, 1]]], default = []
    Dependent variables, with one vector per observation (empty if none are provided).
dependent_variable_bookkeeping : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableBookkeeping, default = None
    Object defining which dependent variables are stored in ``dependent_variables``.
ancillary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings, default = None
    Ancillary settings of the observations.
weights : list[numpy.ndarray[numpy.float64[m, 1]]], default = []
    Weights, with one vector per observation, containing the weight of each observation component. If empty, all
    weights are set to 1.
residuals : list[numpy.ndarray[numpy.float64[m, 1]]], default = []
    Residuals, with one vector per observation. If empty, all residuals are set to 0.
sort_observations : bool, default = False
    Whether to sort the observations by time.
erase_duplicate_observations : bool, default = False
    Whether to remove observations with the same time as a preceding observation in the set.

Returns
-------
int
    Id of the new observation set.
)doc" )
                .def(
                        "add_observation_set",
                        []( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& self,
                            const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >& set ) {
                            return tom::addSingleObservationSetToDataset( self, set );
                        },
                        py::arg( "single_observation_set" ).none( false ),
                        R"doc(
Function to add a copy of a (deprecated) :class:`~tudatpy.estimation.observations.SingleObservationSet` to the dataset,
as a new observation set.

Parameters
----------
single_observation_set : tudatpy.estimation.observations.SingleObservationSet
    Observation set to add.

Returns
-------
int
    Id of the new observation set.
)doc" )
                .def( "add_observation_set_from_dataset",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSetFromDataset,
                      py::arg( "source_dataset" ),
                      py::arg( "source_set_id" ),
                      R"doc(
Function to add a copy of an observation set of another dataset to this dataset.

The observations are copied together with their weights (including correlations within the set), residuals, dependent
variables, rejection status and numerical metadata, and receive new observation ids. Correlations with observations in
other sets are not copied.

Parameters
----------
source_dataset : tudatpy.estimation.observations.ObservationDataset
    Dataset from which the observation set is copied.
source_set_id : int
    Id of the observation set in ``source_dataset``.

Returns
-------
int
    Id of the new observation set in this dataset.
)doc" )
                .def(
                        "add_observation_set_with_weights",
                        []( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& self,
                            const tom::ObservableType observableType,
                            const tom::LinkDefinition& linkDefinition,
                            py::object observations,
                            py::object times,
                            const tom::LinkEndType referenceLinkEnd,
                            const tom::ObservationWeightSettings& weightSettings,
                            py::object dependentVariables,
                            std::shared_ptr< tss::ObservationDependentVariableBookkeeping > dependentVariableBookkeeping,
                            std::shared_ptr< tom::ObservationAncillarySimulationSettings > ancillarySettings,
                            py::object residuals ) {
                            const std::vector< TIME_TYPE > convertedTimes = castObservationTimeList< TIME_TYPE >( times, "times" );
                            const int observationCount = static_cast< int >( convertedTimes.size( ) );
                            return self.addObservationSetWithWeights(
                                    observableType,
                                    linkDefinition,
                                    castObservationVectorList< STATE_SCALAR_TYPE >( observations, "observations", observationCount ),
                                    convertedTimes,
                                    referenceLinkEnd,
                                    weightSettings,
                                    castObservationVectorList< double >( dependentVariables, "dependent_variables", observationCount ),
                                    dependentVariableBookkeeping,
                                    ancillarySettings,
                                    castObservationVectorList< STATE_SCALAR_TYPE >( residuals, "residuals", observationCount ) );
                        },
                        py::arg( "observable_type" ),
                        py::arg( "link_definition" ),
                        py::arg( "observations" ),
                        py::arg( "times" ),
                        py::arg( "reference_link_end" ),
                        py::arg( "weight_settings" ),
                        py::arg( "dependent_variables" ) = py::list( ),
                        py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                        py::arg( "ancillary_settings" ) = nullptr,
                        py::arg( "residuals" ) = py::list( ),
                        R"doc(
Function to add an observation set to the dataset, with weights defined by weight settings.

Same as :meth:`add_observation_set`, but with the weights (which may include correlations between the components of an
observation, or between the observations in the set) defined by an
:class:`~tudatpy.data_input.tracking_data.ObservationWeightSettings` object.

Parameters
----------
observable_type : tudatpy.estimation.observable_models_setup.model_settings.ObservableType
    Observable type of the observations.
link_definition : tudatpy.estimation.observable_models_setup.links.LinkDefinition
    Link ends of the observations.
observations : list[numpy.ndarray[numpy.float64[m, 1]]] or numpy.ndarray[numpy.float64[n, m]]
    Values of the observations (see :meth:`add_observation_set`).
times : list[tudatpy.astro.time_representation.Time] or numpy.ndarray[numpy.float64[n]]
    Times of the observations, in seconds since J2000 (TDB).
reference_link_end : tudatpy.estimation.observable_models_setup.links.LinkEndType
    Link end at which the observation times are defined.
weight_settings : tudatpy.data_input.tracking_data.ObservationWeightSettings
    Settings defining the weights of the observations.
dependent_variables : list[numpy.ndarray[numpy.float64[k, 1]]], default = []
    Dependent variables, with one vector per observation (empty if none are provided).
dependent_variable_bookkeeping : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableBookkeeping, default = None
    Object defining which dependent variables are stored in ``dependent_variables``.
ancillary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings, default = None
    Ancillary settings of the observations.
residuals : list[numpy.ndarray[numpy.float64[m, 1]]], default = []
    Residuals, with one vector per observation. If empty, all residuals are set to 0.

Returns
-------
int
    Id of the new observation set.
)doc" )
                .def(
                        "set_observations_for_set",
                        []( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& self,
                            const unsigned int setId,
                            py::object observations ) {
                            self.setObservationsForSet(
                                    setId,
                                    castObservationVectorList< STATE_SCALAR_TYPE >(
                                            observations, "observations", self.getNumberOfObservationsForSet( setId ) ) );
                        },
                        py::arg( "set_id" ),
                        py::arg( "observations" ),
                        R"doc(
Function to replace the values of all observations in an observation set.

Parameters
----------
set_id : int
    Id of the observation set.
observations : list[numpy.ndarray[numpy.float64[m, 1]]] or numpy.ndarray[numpy.float64[n, m]]
    New values of the observations, in their order within the set (see :meth:`add_observation_set`).
)doc" )
                .def(
                        "set_residuals_for_set",
                        []( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& self,
                            const unsigned int setId,
                            py::object residuals ) {
                            self.setResidualsForSet( setId,
                                                     castObservationVectorList< STATE_SCALAR_TYPE >(
                                                             residuals, "residuals", self.getNumberOfObservationsForSet( setId ) ) );
                        },
                        py::arg( "set_id" ),
                        py::arg( "residuals" ),
                        R"doc(
Function to replace the residuals of all observations in an observation set.

Parameters
----------
set_id : int
    Id of the observation set.
residuals : list[numpy.ndarray[numpy.float64[m, 1]]] or numpy.ndarray[numpy.float64[n, m]]
    New residuals of the observations, in their order within the set.
)doc" )
                .def( "set_constant_single_observation_scalar_weight",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationScalarWeight,
                      py::arg( "condition" ),
                      py::arg( "weight" ),
                      R"doc(
Function to set the same weight for all components of the selected observations.

Correlations between the components of a selected observation are removed; correlations between different
observations are retained.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition
    Condition selecting the observations.
weight : float
    Weight of each observation component.
)doc" )
                .def( "set_constant_single_observation_diagonal_weight",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationDiagonalWeight,
                      py::arg( "condition" ),
                      py::arg( "weight" ),
                      R"doc(
Function to set the same diagonal weights for all selected observations.

Correlations between the components of a selected observation are removed; correlations between different
observations are retained.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition
    Condition selecting the observations.
weight : numpy.ndarray[numpy.float64[m, 1]]
    Weight of each observation component; its size must be equal to the number of components of the selected
    observations.
)doc" )
                .def( "set_constant_single_observation_matrix_weight",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantSingleObservationMatrixWeight,
                      py::arg( "condition" ),
                      py::arg( "weight" ),
                      R"doc(
Function to set the same weight matrix for all selected observations.

The weight matrix defines the weights of, and the correlations between, the components of each selected
observation. Correlations between different observations are retained.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition
    Condition selecting the observations.
weight : numpy.ndarray[numpy.float64[m, m]]
    Symmetric weight matrix, with one row and column per observation component.
)doc" )
                .def( "set_weight_vector_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightVectorForSet,
                      py::arg( "set_id" ),
                      py::arg( "weight_vector" ),
                      R"doc(
Function to set diagonal weights for all observations in an observation set.

Correlations between the observations in the set are removed; correlations with observations in other sets are
retained.

Parameters
----------
set_id : int
    Id of the observation set.
weight_vector : numpy.ndarray[numpy.float64[m, 1]]
    Weight of each component of each observation in the set, concatenated in the order of the observations within
    the set.
)doc" )
                .def( "set_weight_matrix_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightMatrixForSet,
                      py::arg( "set_id" ),
                      py::arg( "weight_matrix" ),
                      R"doc(
Function to set the weight matrix of an observation set, including the correlations between its observations.

Parameters
----------
set_id : int
    Id of the observation set.
weight_matrix : numpy.ndarray[numpy.float64[m, m]]
    Symmetric weight matrix, with one row and column per component of each observation in the set (in the order of
    the observations within the set).
)doc" )
                .def( "has_weight_matrix_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::hasWeightMatrixForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to check whether the weight matrix of an observation set contains correlations (off-diagonal entries).

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
bool
    Whether the weight matrix of the set contains off-diagonal entries.
)doc" )
                .def( "set_weight_matrix_for_observation",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightMatrixForObservation,
                      py::arg( "observation_id" ),
                      py::arg( "weight_matrix" ),
                      R"doc(
Function to set the weight matrix of a single observation, including the correlations between its components.

Parameters
----------
observation_id : int
    Id of the observation.
weight_matrix : numpy.ndarray[numpy.float64[m, m]]
    Symmetric weight matrix, with one row and column per observation component.
)doc" )
                .def( "has_weight_matrix_for_observation",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::hasWeightMatrixForObservation,
                      py::arg( "observation_id" ),
                      R"doc(
Function to check whether the weight matrix of a single observation contains correlations (off-diagonal entries).

Parameters
----------
observation_id : int
    Id of the observation.

Returns
-------
bool
    Whether the weight matrix of the observation contains off-diagonal entries.
)doc" )
                .def( "set_weight_block",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightBlock,
                      py::arg( "row_observation_ids" ),
                      py::arg( "column_observation_ids" ),
                      py::arg( "weight_block" ),
                      py::arg( "row_components" ) = std::vector< unsigned int >( ),
                      py::arg( "column_components" ) = std::vector< unsigned int >( ),
                      R"doc(
Function to set a block of the weight matrix, for instance to define correlations between observations in
different sets.

The block defines the weights between the components of the observations in ``row_observation_ids`` (rows) and
those of the observations in ``column_observation_ids`` (columns). The transposed block is set as well, so that the
weight matrix remains symmetric; where the rows and columns overlap, the block must be consistent with this (to a
relative accuracy of :math:`10^{-12}` of its largest entry). Each observation id may be given only once in each list.
Weights that correlate observations in different sets cannot be used with the Carpino outlier rejection or the
rejection of observations based on a sigma limit, which require the covariance of single observations.

Parameters
----------
row_observation_ids : list[int]
    Ids of the observations corresponding to the rows of the block.
column_observation_ids : list[int]
    Ids of the observations corresponding to the columns of the block.
weight_block : numpy.ndarray[numpy.float64[m, n]]
    Block of the weight matrix, with one row (column) per selected component of each row (column) observation.
row_components : list[int], default = []
    Indices of the components of each row observation that are included in the block. If empty, all components
    are included.
column_components : list[int], default = []
    Same as ``row_components``, for the column observations.
)doc" )
                .def_property_readonly( "has_extra_weight_blocks",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::hasExtraWeightBlocks,
                                        R"doc(
**read-only**

Whether the weight matrix of the dataset contains any off-diagonal entries (correlations).

:type: bool
)doc" )
                .def(
                        "add_observations_to_set",
                        []( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& self,
                            const unsigned int setId,
                            py::object observations,
                            py::object times,
                            py::object dependentVariables,
                            py::object weights,
                            py::object residuals,
                            const bool sortObservations ) {
                            const std::vector< TIME_TYPE > convertedTimes = castObservationTimeList< TIME_TYPE >( times, "times" );
                            const int observationCount = static_cast< int >( convertedTimes.size( ) );
                            self.addObservationsToSet(
                                    setId,
                                    castObservationVectorList< STATE_SCALAR_TYPE >( observations, "observations", observationCount ),
                                    convertedTimes,
                                    castObservationVectorList< double >( dependentVariables, "dependent_variables", observationCount ),
                                    castObservationVectorList< double >( weights, "weights", observationCount ),
                                    castObservationVectorList< STATE_SCALAR_TYPE >( residuals, "residuals", observationCount ),
                                    sortObservations );
                        },
                        py::arg( "set_id" ),
                        py::arg( "observations" ),
                        py::arg( "times" ),
                        py::arg( "dependent_variables" ) = py::list( ),
                        py::arg( "weights" ) = py::list( ),
                        py::arg( "residuals" ) = py::list( ),
                        py::arg( "sort_observations" ) = true,
                        R"doc(
Function to add observations to an existing observation set.

Parameters
----------
set_id : int
    Id of the observation set.
observations : list[numpy.ndarray[numpy.float64[m, 1]]] or numpy.ndarray[numpy.float64[n, m]]
    Values of the new observations (see :meth:`add_observation_set`).
times : list[tudatpy.astro.time_representation.Time] or numpy.ndarray[numpy.float64[n]]
    Times of the new observations, in seconds since J2000 (TDB).
dependent_variables : list[numpy.ndarray[numpy.float64[k, 1]]], default = []
    Dependent variables of the new observations. These must be provided if, and only if, dependent variables are
    stored for the observations already in the set.
weights : list[numpy.ndarray[numpy.float64[m, 1]]], default = []
    Weights of the new observations, with one vector per observation. If empty, all weights are set to 1.
residuals : list[numpy.ndarray[numpy.float64[m, 1]]], default = []
    Residuals of the new observations. If empty, all residuals are set to 0.
sort_observations : bool, default = True
    Whether to sort the observations in the set by time after adding the new observations.
)doc" )
                .def( "remove_observations",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::removeObservations,
                      py::arg( "condition" ),
                      R"doc(
Function to remove the selected observations from the dataset.

Unlike :meth:`reject_observations`, this permanently removes the observations, together with their weights and
numerical metadata. The ids of the remaining observations do not change.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition
    Condition selecting the observations to remove.
)doc" )
                .def( "delete_rejected_observations",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::deleteRejectedObservations,
                      R"doc(
Function to permanently remove all rejected observations from the dataset.

The ids and weights of the remaining observations do not change.
)doc" )
                .def_property_readonly( "observation_time_bounds",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBounds,
                                        R"doc(
**read-only**

Earliest and latest observation time in the dataset, in seconds since J2000 (TDB). Both are NaN if the dataset
contains no observations.

:type: tuple[tudatpy.astro.time_representation.Time, tudatpy.astro.time_representation.Time]
)doc" )
                .def_property_readonly( "number_of_observation_sets",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservationSets,
                                        R"doc(
**read-only**

Number of observation sets in the dataset.

:type: int
)doc" )
                .def_property_readonly( "number_of_observations",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservations,
                                        R"doc(
**read-only**

Number of observations in the dataset, including rejected observations.

:type: int
)doc" )
                .def_property_readonly( "total_scalar_size",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalScalarSize,
                                        R"doc(
**read-only**

Total number of observation components in the dataset, including those of rejected observations.

:type: int
)doc" )
                .def_property_readonly(
                        "observation_set_metadata",
                        py::overload_cast<>( &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSetMetadata,
                                             py::const_ ),
                        py::return_value_policy::copy,
                        R"doc(
**read-only**

Properties of each observation set in the dataset, ordered by set id.

:type: list[tudatpy.estimation.observations.ObservationSetMetadata]
)doc" )
                .def( "get_observation_set_metadata",
                      py::overload_cast< const unsigned int >(
                              &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSetMetadata, py::const_ ),
                      py::return_value_policy::copy,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the properties of an observation set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
tudatpy.estimation.observations.ObservationSetMetadata
    Properties of the observation set.
)doc" )
                .def_property_readonly( "observation_rows",
                                        &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationRows,
                                        py::return_value_policy::copy,
                                        R"doc(
**read-only**

Properties of each observation in the dataset, in the ``"internal"`` order (see :class:`ObservationDataset`).

:type: list[tudatpy.estimation.observations.ObservationDatasetRow]
)doc" )
                .def( "observation_row",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationRow,
                      py::return_value_policy::copy,
                      py::arg( "observation_id" ),
                      R"doc(
Function to retrieve the properties of a single observation.

Parameters
----------
observation_id : int
    Id of the observation.

Returns
-------
tudatpy.estimation.observations.ObservationDatasetRow
    Properties of the observation.
)doc" )
                .def( "scalar_component_row",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getScalarComponentRow,
                      py::arg( "scalar_component_id" ),
                      R"doc(
Function to retrieve the observation and component corresponding to an index in the storage of all observation
components in the dataset.

Parameters
----------
scalar_component_id : int
    Index of the component in the storage of the dataset (see
    :attr:`ObservationVectorData.scalar_component_ids` and :attr:`ObservationDatasetRow.first_scalar_component`).

Returns
-------
tudatpy.estimation.observations.ObservationScalarComponentRow
    Observation id and component index.
)doc" )
                .def( "observation_ids_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationIdsForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the ids of the observations in an observation set, in their order within the set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
list[int]
    Ids of the observations in the set.
)doc" )
                .def( "add_numerical_observation_metadata",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addNumericalObservationMetadata,
                      py::arg( "key" ),
                      py::arg( "observation_ids" ),
                      py::arg( "values" ),
                      R"doc(
Function to store a user-defined numerical value for a number of observations, under a given name.

These values (e.g. the scan angle of Gaia observations) are stored with the observations, but are not used by the
observation models or the estimation. If values are already stored under the same name, the values of the given
observations are replaced, and those of the other observations are retained.

Parameters
----------
key : str
    Name under which the values are stored.
observation_ids : list[int]
    Ids of the observations.
values : list[float]
    Value for each observation in ``observation_ids``.
)doc" )
                .def( "get_numerical_observation_metadata",
                      py::overload_cast< const std::string&, const std::vector< unsigned int >& >(
                              &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumericalObservationMetadata, py::const_ ),
                      py::arg( "key" ),
                      py::arg( "observation_ids" ),
                      R"doc(
Function to retrieve user-defined numerical values (see :meth:`add_numerical_observation_metadata`) for a number of
observations.

Parameters
----------
key : str
    Name under which the values are stored.
observation_ids : list[int]
    Ids of the observations, for instance retrieved with :meth:`get_observation_ids`.

Returns
-------
list[float]
    Value for each observation in ``observation_ids``. An ``IndexError`` is raised if no value is stored for one of
    the observations.
)doc" )
                .def_property_readonly(
                        "numerical_observation_metadata",
                        py::overload_cast<>( &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumericalObservationMetadata,
                                             py::const_ ),
                        py::return_value_policy::copy,
                        R"doc(
**read-only**

All user-defined numerical values stored in the dataset (see :meth:`add_numerical_observation_metadata`), with the
name of the values as key and, as value, a dictionary with the observation id as key and the numerical value as value.

:type: dict[str, dict[int, float]]
)doc" )
                .def( "observations_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the values of all observations in an observation set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
list[numpy.ndarray[numpy.float64[m, 1]]]
    Value of each observation, in their order within the set.
)doc" )
                .def( "observation_vector_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationVectorForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the values of all observations in an observation set, concatenated into a single vector.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
numpy.ndarray[numpy.float64[m, 1]]
    Concatenated values of the observations, in their order within the set.
)doc" )
                .def( "observation_value",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationValue,
                      py::arg( "observation_id" ),
                      R"doc(
Function to retrieve the value of a single observation.

Parameters
----------
observation_id : int
    Id of the observation.

Returns
-------
numpy.ndarray[numpy.float64[m, 1]]
    Value of the observation, with one entry per observation component.
)doc" )
                .def( "observation_times_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimesForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the times of all observations in an observation set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
list[tudatpy.astro.time_representation.Time]
    Time of each observation (in seconds since J2000, TDB), in their order within the set.
)doc" )
                .def( "observation_time",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTime,
                      py::arg( "observation_id" ),
                      R"doc(
Function to retrieve the time of a single observation.

Parameters
----------
observation_id : int
    Id of the observation.

Returns
-------
tudatpy.astro.time_representation.Time
    Time of the observation, in seconds since J2000 (TDB).
)doc" )
                .def( "weights_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the diagonal weights of all observations in an observation set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
list[numpy.ndarray[numpy.float64[m, 1]]]
    Diagonal of the weight matrix of each observation, in their order within the set.
)doc" )
                .def( "weight_matrix_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightMatrixForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the weight matrix of an observation set, including the correlations between its observations.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
numpy.ndarray[numpy.float64[m, m]]
    Weight matrix, with one row and column per component of each observation in the set.
)doc" )
                .def( "weight_value",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightValue,
                      py::arg( "observation_id" ),
                      R"doc(
Function to retrieve the diagonal weights of a single observation.

Parameters
----------
observation_id : int
    Id of the observation.

Returns
-------
numpy.ndarray[numpy.float64[m, 1]]
    Diagonal of the weight matrix of the observation.
)doc" )
                .def( "weight_matrix_for_observation",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightMatrixForObservation,
                      py::arg( "observation_id" ),
                      R"doc(
Function to retrieve the weight matrix of a single observation, including the correlations between its components.

Parameters
----------
observation_id : int
    Id of the observation.

Returns
-------
numpy.ndarray[numpy.float64[m, m]]
    Weight matrix, with one row and column per observation component.
)doc" )
                .def( "residuals_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualsForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the residuals of all observations in an observation set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
list[numpy.ndarray[numpy.float64[m, 1]]]
    Residual of each observation, in their order within the set.
)doc" )
                .def( "residual_vector_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualVectorForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the residuals of all observations in an observation set, concatenated into a single vector.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
numpy.ndarray[numpy.float64[m, 1]]
    Concatenated residuals of the observations, in their order within the set.
)doc" )
                .def( "residual_value",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualValue,
                      py::arg( "observation_id" ),
                      R"doc(
Function to retrieve the residual of a single observation.

Parameters
----------
observation_id : int
    Id of the observation.

Returns
-------
numpy.ndarray[numpy.float64[m, 1]]
    Residual of the observation, with one entry per observation component.
)doc" )
                .def( "dependent_variables_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariablesForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the dependent variables of all observations in an observation set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
list[numpy.ndarray[numpy.float64[k, 1]]]
    Dependent variables of each observation, in their order within the set (empty if no dependent variables are
    stored).
)doc" )
                .def( "dependent_variables",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariables,
                      py::arg( "observation_id" ),
                      R"doc(
Function to retrieve the dependent variables of a single observation.

Parameters
----------
observation_id : int
    Id of the observation.

Returns
-------
numpy.ndarray[numpy.float64[k, 1]]
    Dependent variables of the observation (empty if none are stored).
)doc" )
                .def( "set_dependent_variables_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setDependentVariablesForSet,
                      py::arg( "set_id" ),
                      py::arg( "dependent_variables" ),
                      R"doc(
Function to replace the dependent variables of all observations in an observation set.

Parameters
----------
set_id : int
    Id of the observation set.
dependent_variables : list[numpy.ndarray[numpy.float64[k, 1]]]
    Dependent variables of each observation, in their order within the set. An empty list removes the dependent
    variables of the set.
)doc" )
                .def( "add_dependent_variable",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addDependentVariableToSets,
                      py::arg( "dependent_variable_settings" ),
                      py::arg_v( "condition", tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::all( ), "observation_query.all" ),
                      R"doc(
Function to add a dependent variable to the observation sets that contain the selected observations.

The dependent variable is added to each of these sets for which it is compatible with the observable type and link
ends. Its values are not computed by this function, but by
:func:`~tudatpy.estimation.observations.compute_residuals_and_dependent_variables`. A dependent variable cannot be
added to a set for which dependent variable values are already stored. The estimation does not compute dependent
variables: values stored before an estimation are kept, and dependent variables added afterwards are computed by
:func:`~tudatpy.estimation.observations.compute_residuals_and_dependent_variables`, which also replaces the residuals.

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    Settings of the dependent variable to add.
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations; the dependent variable is added to each set that contains at least one
    selected observation.
)doc" )
                .def( "number_of_observations_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservationsForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the number of observations in an observation set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
int
    Number of observations in the set, including rejected observations.
)doc" )
                .def( "total_scalar_size_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalScalarSizeForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the total number of observation components in an observation set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
int
    Number of observations in the set, multiplied by the number of components per observation.
)doc" )
                .def( "link_definition",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinition,
                      py::arg( "link_definition_id" ),
                      R"doc(
Function to retrieve a link definition from its id (see :attr:`ObservationSetMetadata.link_definition_id`).

Parameters
----------
link_definition_id : int
    Id of the link definition.

Returns
-------
tudatpy.estimation.observable_models_setup.links.LinkDefinition
    Link definition.
)doc" )
                .def( "set_link_end_reference_point",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setLinkEndReferencePoint,
                      py::arg( "body_name" ),
                      py::arg( "reference_point_name" ),
                      py::arg( "link_end_type" ),
                      py::arg_v( "condition", tom::ObservationSelectionCondition< STATE_SCALAR_TYPE, TIME_TYPE >::all( ), "observation_query.all" ),
                      R"doc(
Function to change the reference point of a link end, in the observation sets that contain the selected
observations.

In each of these sets, the reference point of the link end of type ``link_end_type`` is set to
``reference_point_name``, if this link end is on the body ``body_name``. Only the link definition stored in the
dataset is changed: the reference point must be created separately in the system of bodies (e.g. with
:func:`~tudatpy.dynamics.environment_setup.add_ground_station` for a ground station, or with
:meth:`~tudatpy.dynamics.environment.VehicleSystems.set_reference_point` for a spacecraft antenna).
The dependent variables of these sets are defined for the new link ends. Their stored values and the residuals are
not changed: they are recomputed for the new link ends by
:func:`~tudatpy.estimation.observations.compute_residuals_and_dependent_variables` (or, for the residuals, by an
estimation).

Parameters
----------
body_name : str
    Name of the body of the link end.
reference_point_name : str
    Name of the new reference point.
link_end_type : tudatpy.estimation.observable_models_setup.links.LinkEndType
    Role of the link end in the observations (e.g. receiver or transmitter).
condition : tudatpy.estimation.observations.ObservationSelectionCondition, default = observation_query.all
    Condition selecting the observations; the link end is changed in each set that contains selected observations. An
    error is raised if the condition selects only part of a set in which the link end would be changed, since the
    reference point can only be changed for complete sets.
)doc" )
                .def(
                        "link_definitions_for_observable",
                        []( const InspectionDataset& dataset, const tom::ObservableType observableType ) {
                            std::map< unsigned int, tom::LinkDefinition > definitions;
                            for( unsigned int setId = 0; setId < dataset.getNumberOfObservationSets( ); ++setId )
                            {
                                const auto& metadata = dataset.getObservationSetMetadata( setId );
                                if( metadata.observableType_ == observableType && dataset.getNumberOfObservationsForSet( setId ) > 0 )
                                {
                                    definitions.emplace( metadata.linkDefinitionId_,
                                                         dataset.getLinkDefinition( metadata.linkDefinitionId_ ) );
                                }
                            }
                            std::vector< tom::LinkDefinition > result;
                            result.reserve( definitions.size( ) );
                            for( const auto& definition : definitions )
                            {
                                result.push_back( definition.second );
                            }
                            return result;
                        },
                        py::arg( "observable_type" ),
                        R"doc(
Function to retrieve the link definitions of the observation sets of a given observable type.

Parameters
----------
observable_type : tudatpy.estimation.observable_models_setup.model_settings.ObservableType
    Observable type of the observation sets.

Returns
-------
list[tudatpy.estimation.observable_models_setup.links.LinkDefinition]
    Link definitions of the (non-empty) observation sets of the given observable type, each given once.
)doc" )
                .def( "ancillary_settings",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettings,
                      py::arg( "ancillary_settings_id" ),
                      R"doc(
Function to retrieve ancillary settings from their id (see :attr:`ObservationSetMetadata.ancillary_settings_id`).

Parameters
----------
ancillary_settings_id : int
    Id of the ancillary settings.

Returns
-------
tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings
    Ancillary settings (``None`` if the id corresponds to sets without ancillary settings).
)doc" )
                .def( "ancillary_settings_for_set",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsForSet,
                      py::arg( "set_id" ),
                      R"doc(
Function to retrieve the ancillary settings of an observation set.

Parameters
----------
set_id : int
    Id of the observation set.

Returns
-------
tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings
    Ancillary settings of the set (``None`` if the set has none).
)doc" )
                .def( "dependent_variable_bookkeeping",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableBookkeeping,
                      py::arg( "dependent_variable_layout_id" ),
                      R"doc(
Function to retrieve the object defining which dependent variables are stored for an observation set, from its id
(see :attr:`ObservationSetMetadata.dependent_variable_layout_id`).

Parameters
----------
dependent_variable_layout_id : int
    Id of the object.

Returns
-------
tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableBookkeeping
    Object defining the dependent variables (``None`` if the id corresponds to sets without dependent variables).
)doc" )
                .def(
                        "__copy__",
                        []( const tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& dataset ) {
                            return std::make_shared< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >( dataset );
                        },
                        R"doc(Returns an independent copy of the dataset (used through ``copy.copy``).)doc" )
                .def(
                        "__deepcopy__",
                        []( const tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& dataset, const py::dict& ) {
                            return std::make_shared< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >( dataset );
                        },
                        py::arg( "memo" ),
                        R"doc(Returns an independent copy of the dataset (used through ``copy.deepcopy``).)doc" )
                .def( "create_new_and_keep",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createNewAndKeep,
                      py::arg( "condition" ),
                      R"doc(
Function to create a new dataset that contains only the selected observations.

The selected observations keep their ids, set ids, weights and rejection status. This dataset is not modified.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition
    Condition selecting the observations to keep.

Returns
-------
tudatpy.estimation.observations.ObservationDataset
    New dataset containing the selected observations.
)doc" )
                .def( "create_new_and_drop",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createNewAndDrop,
                      py::arg( "condition" ),
                      R"doc(
Function to create a new dataset that contains all observations except the selected ones.

The remaining observations keep their ids, set ids, weights and rejection status. This dataset is not modified.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition
    Condition selecting the observations to leave out.

Returns
-------
tudatpy.estimation.observations.ObservationDataset
    New dataset containing the observations that are not selected.
)doc" )
                .def( "reject_observations",
                      py::overload_cast< const InspectionCondition&, const std::string& >(
                              &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::rejectObservations ),
                      py::arg( "condition" ),
                      py::arg( "reason" ) = "",
                      R"doc(
Function to reject the selected observations.

Rejected observations remain in the dataset, but are not used in the estimation. They can be used again after
calling :meth:`restore_observations`.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition
    Condition selecting the observations to reject.
reason : str, default = ""
    Reason for the rejection, stored with each rejected observation (see
    :attr:`ObservationDatasetRow.rejection_reason`). Observations that are already rejected keep their reason.
)doc" )
                .def( "reject_observations",
                      py::overload_cast< const InspectionCondition&, const double, const std::string& >(
                              &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::rejectObservations ),
                      py::arg( "condition" ),
                      py::arg( "sigma_limit" ),
                      py::arg( "reason" ) = "",
                      R"doc(
Function to reject the selected observations for which a residual is larger than a multiple of the observation
uncertainty.

A selected observation is rejected if the absolute value of any component of its residual is larger than
``sigma_limit`` times the standard deviation of that component. The standard deviations are computed from the
weights, as described for :meth:`ObservationVectorData.inverse_weight_matrix_for_observation`. The residuals must
be computed beforehand, for instance with :func:`~tudatpy.estimation.observations.compute_residuals_and_dependent_variables`
or by an estimation: residuals that have not been computed are zero, so that no observations are rejected.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition
    Condition selecting the observations to check. Observations that are already rejected remain rejected.
sigma_limit : float
    Maximum ratio of the absolute value of a residual component and its standard deviation (must be positive).
reason : str, default = ""
    Reason for the rejection, stored with each rejected observation.
)doc" )
                .def( "restore_observations",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::restoreObservations,
                      py::arg( "condition" ),
                      R"doc(
Function to restore the selected (rejected) observations, so that they are used in the estimation again.

The rejection reason of the restored observations is cleared.

Parameters
----------
condition : tudatpy.estimation.observations.ObservationSelectionCondition
    Condition selecting the observations to restore.
)doc" )
                .def( "observation_vector_data",
                      &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createObservationVectorData,
                      py::arg( "include_rejected" ) = false,
                      R"doc(
Function to create an object holding the concatenated observations, residuals and weights of the dataset, as used
in the estimation.

Parameters
----------
include_rejected : bool, default = False
    Whether to include the rejected observations.

Returns
-------
tudatpy.estimation.observations.ObservationVectorData
    Concatenated data, in the ``"estimation"`` order (see :class:`ObservationDataset`).
)doc" ) TUDATPY_DEF_PICKLE( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > )
                        TUDATPY_DEF_BINARY_IO( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > );
    }

    m.def( "create_observation_dataset_from_tracking_data",
           &createObservationDatasetFromTrackingData,
           py::arg( "tracking_data" ),
           py::arg( "bodies" ),
           py::arg( "apply_corrections" ) = false,
           R"doc(
Function to create an observation dataset from a list of tracking data objects.

Each entry of ``tracking_data`` (as produced e.g. by :func:`~tudatpy.data_input.tracking_data.ifms.read_ifms_data`,
:func:`~tudatpy.data_input.tracking_data.fdets.read_fdets_data`, or a custom reader) is converted into an observation
set, with the observations sorted by time. The observation times are converted to TDB, using the position of the
ground station in ``bodies`` (if available) for the time scale conversion. The weights and numerical metadata stored
in the tracking data are transferred to the dataset.

Parameters
----------
tracking_data : list[:class:`~tudatpy.data_input.tracking_data.TrackingData`]
    List of tracking data objects to convert.
bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
    System of bodies, used to retrieve the ground station positions for the time scale conversions, and to compute the
    angular corrections (see ``apply_corrections``).
apply_corrections : bool, default = False
    Whether to apply the corrections stored in the tracking data to the observation values. If ``True``, the
    corrections requested with
    :meth:`~tudatpy.data_input.tracking_data.TrackingData.set_observation_correction_settings` (light deflection and
    photocenter offset) are also computed and applied, using the environment in ``bodies``. The tracking data objects
    are not modified.

Returns
-------
tudatpy.estimation.observations.ObservationDataset
    Dataset containing one observation set per tracking data object.
)doc" );

    // SINGLE OBSERVATION SET

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
            .def( py::init< const tom::ObservableType,
                            const tom::LinkDefinition,
                            const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >,
                            const std::vector< TIME_TYPE >,
                            const tom::LinkEndType,
                            const std::vector< Eigen::VectorXd >,
                            const std::shared_ptr< tss::ObservationDependentVariableBookkeeping >,
                            const std::shared_ptr< tom::ObservationAncillarySimulationSettings > >( ),
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "observations" ),
                  py::arg( "observation_epochs" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "observation_dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg_v( "dependent_variable_bookkeeping", std::shared_ptr< tss::ObservationDependentVariableBookkeeping >( ), "None" ),
                  py::arg_v( "ancillary_settings", std::shared_ptr< tom::ObservationAncillarySimulationSettings >( ), "None" ),
                  R"doc(

         Create an observation set from its observable type, link, measurement vectors, epochs and reference link end,
         with optional dependent variables and ancillary data.

         Parameters
         ----------
         observable_type : ObservableType
             Type of observable to which the settings or measurements apply.
         link_ends : LinkDefinition
             Definition of the bodies and reference points participating in the observation link.
         observations : list[numpy.ndarray[numpy.float64[m, 1]]]
             Measurement vectors for the observation set.
         observation_epochs : list[Time]
             Epochs of the measurements, in seconds since J2000.
         reference_link_end : LinkEndType
             Link end at which the observation epochs are defined.
         observation_dependent_variables : list[numpy.ndarray[numpy.float64[m, 1]]], default = []
             Dependent-variable values associated with each measurement.
         dependent_variable_bookkeeping : ObservationDependentVariableBookkeeping, default = None
             Optional object defining the dependent-variable layout.
         ancillary_settings : ObservationAncillarySimulationSettings, default = None
             Optional observable-specific ancillary data.

      )doc" )
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
            .def( "filter_observations",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations,
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
dict[astro.time_representation.Time, numpy.ndarray]
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
)doc" ) TUDATPY_DEF_PICKLE( tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > )
                    TUDATPY_DEF_EQ_NE( tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > )
                            TUDATPY_DEF_BINARY_IO( tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > );

    m.def( "create_single_observation_set",
           py::overload_cast< const tom::ObservableType,
                              const tom::LinkEnds&,
                              const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >&,
                              const std::vector< TIME_TYPE >,
                              const tom::LinkEndType,
                              const std::shared_ptr< tom::ObservationAncillarySimulationSettings > >(
                   &tom::createSingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "observations" ),
           py::arg( "observation_times" ),
           py::arg( "reference_link_end" ),
           py::arg_v( "ancillary_settings", std::shared_ptr< tom::ObservationAncillarySimulationSettings >( ), "None" ),
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
        observation_times : list[astro.time_representation.Time]
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

    // OBSERVATION COLLECTION

    py::class_< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                 "ObservationCollection",
                                                                                                 R"doc(

         Class collecting all observations and associated data for use in an estimation.

         Class containing the full set of observations and associated data, typically for input into the estimation. When using simulated data,
         this class is instantiated via a call to the :func:`~tudatpy.estimation.observations.simulate_observations` function. More information is provided
         on the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__





      )doc" )
            .def( py::init< std::vector< std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > > >( ),
                  py::arg( "observation_sets" ),
                  R"doc(
Constructor for the ObservationCollection class.

Parameters
----------
observation_sets : list[tudatpy.estimation.observations.SingleObservationSet]
    List of single observation sets to be included in the collection.
)doc" )
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the observable types for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[ObservableType]
             A list of observable types present in the selected subset.
     )doc" )
            .def( "get_bodies_in_link_ends",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getBodiesInLinkEnds,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the names of bodies present in the link ends of a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[str]
             A list of body names.
     )doc" )
            .def( "get_reference_points_in_link_ends",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getReferencePointsInLinkEnds,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the names of reference points (e.g., ground stations) in the link ends of a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[str]
             A list of reference point names.
     )doc" )
            .def( "get_time_bounds_list",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsListDouble,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the time bounds for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[float, float]]
             A list of time bound pairs (start_time, end_time).
     )doc" )
            .def( "get_time_bounds_list_time_object",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsList,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the time bounds for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[Time, Time]]
             A list of time bound pairs (start_time, end_time) as Time objects.
     )doc" )
            .def( "get_time_bounds_per_set",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsPerSetDouble,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the time bounds for each set in a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[float, float]]
             A list of time bound pairs for each observation set.
     )doc" )
            .def( "get_time_bounds_per_set_time_object",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsPerSet,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the time bounds for each set in a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[Time, Time]]
             A list of time bound pairs for each observation set as Time objects.
     )doc" )
            .def( "get_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of observation vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservations,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the concatenated observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of observation values.
     )doc" )
            .def( "get_observation_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimesDouble,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the observation times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[list[float]]
             A list of lists of observation times, one for each selected observation set.
     )doc" )
            .def( "get_observation_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimes,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the observation times for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[list[Time]]
             A list of lists of observation times as Time objects.
     )doc" )
            .def( "get_concatenated_observation_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationTimes,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the concatenated observation times for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[Time]
             A concatenated list of observation times as Time objects.
     )doc" )
            .def( "get_concatenated_observation_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDoubleObservationTimes,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the concatenated observation times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[float]
             A concatenated list of observation times.
     )doc" )
            .def( "get_observations_and_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsAndTimesDouble,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the observations and times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], list[list[float]]]
             A pair containing a list of observation vectors and a list of lists of observation times.
     )doc" )
            .def( "get_observations_and_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsAndTimes,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the observations and times for a subset of observation sets, with times as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], list[list[Time]]]
             A pair containing a list of observation vectors and a list of lists of observation times as Time objects.
     )doc" )
            .def( "get_concatenated_observations_and_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationsAndTimesDouble,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the concatenated observations and times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[numpy.ndarray, list[float]]
             A pair containing the concatenated observation vector and list of times.
     )doc" )
            .def( "get_concatenated_observations_and_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationsAndTimes,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the concatenated observations and times for a subset of observation sets, with times as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the weights for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of weight vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_weights",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedWeights,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the concatenated weights for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of weights.
     )doc" )
            .def( "get_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getResiduals,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of residual vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedResiduals,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the concatenated residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of residuals.
     )doc" )
            .def( "get_rms_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResiduals,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the RMS of residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of RMS residual vectors, one for each selected observation set.
     )doc" )
            .def( "get_mean_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResiduals,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the mean of residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of mean residual vectors, one for each selected observation set.
     )doc" )
            .def( "get_computed_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservations,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the computed observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of computed observation vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_computed_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedComputedObservations,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get the concatenated computed observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Set a constant weight for a subset of observation sets.

         Parameters
         ----------
         weight : float
             The constant weight to set.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const Eigen::VectorXd, const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Set a constant weight vector for a subset of observation sets.

         Parameters
         ----------
         weight : numpy.ndarray
             The constant weight vector to set.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Set tabulated weights for a subset of observation sets.

         Parameters
         ----------
         tabulated_weights : numpy.ndarray
             The vector of tabulated weights to set.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg( "observation_collection_to_append" ),
                  R"doc(

         Append all observation sets from ``observation_collection_to_append`` to this collection and refresh its
         observation bookkeeping.

         Parameters
         ----------
         observation_collection_to_append : ObservationCollection
             Collection whose observation sets are appended to this collection.

      )doc" )
            .def( "filter_observations",
                  py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                     std::shared_ptr< tom::ObservationFilterBase > >&,
                                     const bool >( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations ),
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
            .def( "filter_observations",
                  py::overload_cast< std::shared_ptr< tom::ObservationFilterBase >,
                                     std::shared_ptr< tom::ObservationCollectionParser >,
                                     const bool >( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  py::arg( "save_filtered_observations" ) = true,
                  R"doc(
         Filter observations using a single filter.

         This function filters a subset of observations (or all) using a single observation filter.

         Parameters
         ----------
         observation_filter : tudatpy.estimation.observations.observations_processing.ObservationFilterBase
             The observation filter to apply.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
         save_filtered_observations : bool, optional
             If true, the filtered-out observations are saved within each observation set, by default True.
     )doc" )
            .def( "split_observation_sets",
                  py::overload_cast< std::shared_ptr< tom::ObservationSetSplitterBase >,
                                     std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::splitObservationSets ),
                  py::arg( "observation_set_splitter" ),
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Split observation sets based on a splitter.

         This function splits a subset of observation sets (or all) into smaller sets based on the criteria defined by the splitter.

         Parameters
         ----------
         observation_set_splitter : tudatpy.estimation.observations.observations_processing.ObservationSetSplitterBase
             The splitter to use for splitting the observation sets.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
            .def( "get_single_observation_sets",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleObservationSets,
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get a subset of single observation sets.

         This function returns a list of the single observation sets that are selected by the provided observation parser.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
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
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
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
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
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
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Parser to select the observation sets.
     )doc" )
            .def( "set_transponder_delay",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTransponderDelay,
                  py::arg( "spacecraft_name" ),
                  py::arg( "transponder_delay" ),
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Deprecated: set the transponder delay for a subset of observations by modifying the
         retransmission delay in their ancillary settings.

         For new simulations, set the default transponder delay on the spacecraft vehicle systems
         before creating the observation model:
         ``bodies.get(spacecraft_name).system_models.transponder_delay = transponder_delay``.

         Parameters
         ----------
         spacecraft_name : str
             Name of the spacecraft with the transponder.
         transponder_delay : float
             The transponder delay in seconds.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Parser to select the observation sets.
     )doc" )
            .def( "remove_empty_observation_sets",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeEmptySingleObservationSets,
                  R"doc(Remove all single observation sets that contain no observations.)doc" )
            .def( "add_dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::addDependentVariable,
                  py::arg( "dependent_variable_settings" ),
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(

         Add an observation dependent variable to a subset of the single observation sets.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to add.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Retrieve the values of a given dependent variable, sorted per single observation set.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Retrieve the concatenated values of a given dependent variable.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[numpy.ndarray, tudatpy.estimation.observations.observations_processing.ObservationCollectionParser]
             A pair containing a matrix with the concatenated dependent variable values and the parser used.
     )doc" )
            .def( "compatible_dependent_variable_settings",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                  py::arg( "dependent_variable_settings" ),
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get a list of all dependent variable settings compatible with the input settings.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[list[list[tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings]], tudatpy.estimation.observations.observations_processing.ObservationCollectionParser]
             A pair containing a list of lists of compatible settings and the parser used.
     )doc" )
            .def( "compatible_dependent_variables_list",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Get all dependent variables compatible with the input settings.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Retrieve the time history of a given dependent variable, sorted per observation set.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Retrieve the time history of a given dependent variable, sorted per observation set, with times as Time objects.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Retrieve the concatenated time history of a given dependent variable.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
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
                  py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
                  R"doc(
         Retrieve the concatenated time history of a given dependent variable, with times as Time objects.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, default = observations_processing.observation_parser()
             Parser to select a subset of observation sets.

         Returns
         -------
         dict[Time, numpy.ndarray]
             A map from time to dependent variable value, with times as Time objects.
     )doc" ) TUDATPY_DEF_PICKLE( tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > )
                    TUDATPY_DEF_EQ_NE( tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > )
                            TUDATPY_DEF_BINARY_IO( tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > );

    tudatpy::estimation::observations::expose_observations_io_bindings( m );
    tudatpy::estimation::observations::expose_observations_simulation_bindings( m );

    m.def(
            "create_observation_dataset_from_single_observation_set",
            []( const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >& observationSet ) {
                return tom::createObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >( observationSet );
            },
            py::arg( "observation_set" ),
            R"doc(
Function to create an observation dataset from a (deprecated) single observation set.

The :class:`~tudatpy.estimation.observations.SingleObservationSet` and
:class:`~tudatpy.estimation.observations.ObservationCollection` classes are deprecated, and replaced by
:class:`~tudatpy.estimation.observations.ObservationDataset`. This function copies the observations, times, weights, residuals, dependent variables and
ancillary settings of the observation set into a new dataset. Observations that were filtered out of the set are not
included.

Parameters
----------
observation_set : tudatpy.estimation.observations.SingleObservationSet
    Observation set to convert.

Returns
-------
tudatpy.estimation.observations.ObservationDataset
    Dataset containing a single observation set.
)doc" );

    m.def(
            "create_observation_dataset_from_collection",
            []( const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >& observationCollection ) {
                return tom::createObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >( observationCollection );
            },
            py::arg( "observation_collection" ),
            R"doc(
Function to create an observation dataset from a (deprecated) observation collection.

The :class:`~tudatpy.estimation.observations.SingleObservationSet` and
:class:`~tudatpy.estimation.observations.ObservationCollection` classes are deprecated, and replaced by
:class:`~tudatpy.estimation.observations.ObservationDataset`. This function copies the observations, times, weights, residuals, dependent variables and
ancillary settings of each single observation set in the collection (in the order of
:meth:`~tudatpy.estimation.observations.ObservationCollection.get_single_observation_sets`) into an observation set
of a new dataset. Observations that were filtered out of the collection are not included.

Parameters
----------
observation_collection : tudatpy.estimation.observations.ObservationCollection
    Observation collection to convert.

Returns
-------
tudatpy.estimation.observations.ObservationDataset
    Dataset containing one observation set per single observation set of the collection.
)doc" );

    m.def(
            "create_single_observation_set_from_dataset",
            []( const std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >& observationDataset,
                const unsigned int setId ) {
                return tom::createSingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >( observationDataset, setId );
            },
            py::arg( "observation_dataset" ),
            py::arg( "set_id" ) = 0,
            R"doc(
Function to create a (deprecated) single observation set from an observation set in an observation dataset.

The :class:`~tudatpy.estimation.observations.SingleObservationSet` and
:class:`~tudatpy.estimation.observations.ObservationCollection` classes are deprecated, and replaced by
:class:`~tudatpy.estimation.observations.ObservationDataset`. Only the diagonal of the weight matrix is transferred, so that any correlations are lost.
The conversion is not possible if the observation set contains rejected observations; these can be removed with
:meth:`~tudatpy.estimation.observations.ObservationDataset.delete_rejected_observations`.

Parameters
----------
observation_dataset : tudatpy.estimation.observations.ObservationDataset
    Dataset containing the observation set.
set_id : int, default = 0
    Id of the observation set to convert.

Returns
-------
tudatpy.estimation.observations.SingleObservationSet
    Single observation set containing a copy of the observations.
)doc" );

    m.def(
            "create_observation_collection_from_dataset",
            []( const std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >& observationDataset ) {
                return tom::createObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >( observationDataset );
            },
            py::arg( "observation_dataset" ),
            R"doc(
Function to create a (deprecated) observation collection from an observation dataset.

The :class:`~tudatpy.estimation.observations.SingleObservationSet` and
:class:`~tudatpy.estimation.observations.ObservationCollection` classes are deprecated, and replaced by
:class:`~tudatpy.estimation.observations.ObservationDataset`. Only the diagonal of the weight matrix is transferred, so that any correlations are lost.
The conversion is not possible if the dataset contains rejected observations; these can be removed with
:meth:`~tudatpy.estimation.observations.ObservationDataset.delete_rejected_observations`.

Parameters
----------
observation_dataset : tudatpy.estimation.observations.ObservationDataset
    Dataset to convert.

Returns
-------
tudatpy.estimation.observations.ObservationCollection
    Observation collection containing one single observation set per observation set of the dataset.
)doc" );

    m.def(
            "create_observation_collection_from_tracking_data",
            []( const std::vector< std::shared_ptr< tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE > > >& trackingData,
                tss::SystemOfBodies& bodies,
                const bool applyCorrections ) {
                return tom::createObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >(
                        createObservationDatasetFromTrackingData( trackingData, bodies, applyCorrections ) );
            },
            py::arg( "tracking_data" ),
            py::arg( "bodies" ),
            py::arg( "apply_corrections" ) = false,
            R"doc(

        Factory function to create an `ObservationCollection` from a list of `TrackingData` objects.

        This function converts each entry of `tracking_data` (as produced e.g. by
        :func:`~tudatpy.data_input.tracking_data.ifms.read_ifms_data`, :func:`~tudatpy.data_input.tracking_data.fdets.read_fdets_data`, or a custom
        reader) into a `SingleObservationSet`, and collects the results into an `ObservationCollection`.

        Parameters
        ----------
        tracking_data : list[:class:`~tudatpy.data_input.tracking_data.TrackingData`]
            List of tracking data objects to convert.
        bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
            System of bodies, used to resolve ground station positions for time scale conversions.
        apply_corrections : bool, optional
            Whether corrections stored in the tracking data should be applied to
            the observation values during conversion. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            An `ObservationCollection` containing one observation set per input `TrackingData` object.
        )doc" );

    m.def( "set_tracking_supplementary_data_in_bodies",
           py::overload_cast< tss::SystemOfBodies&, const std::vector< std::shared_ptr< tdat::TrackingSupplementaryData > >& >(
                   &tom::setTrackingSupplementaryDataInBodies ),
           py::arg( "bodies" ),
           py::arg( "supplementary_data" ),
           R"doc(

        Sets tracking supplementary data (e.g. ground station frequency ramps) in the bodies.

        Parameters
        ----------
        bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
            System of bodies in which the supplementary data is set.
        supplementary_data : list[:class:`~tudatpy.data_input.tracking_data.TrackingSupplementaryData`]
            List of tracking supplementary data objects, as produced e.g. by
            :func:`~tudatpy.data_input.tracking_data.ifms.read_ifms_data` or a custom reader.
        )doc" );

    m.def(
            "compute_residuals_and_dependent_variables",
            []( const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > observationCollection,
                const std::vector< std::shared_ptr< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                        observationSimulators,
                const tss::SystemOfBodies& bodies ) {
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

    m.def( "observation_simulation_settings_from_dataset",
           &tss::getObservationSimulationSettingsFromObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_dataset" ),
           py::arg( "bodies" ),
           R"doc(
Function to create observation simulation settings matching the observable types, link ends and times of the
observations in an observation dataset.

One set of (tabulated) observation simulation settings is created per observation set in the dataset, with the
same observation times, reference link end, ancillary settings and dependent variables. These settings can be used
to simulate observations at the times of (real) observations, for instance with
:func:`~tudatpy.estimation.observations.simulate_observation_dataset`. The weights are not part of these settings: the
observations simulated with them have unit weights.

Parameters
----------
observation_dataset : tudatpy.estimation.observations.ObservationDataset
    Dataset defining the observable types, link ends and observation times.
bodies : tudatpy.dynamics.environment.SystemOfBodies
    System of bodies defining the physical environment.

Returns
-------
list[ObservationSimulationSettings]
    Observation simulation settings, one per observation set in the dataset.
)doc" );

    m.def( "compute_residuals_and_dependent_variables",
           &computeDatasetResiduals,
           py::arg( "observation_dataset" ),
           py::arg( "observation_simulators" ),
           py::arg( "bodies" ),
           R"doc(
Computes residuals and dependent variables for a given observation dataset.

This function simulates observations based on the settings of the input ``observation_dataset`` (which typically
contains real data), see :func:`~tudatpy.estimation.observations.observation_simulation_settings_from_dataset`. It
then computes the residuals by subtracting the simulated observations from the original observations. The computed
residuals, and the dependent variables defined for each observation set (see
:meth:`~tudatpy.estimation.observations.ObservationDataset.add_dependent_variable`), are then stored in the input
``observation_dataset``. Residuals are computed for all observations, including rejected observations.

The observations are simulated with the current environment. After an estimation, the ephemerides of the propagated
bodies correspond to the last propagation of the estimation (also if the parameters were reset to those of the best
iteration), so the residuals stored in the dataset by the estimation (those of the best iteration) can be replaced
by different values.

Parameters
----------
observation_dataset : tudatpy.estimation.observations.ObservationDataset
    The dataset of observations for which to compute residuals and dependent variables. This object is modified
    in-place.
observation_simulators : list[tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator]
    List of observation simulators to be used for computing the observations.
bodies : tudatpy.dynamics.environment.SystemOfBodies
    The system of bodies required for the observation simulation.
)doc" );

    m.def( "create_filtered_observation_set",
           py::overload_cast< const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationFilterBase >,
                              const bool >( &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
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
            "filter_observations",
            []( const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > originalObservationSet,
                const std::shared_ptr< tom::ObservationFilterBase > observationFilter,
                const bool saveFilteredObservations ) {
                tudat::utilities::printDeprecationWarning( "tudatpy.estimation.observations.filter_observations",
                                                           "tudatpy.estimation.observations.create_filtered_observation_set",
                                                           "Deprecated as of v1.1." );
                return tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE >(
                        originalObservationSet, observationFilter, saveFilteredObservations );
            },
            py::arg( "original_observation_set" ),
            py::arg( "observation_filter" ),
            py::arg( "save_filtered_observations" ) = false,
            R"doc(Deprecated as of v1.1. Use :func:`~create_filtered_observation_set` instead.)doc" );

    m.def( "split_observation_set",
           py::overload_cast< const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationSetSplitterBase >,
                              const bool >( &tom::splitObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > ),
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

    m.def( "merge_observation_collections",
           &tss::mergeObservationCollections< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection_list" ),
           R"doc(

         Return an ObservationCollection containing the observation sets from every collection in
         ``observation_collection_list``.

         Parameters
         ----------
         observation_collection_list : list[ObservationCollection]
             Observation collections to merge, in the supplied order.

         Returns
         -------
         ObservationCollection
             Collection containing the merged observation sets.

      )doc" );

    // The following functions create a new ObservationCollection object from an existing one

    m.def( "create_filtered_observation_collection",
           py::overload_cast<
                   const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                   const std::map< std::shared_ptr< tom::ObservationCollectionParser >, std::shared_ptr< tom::ObservationFilterBase > >& >(
                   &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
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

    m.def( "create_filtered_observation_collection",
           py::overload_cast< const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationFilterBase >,
                              const std::shared_ptr< tom::ObservationCollectionParser > >(
                   &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_collection" ),
           py::arg( "observation_filter" ),
           py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
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
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, default = observations_processing.observation_parser()
            Parser to select the subset of observations to filter. Defaults to an empty parser (all observations).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the filtered observations.

        )doc" );

    m.def( "split_observation_collection",
           &tom::splitObservationSets< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "observation_set_splitter" ),
           py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
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
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, default = observations_processing.observation_parser()
            Parser to select which observation sets to split. Defaults to an empty parser (all sets).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection with the split observation sets.
        )doc" );

    m.def( "create_new_observation_collection",
           &tom::createNewObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg_v( "observation_parser", tom::observationParser( ), "..." ),
           R"doc(

        Creates a new observation collection containing a subset of an existing collection.

        This function selects a subset of observation sets from an original collection using a parser
        and creates a new collection containing only those sets.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The collection from which to extract a subset.
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, default = observations_processing.observation_parser()
            Parser to select the observation sets to include in the new collection. Defaults to an empty parser (all sets).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the selected subset of observation sets.
        )doc" );

    m.def( "create_compressed_doppler_collection",
           &tom::createCompressedDopplerCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "compression_ratio" ),
           py::arg( "minimum_number_of_observations" ) = 10,
           py::arg( "max_arc_gap" ) = 300.0,
           py::arg_v( "earth_fixed_ground_station_positions",
                      tss::getApproximateDsnGroundStationPositions( ),
                      "tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions()" ),
           R"doc(
        Create a compressed Doppler observation collection.

        This function takes a collection of Doppler observations, splits them into arcs, and compresses each arc by averaging observations.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The original observation collection containing Doppler data.
        compression_ratio : int
            The number of observations to average into a single compressed observation.
        minimum_number_of_observations : int, optional
            The minimum number of observations required in an arc to be considered for compression, by default 10.
        max_arc_gap : float, optional
            Maximum time gap (in seconds) between consecutive observations before splitting into a new arc, by default 300.0.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection with compressed Doppler data.
        )doc" );

    m.def( "create_count_time_compressed_doppler_collection",
           &tom::createCountTimeCompressedDopplerCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "bodies" ),
           py::arg( "target_count_time" ),
           py::arg( "minimum_number_of_observations" ) = 10,
           py::arg( "max_arc_gap" ) = 300.0,
           R"doc(
        Create a Doppler observation collection compressed to a target count time.

        Same as :func:`create_compressed_doppler_collection`, but with the number of observations to average determined
        separately for each observation set, from its count time and the ``target_count_time``. The ground station
        positions used for the time scale conversions are taken from the ground stations on Earth in ``bodies``.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The original observation collection containing Doppler data.
        bodies : tudatpy.dynamics.environment.SystemOfBodies
            System of bodies, containing the Earth with the ground stations of the Doppler observations.
        target_count_time : float
            Count time (in seconds) of the compressed observations, which must be an integer multiple of the count time of each
            observation set.
        minimum_number_of_observations : int, optional
            The minimum number of observations required in an arc to be considered for compression, by default 10.
        max_arc_gap : float, optional
            Maximum time gap (in seconds) between consecutive observations before splitting into a new arc, by default 300.0.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection with compressed Doppler data.
        )doc" );
}

}  // namespace observations
}  // namespace estimation
}  // namespace tudatpy
