/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_DATASET_ROWS_H
#define TUDAT_OBSERVATION_DATASET_ROWS_H

#include <string>

#include <Eigen/Core>

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/tudatTypeTraits.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class ObservationDataset;

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class ObservationDatasetViewer;

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class ObservationSelectionCondition;

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class FlattenedObservationData;

//! Metadata shared by all observations in one logical observation set.
/*!
 * An observation set is the dataset-backed replacement for the common metadata
 * formerly owned by a SingleObservationSet. Individual observation rows only
 * store a set id; this struct stores the observable type, link definition,
 * reference link end, observable size, ancillary-settings id and dependent-
 * variable bookkeeping id that apply to all observations in the set.
 */
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
struct ObservationSetMetadata {
    //! Observable type of the set, e.g. one-way range or angular position.
    ObservableType observableType_;

    //! Registry id of the full link definition shared by the set.
    unsigned int linkDefinitionId_;

    //! Link end whose time tag defines the observation time.
    LinkEndType referenceLinkEnd_;

    //! Observable size for one observation in this set, e.g. 1 for range and 2 for angular position.
    unsigned int observableSize_;

    //! Registry id of ancillary settings shared by the set; may point to nullptr.
    unsigned int ancillarySettingsId_;

    //! Registry id of dependent-variable layout/bookkeeping; may point to nullptr.
    unsigned int dependentVariableLayoutId_;
};

//! One row per observation event, independent of observable dimension.
/*!
 * A row represents one measurement epoch/event. Vector observables are not split
 * into multiple rows. Instead, scalarSize_ gives the number of scalar values in
 * this observation, and firstScalarComponent_ points to the first of those
 * values in the dataset-wide observed-value, residual-value and weight storage.
 * This is the distinction between the scientific observation event and the
 * scalar entries used by estimation vectors.
 */
template< typename TimeType = double >
struct ObservationDatasetRow {
    //! Observation time at the row's reference link end.
    TimeType time_;

    //! Metadata set to which this observation belongs.
    unsigned int setId_;

    //! Index of this observation's first scalar value in the dataset-wide scalar-value storage.
    unsigned int firstScalarComponent_;

    //! Observable size of this observation: the number of scalar values it contributes.
    unsigned int scalarSize_;

    //! Zero-based index of this observation within its set.
    unsigned int indexInSet_;

    //! Dependent-variable values computed for this observation event.
    Eigen::VectorXd dependentVariableValues_;

    //! Status flag used by flattened data objects that exclude inactive observations.
    bool isActive_;

    //! Optional human-readable reason for rejection or deactivation.
    std::string rejectionReason_;
};

//! Reverse mapping from scalar component storage to its observation event.
/*!
 * Each scalar component has one row in scalarComponentRows_. It maps a scalar
 * index in observedValues_/residualValues_ back to the owning observation event
 * and its component number within that event.
 */
struct ObservationScalarComponentRow {
    //! Observation event that owns this scalar component.
    unsigned int observationId_;

    //! Component number inside the owning observation event.
    unsigned int componentIndex_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_DATASET_ROWS_H
