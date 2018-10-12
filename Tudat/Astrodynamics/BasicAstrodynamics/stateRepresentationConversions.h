/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_STATEREPRESENTATIONCONVERSIONS_H
#define TUDAT_STATEREPRESENTATIONCONVERSIONS_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelQuaternionElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelModifiedRodriguesParameterElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelExponentialMapElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/attitudeElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"

namespace tudat
{

namespace coordinate_conversions
{

//! Enum defining available types of position representations
enum StateElementTypes
{
    cartesian_state,
    keplerian_state,
    modified_equinoctial_state,
    unified_state_model_quaternions_state,
    unified_state_model_modified_rodrigues_parameters_state,
    unified_state_model_exponential_map_state
};

//! Enum defining available types of position representations
enum PositionElementTypes
{
    cartesian_position,
    spherical_position,
    geodetic_position
};

//! Function to convert a position from one representation to another
/*!
 * Function to convert a position from one representation to another
 * \param originalElements Position in element type given by originalElementTypes
 * \param originalElementTypes Element type used for input.
 * \param convertedElementTypes Element type to which originalElements is to be converted.
 * \param shapeModel Shape model associated with position (only required for specific element types, e.g. geodetic)
 * default nullptr.
 * \param tolerance Tolerance used for conversion (only required for specific element types, e.g. geodetic), default 0.1 mm.
 * \return Position in requested element type.
 */
Eigen::Vector3d convertPositionElements(
        const Eigen::Vector3d& originalElements,
        const PositionElementTypes originalElementTypes,
        const PositionElementTypes convertedElementTypes,
        const std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel = nullptr,
        const double tolerance = 1.0E-4 );

}

}

#endif // TUDAT_STATEREPRESENTATIONCONVERSIONS_H
