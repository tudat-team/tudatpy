#ifndef TUDAT_STATEREPRESENTATIONCONVERSIONS_H
#define TUDAT_STATEREPRESENTATIONCONVERSIONS_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"

namespace tudat
{

namespace coordinate_conversions
{

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
 * default NULL.
 * \param tolerance Tolerance used for conversion (only required for specific element types, e.g. geodetic), default 0.1 mm.
 * \return Position in requested element type.
 */
Eigen::Vector3d convertPositionElements(
        const Eigen::Vector3d& originalElements,
        const PositionElementTypes originalElementTypes,
        const PositionElementTypes convertedElementTypes,
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel = NULL,
        const double tolerance = 1.0E-4 );

}

}

#endif // TUDAT_STATEREPRESENTATIONCONVERSIONS_H
