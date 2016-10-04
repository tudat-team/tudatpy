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

enum PositionElementTypes
{
    cartesian_position,
    spherical_position,
    geodetic_position
};

Eigen::Vector3d convertPositionElements(
        const Eigen::Vector3d& originalElements,
        const PositionElementTypes originalElementTypes,
        const PositionElementTypes convertedElementTypes,
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel = NULL,
        const double tolerance = 1.0E-4 );

}

}

#endif // TUDAT_STATEREPRESENTATIONCONVERSIONS_H
