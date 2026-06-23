/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PIXEL_COORDINATES_POINTING_PARTIAL_H
#define TUDAT_PIXEL_COORDINATES_POINTING_PARTIAL_H

#include <memory>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/observation_partials/pixelCoordinatesPartial.h"

namespace tudat
{

namespace observation_partials
{

class PixelCoordinatesPointingPartial : public ObservationPartial< 2 >
{
public:
    PixelCoordinatesPointingPartial(
            const std::shared_ptr< PixelCoordinatesScaling > pixelCoordinatesScaling,
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > pointingCorrectionParameter );

    ~PixelCoordinatesPointingPartial( ) {}

    std::vector< std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, 2, 1 >& currentObservation = Eigen::Matrix< double, 2, 1 >::Constant( TUDAT_NAN ) ) override;

private:
    std::shared_ptr< PixelCoordinatesScaling > pixelCoordinatesScaling_;
    std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > pointingCorrectionParameter_;
    int receiverIndex_;
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_PIXEL_COORDINATES_POINTING_PARTIAL_H
