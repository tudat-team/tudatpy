/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */

#ifndef TUDAT_BODYDEFORMATIONMODEL_H
#define TUDAT_BODYDEFORMATIONMODEL_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <functional>
#include <memory>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/ground_stations/groundStationState.h"

namespace tudat
{

namespace basic_astrodynamics
{

class BodyDeformationModel
{
public:

    virtual ~BodyDeformationModel( ){ }

    virtual Eigen::Vector3d calculateDisplacement(
            const double time,
            const Eigen::Vector3d& bodyFixedPosition ) = 0;

    virtual Eigen::Vector3d calculateDisplacement(
            const double time,
            const std::shared_ptr< ground_stations::GroundStationState > stationState ) = 0;

};

} // namespace basic_astrodynamics

} // namespace tudat

#endif // TUDAT_BODYDEFORMATIONMODEL_H
