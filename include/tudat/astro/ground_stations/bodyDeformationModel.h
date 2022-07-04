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

namespace ground_stations
{

struct BodyDeformationStationMotionModel: public StationMotionModel
{
public:
    BodyDeformationStationMotionModel(
            std::function< std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >& ( ) > modelList ):
    modelList_( modelList ){ }

    Eigen::Vector6d getBodyFixedStationMotion( const double time )
    {
        Eigen::Vector6d motion = Eigen::Vector6d::Zero( );
        std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >& currentModels = modelList_( );
        for( unsigned int i = 0; i < currentModels.size( ); i++ )
        {
            motion.segment( 0, 3 ) += currentModels.at( i )->calculateDisplacement(
                        time, groundStationState_ );
        }
        return motion;
    }

    void setNominalStationState(
                const std::shared_ptr< ground_stations::GroundStationState > groundStationState )
    {
        groundStationState_ = groundStationState;
    }

protected:

    std::function< std::vector< std::shared_ptr< basic_astrodynamics::BodyDeformationModel > >& ( ) > modelList_;

    std::shared_ptr< ground_stations::GroundStationState > groundStationState_;
};

}

} // namespace tudat

#endif // TUDAT_BODYDEFORMATIONMODEL_H
