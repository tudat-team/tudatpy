/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      Note that the exact implementation of Newton-Raphson as root finder should be updated if
 *      someone would want to use a different root-finding technique.
 *
 *      By default the eccentricity is used as the iteration procedure. This is because in
 *      optimizing a Cassini-like trajectory, the pericenter radius had about 2-4 NaN values in
 *      100000 times the gravity assist calculation. The eccentricity iteration had no NaN values
 *      for a similar run in which 100000 gravity assist calculations were done. Also the
 *      eccentricity seemed to require slightly less iterations (does not necessarily mean it is
 *      faster or more accurate).
 *
 */

#ifndef TUDAT_TRANSFER_TRAJECTORY_H
#define TUDAT_TRANSFER_TRAJECTORY_H



#include <Eigen/Core>

#include "tudat/astro/mission_segments/transferLeg.h"
#include "tudat/astro/mission_segments/transferNode.h"

#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

namespace tudat
{

namespace mission_segments
{

class TransferTrajectory
{
public:

    //! Class constructor
    TransferTrajectory(
            const std::vector< std::shared_ptr< TransferLeg > > legs,
            const std::vector< std::shared_ptr< TransferNode > > nodes ):
        legs_( legs ), nodes_( nodes ), isComputed_( false )
    {
        totalDeltaV_ = 0.0;
    }

    //! Update trajectory with new independent variables
    void evaluateTrajectory(
            const std::vector< double >& nodeTimes,
            const std::vector< Eigen::VectorXd >& legFreeParameters,
            const std::vector< Eigen::VectorXd >& nodeFreeParameters );

    //! Retrieve total trajectory Delta V
    double getTotalDeltaV( );

    double getNodeDeltaV( const int nodeIndex );

    double getLegDeltaV( const int legIndex );

    std::vector< double > getDeltaVPerNode( );

    std::vector< double > getDeltaVPerLeg( );

    double getTotalTimeOfFlight ( );

    int getNumberOfNodes( )
    {
        return static_cast< int >( nodes_.size( ) );
    }

    int getNumberOfLegs( )
    {
        return static_cast< int >( legs_.size( ) );
    }

    std::vector< std::shared_ptr< TransferLeg > > getLegs( )
    {
        return legs_;
    }

    std::vector< std::shared_ptr< TransferNode > > getNodes( )
    {
        return nodes_;
    }

    //! Get inertial cartesian position and velocity along full trajectory.
    void getStatesAlongTrajectoryPerLeg( std::vector< std::map< double, Eigen::Vector6d > >& statesAlongTrajectoryPerLeg,
                                        const int numberOfDataPointsPerLeg );

    std::vector< std::map< double, Eigen::Vector6d > > getStatesAlongTrajectoryPerLeg( const int numberOfDataPointsPerLeg );

    void getStatesAlongTrajectory( std::map< double, Eigen::Vector6d >& statesAlongTrajectory, const int numberOfDataPointsPerLeg );

    std::map< double, Eigen::Vector6d > getStatesAlongTrajectory( const int numberOfDataPointsPerLeg );

    //! Get inertial cartesian thrust acceleration along full trajectory
    void getInertialThrustAccelerationsAlongTrajectoryPerLeg(
            std::vector< std::map< double, Eigen::Vector3d > >& thrustAccelerationsAlongTrajectoryPerLeg,
            const int numberOfDataPointsPerLeg );

    std::vector< std::map< double, Eigen::Vector3d > > getInertialThrustAccelerationAlongTrajectoryPerLeg(const int numberOfDataPointsPerLeg );

    void getInertialThrustAccelerationsAlongTrajectory(std::map< double, Eigen::Vector3d >& thrustAccelerationsAlongTrajectory,
                                                       const int numberOfDataPointsPerLeg );

    std::map< double, Eigen::Vector3d > getInertialThrustAccelerationsAlongTrajectory(const int numberOfDataPointsPerLeg );

    //! Get RSW thrust acceleration along full trajectory
    std::map< double, Eigen::Vector3d > getRswThrustAccelerationsAlongTrajectory( const int numberOfDataPointsPerLeg );

    //! Get TNW thrust acceleration along full trajectory
    std::map< double, Eigen::Vector3d > getTnwThrustAccelerationsAlongTrajectory( const int numberOfDataPointsPerLeg );


private:

    //! Retrieve full set of parameters for single leg
    void getLegTotalParameters(
            const std::vector< double >& nodeTimes,
            const Eigen::VectorXd& legFreeParameters,
            const int legIndex,
            Eigen::VectorXd& legTotalParameters );

    //! Retrieve full set of parameters for single node
    void getNodeTotalParameters(
            const std::vector< double >& nodeTimes,
            const Eigen::VectorXd& nodeFreeParameters,
            const int nodeIndex,
            Eigen::VectorXd& nodeTotalParameters );

    //! List of legs
    const std::vector< std::shared_ptr< TransferLeg > > legs_;

    //! List of nodes
    const std::vector< std::shared_ptr< TransferNode > > nodes_;

    //! Total Delta V over trajectory
    double totalDeltaV_;

    //! Total time of flight over trajectory
    double totalTimeOfFlight_;

    //! Boolean defining whether the object is in a valid state (trajectory parameters have been set)
    bool isComputed_;
};


} // namespace mission_segments

} // namespace tudat

#endif // TUDAT_TRANSFER_TRAJECTORY_H
