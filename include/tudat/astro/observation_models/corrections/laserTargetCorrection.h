/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_LASER_TARGET_CORRECTION_H
#define TUDAT_LASER_TARGET_CORRECTION_H

#include <functional>
#include <memory>
#include <vector>
#include <Eigen/Geometry>
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"

namespace tudat
{
namespace observation_models
{

//! Mean optical path changes and the width of the equivalent one-way return, all in metres.
class LaserTargetResponse
{
public:
    //! Mean path changes for uplink and downlink, respectively.
    Eigen::Vector2d meanPathCorrections;

    //! Weighted standard deviation of half the round-trip path changes, not a range bias.
    double rangeSpread;

    //! Equivalent one-way centroid shift in metres.
    double getMeanRangeCorrection( ) const
    {
        return 0.5 * meanPathCorrections.sum( );
    }
};

//! Scalar spherical-target CoM correction or explicitly weighted optical points relative to the target endpoint.
class LaserTargetSignature
{
public:
    //! Return weights from reflection TDB and the two unit directions from target to observer, in target-fixed axes.
    using ReturnWeightFunction = std::function< std::vector< double >( double, const Eigen::Vector3d&, const Eigen::Vector3d& ) >;

    //! Positive published CoM distance is added to observed range, hence subtracted from computed surface range.
    explicit LaserTargetSignature( double centerOfMassDistance );

    //! Optical point coordinates are metres from the nominated endpoint, in the target rotation model's fixed axes.
    //! Weights describe the selected return centroid; an optional callback replaces them for the two viewing directions.
    //! A supplied internal optical path (metres, nonnegative) adds equally to each leg, without changing return width.
    LaserTargetSignature( std::vector< Eigen::Vector3d > opticalPointOffsets,
                          std::vector< double > returnWeights,
                          double internalOpticalPathLength = 0.0,
                          ReturnWeightFunction returnWeightFunction = nullptr );

    //! Only optical arrays require a target orientation; a spherical CoM correction is direction independent.
    bool requiresTargetRotation( ) const
    {
        return !opticalPointOffsets_.empty( );
    }

    //! Compute both mean leg changes and one-way spread from target-to-observer vectors in target-fixed axes (metres).
    //! No glass delay, incidence law, or detector response is inferred beyond the explicitly supplied optical path and weights.
    LaserTargetResponse calculateResponse( const Eigen::Vector3d& targetToTransmitter,
                                           const Eigen::Vector3d& targetToReceiver,
                                           double reflectionTime ) const;

private:
    //! Positive scalar CoM distance (metres); an empty point list selects this model.
    double centerOfMassDistance_;

    //! Common internal optical path added to each array leg, in metres, relative to the supplied optical points.
    double internalOpticalPathLength_;

    //! Optical points (metres) and nonnegative relative return weights in matching order.
    std::vector< Eigen::Vector3d > opticalPointOffsets_;

    std::vector< double > returnWeights_;

    //! Optional measured or otherwise explicitly specified directional return response.
    ReturnWeightFunction returnWeightFunction_;
};

//! Apply each leg's share of one common laser target response during the coupled light-time solution.
class LaserTargetLightTimeCorrection : public LightTimeCorrection
{
public:
    //! Share the physical model and use the target's rotation from inertial to target-fixed axes.
    LaserTargetLightTimeCorrection( std::shared_ptr< LaserTargetSignature > signature,
                                    std::function< Eigen::Quaterniond( double ) > rotationToTargetFrame );

    //! Return the current leg's mean optical path change divided by c, using all four retarded link events.
    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            unsigned int currentMultiLegTransmitterIndex = 0,
            std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings = nullptr ) override;

private:
    //! Shared immutable geometry/response definition and target orientation at each reflection event.
    std::shared_ptr< LaserTargetSignature > signature_;

    std::function< Eigen::Quaterniond( double ) > rotationToTargetFrame_;
};

}  // namespace observation_models
}  // namespace tudat
#endif
