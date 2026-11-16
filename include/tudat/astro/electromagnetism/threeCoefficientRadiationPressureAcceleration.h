/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved.
 *
 *    This file is part of Tudat. Redistribution and use in source and binary
 *    forms, with or without modification, are permitted exclusively under the
 *    terms of the Modified BSD license.
 */

#ifndef TUDAT_THREECOEFFICIENTRADIATIONPRESSUREACCELERATION_H
#define TUDAT_THREECOEFFICIENTRADIATIONPRESSUREACCELERATION_H

#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"

namespace tudat
{
namespace electromagnetism
{

//! Radiation-pressure acceleration resolved in a frame tied to a reference body's source-centered orbit.
/*!
 * This class implements the three-coefficient model of McMahon and Scheeres (2015),
 * https://doi.org/10.2514/1.G000666. The basis is
 * \f$\left(\hat{\mathbf{U}},\hat{\mathbf{V}},\hat{\mathbf{W}}\right)\f$, with
 *
 * \f[
 * \hat{\mathbf{U}}=-\frac{\mathbf{r}_{R/S}}{\|\mathbf{r}_{R/S}\|},\qquad
 * \hat{\mathbf{W}}=\cos\phi\,\hat{\mathbf{Z}}-
 *     \sin\phi\left(\hat{\mathbf{Z}}\times\hat{\mathbf{U}}\right),\qquad
 * \hat{\mathbf{V}}=\hat{\mathbf{W}}\times\hat{\mathbf{U}},
 * \f]
 *
 * where \f$\hat{\mathbf{Z}}\f$ is the reference body's source-centered orbit-normal direction and
 * \f$\phi=23.4^\circ\f$, as used in the source model. For the published Earth/Sun model, the reference body is Earth
 * and the source is the Sun. All input states must be expressed in a common frame with inertial orientation. The
 * origin may be chosen freely because only relative states are used.
 */
class ThreeCoefficientRadiationPressureAcceleration : public RadiationPressureAcceleration
{
public:
    ThreeCoefficientRadiationPressureAcceleration( const std::shared_ptr< IsotropicPointRadiationSourceModel >& sourceModel,
                                                   const std::shared_ptr< basic_astrodynamics::BodyShapeModel >& sourceBodyShapeModel,
                                                   const std::function< Eigen::Vector3d( ) >& sourcePositionFunction,
                                                   const std::function< Eigen::Vector3d( ) >& sourceVelocityFunction,
                                                   const std::shared_ptr< RadiationPressureTargetModel >& targetModel,
                                                   const std::function< Eigen::Vector3d( ) >& targetPositionFunction,
                                                   const std::function< double( ) >& targetMassFunction,
                                                   const std::shared_ptr< OccultationModel >& sourceToTargetOccultationModel,
                                                   const std::function< Eigen::Vector3d( ) >& referenceBodyPositionFunction,
                                                   const std::function< Eigen::Vector3d( ) >& referenceBodyVelocityFunction,
                                                   const Eigen::Vector3d& coefficients,
                                                   const std::string& referenceBodyName );

    ~ThreeCoefficientRadiationPressureAcceleration( ) override = default;

    std::shared_ptr< RadiationSourceModel > getSourceModel( ) const override
    {
        return sourceModel_;
    }

    Eigen::Vector3d getCoefficients( ) const
    {
        return coefficients_;
    }

    void resetCoefficients( const Eigen::Vector3d& coefficients );

    Eigen::Matrix3d getCurrentBasis( ) const
    {
        return currentBasis_;
    }

    double getCurrentAccelerationScalingFactor( ) const
    {
        return currentAccelerationScalingFactor_;
    }

    Eigen::Vector3d getCurrentSourceToTargetPosition( ) const
    {
        return targetCenterPositionInSourceFrame_;
    }

    Eigen::Vector3d getCurrentSourceToReferencePosition( ) const
    {
        return currentSourceToReferencePosition_;
    }

    Eigen::Vector3d getCurrentSourceToReferenceVelocity( ) const
    {
        return currentSourceToReferenceVelocity_;
    }

    double getCurrentTargetMass( ) const
    {
        return currentTargetMass_;
    }

    double getSourceToTargetReceivedFraction( ) const
    {
        return sourceToTargetReceivedFraction_;
    }

    std::string getReferenceBodyName( ) const
    {
        return referenceBodyName_;
    }

private:
    void computeAcceleration( ) override;

    std::shared_ptr< IsotropicPointRadiationSourceModel > sourceModel_;
    std::shared_ptr< basic_astrodynamics::BodyShapeModel > sourceBodyShapeModel_;
    std::function< Eigen::Vector3d( ) > sourceVelocityFunction_;
    std::function< Eigen::Vector3d( ) > referenceBodyPositionFunction_;
    std::function< Eigen::Vector3d( ) > referenceBodyVelocityFunction_;

    Eigen::Vector3d coefficients_;
    const std::string referenceBodyName_;

    Eigen::Matrix3d currentBasis_;
    Eigen::Vector3d currentSourceToReferencePosition_;
    Eigen::Vector3d currentSourceToReferenceVelocity_;
    double currentAccelerationScalingFactor_;
    double currentTargetMass_;
    double sourceToTargetReceivedFraction_;
};

}  // namespace electromagnetism
}  // namespace tudat

#endif  // TUDAT_THREECOEFFICIENTRADIATIONPRESSUREACCELERATION_H
