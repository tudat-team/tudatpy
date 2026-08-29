/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RIGIDBODYPROPERTIES_H
#define TUDAT_RIGIDBODYPROPERTIES_H

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace gravitation
{

class GravityFieldModel;

}  // namespace gravitation

namespace simulation_setup
{

class RigidBodyProperties
{
public:
    RigidBodyProperties( );

    virtual ~RigidBodyProperties( );

    void update( const double currentTime );

    virtual void updateMass( const double currentTime ) = 0;

    virtual void updateMassDistribution( const double currentTime ) = 0;

    virtual void updateInertiaTensorDerivative( const Eigen::Vector5d& ) {}

    virtual void resetCurrentTime( );

    double getCurrentMass( );

    Eigen::Vector3d getCurrentCenterOfMass( );

    Eigen::Matrix3d getCurrentInertiaTensor( );

    Eigen::Matrix3d getCurrentDerivativeInertiaTensor( );

    bool isInertiaTensorAvailable( ) const;

    bool isInertiaTensorDerivativeAvailable( ) const;

    virtual void setCurrentMass( const double currentMass ) = 0;

    virtual void setIsBodyInPropagation( const bool isBodyInPropagation );

protected:
    double currentMass_;

    Eigen::Vector3d currentCenterOfMass_;

    Eigen::Matrix3d currentInertiaTensor_;

    Eigen::Matrix3d currentDerivativeInertiaTensor_;

    bool isBodyInPropagation_;

    bool isMassComputed_;

    bool isComComputed_;

    bool isInertiaTensorComputed_;

    bool isDerivativeInertiaTensorComputed_;

    bool isInertiaTensorAvailable_;

    bool isDerivativeInertiaTensorAvailable_;
};

class TimeDependentRigidBodyProperties : public RigidBodyProperties
{
public:
    TimeDependentRigidBodyProperties( const std::function< double( const double ) > massFunction,
                                      const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction = nullptr,
                                      const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction = nullptr );

    TimeDependentRigidBodyProperties( const double mass,
                                      const Eigen::Vector3d& centerOfMass = Eigen::Vector3d::Constant( TUDAT_NAN ),
                                      const Eigen::Matrix3d& inertiaTensor = Eigen::Matrix3d::Constant( TUDAT_NAN ) );

    virtual ~TimeDependentRigidBodyProperties( );

    virtual void resetCurrentTime( );

    std::function< double( const double ) > getMassFunction( );

    void setMassFunction( const std::function< double( const double ) > massFunction );

    virtual void setCurrentMass( const double currentMass );

    void setInertiaTensor( const Eigen::Matrix3d& inertiaTensor );

    virtual void updateMass( const double currentTime );

    virtual void updateMassDistribution( const double currentTime );

protected:
    std::function< double( const double ) > massFunction_;

    std::function< Eigen::Vector3d( const double ) > centerOfMassFunction_;

    std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction_;
};

class MassDependentRigidBodyProperties : public RigidBodyProperties
{
public:
    MassDependentRigidBodyProperties( const double currentMass,
                                      const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction,
                                      const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction );

    virtual ~MassDependentRigidBodyProperties( );

    virtual void updateMass( const double currentTime );

    virtual void updateMassDistribution( const double currentTime );

    virtual void setCurrentMass( const double currentMass );

protected:
    std::function< Eigen::Vector3d( const double ) > centerOfMassFunction_;

    std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction_;
};

//! Rigid-body properties whose mass distribution is derived from an associated gravity model.
/*!
 * This class is the runtime owner of the scaled mean moment of inertia and of all derived
 * inertia state. The gravity model supplies only gravitational data and keeps a non-owning
 * reverse link, installed by Body, so changes to that data can trigger synchronization.
 */
class FromGravityFieldRigidBodyProperties : public RigidBodyProperties
{
public:
    //! Create properties from gravity data, optionally enabling spherical-harmonic inertia.
    /*!
     * Direct C++ construction of a spherical-harmonic gravity field does not accept or retain a
     * scaled mean moment. Supply it here and install these properties on the same Body.
     */
    FromGravityFieldRigidBodyProperties( const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel,
                                         const double scaledMeanMomentOfInertia = TUDAT_NAN );

    virtual ~FromGravityFieldRigidBodyProperties( );

    virtual void resetCurrentTime( );

    virtual void updateMass( const double currentTime );

    virtual void updateMassDistribution( const double currentTime );

    virtual void updateInertiaTensorDerivative( const Eigen::Vector5d& derivativeDegreeTwoCoefficients );

    virtual void setCurrentMass( const double currentMass );

    virtual void setIsBodyInPropagation( const bool isBodyInPropagation );

    //! Reset the gravity field from which these properties are derived, retaining their owned configuration.
    void resetGravityFieldModel( const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel );

    //! Immediately synchronize mass after the linked gravity field's gravitational parameter changes.
    void synchronizeMassFromGravityField( );

    //! Immediately synchronize center of mass and inertia after linked gravity data change.
    void synchronizeMassDistributionFromGravityField( );

    //! Return the mean principal moment divided by mass times squared gravity reference radius.
    double getScaledMeanMomentOfInertia( ) const;

    //! Reset the owned scaled mean moment and immediately refresh gravity-derived inertia.
    void setScaledMeanMomentOfInertia( const double scaledMeanMomentOfInertia );

protected:
    //! Gravity data source retained by these derived properties.
    std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel_;

    //! Canonical runtime value; the gravity model does not retain a copy.
    double scaledMeanMomentOfInertia_;

    bool modelIsTimeDependent_;
};

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_RIGIDBODYPROPERTIES_H
