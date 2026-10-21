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

    virtual void resetCurrentTime( );

    double getCurrentMass( );

    Eigen::Vector3d getCurrentCenterOfMass( );

    Eigen::Matrix3d getCurrentInertiaTensor( );

    //! Function to retrieve whether these properties define an inertia tensor
    bool isInertiaTensorAvailable( ) const;

    virtual void setCurrentMass( const double currentMass ) = 0;

    virtual void setIsBodyInPropagation( const bool isBodyInPropagation );

protected:
    double currentMass_;

    Eigen::Vector3d currentCenterOfMass_;

    Eigen::Matrix3d currentInertiaTensor_;

    bool isBodyInPropagation_;

    bool isMassComputed_;

    bool isComComputed_;

    bool isInertiaTensorComputed_;

    //! Boolean denoting whether these properties define an inertia tensor
    bool isInertiaTensorAvailable_;
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

//! Rigid-body properties (mass, center of mass and inertia tensor) derived from a body's gravity field.
/*!
 * Rigid-body properties derived from a body's gravity field. The properties are recomputed from the gravity field when updated
 * (for a time-dependent gravity field, once per time step during propagation).
 */
class FromGravityFieldRigidBodyProperties : public RigidBodyProperties
{
public:
    //! Constructor, with scaled mean moment of inertia required to compute an inertia tensor from a spherical harmonic gravity field.
    FromGravityFieldRigidBodyProperties( const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel,
                                         const double scaledMeanMomentOfInertia = TUDAT_NAN );

    virtual ~FromGravityFieldRigidBodyProperties( );

    void resetCurrentTime( ) override;

    void updateMass( const double currentTime ) override;

    void updateMassDistribution( const double currentTime ) override;

    void setCurrentMass( const double currentMass ) override;

    void setIsBodyInPropagation( const bool isBodyInPropagation ) override;

    //! Reset the gravity field from which these properties are derived
    void resetGravityFieldModel( const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel );

    //! Return the scaled mean moment of inertia (I_xx + I_yy + I_zz) / ( 3 M R^2 )
    double getScaledMeanMomentOfInertia( ) const;

    //! Reset the scaled mean moment of inertia, and recompute the inertia tensor
    void setScaledMeanMomentOfInertia( const double scaledMeanMomentOfInertia );

protected:
    //! Function to compute the mass from the gravity field
    void computeMassFromGravityField( );

    //! Function to compute the center of mass and inertia tensor from the gravity field
    void computeMassDistributionFromGravityField( );

    //! Gravity field from which properties are derived.
    std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel_;

    //! Scaled mean moment of inertia, used with the degree-two coefficients to compute the inertia tensor.
    double scaledMeanMomentOfInertia_;

    bool modelIsTimeDependent_;

    //! Inertia tensor per unit mass of (homogeneous) polyhedron
    Eigen::Matrix3d polyhedronInertiaTensorPerUnitMass_;
};

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_RIGIDBODYPROPERTIES_H
