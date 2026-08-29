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

class FromGravityFieldRigidBodyProperties : public RigidBodyProperties
{
public:
    FromGravityFieldRigidBodyProperties( const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel );

    virtual ~FromGravityFieldRigidBodyProperties( );

    virtual void resetCurrentTime( );

    virtual void updateMass( const double currentTime );

    virtual void updateMassDistribution( const double currentTime );

    virtual void updateInertiaTensorDerivative( const Eigen::Vector5d& derivativeDegreeTwoCoefficients );

    virtual void setCurrentMass( const double currentMass );

    virtual void setIsBodyInPropagation( const bool isBodyInPropagation );

protected:
    const std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel_;

    bool modelIsTimeDependent_;
};

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_RIGIDBODYPROPERTIES_H
