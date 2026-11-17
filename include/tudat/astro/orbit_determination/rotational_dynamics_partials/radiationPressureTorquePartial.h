/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RADIATIONPRESSURETORQUEPARTIAL_H
#define TUDAT_RADIATIONPRESSURETORQUEPARTIAL_H

#include "tudat/astro/electromagnetism/radiationPressureTorque.h"
#include "tudat/astro/orbit_determination/rotational_dynamics_partials/torquePartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class to compute the partial derivatives of the radiation pressure torque of an isotropic point source.
/*!
 *  Class to compute the partial derivatives of the radiation pressure torque of an isotropic point source, for a cannonball or
 *  panelled target. The torque depends on the orientation and position of the target only through the radiation pressure and the
 *  source-target direction in the body-fixed frame. The partials are analytical, using the derivatives of the panel reflection laws.
 *  The dependency of self-shadowing illumination fractions and of occultation on the states is not included.
 */
class RadiationPressureTorquePartial : public TorquePartial
{
public:
    RadiationPressureTorquePartial( const std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureTorque > torqueModel,
                                    const std::string& acceleratedBody,
                                    const std::string& acceleratingBody );

    ~RadiationPressureTorquePartial( ) {}

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    void wrtOrientationOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                          const bool addContribution = 1,
                                          const int startRow = 0,
                                          const int startColumn = 0 )
    {
        partialMatrix.block( startRow, startColumn, 3, 4 ) += ( addContribution ? 1.0 : -1.0 ) * currentPartialWrtQuaternion_;
    }

    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes( const std::pair< std::string, std::string >& stateReferencePoint,
                                                                     const propagators::IntegratedStateType integratedStateType )
    {
        return ( stateReferencePoint.first == bodyUndergoingTorque_ || stateReferencePoint.first == bodyExertingTorque_ ) &&
                integratedStateType == propagators::translational_state;
    }

    void wrtNonRotationalStateOfAdditionalBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                const std::pair< std::string, std::string >& stateReferencePoint,
                                                const propagators::IntegratedStateType integratedStateType );

    void update( const double currentTime = TUDAT_NAN );

private:
    void wrtSourceDirectionScaling( Eigen::MatrixXd& partial )
    {
        partial = currentPartialWrtSourceDirectionScaling_;
    }

    void wrtPerpendicularSourceDirectionScaling( Eigen::MatrixXd& partial )
    {
        partial = currentPartialWrtPerpendicularScaling_;
    }

    void wrtRadiationPressureCoefficient( Eigen::MatrixXd& partial );

    std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureTorque > torqueModel_;

    std::shared_ptr< electromagnetism::RadiationPressureTargetModel > targetModel_;

    std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > paneledTargetModel_;

    std::vector< Eigen::Matrix3d > currentRotationMatrixDerivativesWrtQuaternion_;

    Eigen::Matrix< double, 3, 4 > currentPartialWrtQuaternion_;

    Eigen::Matrix3d currentPartialWrtTargetPosition_;

    Eigen::Vector3d currentPartialWrtSourceDirectionScaling_;

    Eigen::Vector3d currentPartialWrtPerpendicularScaling_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_RADIATIONPRESSURETORQUEPARTIAL_H
