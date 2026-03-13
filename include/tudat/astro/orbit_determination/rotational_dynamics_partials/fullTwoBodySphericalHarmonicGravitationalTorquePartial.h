/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FULLTWOBODYSPHERICALHARMONICGRAVITATIONALTORQUEPARTIAL_H
#define TUDAT_FULLTWOBODYSPHERICALHARMONICGRAVITATIONALTORQUEPARTIAL_H

#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicTorque.h"
#include "tudat/astro/orbit_determination/rotational_dynamics_partials/torquePartial.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class to calculate partials of full two-body spherical harmonic gravitational torque.
class FullTwoBodySphericalHarmonicGravitationalTorquePartial : public TorquePartial
{
public:
    FullTwoBodySphericalHarmonicGravitationalTorquePartial(
            const std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicTorque > torqueModel,
            const std::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
            const std::shared_ptr< simulation_setup::Body > bodyExertingTorque,
            const std::string& acceleratedBody,
            const std::string& acceleratingBody );

    ~FullTwoBodySphericalHarmonicGravitationalTorquePartial( ) { }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

    void wrtOrientationOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                          const bool addContribution = 1,
                                          const int startRow = 0,
                                          const int startColumn = 0 ) override;

    void wrtOrientationOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                           const bool addContribution = 1,
                                           const int startRow = 0,
                                           const int startColumn = 0 ) override;

    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) override;

    void wrtNonRotationalStateOfAdditionalBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                const std::pair< std::string, std::string >& stateReferencePoint,
                                                const propagators::IntegratedStateType integratedStateType ) override;

    void update( const double currentTime = TUDAT_NAN ) override;

private:
    void wrtSphericalHarmonicCoefficientParameter(
            Eigen::MatrixXd& partialMatrix,
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >& parameter,
            const Eigen::VectorXd& perturbation );

    std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicTorque > torqueModel_;

    std::shared_ptr< simulation_setup::Body > bodyUndergoingTorqueObject_;
    std::shared_ptr< simulation_setup::Body > bodyExertingTorqueObject_;

    Eigen::Matrix< double, 3, 4 > currentPartialWrtQuaternionOfBodyUndergoingTorque_;
    Eigen::Matrix< double, 3, 4 > currentPartialWrtQuaternionOfBodyExertingTorque_;
    Eigen::Matrix3d currentPartialWrtPositionOfBodyUndergoingTorque_;
    Eigen::Matrix3d currentPartialWrtPositionOfBodyExertingTorque_;

    Eigen::Vector4d quaternionPerturbation_;
    Eigen::Vector3d positionPerturbation_;
    double sphericalHarmonicCoefficientPerturbation_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_FULLTWOBODYSPHERICALHARMONICGRAVITATIONALTORQUEPARTIAL_H
