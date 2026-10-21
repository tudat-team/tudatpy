/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_THREECOEFFICIENTRADIATIONPRESSUREACCELERATIONPARTIAL_H
#define TUDAT_THREECOEFFICIENTRADIATIONPRESSUREACCELERATIONPARTIAL_H

#include "tudat/astro/electromagnetism/threeCoefficientRadiationPressureAcceleration.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{
namespace acceleration_partials
{

//! Analytical state and coefficient partials of the three-coefficient radiation-pressure acceleration.
class ThreeCoefficientRadiationPressureAccelerationPartial : public AccelerationPartial
{
public:
    using AccelerationPartial::getParameterPartialFunctionDerivedAcceleration;

    ThreeCoefficientRadiationPressureAccelerationPartial(
            const std::shared_ptr< electromagnetism::ThreeCoefficientRadiationPressureAcceleration >& accelerationModel,
            const std::string& acceleratedBody,
            const std::string& acceleratingBody );

    ~ThreeCoefficientRadiationPressureAccelerationPartial( ) override = default;

    void wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       bool addContribution = true,
                                       int startRow = 0,
                                       int startColumn = 0 ) override;

    void wrtVelocityOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       bool addContribution = true,
                                       int startRow = 0,
                                       int startColumn = 3 ) override;

    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        bool addContribution = true,
                                        int startRow = 0,
                                        int startColumn = 0 ) override;

    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        bool addContribution = true,
                                        int startRow = 0,
                                        int startColumn = 3 ) override;

    void wrtPositionOfAdditionalBody( const std::string& bodyName,
                                      Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                      bool addContribution = true,
                                      int startRow = 0,
                                      int startColumn = 0 ) override;

    void wrtVelocityOfAdditionalBody( const std::string& bodyName,
                                      Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                      bool addContribution = true,
                                      int startRow = 0,
                                      int startColumn = 3 ) override;

    bool isAccelerationPartialWrtAdditionalBodyNonnullptr( const std::string& bodyName ) override;

    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes( const std::pair< std::string, std::string >& stateReferencePoint,
                                                                     propagators::IntegratedStateType integratedStateType ) override;

    void wrtNonTranslationalStateOfAdditionalBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                   const std::pair< std::string, std::string >& stateReferencePoint,
                                                   propagators::IntegratedStateType integratedStateType,
                                                   bool addContribution = true ) override;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunctionDerivedAcceleration(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override;

    void update( double currentTime = 0.0 ) override;

private:
    void wrtThreeCoefficientRadiationPressureCoefficients( Eigen::MatrixXd& partial );

    void addPartial( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                     const Eigen::Matrix3d& partial,
                     bool addContribution,
                     int startRow,
                     int startColumn ) const;

    std::shared_ptr< electromagnetism::ThreeCoefficientRadiationPressureAcceleration > accelerationModel_;
    const std::string referenceBody_;

    Eigen::Matrix3d currentPartialWrtTargetPosition_;
    Eigen::Matrix3d currentPartialWrtSourcePosition_;
    Eigen::Matrix3d currentPartialWrtSourceVelocity_;
    Eigen::Matrix3d currentPartialWrtReferencePosition_;
    Eigen::Matrix3d currentPartialWrtReferenceVelocity_;
};

}  // namespace acceleration_partials
}  // namespace tudat

#endif  // TUDAT_THREECOEFFICIENTRADIATIONPRESSUREACCELERATIONPARTIAL_H
