/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_HODOGRAPHIC_SHAPING_LEG_H
#define TUDAT_HODOGRAPHIC_SHAPING_LEG_H

#include "tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h"
#include "tudat/astro/low_thrust/shape_based/compositeFunctionHodographicShaping.h"
#include "tudat/astro/mission_segments/transferLeg.h"
#include "tudat/math/quadrature/createNumericalQuadrature.h"

namespace tudat
{
namespace shape_based_methods
{


class HodographicShapingLeg : public mission_segments::TransferLeg
{
public:

    //! Constructor which sets radial, normal and axial velocity functions and boundary conditions.
    HodographicShapingLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const double centralBodyGravitationalParameter,
            const int numberOfRevolutions,
            const HodographicBasisFunctionList& radialVelocityFunctionComponents,
            const HodographicBasisFunctionList& normalVelocityFunctionComponents,
            const HodographicBasisFunctionList& axialVelocityFunctionComponents );

    //! Default destructor.
    ~HodographicShapingLeg( ) { }

    //! Compute radial distance from the central body.
    double computeCurrentRadialDistance( const double timeSinceDeparture );

    //! Compute polar angle.
    double computeCurrentPolarAngle( double timeSinceDeparture );

    //! Compute axial distance from central body.
    double computeCurrentAxialDistance( const double timeSinceDeparture );

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentStateVector( const double timeSinceDeparture );

    virtual void getStateAlongTrajectory( Eigen::Vector6d& stateAlongTrajectory,
                                          const double time )
    {
        stateAlongTrajectory = computeCurrentStateVector( time - departureTime_ );
    }

    //! Compute DeltaV.
    double computeDeltaV( );

    //! Return thrust acceleration profile.
    void getCylindricalThrustAccelerationProfile(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::VectorXd >& thrustAccelerationProfile )
    {
        thrustAccelerationProfile.clear();

        for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
        {
            thrustAccelerationProfile[ epochsVector.at( i ) ] = computeThrustAccelerationInCylindricalCoordinates( epochsVector.at( i ) );
        }

    }

    //! Compute magnitude thrust acceleration.
    double computeCurrentThrustAccelerationMagnitude(
            const double timeSinceDeparture );

    //! Compute direction thrust acceleration in cartesian coordinates.
    Eigen::Vector3d computeCurrentThrustAccelerationDirection(
            double timeSinceDeparture );

    Eigen::Vector3d computeCurrentThrustAcceleration( double timeSinceDeparture );

    Eigen::Vector3d computeCurrentThrustAcceleration(
                const double currentTime,
                const double timeOffset  );



protected:

    void computeTransfer(  );

    void updateFreeCoefficients( );

    void satisfyConstraints( );


    Eigen::Matrix2d computeInverseMatrixNormalBoundaries( std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction );

    Eigen::Matrix3d computeInverseMatrixRadialOrAxialBoundaries( std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction );

    //! Satisfy boundary conditions in radial direction.
    void satisfyRadialBoundaryConditions( Eigen::VectorXd freeCoefficients );

    //! Satisfy boundary conditions in axial direction.
    void satisfyAxialBoundaryConditions( Eigen::VectorXd freeCoefficients );

    //! Satisfy boundary conditions in normal direction.
    void satisfyNormalBoundaryConditions( Eigen::VectorXd freeCoefficients );

    //! Define angular velocity due to the third component of the composite function only.
    double computeDerivativePolarAngleDueToThirdComponent(
            const double timeSinceDeparture, const Eigen::Vector2d& matrixK );

    //! Define the angular velocity due to all the other components of the composite function, once combined.
    double computeDerivativePolarAngleDueToOtherComponents(
            const double timeSinceDeparture, const Eigen::Vector2d& matrixL, const Eigen::VectorXd& freeCoefficients );

    //! Compute third fixed coefficient of the normal velocity composite function, so that the condition on the final polar angle
    //! is fulfilled.
    double computeThirdFixedCoefficientAxialVelocityFromFinalPolarAngle( Eigen::VectorXd freeCoefficients );

    //! Compute velocity vector in cylindrical coordinates.
    Eigen::Vector3d computeVelocityVectorInCylindricalCoordinates( double timeSinceDeparture );

    //! Compute angular velocity.
    double evaluateDerivativePolarAngleWrtTime( const double timeSinceDeparture );

    //! Compute state vector in cylindrical coordinates.
    Eigen::Vector6d computeStateVectorInCylindricalCoordinates( const double timeSinceDeparture );

    //! Compute thrust acceleration vector in cylindrical coordinates.
    Eigen::Vector3d computeThrustAccelerationInCylindricalCoordinates( double timeSinceDeparture );


private:

    //! Central body gravitational parameter.
    double centralBodyGravitationalParameter_;

    //! Number of revolutions.
    int numberOfRevolutions_;

    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_;

    //! Radial velocity function.
    std::shared_ptr< CompositeFunctionHodographicShaping > radialVelocityFunction_;

    //! Normal velocity function.
    std::shared_ptr< CompositeFunctionHodographicShaping > normalVelocityFunction_;

    //! Axial velocity function.
    std::shared_ptr< CompositeFunctionHodographicShaping > axialVelocityFunction_;

    //! Vector containing the coefficients of the radial function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction_;

    //! Vector containing the coefficients of the normal function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction_;

    //! Vector containing the coefficients of the axial function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction_;

    Eigen::VectorXd fullCoefficientsRadialVelocityFunction_;

    Eigen::VectorXd fullCoefficientsNormalVelocityFunction_;

    Eigen::VectorXd fullCoefficientsAxialVelocityFunction_;

    //! Number of free coefficients for the radial velocity function.
    int numberOfFreeRadialCoefficients;

    //! Number of free coefficients for the normal velocity function.
    int numberOfFreeNormalCoefficients;

    //! Number of free coefficients for the axial velocity function.
    int numberOfFreeAxialCoefficients;

    //! Radial boundary conditions.
    std::vector< double > radialBoundaryConditions_;

    //! Normal boundary conditions.
    std::vector< double > normalBoundaryConditions_;

    //! Axial boundary conditions.
    std::vector< double > axialBoundaryConditions_;

    /*! Inverse of matrix containing the boundary values of the terms in the radial velocity
     *  function which are used to satisfy the radial boundary conditions.
     */
    Eigen::MatrixXd inverseMatrixRadialBoundaryValues_;

    /*! Inverse of matrix containing the boundary values of the terms in the normal velocity
     *  function which are used to satisfy the normal boundary conditions.
     */
    Eigen::MatrixXd inverseMatrixNormalBoundaryValues_;

    /*! Inverse of matrix containing the boundary values of the terms in the axial velocity
     *  function which are used to satisfy the normal boundary conditions.
     */
    Eigen::MatrixXd inverseAxialMatrixBoundaryValues_;

    std::map< double, Eigen::Vector3d > thrustAccelerationVectorCache_;
};


} // namespace shape_based_methods
} // namespace tudat

#endif
