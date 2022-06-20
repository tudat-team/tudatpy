/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *  References:
 *      "Hodographic-shaping method for low-thrust interplanetary trajectory design.", Gondelach, D. and Noomen, R.,
 *      Journal of Spacecraft and Rockets 52.3 (2015): 728-738
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

    //! Constructor which sets radial, normal and axial velocity functions and boundary conditions. Also sets quadrature
    //! settings.
    HodographicShapingLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const double centralBodyGravitationalParameter,
            const std::function< Eigen::Vector3d( ) > departureVelocityFunction,
            const std::function< Eigen::Vector3d( ) > arrivalVelocityFunction,
            const HodographicBasisFunctionList& radialVelocityFunctionComponents,
            const HodographicBasisFunctionList& normalVelocityFunctionComponents,
            const HodographicBasisFunctionList& axialVelocityFunctionComponents );

    //! Default destructor.
    ~HodographicShapingLeg( ) { }

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentCartesianState(const double timeSinceDeparture );

    //! Get single value of state
    virtual void getStateAlongTrajectory( Eigen::Vector6d& stateAlongTrajectory,
                                          const double time )
    {
        stateAlongTrajectory = computeCurrentCartesianState(time - departureTime_);
    }

    //! Compute magnitude thrust acceleration.
    double computeThrustAccelerationMagnitude(const double timeSinceDeparture );

    //! Compute direction thrust acceleration in cartesian coordinates.
    Eigen::Vector3d computeThrustAccelerationDirection(const double timeSinceDeparture );

    //! Compute the thrust acceleration in cartesian coordinates at the specified time.
    Eigen::Vector3d computeThrustAcceleration( const double timeSinceDeparture );

    //! Get single value of thrust acceleration.
    virtual void getThrustAccelerationAlongTrajectory(Eigen::Vector3d& thrustAccelerationAlongTrajectory,
                                                      const double time )
    {
        thrustAccelerationAlongTrajectory = computeThrustAcceleration(time - departureTime_);
    }

    // Return number of free radial coefficients
    int getNumberOfFreeRadialCoefficients( )
    {
        return numberOfFreeRadialCoefficients_;
    }

    // Return number of free normal coefficients
    int getNumberOfFreeNormalCoefficients( )
    {
        return numberOfFreeNormalCoefficients_;
    }

    // Return number of free axial coefficients
    int getNumberOfFreeAxialCoefficients( )
    {
        return numberOfFreeAxialCoefficients_;
    }

    // Return total number of free coefficients
    int getNumberOfFreeCoefficients( )
    {
        return numberOfFreeRadialCoefficients_ + numberOfFreeNormalCoefficients_ + numberOfFreeAxialCoefficients_;
    }

protected:

    //! Evaluate the transfer trajectory (i.e. do all the operations necessary to compute the DeltaV and compute it)
    void computeTransfer( );

private:

    //! Compute DeltaV.
    double computeDeltaV( );

    //! Update the value of the coefficients being used in the velocity functions, according to the provided leg parameters
    void updateFreeCoefficients( );

    //! Select value of first three coefficient, in order to meet the boundary conditions
    void satisfyBoundaryConditions( );

    //! Compute radial distance from the central body.
    double computeCurrentRadialDistance( const double timeSinceDeparture );

    //! Compute polar angle.
    double computeCurrentPolarAngle( const double timeSinceDeparture );

    //! Compute axial distance from central body.
    double computeCurrentAxialDistance( const double timeSinceDeparture );

    //! Compute inverse of matrix used to satisfy normal boundary conditions
    Eigen::Matrix2d computeInverseMatrixNormalBoundaries( std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction );

    //! Compute inverse of matrix used to satisfy radial or axial boundary conditions
    Eigen::Matrix3d computeInverseMatrixRadialOrAxialBoundaries( std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction );

    //! Satisfy boundary conditions in radial direction.
    void satisfyRadialBoundaryConditions( const Eigen::VectorXd& freeCoefficients );

    //! Satisfy boundary conditions in axial direction.
    void satisfyAxialBoundaryConditions( const Eigen::VectorXd& freeCoefficients );

    //! Satisfy boundary conditions in normal direction.
    void satisfyNormalBoundaryConditions( const Eigen::VectorXd& freeCoefficients );

    //! Define angular velocity due to the third component of the composite function only.
    double computeDerivativePolarAngleDueToThirdComponent(
            const double timeSinceDeparture, const Eigen::Vector2d& matrixK );

    //! Define the angular velocity due to all the other components of the composite function, once combined.
    double computeDerivativePolarAngleDueToOtherComponents(
            const double timeSinceDeparture, const Eigen::Vector2d& matrixL, const Eigen::VectorXd& freeCoefficients );

    //! Compute third fixed coefficient of the normal velocity composite function, so that the condition on the final polar angle
    //! is fulfilled.
    double computeThirdFixedCoefficientAxialVelocity ( const Eigen::VectorXd& freeCoefficients );

    //! Compute velocity vector in cylindrical coordinates.
    Eigen::Vector3d computeVelocityVectorInCylindricalCoordinates( const double timeSinceDeparture );

    //! Compute angular velocity.
    double evaluateDerivativePolarAngleWrtTime( const double timeSinceDeparture );

    //! Compute state vector in cylindrical coordinates.
    Eigen::Vector6d computeStateVectorInCylindricalCoordinates( const double timeSinceDeparture );

    //! Compute thrust acceleration vector in cylindrical coordinates.
    Eigen::Vector3d computeThrustAccelerationInCylindricalCoordinates( const double timeSinceDeparture );


    //! Central body gravitational parameter.
    const double centralBodyGravitationalParameter_;

    //! Function that outputs the departure velocity
    std::function< Eigen::Vector3d( ) > departureVelocityFunction_;

    //! Function that outputs the arrival velocity
    std::function< Eigen::Vector3d( ) > arrivalVelocityFunction_;

    //! Number of revolutions.
    int numberOfRevolutions_;

    //! Quadrature settings (used when computing multiple things)
    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_;

    //! Velocity functions.
    std::shared_ptr< CompositeFunctionHodographicShaping > radialVelocityFunction_;
    std::shared_ptr< CompositeFunctionHodographicShaping > normalVelocityFunction_;
    std::shared_ptr< CompositeFunctionHodographicShaping > axialVelocityFunction_;

    //! Coefficients of velocity functions (including
    //! fixed coefficients and free coefficients)
    Eigen::VectorXd fullCoefficientsRadialVelocityFunction_;
    Eigen::VectorXd fullCoefficientsNormalVelocityFunction_;
    Eigen::VectorXd fullCoefficientsAxialVelocityFunction_;

    //! Number of free coefficients
    const int numberOfFreeRadialCoefficients_;
    const int numberOfFreeNormalCoefficients_;
    const int numberOfFreeAxialCoefficients_;

    //! Boundary conditions.
    std::vector< double > radialBoundaryConditions_;
    std::vector< double > normalBoundaryConditions_;
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
    Eigen::MatrixXd inverseMatrixAxialBoundaryValues_;

    //! Previously computed thrust acceleration values
    std::map< double, Eigen::Vector3d > thrustAccelerationVectorCache_;
};


} // namespace shape_based_methods
} // namespace tudat

#endif
