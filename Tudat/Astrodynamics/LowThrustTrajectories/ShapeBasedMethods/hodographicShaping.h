/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_HODOGRAPHIC_SHAPING_H
#define TUDAT_HODOGRAPHIC_SHAPING_H

#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/shapeBasedMethod.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/compositeFunctionHodographicShaping.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Mathematics/NumericalQuadrature/createNumericalQuadrature.h"
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <map>

namespace tudat
{
namespace shape_based_methods
{


class HodographicShaping : public ShapeBasedMethod
{
public:

    //! Constructor which sets radial, normal and axial velocity functions and boundary conditions.
    HodographicShaping(
            const Eigen::Vector6d& initialState,
            const Eigen::Vector6d& finalState,
            const double timeOfFlight,
            const double centralBodyGravitationalParameter,
            const int numberOfRevolutions,
            const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
            const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
            const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
            const Eigen::VectorXd& freeCoefficientsRadialVelocityFunction,
            const Eigen::VectorXd& freeCoefficientsNormalVelocityFunction,
            const Eigen::VectorXd& freeCoefficientsAxialVelocityFunction,
            const double initialBodyMass = TUDAT_NAN );

    //! Default destructor.
    ~HodographicShaping( ) { }

    //! Convert time to independent variable.
    double convertTimeToIndependentVariable( const double time )
    {
        return time;
    }

    //! Convert independent variable to time.
    double convertIndependentVariableToTime( const double independentVariable )
    {
        return independentVariable;
    }

    //! Returns initial value of the independent variable.
    double getInitialValueInpendentVariable( )
    {
        return 0.0;
    }

    //! Returns final value of the independent variable.
    double getFinalValueInpendentVariable( )
    {
        return timeOfFlight_;
    }

    //! Compute radial distance from the central body.
    double computeCurrentRadialDistance( const double currentTime );

    //! Compute polar angle.
    double computeCurrentPolarAngle( double currentTime );

    //! Compute axial distance from central body.
    double computeCurrentAxialDistance( const double currentTime );

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentStateVector( const double currentTime );

    //! Compute thrust acceleration in cartesian coordinates.
    Eigen::Vector3d computeThrustAccelerationVector( double currentTime );

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
            const double currentTime, std::function< double ( const double ) > specificImpulseFunction =
            [ ]( const double ){ return TUDAT_NAN; },
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings = nullptr );

protected:

    Eigen::Matrix2d computeInverseMatrixNormalBoundaries( std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction );

    Eigen::Matrix3d computeInverseMatrixRadialOrAxialBoundaries( std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction );

    //! Satisfy boundary conditions in radial direction.
    void satisfyRadialBoundaryConditions( Eigen::VectorXd freeCoefficients );

    //! Satisfy boundary conditions in axial direction.
    void satisfyAxialBoundaryConditions( Eigen::VectorXd freeCoefficients );

    //! Satisfy boundary conditions in normal direction.
    void satisfyNormalBoundaryConditions( Eigen::VectorXd freeCoefficients );

    //! Compute third fixed coefficient of the normal velocity composite function, so that the condition on the final polar angle
    //! is fulfilled.
    double computeThirdFixedCoefficientAxialVelocityFromFinalPolarAngle( Eigen::VectorXd freeCoefficients );

    //! Compute velocity vector in cylindrical coordinates.
    Eigen::Vector3d computeVelocityVectorInCylindricalCoordinates( double currentTime );

    //! Compute angular velocity.
    double evaluateDerivativePolarAngleWrtTime( const double currentTime );

    //! Compute state vector in cylindrical coordinates.
    Eigen::Vector6d computeStateVectorInCylindricalCoordinates( const double currentTime );

    //! Compute thrust acceleration vector in cylindrical coordinates.
    Eigen::Vector3d computeThrustAccelerationInCylindricalCoordinates( double currentTime );


    //! Compute direction thrust acceleration in cartesian coordinates.
    Eigen::Vector3d computeCurrentThrustAccelerationDirection(
            double currentTime, std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );


private:

    //! Central body gravitational parameter.
    double centralBodyGravitationalParameter_;

    //! Number of revolutions.
    int numberOfRevolutions_;

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

    //! Number of free coefficients for the radial velocity function.
    int numberOfFreeCoefficientsRadialVelocityFunction_;

    //! Number of free coefficients for the normal velocity function.
    int numberOfFreeCoefficientsNormalVelocityFunction_;

    //! Number of free coefficients for the axial velocity function.
    int numberOfFreeCoefficientsAxialVelocityFunction;

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

#endif // TUDAT_HODOGRAPHIC_SHAPING_H
