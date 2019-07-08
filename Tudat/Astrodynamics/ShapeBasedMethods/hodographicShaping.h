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

#ifndef HODOGRAPHICSHAPING_H
#define HODOGRAPHICSHAPING_H

#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionHodographicShaping.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Mathematics/NumericalQuadrature/createNumericalQuadrature.h"
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <map>

namespace tudat
{
namespace shape_based_methods
{


class HodographicShaping
{
public:

//    //! Constructor which sets radial, normal and axial velocity functions.
//    HodographicShaping( CompositeFunction& radialVelocityFunction,
//                        CompositeFunction& normalVelocityFunction,
//                        CompositeFunction& axialVelocityFunction ):
//        radialVelocityFunction_( radialVelocityFunction ),
//        normalVelocityFunction_( normalVelocityFunction ),
//        axialVelocityFunction_( axialVelocityFunction ),
//        numberOfIntegrationSteps_( 25 ) { }

    //! Constructor which sets radial, normal and axial velocity functions and boundary conditions.
    HodographicShaping( CompositeFunction& radialVelocityFunction,
                        CompositeFunction& normalVelocityFunction,
                        CompositeFunction& axialVelocityFunction,
                        Eigen::Vector6d initialState,
                        std::vector< double > boundaryConditionsRadial,
                        std::vector< double > boundaryConditionsNormal,
                        std::vector< double > boundaryConditionsAxial,
                        std::vector< double > boundaryValuesTime,
                        double initialTime,
                        double finalTime );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~HodographicShaping( ) { }

    //! Set radial, normal and axial velocity functions.
    void setVelocityFunctions( CompositeFunction& radialVelocityFunction,
                               CompositeFunction& normalVelocityFunction,
                               CompositeFunction& axialVelocityFunction )
    {
        radialVelocityFunction_ = radialVelocityFunction;
        normalVelocityFunction_ = normalVelocityFunction;
        axialVelocityFunction_  = axialVelocityFunction;
    }

    //! Set boundary conditions.
    void setBoundaryConditions( std::vector< double > boundaryConditionsRadial,
                                std::vector< double > boundaryConditionsNormal,
                                std::vector< double > boundaryConditionsAxial,
                                std::vector< double > boundaryValuesTime );

    //! Satisfy boundary conditions in radial direction.
    void satisfyRadialBoundaryConditions( Eigen::VectorXd freeCoefficients );

    //! Satisfy boundary conditions in normal direction.
    void satisfyNormalBoundaryConditions( Eigen::VectorXd freeCoefficients );

    //! Satisfy boundary conditions in axial direction.
    void satisfyAxialBoundaryConditions( Eigen::VectorXd freeCoefficients );

    //! Satisfy boundary conditions in normal direction with fixed final polar angle.
    void satisfyNormalBoundaryConditionsWithFinalPolarAngle( Eigen::VectorXd freeCoefficients );

    //! Compute thrust acceleration.
    double computeThrustAccelerationCurrentTime( const double currentTime );

    //! Compute cartesian acceleration.
    Eigen::Vector3d computeCartesianAcceleration( double currentTime );

    //! Convert cartesian acceleration into cylindrical one.
    Eigen::Vector3d convertCartesianToCylindricalAcceleration( Eigen::Vector3d cartesianAcceleration, Eigen::Vector3d cartesianState );

    //! Compute magnitude cartesian acceleration.
    double computeMagnitudeCartesianAcceleration( double currentTime );

    //! Compute direction cartesian acceleration.
    Eigen::Vector3d computeDirectionCartesianAcceleration( double currentTime );

    //! Compute DeltaV.
    /*!
     * Computes the required DeltaV by integrating the required thrust acceleration over time using a fast RK4-integrator.
     */
    double computeDeltaV( );

    //! Compute angular velocity.
    double computeAngularVelocityCurrentTime( const double currentTime );

    //! Compute final polar angle.
    /*!
     * Computes the final polar angle by integrating the angular velocity over time
       using a fast RK4-integrator.
     */
    double computeFinalPolarAngle( );

    //! Compute polar angle.
    double computePolarAngle( double currentTime );

    //! Compute components of the thrust acceleration.
    Eigen::Vector3d computeThrustAccelerationComponents( double currentTime );

    //! Compute thrust acceleration direction.
    Eigen::Vector3d computeThrustAccelerationDirection( double currentTime );

    //! Compute components of the velocity.
    std::vector< double > computeVelocityComponents( double currentTime );


    //! Compute radial distance from the Sun.
    double computeRadialDistanceCurrentTime( const double currentTime );

    //! Compute axial distance from central body.
    double computeAxialDistanceCurrentTime( const double currentTime );

    //! Compute current cylindrical state.
    Eigen::Vector6d computeCurrentCylindricalState( const double currentTime );

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentCartesianState( const double currentTime );

    //! Set velocity function for the radial direction.
    void setRadialVelocityFunction( CompositeFunction& radialVelocityFunction )
    {
        radialVelocityFunction_ = radialVelocityFunction;
    }

    //! Set velocity function for the normal direction.
    void setNormalVelocityFunction( CompositeFunction& normalVelocityFunction )
    {
        normalVelocityFunction_ = normalVelocityFunction;
    }

    //! Set velocity function for the axial direction.
    void setAxialVelocityFunction( CompositeFunction& axialVelocityFunction )
    {
        axialVelocityFunction_ = axialVelocityFunction;
    }

    //! Set boundary conditions for the radial direction.
    void setBoundaryConditionsRadial( std::vector< double > boundaryConditionsRadial )
    {
        boundaryConditionsRadial_ = boundaryConditionsRadial;
        numberOfFreeCoefficientsRadial = radialVelocityFunction_.getNumberOfCompositeFunctionComponents()
                                        - std::min( 3, static_cast< int >( boundaryConditionsRadial.size() ) );
    }

    //! Set boundary conditions for the normal direction.
    void setBoundaryConditionsNormal( std::vector< double > boundaryConditionsNormal )
    {
        boundaryConditionsNormal_ = boundaryConditionsNormal;
        numberOfFreeCoefficientsNormal = normalVelocityFunction_.getNumberOfCompositeFunctionComponents()
                                        - std::min( 3, static_cast< int >( boundaryConditionsNormal.size() ) );
    }

    //! Set boundary conditions for the axial direction.
    void setBoundaryConditionsAxial( std::vector< double > boundaryConditionsAxial )
    {
        boundaryConditionsAxial_ = boundaryConditionsAxial;
        numberOfFreeCoefficientsAxial = axialVelocityFunction_.getNumberOfCompositeFunctionComponents()
                                        - std::min( 3, static_cast< int > ( boundaryConditionsAxial.size() ) );
    }

    //! Get low-thrust acceleration model from shaping method.
    std::shared_ptr< propulsion::ThrustAcceleration > getLowThrustAccelerationModel(
            simulation_setup::NamedBodyMap& bodyMap,
            const std::string& bodyToPropagate,
            std::function< double( const double ) > specificImpulseFunction );

    //! Function to compute the shaped trajectory and the propagation fo the full problem.
    void computeShapingTrajectoryAndFullPropagation(simulation_setup::NamedBodyMap& bodyMap,
            basic_astrodynamics::AccelerationMap& accelerationMap,
            const std::string& centralBody,
            const std::string& bodyToPropagate,
            std::function< double ( const double ) > specificImpulseFunction,
            const std::shared_ptr<numerical_integrators::IntegratorSettings<double> > integratorSettings,
            std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
                                                    std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettings,
            std::map< double, Eigen::VectorXd >& fullPropagationResults,
            std::map< double, Eigen::VectorXd >& shapingMethodResults,
            std::map< double, Eigen::VectorXd>& dependentVariablesHistory,
            propagators::TranslationalPropagatorType propagatorType = propagators::cowell,
            const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< propagators::DependentVariableSaveSettings >( ) );



protected:

    Eigen::Matrix2d computeInverseMatrixNormalBoundaries( CompositeFunction& velocityFunction );

    Eigen::Matrix3d computeInverseMatrixRadialOrAxialBoundaries( CompositeFunction& velocityFunction );

    //! Compute first part of final polar angle integral.
    double computeFirstPartialFinalPolarAngle( );

    //! Compute second part of final polar angle integral.
    double computeSecondPartialFinalPolarAngle( Eigen::VectorXd freeCoefficients );

    //! Compute third fixed coefficient of the normal velocity composite function, so that the condition on the final polar angle
    //! is fulfilled.
    double computeThirdFixedCoefficientAxialVelocityFromFinalPolarAngle( Eigen::VectorXd freeCoefficients );

private:

    //! Initial state in cylindrical coordinates.
    Eigen::Vector6d initialState_;

    //! Radial velocity function.
    CompositeFunction& radialVelocityFunction_;

    //! Normal velocity function.
    CompositeFunction& normalVelocityFunction_;

    //! Axial velocity function.
    CompositeFunction& axialVelocityFunction_;

    //! Radial boundary conditions.
    std::vector< double > boundaryConditionsRadial_;

    //! Normal boundary conditions.
    std::vector< double > boundaryConditionsNormal_;

    //! Axial boundary conditions.
    std::vector< double > boundaryConditionsAxial_;

    //! Time boundary values.
    std::vector< double > boundaryValuesTime_;

    //! Initial time.
    double initialTime_;

    //! Final time.
    double finalTime_;

    //! Numerical quadrature settings (required to compute the final deltaV and polar angle values).
    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_;

    /*! Inverse of matrix containing the boundary values of the terms in the radial velocity
     *  function which are used to satisfy the radial boundary conditions.
     */
    Eigen::MatrixXd inverseMatrixBoundaryValuesRadial;

    /*! Inverse of matrix containing the boundary values of the terms in the normal velocity
     *  function which are used to satisfy the normal boundary conditions.
     */
    Eigen::MatrixXd inverseMatrixBoundaryValuesNormal;

    /*! Inverse of matrix containing the boundary values of the terms in the axial velocity
     *  function which are used to satisfy the normal boundary conditions.
     */
    Eigen::MatrixXd inverseMatrixBoundaryValuesAxial;

    //! Number of free coefficients for the radial velocity function.
    int numberOfFreeCoefficientsRadial;

    //! Number of free coefficients for the normal velocity function.
    int numberOfFreeCoefficientsNormal;

    //! Number of free coefficients for the axial velocity function.
    int numberOfFreeCoefficientsAxial;

};


} // namespace shape_based_methods
} // namespace tudat

#endif // HODOGRAPHICSHAPING_H
