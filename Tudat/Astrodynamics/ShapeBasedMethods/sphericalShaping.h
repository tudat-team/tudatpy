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

#ifndef SPHERICALSHAPING_H
#define SPHERICALSHAPING_H

#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionSphericalShaping.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <map>

namespace tudat
{
namespace shape_based_methods
{


class SphericalShaping
{
public:

    //! Constructor for spherical shaping.
    SphericalShaping(Eigen::Vector6d initialState,
                     Eigen::Vector6d finalState,
                     double requiredTimeOfFlight,
                     double initialValueCoefficientRadialInversePolynomial,
                     Eigen::VectorXd freeCoefficientsRadialFunction,
                     Eigen::VectorXd freeCoefficientsElevationFunction,
                     double centralBodyGravitationalParameter ,
                     root_finders::RootFinderType rootFinderType = root_finders::bisection_root_finder,
                     const double lowerBoundFreeCoefficient = TUDAT_NAN,
                     const double upperBoundFreeCoefficient = TUDAT_NAN,
                     const double initialGuessForFreeCoefficient = TUDAT_NAN,
                     const int maxNumberOfIterations = 30,
                     const double requiredToleranceForTimeOfFlight = 1.0e-6 );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SphericalShaping( ) { }

    //! Return the coefficients of the radial composite function.
    Eigen::VectorXd getRadialCompositionFunctionCoefficients( )
    {
        return freeCoefficientsRadialFunction_;
    }

    //! Return the coefficients of the elevation composite function.
    Eigen::VectorXd getElevationCompositionFunctionCoefficients( )
    {
        return freeCoefficientsElevationFunction_;
    }

    //! Return initial azimuthal angle.
    double getInitialAzimuthalAngle( )
    {
        return initialAzimuthalAngle_;
    }

    //! Return final azimuthal angle.
    double getFinalAzimuthalAngle( )
    {
        return finalAzimuthalAngle_;
    }

    //! Compute current spherical position.
    Eigen::Vector3d computeCurrentSphericalPosition( const double currentAzimuthalAngle );

    //! Compute current cartesian position.
    Eigen::Vector6d computeCurrentCartesianPosition( const double currentAzimuthalAngle );

    //! Compute first derivative of the azimuthal angle.
    double computeFirstDerivativeAzimuthalAngle( const double currentAzimuthalAngle );

    //! Compute second derivative of the azimuthal angle.
    double computeSecondDerivativeAzimuthalAngle( const double currentAzimuthalAngle );

    //! Compute current velocity in spherical coordinates.
    Eigen::Vector3d computeCurrentSphericalVelocity( const double currentAzimuthalAngle );

    //! Compute current velocity parametrized by azimuthal angle theta.
    Eigen::Vector3d computeCurrentVelocityParametrizedByAzimuthAngle( const double currentAzimuthalAngle );

    //! Compute current acceleration in spherical coordinates.
    Eigen::Vector3d computeCurrentSphericalAcceleration( const double currentAzimuthalAngle );

    //! Compute current acceleration parametrized by azimuthal angle theta.
    Eigen::Vector3d computeCurrentAccelerationParametrizedByAzimuthAngle( const double currentAzimuthalAngle );

    //! Compute control acceleration vector in spherical coordinates.
    Eigen::Vector3d computeSphericalControlAccelerationVector( const double currentAzimuthalAngle );

    //! Compute final deltaV.
    double computeDeltav( );

    //! Compute current spherical state.
    Eigen::Vector6d computeCurrentSphericalState( const double currentAzimuthalAngle );

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentCartesianState( const double currentAzimuthalAngle );

    void satisfyBoundaryConditions( );

    void resetValueFreeCoefficient( const double freeCoefficient )
    {
        initialValueCoefficientRadialInversePolynomial_ = freeCoefficient;
        freeCoefficientsRadialFunction_[ 2 ] = freeCoefficient;
    }

    //! Returns the required TOF
    double getRequiredTimeOfFlight( )
    {
        return requiredTimeOfFlight_;
    }

    //! Time of flight function for the root-finders.
    struct TOFfunction : public basic_mathematics::BasicFunction< double, double >
    {

        //! Create a function, where aMaximumDerivativeOrder is the maximum order of the derivative.
        TOFfunction( std::function< void ( const double ) > resetFreeCoefficientFunction,
                     std::function< void( ) > satisfyBoundaryConditionsFunction,
                     std::function< double ( ) >  computeTOFfunction,
                     std::function< double( ) > getRequiredTOFfunction ):
        resetFreeCoefficientFunction_( resetFreeCoefficientFunction ),
        satisfyBoundaryConditionsFunction_( satisfyBoundaryConditionsFunction ),
        computeTOFfunction_( computeTOFfunction ),
        getRequiredTOFfunction_( getRequiredTOFfunction ){ }

        //! Mathematical test function.
        double evaluate( const double inputValue )
        {
            resetFreeCoefficientFunction_( inputValue );
            satisfyBoundaryConditionsFunction_( );
            double currentTimeOfFlight = computeTOFfunction_( );

            return getRequiredTOFfunction_( ) - currentTimeOfFlight;
        }

        //! Derivatives of mathematical test function.
        double computeDerivative( const unsigned int order, const double inputValue )
        {
            throw std::runtime_error( "The rootfinder for TOF should not evaluate derivatives!" );
        }

        //! Crash on integration as root_finders should not execute these.
        double computeDefiniteIntegral( unsigned int order, double lowerBound, double upperbound )
        {
            throw std::runtime_error( "The rootfinder for TOF should not evaluate integrals!" );
        }

        //! Get the expected true location of the root.
        /*!
         * Not implemented.
         */
        double getTrueRootLocation( ) { return TUDAT_NAN; }

        //! Get the accuracy of the true location of the root.
        /*!
         * Not implemented.
         */
        double getTrueRootAccuracy( ) { return TUDAT_NAN; }

        //! Get a reasonable initial guess of the root location.
        /*!
         * Not implemented.
         */
        double getInitialGuess( ) { return TUDAT_NAN; }

        //! Get a reasonable lower boundary for the root location.
        /*!
         * Not implemented.
         */
        double getLowerBound( ) { return TUDAT_NAN; }

        //! Get a reasonable upper boundary for the root location.
        /*!
         * Not implemented.
         */
        double getUpperBound( ) { return TUDAT_NAN; }

    protected:

    private:            
        std::function< void ( const double ) > resetFreeCoefficientFunction_;
        std::function< void ( ) > satisfyBoundaryConditionsFunction_;
        std::function< double ( ) >  computeTOFfunction_;
        std::function< double ( ) > getRequiredTOFfunction_;
    };

    //! Compute time of flight.
    double computeTimeOfFlight();

    //! Iterate to match the required time of flight
    void iterateToMatchRequiredTimeOfFlight2( int maximumNumberOfIterations );

    //! Iterate to match the required time of flight
    void iterateToMatchRequiredTimeOfFlight( std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings,
                                             const double lowerBound = TUDAT_NAN,
                                             const double upperBound = TUDAT_NAN,
                                             const double initialGuess = TUDAT_NAN );


    void computeShapingTrajectoryAndFullPropagation(simulation_setup::NamedBodyMap& bodyMap,
            basic_astrodynamics::AccelerationMap& accelerationMap,
            const Eigen::Vector6d initialCartesianState,
            const std::string& centralBody,
            const std::string& bodyToPropagate,
            const propagators::TranslationalPropagatorType propagator/*propagators::cowell*/,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
            std::map< double, Eigen::VectorXd >& fullPropagationResults,
            std::map< double, Eigen::VectorXd >& shapingMethodResults,
            std::map<double, Eigen::VectorXd>& dependentVariables );

    double computeScalarFunctionTimeEquation( double currentAzimuthalAngle );
    double computeDerivativeScalarFunctionTimeEquation( double currentAzimuthalAngle );



protected:

    Eigen::MatrixXd computeInverseMatrixBoundaryConditions( );

    double computeInitialAlphaValue( );

    double computeFinalAlphaValue( );

    double computeInitialValueBoundaryConstant( );

    double computeFinalValueBoundaryConstant( );



private:

    //! Initial state in cartesian coordinates.
    Eigen::Vector6d initialState_;

    //! Final state in cartesian coordinates.
    Eigen::Vector6d finalState_;

    //! Initial state in spherical coordinates.
    Eigen::Vector6d initialStateSphericalCoordinates_;

    //! Final state in spherical coordinates.
    Eigen::Vector6d finalStateSphericalCoordinates_;

    //! Initial state parametrised w.r.t. the elevation angle theta.
    Eigen::Vector6d initialStateThetaParametrized_;

    //! Final state parametrised w.r.t. the elevation angle theta.
    Eigen::Vector6d finalStateThetaParametrized_;

    //! Pointer to the spherical shaping radial composite function.
    std::shared_ptr< CompositeRadialFunctionSphericalShaping > compositeRadialFunction_;

    //! Pointer to the spherical shaping elevation composite function.
    std::shared_ptr< CompositeElevationFunctionSphericalShaping > compositeElevationFunction_;

    //! Free coefficients for the radial composite function.
    Eigen::VectorXd freeCoefficientsRadialFunction_;

    //! Free coefficients for the elevation composite function.
    Eigen::VectorXd freeCoefficientsElevationFunction_;

    //! Targeted value for the time of flight.
    double requiredTimeOfFlight_;

    //! Initial guess for the coefficient of the second order component of the radial inverse polynomial.
    double initialValueCoefficientRadialInversePolynomial_;

    //! Initial value azimuthal angle.
    double initialAzimuthalAngle_;

    //! Final value azimuthal angle.
    double finalAzimuthalAngle_;

    //! Central body gravitational parameter.
    double centralBodyGravitationalParameter_;

    //! Root finder type.
    root_finders::RootFinderType rootFinderType_;

    //! Lower bound for the free coefficient, to be used when trying to match the required time of flight.
    const double lowerBoundFreeCoefficient_;

    //! Upper bound for the free coefficient, to be used when trying to match the required time of flight.
    const double upperBoundFreeCoefficient_;

    //! Initial guess for the value of the free coefficient.
    const double initialGuessForFreeCoefficient_;

    //! Maximum number of iterations for the root finder to find the free coefficient value that matches the required time of flight.
    const int maxNumberOfIterations_;

    //! Required tolerance between the actual time of flight and the required one.
    const double requiredToleranceForTimeOfFlight_;


    /*! Inverse of matrix containing the boundary conditions
     */
    Eigen::MatrixXd inverseMatrixBoundaryConditions;
};


} // namespace shape_based_methods
} // namespace tudat

#endif // SPHERICALSHAPING_H
