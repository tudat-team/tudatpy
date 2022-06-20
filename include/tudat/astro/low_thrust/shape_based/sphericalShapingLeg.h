/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Roegiers, T., Application of the Spherical Shaping Method to a Low-Thrust Multiple Asteroid Rendezvous Mission,
 *          TU Delft (MSc thesis), 2014
 *
 */

#ifndef TUDAT_SPHERICAL_SHAPING_LEG_H
#define TUDAT_SPHERICAL_SHAPING_LEG_H

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/mission_segments/transferLeg.h"
#include "tudat/astro/low_thrust/shape_based/baseFunctionsSphericalShaping.h"
#include "tudat/astro/low_thrust/shape_based/compositeFunctionSphericalShaping.h"
#include "tudat/math/basic/basicFunction.h"
#include "tudat/math/root_finders/createRootFinder.h"
#include "tudat/math/quadrature/createNumericalQuadrature.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <map>

namespace tudat
{
namespace shape_based_methods
{


class SphericalShapingLeg : public mission_segments::TransferLeg {
public:

    //! Constructor for spherical shaping, with initial and final velocity given by departure and arrival bodies' ephemeris.
    SphericalShapingLeg(const std::shared_ptr<ephemerides::Ephemeris> departureBodyEphemeris,
                        const std::shared_ptr<ephemerides::Ephemeris> arrivalBodyEphemeris,
                        const double centralBodyGravitationalParameter,
                        const std::shared_ptr<root_finders::RootFinderSettings> rootFinderSettings,
                        const double lowerBoundFreeCoefficient = TUDAT_NAN,
                        const double upperBoundFreeCoefficient = TUDAT_NAN,
                        const double initialValueFreeCoefficient = TUDAT_NAN,
                        const double timeToAzimuthInterpolatorStepSize = physical_constants::JULIAN_DAY);

    //! Constructor for spherical shaping, with initial and final velocity given as arguments.
    SphericalShapingLeg(const std::shared_ptr<ephemerides::Ephemeris> departureBodyEphemeris,
                        const std::shared_ptr<ephemerides::Ephemeris> arrivalBodyEphemeris,
                        const double centralBodyGravitationalParameter,
                        const std::function< Eigen::Vector3d( ) > departureVelocityFunction,
                        const std::function< Eigen::Vector3d( ) > arrivalVelocityFunction,
                        const std::shared_ptr<root_finders::RootFinderSettings> rootFinderSettings,
                        const double lowerBoundFreeCoefficient = TUDAT_NAN,
                        const double upperBoundFreeCoefficient = TUDAT_NAN,
                        const double initialValueFreeCoefficient = TUDAT_NAN,
                        const double timeToAzimuthInterpolatorStepSize = physical_constants::JULIAN_DAY);

    //! Default destructor.
    ~SphericalShapingLeg() {}

    //! Compute dimensional current cartesian state.
    Eigen::Vector6d computeStateFromAzimuth(const double currentAzimuthAngle);

    Eigen::Vector6d computeState (const double timeSinceDeparture);

    virtual void getStateAlongTrajectory( Eigen::Vector6d& stateAlongTrajectory,
                                          const double time )
    {
        stateAlongTrajectory = computeState(time - departureTime_);
    }

    //! Compute thrust acceleration in cartesian coordinates at the current azimuth.
    Eigen::Vector3d computeThrustAccelerationFromAzimuth (const double currentAzimuthAngle );

    //! Compute thrust acceleration in cartesian coordinates at the current time.
    Eigen::Vector3d computeThrustAcceleration (const double timeSinceDeparture );

    //! Compute magnitude of thrust acceleration at the current azimuth.
    double computeThrustAccelerationMagnitudeFromAzimuth (const double currentAzimuthAngle );

    //! Compute magnitude of thrust acceleration at the current time.
    double computeThrustAccelerationMagnitude (const double timeSinceDeparture);

    //! Compute direction of thrust acceleration in cartesian coordinates at the current azimuth.
    Eigen::Vector3d computeThrustAccelerationDirectionFromAzimuth (const double currentAzimuthAngle );

    //! Compute direction of thrust acceleration in cartesian coordinates at the current time.
    Eigen::Vector3d computeThrustAccelerationDirection (const double timeSinceDeparture);

    //! Get single value of thrust acceleration.
    virtual void getThrustAccelerationAlongTrajectory(Eigen::Vector3d& thrustAccelerationAlongTrajectory,
                                                      const double time )
    {
        thrustAccelerationAlongTrajectory = computeThrustAcceleration(time - departureTime_);
    }

    //! Convert time to azimuth (azimuth is the independent variable used in the spherical shaping method).
    double convertTimeToAzimuth ( const double timeSinceDeparture );

    //! Convert azimuth to time (azimuth is the independent variable used in the spherical shaping method).
    double convertAzimuthToTime ( const double currentAzimuthAngle );

    //! Returns initial value of the azimuth (i.e. the independent variable).
    double getInitialValueAzimuth( )
    {
        return initialAzimuthAngle_;
    }

    //! Returns final value of the (i.e. the independent variable).
    double getFinalValueAzimuth( )
    {
        return finalAzimuthAngle_;
    }


protected:

    void computeTransfer( );

    //! Compute deltaV.
    double computeDeltaV( );

    //! Compute the inverse of the boundary conditions matrix.
    Eigen::MatrixXd computeInverseMatrixBoundaryConditions( );

    //! Compute the initial value of the constant alpha, as defined in Eq. 7.16 of Roegiers (2014) to express the boundary conditions.
    double computeValueConstantAlpha ( Eigen::Vector6d stateParametrizedByAzimuthAngle );

    //! Compute the initial value of the constant C, as defined in section 7.4 of Roegiers (2014) to express the boundary conditions.
    double computeValueConstantC ( Eigen::Vector6d stateParametrizedByAzimuthAngle, Eigen::Vector6d stateSphericalCoordinates );

    //! Ensure that the boundary conditions are respected.
    void satisfyBoundaryConditions( double freeCoefficient );

    //! Iterate to match the required time of flight, by updating the value of the free coefficient.
    void iterateToMatchRequiredTimeOfFlight( std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings,
                                             const double lowerBound = TUDAT_NAN,
                                             const double upperBound = TUDAT_NAN,
                                             const double initialGuess = TUDAT_NAN );

    //! Compute scalar function D (Eq. 7.62 of Roegiers (2014)) w.r.t. azimuth angle; used in the time equation.
    double computeScalarFunctionD ( double currentAzimuthAngle );

    //! Compute derivative of the scalar function D (Eq. 7.63 of Roegiers (2014)) w.r.t. azimuth angle; used in the time equation.
    double computeDerivativeScalarFunctionD ( double currentAzimuthAngle );

    //! Return the required value for the time of flight (normalized w.r.t. julian year).
    double getNormalizedRequiredTimeOfFlight( )
    {
        return timeOfFlight_ / physical_constants::JULIAN_YEAR;
    }

    //! Compute normalized time of flight.
    double computeNormalizedTimeOfFlight();

    //! Compute first derivative of the azimuth angle w.r.t. time.
    double computeFirstDerivativeAzimuthAngleWrtTime ( const double currentAzimuthAngle );

    //! Compute second derivative of the azimuth angle w.r.t. time.
    double computeSecondDerivativeAzimuthAngleWrtTime ( const double currentAzimuthAngle );

    //! Compute current normalized velocity in spherical coordinates parametrized by azimuth angle theta.
    Eigen::Vector3d computeNormalizedVelocityParametrizedByAzimuthAngle (const double currentAzimuthAngle );

    //! Compute normalized state vector in spherical coordinates.
    Eigen::Vector6d computeNormalizedStateInSphericalCoordinates (const double currentAzimuthAngle );

    //! Compute current normalized acceleration in spherical coordinates parametrized by azimuth angle theta.
    Eigen::Vector3d computeNormalizedThrustAccelerationParametrizedByAzimuthAngle (const double currentAzimuthAngle );

    //! Compute normalized thrust acceleration vector in spherical coordinates.
    Eigen::Vector3d computeNormalizedThrustAccelerationInSphericalCoordinates (const double currentAzimuthAngle );

    //! Time of flight function for the root-finder.
    struct TimeOfFlightFunction : public basic_mathematics::BasicFunction< double, double >
    {

        //! Create a function that returns the time of flight associated with the current trajectory.
        TimeOfFlightFunction( std::function< void( double ) > satisfyBoundaryConditionsFunction,
                              std::function< double ( ) >  computeTimeOfFlightFunction,
                              std::function< double( ) > getRequiredTimeOfFlightfunction ):
        satisfyBoundaryConditionsFunction_( satisfyBoundaryConditionsFunction ),
        computeTimeOfFlightFunction_( computeTimeOfFlightFunction ),
        getRequiredTimeOfFlightfunction_( getRequiredTimeOfFlightfunction ){ }

        //! Evaluate the difference between the current and required time of flight values, from the current value of the free coefficient.
        double evaluate( const double inputValue )
        {
            satisfyBoundaryConditionsFunction_( inputValue );
            double currentTimeOfFlight = computeTimeOfFlightFunction_( );

            return getRequiredTimeOfFlightfunction_( ) - currentTimeOfFlight;
        }

        //! Derivative function (not provided).
        double computeDerivative( const unsigned int order, const double inputValue )
        {
            throw std::runtime_error( "The rootfinder for TOF should not evaluate derivatives!" );
        }

        //! Integral function (not provided).
        double computeDefiniteIntegral( unsigned int order, double lowerBound, double upperbound )
        {
            throw std::runtime_error( "The rootfinder for TOF should not evaluate integrals!" );
        }

        //! Get the expected true location of the root (not implemented).
        double getTrueRootLocation( ) { return TUDAT_NAN; }

        //! Get the accuracy of the true location of the root (not implemented).
        double getTrueRootAccuracy( ) { return TUDAT_NAN; }

        //! Get a reasonable initial guess of the root location (not implemented).
        double getInitialGuess( ) { return TUDAT_NAN; }

        //! Get a reasonable lower boundary for the root location (not implemented).
        double getLowerBound( ) { return TUDAT_NAN; }

        //! Get a reasonable upper boundary for the root location (not implemented).
        double getUpperBound( ) { return TUDAT_NAN; }

    protected:

    private:

        std::function< void ( double ) > satisfyBoundaryConditionsFunction_;
        std::function< double ( ) >  computeTimeOfFlightFunction_;
        std::function< double ( ) > getRequiredTimeOfFlightfunction_;
    };



private:

    //! Central body gravitational parameter.
    double centralBodyGravitationalParameter_;

    //! Number of revolutions.
    int numberOfRevolutions_;

    //! Function that outputs the departure velocity
    std::function< Eigen::Vector3d( ) > departureVelocityFunction_;

    //! Function that outputs the arrival velocity
    std::function< Eigen::Vector3d( ) > arrivalVelocityFunction_;

    //! Initial state in spherical coordinates.
    Eigen::Vector6d initialStateSphericalCoordinates_;

    //! Final state in spherical coordinates.
    Eigen::Vector6d finalStateSphericalCoordinates_;

    //! Initial value azimuth angle.
    double initialAzimuthAngle_;

    //! Final value azimuth angle.
    double finalAzimuthAngle_;

    //! Initial state parametrised by the azimuth angle theta.
    Eigen::Vector6d initialStateParametrizedByAzimuthAngle_;

    //! Final state parametrised by the azimuth angle theta.
    Eigen::Vector6d finalStateParametrizedByAzimuthAngle_;

    //! Pointer to the spherical shaping radial distance composite function.
    std::shared_ptr< CompositeRadialFunctionSphericalShaping > radialDistanceCompositeFunction_;

    //! Pointer to the spherical shaping elevation angle composite function.
    std::shared_ptr< CompositeElevationFunctionSphericalShaping > elevationAngleCompositeFunction_;

    //! Coefficients for the radial distance composite function.
    Eigen::Vector7d coefficientsRadialDistanceFunction_;

    //! Coefficients for the elevation angle composite function.
    Eigen::Vector4d coefficientsElevationAngleFunction_;

    //! Root finder settings, to be used to find the free coefficient value that ensures the time of flight is correct.
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    //! Initial guess for the free coefficient (i.e. coefficient of the second order component of the radial inverse polynomial).
    const double initialValueFreeCoefficient_;

    //! Lower bound for the free coefficient, to be used when trying to match the required time of flight.
    const double lowerBoundFreeCoefficient_;

    //! Upper bound for the free coefficient, to be used when trying to match the required time of flight.
    const double upperBoundFreeCoefficient_;

    //! Step size used when creating the interpolator to convert time to azimuth
    const double timeToAzimuthInterpolatorStepSize_;

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > interpolator_;

    std::map< double, Eigen::Vector3d > thrustAccelerationVectorCache_;

    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_;

};


} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_SPHERICAL_SHAPING_LEG_H
