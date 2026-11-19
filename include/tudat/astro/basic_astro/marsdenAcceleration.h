/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Marsden, B. G., Sekanina, Z., Yeomans, D. K. (1973), Comets and nongravitational forces. V,
 *          The Astronomical Journal, 78, 211.
 */

// header guard
#ifndef TUDAT_MARSDENACCELERATION_H
#define TUDAT_MARSDENACCELERATION_H

#include <memory>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <functional>
#include <iostream>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

namespace tudat
{
namespace basic_astrodynamics
{

// Marsden-type non-gravitational acceleration (for cometary outgassing / some Yarkovsky cases).
/*
 *  The acceleration is defined in the RSW frame (radial, transverse, normal) of the body w.r.t. the Sun:
 *
 *  a_RSW = [A1, A2, A3] * g(r)
 *  g(r)  = alpha * (r/r0)^(-m) * ( 1 + (r/r0)^n )^(-k)
 *
 *  with r the heliocentric distance in AU, evaluated at t - tau (in this function approximated with a linear motion law).
 *  The result is rotated to the integration (inertial) frame, using the RSW frame at time t.
 *  The model's parameters can be set by the user - default to standard Marsden Acceleration model's parameters
 */
class MarsdenAcceleration : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:
    // Constructor.
    /*
     *  \param A1 Radial acceleration coefficient [m/s^2].
     *  \param A2 Transverse acceleration coefficient [m/s^2].
     *  \param A3 Normal acceleration coefficient [m/s^2].
     *  \param alpha Normalization constant of g(r).
     *  \param r0 Scale distance of g(r) [AU].
     *  \param m Exponent m of g(r).
     *  \param n Exponent n of g(r).
     *  \param k Exponent k of g(r).
     *  \param relativeStateToSunFunction Function returning the Cartesian state (6D) of the body undergoing
     *         acceleration w.r.t. the Sun, in the integration frame [m, m/s].
     *  \param tau Time offset [s] between perihelion and peak activity (default 0). The heliocentric distance in
     *         g(r) is evaluated at t - tau (positive or negative), using a linear motion law:
     *         r(t - tau) = r(t) - v(t) * tau. The RSW rotation matrix is always taken at time t.
     */
    MarsdenAcceleration( const double A1,
                         const double A2,
                         const double A3,
                         const double alpha = 0.1112620426,
                         const double r0 = 2.808,
                         const double m = 2.15,
                         const double n = 5.093,
                         const double k = 4.6142,
                         const std::function< Eigen::Matrix< double, 6, 1 >( ) > relativeStateToSunFunction,
                         const double tau = 0.0 ):
        A1_( A1 ), A2_( A2 ), A3_( A3 ), alpha_( alpha ), r0_( r0 ), m_( m ), n_( n ), k_( k ),
        relativeStateToSunFunction_( relativeStateToSunFunction ), tau_( tau )
    {}

    //! Destructor
    ~MarsdenAcceleration( ) {}

    // Function to update constituent elements of the Marsden acceleration to current time
    /*
     * Function to update constituent elements of the Marsden acceleration to current time
     * \param currentTime Time to which the Marsden acceleration elements are to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )  // This checks that the time is different from the "last time" (so we are
                                                                // not recomputing at the same epoch)
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            // Distance (even when tau is different from 0) at which the rotation matrix gets evaluated
            currentRelativeState_ = relativeStateToSunFunction_( );  // Function defined in createAccelerationModels.cpp

            // Heliocentric position shifted when the tau parameter is != 0.
            // It evaluates the heliocentric position at t - tau, from a linear motion law: r(t - tau) = r(t) - v(t) * tau
            // (works for positive and negative tau; reduces to r(t) for tau = 0)
            // IT IS AN APPROXIMATION!! It would be better to be able to access the past or future state (how??)
            currentDelayedPosition_ = currentRelativeState_.segment( 0, 3 ) - currentRelativeState_.segment( 3, 3 ) * tau_;

            // Heliocentric distance in AU, at t - tau. Distance used in the computation of g(r)
            currentHeliocentricDistanceAu_ = currentDelayedPosition_.norm( ) / physical_constants::ASTRONOMICAL_UNIT;

            // Scaling function g(r)
            const double rRatio = currentHeliocentricDistanceAu_ / r0_;
            // Calculation of the g(r) function
            currentGFunction_ = alpha_ * std::pow( rRatio, -m_ ) * std::pow( 1.0 + std::pow( rRatio, n_ ), -k_ );

            // Rotation RSW -> inertial, from the heliocentric state
            currentRotationRswToInertial_ = reference_frames::getRswSatelliteCenteredToInertialFrameRotationMatrix( currentRelativeState_ );

            currentAccelerationInRswFrame_ = Eigen::Vector3d( A1_ * currentGFunction_, A2_ * currentGFunction_, A3_ * currentGFunction_ );
            currentAcceleration_ = currentRotationRswToInertial_ * currentAccelerationInRswFrame_;
        }
    }

    // Getters
    // the "const" indicates that these methods only read the object, they do not modify it
    double getCurrentHeliocentricDistanceAu( ) const
    {
        return currentHeliocentricDistanceAu_;
    }

    double getCurrentGFunction( ) const
    {
        return currentGFunction_;
    }

    Eigen::Vector3d getCurrentAccelerationInRswFrame( ) const
    {
        return currentAccelerationInRswFrame_;
    }

    Eigen::Matrix3d getCurrentRotationRswToInertial( ) const
    {
        return currentRotationRswToInertial_;
    }

    double getA1( ) const
    {
        return A1_;
    }

    double getA2( ) const
    {
        return A2_;
    }

    double getA3( ) const
    {
        return A3_;
    }
    // Not sure if we need to be able to retrieve this - maybe needed for the partials?
    /*
    double getAlpha( ) const { return alpha_; }
    double getR0( ) const { return r0_; }
    double getM( ) const { return m_; }
    double getN( ) const { return n_; }
    double getK( ) const { return k_; }
    */
    double getTau( ) const
    {
        return tau_;
    }

    // Partials of the inertial acceleration w.r.t. each coefficient (valid after the last updateMembers call).
    Eigen::Vector3d getCurrentPartialWrtA1( ) const
    {
        return currentRotationRswToInertial_.col( 0 ) * currentGFunction_;
    }

    Eigen::Vector3d getCurrentPartialWrtA2( ) const
    {
        return currentRotationRswToInertial_.col( 1 ) * currentGFunction_;
    }

    Eigen::Vector3d getCurrentPartialWrtA3( ) const
    {
        return currentRotationRswToInertial_.col( 2 ) * currentGFunction_;
    }

    // Setters (for the estimatable parameters)
    // Reset A1 [m/s^2].
    void resetA1( const double newA1 )
    {
        A1_ = newA1;
    }

    // Reset A2 [m/s^2].
    void resetA2( const double newA2 )
    {
        A2_ = newA2;
    }

    // Reset A3 [m/s^2].
    void resetA3( const double newA3 )
    {
        A3_ = newA3;
    }

    // Reset the time offset tau [s].
    void resetTau( const double newTau )
    {
        tau_ = newTau;
    }

private:
    //! Acceleration coefficients in RSW frame [m/s^2]
    double A1_;  // radial
    double A2_;  // transverse
    double A3_;  // normal

    //! Parameters of g(r)
    double alpha_;
    double r0_;
    double m_;
    double n_;
    double k_;

    //! Function returning the state of the body w.r.t. the Sun
    std::function< Eigen::Matrix< double, 6, 1 >( ) > relativeStateToSunFunction_;

    //! Current heliocentric state
    Eigen::Matrix< double, 6, 1 > currentRelativeState_;

    //! Time offset between perihelion and peak activity [s]
    double tau_;

    //! Current heliocentric position at t - tau (linear motion law) [m]
    Eigen::Vector3d currentDelayedPosition_;

    //! Current heliocentric distance at t - tau [AU]
    double currentHeliocentricDistanceAu_;

    //! Current value of g(r)
    double currentGFunction_;

    //! Current rotation from RSW frame to inertial frame
    Eigen::Matrix3d currentRotationRswToInertial_;

    //! Current acceleration in RSW frame
    Eigen::Vector3d currentAccelerationInRswFrame_;
};

}  // namespace basic_astrodynamics

}  // namespace tudat

#endif  // TUDAT_MARSDENACCELERATION_H