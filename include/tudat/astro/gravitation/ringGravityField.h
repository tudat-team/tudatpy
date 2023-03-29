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
 *          Precise computation of acceleration due to uniform ring or disk, Toshio Fukushima (2010), Celestial Mechanics
 *          and Dynamical Astronomy, 108:339–356.
 */

#ifndef TUDAT_RINGGRAVITYFIELD_H
#define TUDAT_RINGGRAVITYFIELD_H

#include <memory>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <iostream>

#include "tudat/astro/gravitation/gravityFieldModel.h"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_d.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>

namespace tudat
{

namespace gravitation
{

//! Computes the gravitational potential of a one-dimensional ring.
/*!
 * Computes the gravitational potential of a one-dimensional ring, according to Fukushima (2010), eq. 1. The ring is
 * assumed to be contained in the xy plane, with center at the origin of the reference frame.
 *
 * @param positionOfBodySubjectToAcceleration Position of the body subject to the acceleration wrt the ring's body-fixed frame.
 * @param ringRadius Radius of the ring.
 * @param gravitationalParameter Gravitational parameter of the ring.
 * @param gravitationalConstant Universal gravitational constant.
 * @return Gravitational potential.
 */
double computeRingGravitationalPotential(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double ringRadius,
        const double gravitationalParameter,
        const double ellipticIntegralK );

//! Computes the gravitational acceleration of a one-dimensional ring
/*!
 * Computes the gravitational potential of a one-dimensional ring, according to Fukushima (2010), eqs. 16, 18, 30. The
 * ring is assumed to be contained in the xy plane, with center at the origin of the reference frame.
 *
 * Contrary to other formulations, the one by Fukushima (2010) is not singular in the ring's axis (i.e. x=y=0).
 * The elliptic integral S(m) is calculated from S(m) = ((2-m) K(m) - 2 E(m)) / m^2 if ellipticIntegralSFromDAndB is false,
 * and from S(m) = (D(m) - B(m)) / m if ellipticIntegralSFromDAndB is true. The second formulation is less sensitive to
 * numerical cancellation issues near the ring's axis.
 * A third even more robust formulation is also presented by Fukushima (2010), though that would require implementing
 * the method for integrating the S(m) elliptic integral; hence this formulation is not implemented here.
 *
 * @param positionOfBodySubjectToAcceleration Position of the body subject to the acceleration wrt the ring's body-fixed frame.
 * @param ringRadius Radius of the ring.
 * @param gravitationalParameter Gravitational parameter of the ring.
 * @param gravitationalConstant Universal gravitational constant.
 * @param ellipticIntegralSFromDAndB
 * @return Gravitational acceleration.
 */
Eigen::Vector3d computeRingGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double ringRadius,
        const double gravitationalParameter,
        const double ellipticIntegralB,
        const double ellipticIntegralE,
        const double ellipticIntegralS );

//! Cache object in which variables that are required for the computation of ring gravity field are stored.
class RingGravityCache
{
public:

    /*! Constructor.
     *
     * Constructor.
     * @param ringRadius Radius of the ring.
     * @param ellipticIntegralSFromDAndB Flag indicating whether to compute S(m) from D(m) and B(m) (if true),
     *      or from K(m) and E(m) (if false)
     */
    RingGravityCache(
            const double ringRadius,
            const bool ellipticIntegralSFromDAndB ):
        ringRadius_( ringRadius ),
        ellipticIntegralSFromDAndB_( ellipticIntegralSFromDAndB )
    {
        currentBodyFixedPosition_ = ( Eigen::Vector3d( ) << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN ).finished( );
    }

    /*! Update cached variables to current state.
     *
     * Update cached variables to current state.
     * @param currentBodyFixedPosition Current body fixed position.
     */
    void update( const Eigen::Vector3d& currentBodyFixedPosition );

    // Fucntion to get the elliptic integral K
    double getEllipticIntegralK( )
    {
        return currentEllipticIntegralK_;
    }

    // Fucntion to get the elliptic integral B
    double getEllipticIntegralB( )
    {
        return currentEllipticIntegralB_;
    }

    // Fucntion to get the elliptic integral S
    double getEllipticIntegralS( )
    {
        return currentEllipticIntegralS_;
    }

    // Fucntion to get the elliptic integral E
    double getEllipticIntegralE( )
    { 
        return currentEllipticIntegralE_;
    }

private:

    // Current body fixed position.
    Eigen::Vector3d currentBodyFixedPosition_;

    // Radius of the ring
    const double ringRadius_;

    // Flag indicating whether to compute S(m) from D(m) and B(m) (if true), or from K(m) and E(m) (if false)
    const bool ellipticIntegralSFromDAndB_;

    // Current value of the elliptic integral K(m)
    double currentEllipticIntegralK_;

    // Current value of the elliptic integral E(m)
    double currentEllipticIntegralE_;

    // Current value of the elliptic integral B(m)
    double currentEllipticIntegralB_;

    // Current value of the elliptic integral S(m)
    double currentEllipticIntegralS_;
};

} // namespace gravitation

} // namespace tudat

#endif //TUDAT_RINGGRAVITYFIELD_H
