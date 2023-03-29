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
 * @param ellipticIntegralK Complete elliptic integral K.
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
 * @param ellipticIntegralB Complete elliptic integral B.
 * @param ellipticIntegralE Complete elliptic integral E.
 * @param ellipticIntegralS Complete elliptic integral S.
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
     * Update cached variables to current state. Computes the various elliptic integrals according to Fukushima (2010),
     * Eq. 31, 32, section A.1.
     *
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

//! Class to represent the gravity field of a constant density ring
class RingGravityField: public GravityFieldModel
{
public:

    /*! Constructor.
     *
     * Constructor.
     * @param gravitationalParameter Gravitational parameter of the ring.
     * @param ringRadius Radius of the ring.
     * @param ellipticIntegralSFromDAndB Flag indicating whether to compute S(m) from D(m) and B(m) (if true),
     *      or from K(m) and E(m) (if false). The former has a lower loss of accuracy due to numerical cancellation.
     * @param fixedReferenceFrame Identifier for body-fixed reference frame to which the field is fixed (optional).
     * @param updateInertiaTensor Function that is to be called to update the inertia tensor (typically in Body class;
     *      default empty)
     */
    RingGravityField(
            const double gravitationalParameter,
            const double ringRadius,
            const bool ellipticIntegralSFromDAndB,
            const std::string& fixedReferenceFrame = "",
            const std::function< void( ) > updateInertiaTensor = std::function< void( ) > ( ) ):
        GravityFieldModel(gravitationalParameter, updateInertiaTensor),
        gravitationalParameter_( gravitationalParameter ),
        ringRadius_( ringRadius ),
        ellipticIntegralSFromDAndB_( ellipticIntegralSFromDAndB ),
        fixedReferenceFrame_( fixedReferenceFrame )
    {
        ringGravityCache_ = std::make_shared< RingGravityCache >( ringRadius_, ellipticIntegralSFromDAndB_ );
    }

    /*! Function to calculate the gravitational potential.
     *
     * Function to calculate the gravitational potential.
     * @param bodyFixedPosition Position of point at which potential is to be calculated, in body-fixed frame.
     * @return Gravitational potential.
     */
    virtual double getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        ringGravityCache_->update( bodyFixedPosition );

        return computeRingGravitationalPotential(
                bodyFixedPosition,
                ringRadius_,
                gravitationalParameter_,
                ringGravityCache_->getEllipticIntegralK( ) );
    }

    /*! Function to calculate the gradient of the gravitational potential (i.e. the acceleration).
     *
     * Function to calculate the gradient of the gravitational potential (i.e. the acceleration).
     * @param bodyFixedPosition Position of point at which potential is to be calculated, in body-fixed frame.
     * @return Gradient of the gravitational potential.
     */
    virtual Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        ringGravityCache_->update( bodyFixedPosition );

        return computeRingGravitationalAcceleration(
                bodyFixedPosition,
                ringRadius_,
                gravitationalParameter_,
                ringGravityCache_->getEllipticIntegralB( ),
                ringGravityCache_->getEllipticIntegralE( ),
                ringGravityCache_->getEllipticIntegralS( ) );
    }

    //! Function to retrieve the identifier for the body-fixed reference frame.
    std::string getFixedReferenceFrame( )
    { return fixedReferenceFrame_; }

    //! Function to ring radius.
    double getRingRadius( )
    { return ringRadius_; }

    //! Function to the ellipticIntegralSFromDAndB flag
    bool getEllipticIntegralSFromDAndB( )
    { return ellipticIntegralSFromDAndB_; }

private:

    // Gravitational parameter
    const double gravitationalParameter_;

    // Radius of the ring
    const double ringRadius_;

    // Flag indicating whether to compute S(m) from D(m) and B(m) (if true), or from K(m) and E(m) (if false)
    const bool ellipticIntegralSFromDAndB_;

    //! Ring cache.
    std::shared_ptr< RingGravityCache > ringGravityCache_;

    //! Identifier for body-fixed reference frame
    std::string fixedReferenceFrame_;
};

} // namespace gravitation

} // namespace tudat

#endif //TUDAT_RINGGRAVITYFIELD_H
