/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RELATIVISTIVTIMECONVERSION_H
#define TUDAT_RELATIVISTIVTIMECONVERSION_H

#include <Eigen/Core>
#include <vector>

#include "tudat/basics/basicTypedefs.h"
namespace tudat
{

namespace relativity
{

//! Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, from a speed and scalar potential
/*!
 * Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, from a speed and scalar potential, allowing
 * for the possibility of EP violation (nominally none)
 * \param relativeSpeed Speed of object for which proper time is to be computed, w.r.t. frame of which coordinate time is used
 * \param gravitationalScalarPotential Scalar potential at location of object for which proper time is to be computed
 * \param equivalencePrincipleLpiViolationParameter Violation parameter of equivalence principle (default to GR value of 0.0
 * \return Proper time rate minus one (d tau/dt - 1.0)
 */
double calculateFirstCentralBodyProperTimeRateDifference( const double computationPointSpeed,
                                                          const double gravitationalScalarPotential,
                                                          const double equivalencePrincipleLpiViolationParameter = 0.0 );

//! Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, for a 1/c^2 potential from a static mass monopole
/*!
 * Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, for a 1/c^2 potential from a static mass monopole
 * located at the origin of the reference frame of which the coordinate time is used.,
 * allowing for the possibility of EP violation (nominally none)
 * \param computationPointState Cartesian state of point where proper time rate is to be computed
 * \param centralBodyGravitationalParameter Newtonian gravitational potnetial of central body.
 * \param equivalencePrincipleLpiViolationParameter Violation parameter of equivalence principle (default to GR value of 0.0
 * \return Proper time rate minus one (d tau/dt - 1.0)
 */
double calculateFirstCentralBodyProperTimeRateDifference( const Eigen::Vector6d& computationPointState,
                                                          const std::vector< Eigen::Vector6d >& perturbedStates,
                                                          const std::vector< double >& centralBodyGravitationalParameters,
                                                          const double equivalencePrincipleLpiViolationParameter = 0.0 );

//! Function to calculate the 1st order integrand term in Eq. (58) of Soffel et al. (2003)
/*!
 *  Function to calculate the 1st order integrand term in Eq. (58) of Soffel et al. (2003), for barycentric to bodycentric coordinate time
 *  conversion.
 *  \param velocityVector Velocity vector of body in barycentric frame.
 *  \param gravitationalScalarPotential Gravitational (scalar) potential, evaluated at the center of the body (excluding its own
 * gravitational potential).
 *  \return First-order (1/c^2) contribution to barycentric to bodycentric coordinate time integral
 */
double calculateFirstOrderTcbToTcgIntegrand( const Eigen::Vector3d velocityVector, const double gravitationalScalarPotential );

//! Function to calculate the 1st order integrand term in Eq. (58) of Soffel et al. (2003)
/*!
 *  Function to calculate the 1st order integrand term in Eq. (58) of Soffel et al. (2003), for barycentric to bodycentric coordinate time
 *  conversion.
 *  \param barycentricSpeed Speed (norm of velocity vector) of body in barycentric frame.
 *  \param gravitationalScalarPotential Gravitational (scalar) potential, evaluated at the center of the body (excluding its own
 * gravitational potential).
 *  \return First-order (1/c^2) contribution to barycentric to bodycentric coordinate time integral
 */
double calculateFirstOrderTcbToTcgIntegrand( const double barycentricSpeed, const double gravitationalScalarPotential );

//! Function to calculate the 1st order direct term (i.e. not in integral) in Eq. (58) of Soffel et al. (2003)
double calculateFirstOrderTcbToTcgDirectCorrection( const Eigen::Vector3d& geocentricPositionVector,
                                                    const Eigen::Vector3d& barycentricVelocityVector );

//! Function to calculate the 2nd order integrand term in Eq. (58) of Soffel et al. (2003)
double calculateSecondOrderTcbToTcgIntegrand( const double gravitationalScalarPotential,
                                              const Eigen::Vector3d& barycentricVelocity,
                                              const Eigen::Vector3d& gravitationalVectorPotential,
                                              const double currentSecondOrderExternalPotentialCorrection );

//! Function to calculate the 2nd order integrand term in Eq. (58) of Soffel et al. (2003)
double calculateSecondOrderTcbToTcgIntegrand( const double barycentricSpeed,
                                              const double gravitationalScalarPotential,
                                              const Eigen::Vector3d& barycentricVelocity,
                                              const Eigen::Vector3d& gravitationalVectorPotential,
                                              const double currentSecondOrderExternalPotentialCorrection );

//! Function to calculate the 2nd order direct term (i.e. not in integral) in Eq. (58) of Soffel et al. (2003)
double calculateSecondOrderTcbToTcgDirectCorrection( const Eigen::Vector3d& geocentricPositionVector,
                                                     const Eigen::Vector3d& barycentricVelocityVector,
                                                     const double gravitationalScalarPotential,
                                                     const Eigen::Vector3d& gravitationalVectorPotential );

//! Function to calculate Eq. (22) of Turyshev et al. (2012).
double calculateFirstOrderPlanetocentricToTopocentricConversion(
        const Eigen::Vector3d planetFixedPointPosition,
        const Eigen::Vector3d inertialPointVelocity,
        const double localCentralBodyPotential,
        const Eigen::Vector3d centralBodyAcceleration = Eigen::Vector3d::Zero( ),
        const std::vector< Eigen::Vector3d >& relativeTidalBodyPositions = std::vector< Eigen::Vector3d >( ),
        const std::vector< double >& tidalBodyGravitationalParameters = std::vector< double >( ) );

}  // namespace relativity

}  // namespace tudat

#endif  // TUDAT_RELATIVISTIVTIMECONVERSION_H
