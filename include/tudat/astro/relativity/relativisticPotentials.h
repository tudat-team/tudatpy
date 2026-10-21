/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RELATIVISTICPOTENTIALS_H
#define TUDAT_RELATIVISTICPOTENTIALS_H

#include <vector>
#include <map>
#include <string>

#include <functional>

#include <Eigen/Core>

#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{

namespace relativity
{

//! Function to calculate the leading term, i.e. of O(c^{0}), in Eq. (51) of Soffel et al. (2003).
/*!
 *  \param gravitationalParameters Gravitational parameters of perturbing bodies.
 *  \param distancesFromBody Distances from evaluation point to each perturbing body.
 *  \param sphericalHarmonicPotentials Optional additional scalar potential contributions keyed by body index.
 *  \return First-order scalar potential contribution per body.
 */
std::vector< double > calculateFirstOrderExternalScalarPotentials(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::map< int, double >& sphericalHarmonicPotentials = std::map< int, double >( ) );

//! Function to calculate the leading term, i.e. of O(c^{0}), in Eq. (51) of Soffel et al. (2003).
/*!
 *  \param gravitationalParameters Gravitational parameters of perturbing bodies.
 *  \param distancesFromBody Distances from evaluation point to each perturbing body.
 *  \param sphericalHarmonicPotentials Optional additional scalar potential contributions keyed by body index.
 *  \return Total first-order scalar potential.
 */
double calculateFirstOrderExternalScalarPotential(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::map< int, double >& sphericalHarmonicPotentials = std::map< int, double >( ) );

//! Function to calculate Eq. (51) of Soffel et al. (2003).
/*!
 *  \param gravitationalParameters Gravitational parameters of perturbing bodies.
 *  \param distancesFromBody Distances from evaluation point to each perturbing body.
 *  \param velocities Magnitudes of perturbing-body velocities.
 *  \param stateVectors Perturbing-body states in the global frame.
 *  \param centralBodyVelocity Central-body barycentric velocity.
 *  \param accelerations Optional perturbing-body accelerations.
 *  \return Second-order correction to the scalar potential.
 */
double calculateSecondOrderExternalScalarPotentialCorrection(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::vector< double >& velocities,
        const std::vector< Eigen::Vector6d >& stateVectors,
        const Eigen::Vector3d centralBodyVelocity,
        const std::vector< Eigen::Vector3d >& accelerations = std::vector< Eigen::Vector3d >( ) );

//! Function to calculate Eq. at bottom on left column of p. 2695, Soffel et al. (2003).
/*!
 *  \param gravitationalParameters Gravitational parameters of perturbing bodies.
 *  \param distancesFromBody Distances from evaluation point to each perturbing body.
 *  \param stateVectors Perturbing-body states in the global frame.
 *  \param currentAngularMomenta Optional body angular-momentum terms keyed by body index.
 *  \return Vector-potential contribution per body.
 */
std::vector< Eigen::Vector3d > calculateExternalVectorPotentials(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::vector< Eigen::Vector6d >& stateVectors,
        const std::map< int, Eigen::Vector3d >& currentAngularMomenta = std::map< int, Eigen::Vector3d >( ) );

//! Function to calculate Eq. at bottom on left column of p. 2695, Soffel et al. (2003).
/*!
 *  \param gravitationalParameters Gravitational parameters of perturbing bodies.
 *  \param distancesFromBody Distances from evaluation point to each perturbing body.
 *  \param stateVectors Perturbing-body states in the global frame.
 *  \param currentAngularMomenta Optional body angular-momentum terms keyed by body index.
 *  \return Total external vector potential.
 */
Eigen::Vector3d calculateExternalVectorPotential(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::vector< Eigen::Vector6d >& stateVectors,
        const std::map< int, Eigen::Vector3d >& currentAngularMomenta = std::map< int, Eigen::Vector3d >( ) );

//! Function to calculate vector potential from pre-computed scalar potentials and body states.
/*!
 *  \param externalBodyScalarPotentials Scalar potential contributions per perturbing body.
 *  \param stateVectors Perturbing-body states in the global frame.
 *  \param currentAngularMomenta Optional body angular-momentum terms keyed by body index.
 *  \param distancesFromBody Optional distances from evaluation point to each perturbing body.
 *  \return Total external vector potential.
 */
Eigen::Vector3d calculateExternalVectorPotential(
        const std::vector< double >& externalBodyScalarPotentials,
        const std::vector< Eigen::Vector6d >& stateVectors,
        const std::map< int, Eigen::Vector3d >& currentAngularMomenta = std::map< int, Eigen::Vector3d >( ),
        const std::vector< double >& distancesFromBody = std::vector< double >( ) );

//! Compute anisotropic matrix potential contribution for a single body.
/*!
 *  \param gravitationalParameter Body gravitational parameter.
 *  \param distanceFromBody Distance from evaluation point to body.
 *  \param relativePosition Relative position from body to evaluation point.
 *  \return Matrix potential contribution.
 */
Eigen::Matrix3d calculateExternalMatrixPotential( const double gravitationalParameter,
                                                  const double distanceFromBody,
                                                  const Eigen::Vector3d& relativePosition );

//! Compute anisotropic matrix potential contribution for multiple bodies.
/*!
 *  \param gravitationalParameters Gravitational parameters of perturbing bodies.
 *  \param distancesFromBody Distances from evaluation point to perturbing bodies.
 *  \param relativePositions Relative positions from each perturbing body to evaluation point.
 *  \return Sum of matrix potential contributions.
 */
Eigen::Matrix3d calculateExternalMatrixPotential( const std::vector< double >& gravitationalParameters,
                                                  const std::vector< double >& distancesFromBody,
                                                  const std::vector< Eigen::Vector3d >& relativePositions );
}  // namespace relativity

}  // namespace tudat

#endif  // TUDAT_RELATIVISTICPOTENTIALS_H
