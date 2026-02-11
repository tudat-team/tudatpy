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
std::vector< double > calculateFirstOrderExternalScalarPotentials(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::map< int, double >& sphericalHarmonicPotentials = std::map< int, double >( ) );

//! Function to calculate the leading term, i.e. of O(c^{0}), in Eq. (51) of Soffel et al. (2003).
double calculateFirstOrderExternalScalarPotential(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::map< int, double >& sphericalHarmonicPotentials = std::map< int, double >( ) );

//! Function to calculate Eq. (51) of Soffel et al. (2003).
double calculateSecondOrderExternalScalarPotentialCorrection(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::vector< double >& velocities,
        const std::vector< Eigen::Vector6d >& stateVectors,
        const Eigen::Vector3d centralBodyVelocity,
        const std::vector< Eigen::Vector3d >& accelerations  = std::vector< Eigen::Vector3d >( ) );

//! Function to calculate Eq. at bottom on left column of p. 2695, Soffel et al. (2003).
std::vector< Eigen::Vector3d > calculateExternalVectorPotentials(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::vector< Eigen::Vector6d >& stateVectors,
        const std::map< int, Eigen::Vector3d >& currentAngularMomenta = std::map< int, Eigen::Vector3d >( ) );

//! Function to calculate Eq. at bottom on left column of p. 2695, Soffel et al. (2003).
Eigen::Vector3d calculateExternalVectorPotential(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::vector< Eigen::Vector6d >& stateVectors,
        const std::map< int, Eigen::Vector3d >& currentAngularMomenta = std::map< int, Eigen::Vector3d >( ) );

Eigen::Vector3d calculateExternalVectorPotential(
        const std::vector< double >& externalBodyScalarPotentials,
        const std::vector< Eigen::Vector6d >& stateVectors,
        const std::map< int, Eigen::Vector3d >& currentAngularMomenta = std::map< int, Eigen::Vector3d >( ),
        const std::vector< double >& distancesFromBody = std::vector< double >( ) );

Eigen::Matrix3d calculateExternalMatrixPotential(
        const double gravitationalParameter,
        const double distanceFromBody,
        const Eigen::Vector3d& relativePosition );

Eigen::Matrix3d calculateExternalMatrixPotential(
        const std::vector< double >& gravitationalParameters,
        const std::vector< double >& distancesFromBody,
        const std::vector< Eigen::Vector3d >& relativePositions );
}

}

#endif // TUDAT_RELATIVISTICPOTENTIALS_H
