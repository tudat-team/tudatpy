/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATESTATETRANSITIONANDSENSITIVITYMATRIXINTERPOLATOR_H
#define TUDAT_CREATESTATETRANSITIONANDSENSITIVITYMATRIXINTERPOLATOR_H

#include <map>
#include <memory>

#include <Eigen/Core>

#include "tudat/math/interpolators/oneDimensionalInterpolator.h"

namespace tudat
{

namespace propagators
{

//! Function to create interpolators for state transition and sensitivity matrices from numerical results.
/*!
 * Function to create interpolators for state transition and sensitivity matrices from numerical results.
 * \param stateTransitionMatrixInterpolator Interpolator object for state transition matrix (returned by reference).
 * \param sensitivityMatrixInterpolator Interpolator object for sensitivity matrix (returned by reference).
 * \param variationalEquationsSolution Vector of two matrix histories. First vector entry
 *  is state transition matrix history, second entry is sensitivity matrix history.
 * \param clearRawSolution Boolean denoting whether to clear entries of variationalEquationsSolution after creation
 * of interpolators.
 */
void createStateTransitionAndSensitivityMatrixInterpolator(
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& stateTransitionMatrixInterpolator,
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& sensitivityMatrixInterpolator,
        std::map< double, Eigen::MatrixXd >& stateTransitionSolution,
        std::map< double, Eigen::MatrixXd >& sensitivitySolution,
        const bool clearRawSolution = 1 );

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_CREATESTATETRANSITIONANDSENSITIVITYMATRIXINTERPOLATOR_H
