/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <cereal/types/polymorphic.hpp>

#include "tudat/io/serialization/archives.h"

/**
 * @file registrationsMath.cpp
 * @brief Centralized cereal polymorphic registrations for math-related types.
 *
 * Includes registrations for:
 *   - InterpolatorSettings hierarchy (createInterpolator.h)
 *   - Interpolator / LookUpScheme hierarchies (interpolator.h, oneDimensionalInterpolator.h,
 *     lagrangeInterpolator.h, cubicSplineInterpolator.h, lookupScheme.h)
 *
 * Consumers should include registrations_math.h to retain this translation unit.
 */

#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/math/interpolators/lookupScheme.h"

// =====================================================================
// InterpolatorSettings hierarchy
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::interpolators::InterpolatorSettings )
CEREAL_REGISTER_TYPE( tudat::interpolators::LagrangeInterpolatorSettings )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::interpolators::InterpolatorSettings, tudat::interpolators::LagrangeInterpolatorSettings )

// =====================================================================
// Interpolator hierarchy (double, Eigen::MatrixXd) — used by the state
// transition and sensitivity matrix interface.
// =====================================================================

namespace tudat
{
namespace interpolators
{
using InterpolatorDoubleMatrixXd = Interpolator< double, Eigen::MatrixXd >;
using OneDimensionalInterpolatorDoubleMatrixXd = OneDimensionalInterpolator< double, Eigen::MatrixXd >;
using LagrangeInterpolatorDoubleMatrixXd = LagrangeInterpolator< double, Eigen::MatrixXd, double >;
using CubicSplineInterpolatorDoubleMatrixXd = CubicSplineInterpolator< double, Eigen::MatrixXd, double >;
}  // namespace interpolators
}  // namespace tudat

CEREAL_REGISTER_TYPE( tudat::interpolators::OneDimensionalInterpolatorDoubleMatrixXd )
CEREAL_REGISTER_TYPE( tudat::interpolators::LagrangeInterpolatorDoubleMatrixXd )
CEREAL_REGISTER_TYPE( tudat::interpolators::CubicSplineInterpolatorDoubleMatrixXd )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::interpolators::InterpolatorDoubleMatrixXd,
                                      tudat::interpolators::OneDimensionalInterpolatorDoubleMatrixXd )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::interpolators::OneDimensionalInterpolatorDoubleMatrixXd,
                                      tudat::interpolators::LagrangeInterpolatorDoubleMatrixXd )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::interpolators::OneDimensionalInterpolatorDoubleMatrixXd,
                                      tudat::interpolators::CubicSplineInterpolatorDoubleMatrixXd )

// =====================================================================
// LookUpScheme hierarchy (double)
// =====================================================================

namespace tudat
{
namespace interpolators
{
using LookUpSchemeDouble = LookUpScheme< double >;
using HuntingAlgorithmLookupSchemeDouble = HuntingAlgorithmLookupScheme< double >;
using BinarySearchLookupSchemeDouble = BinarySearchLookupScheme< double >;
}  // namespace interpolators
}  // namespace tudat

CEREAL_REGISTER_TYPE( tudat::interpolators::HuntingAlgorithmLookupSchemeDouble )
CEREAL_REGISTER_TYPE( tudat::interpolators::BinarySearchLookupSchemeDouble )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::interpolators::LookUpSchemeDouble, tudat::interpolators::HuntingAlgorithmLookupSchemeDouble )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::interpolators::LookUpSchemeDouble, tudat::interpolators::BinarySearchLookupSchemeDouble )

CEREAL_REGISTER_DYNAMIC_INIT( tudat_serialization_math )
