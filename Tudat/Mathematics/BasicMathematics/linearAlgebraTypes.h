/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_LINEAR_ALGEBRA_TYPES_H
#define TUDAT_LINEAR_ALGEBRA_TYPES_H

#include <Eigen/Core>

namespace tudat
{
namespace basic_mathematics
{

//! Typedef for Vector6d.
typedef Eigen::Matrix< double, 6, 1 > Vector6d;

//! Typedef for Vector6i.
typedef Eigen::Matrix< int, 6, 1 > Vector6i;

//! Typedef for Vector6f.
typedef Eigen::Matrix< float, 6, 1 > Vector6f;

//! Typedef for Matrix6d.
typedef Eigen::Matrix< double, 6, 6 > Matrix6d;

//! Typedef for Matrix6i.
typedef Eigen::Matrix< int, 6, 6 > Matrix6i;

//! Typedef for Matrix6f.
typedef Eigen::Matrix< float, 6, 6 > Matrix6f;

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_LINEAR_ALGEBRA_TYPES_H
