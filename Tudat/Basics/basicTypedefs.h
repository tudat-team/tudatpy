/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_BASIC_TYPEDEFS_H
#define TUDAT_BASIC_TYPEDEFS_H

#include <Eigen/Core>

namespace Eigen
{

//! Typedef for Vector1d.
typedef Eigen::Matrix< double, 1, 1 > Vector1d;

//! Typedef for Vector5d.
typedef Eigen::Matrix< double, 5, 1 > Vector5d;

//! Typedef for Vector6d.
typedef Eigen::Matrix< double, 6, 1 > Vector6d;

//! Typedef for Vector6d.
typedef Eigen::Matrix< double, 7, 1 > Vector7d;

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

} // namespace Eigen

#endif // TUDAT_BASIC_TYPEDEFS_H
