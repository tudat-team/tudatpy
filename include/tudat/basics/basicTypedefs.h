/*    Copyright (c) 2010-2019, Delft University of Technology
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

//! Typedef for Vector6i.
typedef Eigen::Matrix< int, 6, 1 > Vector6i;

//! Typedef for Vector6f.
typedef Eigen::Matrix< float, 6, 1 > Vector6f;

//! Typedef for Vector6d.
typedef Eigen::Matrix< double, 6, 1 > Vector6d;

//! Typedef for Vector3ld.
typedef Eigen::Matrix< long double, 3, 1 > Vector3ld;

//! Typedef for Vector6ld.
typedef Eigen::Matrix< long double, 6, 1 > Vector6ld;

//! Typedef for Vector7d.
typedef Eigen::Matrix< double, 7, 1 > Vector7d;

//! Typedef for Vector7ld.
typedef Eigen::Matrix< long double, 7, 1 > Vector7ld;

//! Typedef for Matrix6i.
typedef Eigen::Matrix< int, 6, 6 > Matrix6i;

//! Typedef for Matrix6f.
typedef Eigen::Matrix< float, 6, 6 > Matrix6f;


//! Typedef for Matrix6d.
typedef Eigen::Matrix< double, 6, 6 > Matrix6d;

//! Typedef for Matrix7d.
typedef Eigen::Matrix< double, 7, 7 > Matrix7d;

//! Typedef for MatrixXi.
typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > MatrixXi;

//! Typedef for MatrixXl.
typedef Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic > MatrixXl;

//! Typedef for MatrixXll.
typedef Eigen::Matrix< long long, Eigen::Dynamic, Eigen::Dynamic > MatrixXll;

//! Typedef for MatrixXf.
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic > MatrixXf;

//! Typedef for MatrixXld.
typedef Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > MatrixXld;

//! Typedef for VectorXld.
typedef Eigen::Matrix< long double, Eigen::Dynamic, 1 > VectorXld;

} // namespace Eigen

#endif // TUDAT_BASIC_TYPEDEFS_H
