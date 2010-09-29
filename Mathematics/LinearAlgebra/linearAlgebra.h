/*
 *  EigenRoutines.h
 *  Propagator
 *
 *  Created by Melman, J.C.P. (jcpmelman) on 8/7/09.
 *  Copyright 2009 Delft University of Technology. All rights reserved.
 *
 */

#ifndef EIGENROUTINES_H_
#define EIGENROUTINES_H_

// Notice that coefficient access methods in Eigen have assertions
// checking the ranges. So if you do a lot of coefficient access, 
// these assertion can have an important cost. If you want to 
// save cost, define EIGEN_NO_DEBUG, and it won't check assertions.
//#ifndef EIGEN_NO_DEBUG
//#define EIGEN_NO_DEBUG
//#endif

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/QR>
#include <Eigen/Cholesky>

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

typedef VectorXd Vector;

double cosAngle(const Vector3d&, const Vector3d&);
double angle(const Vector3d&, const Vector3d&);
double average(const VectorXd&);
double standDev(const VectorXd&);

#endif // EIGENROUTINES_H
