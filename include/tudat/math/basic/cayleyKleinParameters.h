/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Eberly, D. Spherical Harmonics. Help documentation of Geometric Tools, 2008. Available at
 *        URL http://www.geometrictools.com/Documentation/SphericalHarmonics.pdf. Last access:
 *        09-09-2012.
 *      Heiskanen, W.A., Moritz, H. Physical geodesy. Freeman, 1967.
 *      Holmes, S.A., Featherstone, W.E. A unified approach to the Clenshaw summation and the
 *        recursive computation of very high degree and order normalised associated Legendre
 *        functions. Journal of Geodesy, 76(5):279-299, 2002.
 *      Vallado, D. and McClain, W. Fundamentals of astrodynammics and applications. Microcosm
 *        Press, 2001.
 *      Weisstein, E.W. Associated Legendre Polynomial, 2012.
 *        URL http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html. Last access:
 *        12-09-2012.
 *
 *    Notes
 *      For information on how the caching mechanism works, please contact S. Billemont
 *      (S.Billemont@studelft.tudelft.nl).
 *
 */

#ifndef TUDAT_CALYLEY_KLEIN_PARAMETERS_H
#define TUDAT_CALYLEY_KLEIN_PARAMETERS_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <complex>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace basic_mathematics
{

void convertQuaterionToCayleyKleinParameters(
        const Eigen::Quaterniond quaternion, std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB );

void convert323EulerAnglesToCayleyKleinParameters(
        const double firstZRotation, const double yRotation, const double secondZRotation,
        std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB );

void convert313EulerAnglesToCayleyKleinParameters(
        const double firstZRotation, const double xRotation, const double secondZRotation,
        std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB );

} // namespace basic_mathematics

} // namespace tudat

#endif // TUDAT_CALYLEY_KLEIN_PARAMETERS_H
