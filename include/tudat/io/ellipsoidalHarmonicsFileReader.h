/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Garmier, R., Barriot, J.-P. (2001), Celestial Mechanics and Dynamical Astronomy, 79(4), 235-275.
 *          Hu, X., Jekeli, C. (2015), Journal of Geodesy, 89(2), 159-177.
 *          See tudat/math/basic/ellipsoidalHarmonics.h for the full citations.
 */

#ifndef TUDAT_ELLIPSOIDAL_HARMONICS_FILE_READER_H
#define TUDAT_ELLIPSOIDAL_HARMONICS_FILE_READER_H

/*  Readers for the three text files that define an ellipsoidal harmonic gravity field: the
 *  field coefficients, the Lame polynomials and the transform factors. See
 *  tudat/math/basic/ellipsoidalHarmonics.h for what the model does with them.
 *
 *  All three use Fortran-style exponents ("0.123D+00" as well as "0.123E+00"), which the
 *  standard Tudat text readers do not accept, hence the dedicated readers here.
 */

#include <string>

#include "tudat/math/basic/ellipsoidalHarmonics.h"

namespace tudat
{

namespace input_output
{

//! Reads a number written with either a Fortran "D" exponent or a normal "E" exponent.
double parseFortranDouble( const std::string& text );

/*! Reads the field coefficients of an ellipsoidal harmonic gravity field.
 *
 * The file lists one line per degree and order, "degree order cosineTerm sineTerm", in the same
 * layout a spherical harmonic coefficient file uses. The 2*degree+1 values of each degree are
 * then regrouped into the four term groups the expansion actually uses.
 * Throws if any term up to the maximum degree is missing from the file.
 *
 * @param fileName Full path of the coefficient file.
 * @param maximumDegree Highest degree to read.
 * @return Coefficients of every term, as read, still normalised to the reference ellipsoid.
 */
basic_mathematics::EllipsoidalHarmonicCoefficients readEllipsoidalHarmonicCoefficients( const std::string& fileName,
                                                                                        const int maximumDegree );

/*! Reads the Lame polynomials of an ellipsoidal harmonic gravity field.
 *
 * The file lists one line per term, "degree orderIndex groupLabel coefficient...". The
 * coefficients are stored in powers of (1 - (rho/h)^2) and are converted here to powers of
 * (rho/h)^2, which is the form the Lame functions are evaluated in.
 * Throws if any term up to the maximum degree is missing from the file.
 *
 * @param fileName Full path of the Lame coefficient file.
 * @param maximumDegree Highest degree to read.
 * @return One square matrix per group per degree, where column j holds the polynomial of term
 *      j in ascending powers of (rho/h)^2.
 */
basic_mathematics::EllipsoidalHarmonicMatrices readLameFunctionCoefficients( const std::string& fileName, const int maximumDegree );

/*! Reads the transform factors of an ellipsoidal harmonic gravity field.
 *
 * These mix the terms within each group, and the file has the same layout as the Lame
 * coefficient file, one line per term holding that term's row of its group's matrix.
 * Throws if any term up to the maximum degree is missing from the file.
 *
 * @param fileName Full path of the transform factor file.
 * @param maximumDegree Highest degree to read.
 * @return One square matrix per group per degree.
 */
basic_mathematics::EllipsoidalHarmonicMatrices readEllipsoidalHarmonicTransformFactors( const std::string& fileName,
                                                                                        const int maximumDegree );

/*! Converts polynomial coefficients from powers of (1 - x) to powers of x.
 *
 * The Lame coefficient file stores each polynomial in powers of (1 - (rho/h)^2); evaluating it
 * is cheaper in powers of (rho/h)^2, so the coefficients are converted once on reading.
 *
 * @param coefficients Coefficients in ascending powers of (1 - x).
 * @return Coefficients in ascending powers of x, same length.
 */
Eigen::VectorXd convertPolynomialToPowersOfArgument( const Eigen::VectorXd& coefficients );

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_ELLIPSOIDAL_HARMONICS_FILE_READER_H
