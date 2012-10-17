/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      121017    E. Dekens         Code created.
 *
 *    References
 *      Heiskanen, W.A., Moritz, H. Physical geodesy. Freeman, 1967.
 *
 *    Notes
 *
 */

#ifndef TUDAT_SPHERICAL_HARMONICS_GRAVITY_H
#define TUDAT_SPHERICAL_HARMONICS_GRAVITY_H

#include <Eigen/Core>

namespace tudat
{
namespace gravitation
{

//! Compute gravitational acceleration due to multiple spherical harmonics terms, defined using
//! geodesy-normalization.
/*!
 * This function computes the acceleration caused by gravitational spherical harmonics, with the
 * coefficients expressed using a geodesy-normalization. This acceleration is the summation of all
 * harmonic terms from degree and order zero, up to a user-specified highest degree and order. The
 * harmonic coefficients for the function must be provided in geodesy-normalized format. This
 * geodesy-normalization is defined as:
 * \f{eqnarray*}
 *  {
 *     \bar{ C }_{ n, m } = \Pi_{ n, m } C_{ n, m } \\
 *     \bar{ S }_{ n, m } = \Pi_{ n, m } S_{ n, m }
 * \f}
 * in which \f$ \bar{ C }_{ n, m } \f$ and \f$ \bar{ S }_{ n, m } \f$ are a geodesy-normalized
 * cosine and sine harmonic coefficient respectively (of degree \f$ n \f$ and order \f$ m \f$). The
 * unnormalized harmonic coefficients are represented by \f$ C_{ n, m } \f$ and \f$ S_{ n, m } \f$.
 * The normalization factor \f$ \Pi_{ n, m } \f$ is given by Heiskanen & Moritz [1967] as:
 * \f[
 *     \Pi_{ n, m } = \sqrt{ \frac{ ( n + m )! }{ ( 2 - \delta_{ 0, m } ) ( 2 n + 1 ) ( n - m )! } }
 * \f]
 * in which \f$ n \f$ is the degree, \f$ m \f$ is the order and \f$ \delta_{ 0, m } \f$ is the
 * Kronecker delta.
 * \param position Cartesian position vector with respect to the reference frame that is associated
 *          with the harmonic coefficients.
 *          The order is important!
 *          position( 0 ) = x coordinate [m],
 *          position( 1 ) = y coordinate [m],
 *          position( 2 ) = z coordinate [m].
 * \param highestDegree Highest harmonics degree which is to included.
 * \param highestOrder Highest harmonics order which is to be included.
 * \param cosineHarmonicCoefficients Matrix with <B>geodesy-normalized</B> cosine harmonic
 *          coefficients. The row index indicates the degree and the column index indicates the order
 *          of coefficients. The matrix must have a number of rows equal to 'highestDegree' and a
 *          number of columns equal to 'highestOrder'.
 * \param sineHarmonicCoefficients Matrix with <B>geodesy-normalized</B> sine harmonic coefficients.
 *          The row index indicates the degree and the column index indicates the order of
 *          coefficients. The matrix must have a number of rows equal to 'highestDegree' and a
 *          number of columns equal to 'highestOrder'.
 * \param gravitationalParameter Gravitational parameter associated with the spherical harmonics
 *          [m^3 s^-2].
 * \param planetaryRadius Reference radius of the spherical harmonics [m].
 * \return Cartesian acceleration vector resulting from the summation of all harmonic terms.
 *           The order is important!
 *           acceleration( 0 ) = x acceleration [m s^-2],
 *           acceleration( 1 ) = y acceleration [m s -2],
 *           acceleration( 2 ) = z acceleration [m s^-2].
 */
Eigen::Vector3d computeGeodesyNormalizedGravitationalAccelerationSum(
        const Eigen::Vector3d &position,
        const int highestDegree,
        const int highestOrder,
        const Eigen::MatrixXd& cosineHarmonicCoefficients,
        const Eigen::MatrixXd& sineHarmonicCoefficients,
        const double gravitationalParameter,
        const double planetaryRadius );

//! Compute gravitational acceleration due to single spherical harmonics term.
/*!
 * This function computes the acceleration caused by a single gravitational spherical harmonics
 * term, with the coefficients expressed using a geodesy-normalization. The harmonic coefficients
 * for the function must be provided in geodesy-normalized format. This geodesy-normalization is
 * defined as:
 * \f{eqnarray*}
 *  {
 *     \bar{ C }_{ n, m } = \Pi_{ n, m } C_{ n, m } \\
 *     \bar{ S }_{ n, m } = \Pi_{ n, m } S_{ n, m }
 * \f}
 * in which \f$ \bar{ C }_{ n, m } \f$ and \f$ \bar{ S }_{ n, m } \f$ are a geodesy-normalized
 * cosine and sine harmonic coefficient respectively (of degree \f$ n \f$ and order \f$ m \f$). The
 * unnormalized harmonic coefficients are represented by \f$ C_{ n, m } \f$ and \f$ S_{ n, m } \f$.
 * The normalization factor \f$ \Pi_{ n, m } \f$ is given by Heiskanen & Moritz [1967] as:
 * \f[
 *     \Pi_{ n, m } = \sqrt{ \frac{ ( n + m )! }{ ( 2 - \delta_{ 0, m } ) ( 2 n + 1 ) ( n - m )! } }
 * \f]
 * in which \f$ n \f$ is the degree, \f$ m \f$ is the order and \f$ \delta_{ 0, m } \f$ is the
 * Kronecker delta.
 * \param position Cartesian position vector with respect to the reference frame that is associated
 *          with the harmonic coefficients.
 *          The order is important!
 *          position( 0 ) = x coordinate [m],
 *          position( 1 ) = y coordinate [m],
 *          position( 2 ) = z coordinate [m].
 * \param degree Degree of the harmonic term.
 * \param order Order of the harmonic term.
 * \param cosineHarmonicCoefficient <B>Geodesy-normalized</B> cosine harmonic
 *          coefficient.
 * \param sineHarmonicCoefficient <B>Geodesy-normalized</B> sine harmonic coefficient.
 * \param gravitationalParameter Gravitational parameter associated with the spherical harmonic
 *          [m^3 s^-2].
 * \param planetaryRadius Reference radius of the spherical harmonic [m].
 * \return Cartesian acceleration vector resulting from the spherical harmonic term.
 *           The order is important!
 *           acceleration( 0 ) = x acceleration [m s^-2],
 *           acceleration( 1 ) = y acceleration [m s -2],
 *           acceleration( 2 ) = z acceleration [m s^-2].
 */
Eigen::Vector3d computeSingleGeodesyNormalizedGravitationalAcceleration(
        const Eigen::Vector3d &position,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const double gravitationalParameter,
        const double planetaryRadius );

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_SPHERICAL_HARMONICS_GRAVITY_H
