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
 *          Garmier, R., Barriot, J.-P. (2001), "Ellipsoidal harmonic expansions of the gravitational
 *          potential: theory and application", Celestial Mechanics and Dynamical Astronomy, 79(4),
 *          235-275.
 *
 *          Hu, X., Jekeli, C. (2015), "A numerical comparison of spherical, spheroidal and
 *          ellipsoidal harmonic gravitational field models for small non-spherical bodies: examples
 *          for the Martian moons", Journal of Geodesy, 89(2), 159-177.
 */

#ifndef TUDAT_ELLIPSOIDAL_HARMONICS_H
#define TUDAT_ELLIPSOIDAL_HARMONICS_H

/*  Ellipsoidal harmonics: an expansion of a gravitational potential that is valid outside the
 *  smallest ellipsoid enclosing the body, instead of outside the smallest enclosing sphere as
 *  spherical harmonics are. For elongated bodies the enclosing ellipsoid is much tighter, so the
 *  expansion stays valid far closer to the surface.
 *
 *  Position is described by ellipsoidal coordinates (see convertCartesianToEllipsoidalCoordinates)
 *  and the radial part of each term is a Lame function instead of a power of 1/r. Each degree n
 *  has 2n+1 terms, split into four groups labelled K, L, M and N. The coefficients, the Lame
 *  polynomials and the transform factors are read from files (see
 *  tudat/io/ellipsoidalHarmonicsFileReader.h). The Legendre functions are the geodesy-normalised
 *  ones from LegendreCache.
 *
 *  The expansion should be truncated at degree 9, because the transform factors in the input
 *  files lose accuracy above that degree.
 */

#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/math/basic/legendrePolynomials.h"

namespace tudat
{

namespace basic_mathematics
{

//! The four groups into which the 2n+1 terms of one degree are divided.
enum EllipsoidalHarmonicTermGroup {
    ellipsoidal_harmonic_k_terms = 0,
    ellipsoidal_harmonic_l_terms = 1,
    ellipsoidal_harmonic_m_terms = 2,
    ellipsoidal_harmonic_n_terms = 3
};

//! Number of terms in each group at the given degree. The four counts add up to 2*degree+1.
inline int getNumberOfKTerms( const int degree )
{
    return degree / 2 + 1;
}
inline int getNumberOfLTerms( const int degree )
{
    return degree - degree / 2;
}
inline int getNumberOfMTerms( const int degree )
{
    return degree - degree / 2;
}
inline int getNumberOfNTerms( const int degree )
{
    return degree / 2;
}

//! Number of terms in one group at the given degree.
int getNumberOfTermsInGroup( const int degree, const EllipsoidalHarmonicTermGroup group );

//! Maps a file order index to a term group and a position within that group.
/*!
 * The coefficient files label the 2*degree+1 terms of a degree with an order index running from
 * 1 to 2*degree+1. This splits such an index into the group it belongs to and its position
 * inside that group.
 *
 * \param degree Degree of the term.
 * \param fileOrderIndex Order index as it appears in the file, from 1 to 2*degree+1.
 * \param group Returned term group.
 * \param indexInGroup Returned zero-based position within the group.
 */
void getTermGroupAndIndex( const int degree, const int fileOrderIndex, EllipsoidalHarmonicTermGroup& group, int& indexInGroup );

//! One entry per term, grouped by degree and by term group.
/*!
 * Each of the four members holds one entry per degree, so kTerms.at( n ) describes the K terms
 * of degree n. What a single entry contains depends on EntryType (see the two typedefs below).
 */
template< typename EntryType >
struct EllipsoidalHarmonicTermSet {
    std::vector< EntryType > kTerms;
    std::vector< EntryType > lTerms;
    std::vector< EntryType > mTerms;
    std::vector< EntryType > nTerms;

    //! Highest degree stored.
    int getMaximumDegree( ) const
    {
        return static_cast< int >( kTerms.size( ) ) - 1;
    }

    //! Entries of one group at one degree.
    const EntryType& getGroup( const int degree, const EllipsoidalHarmonicTermGroup group ) const
    {
        switch( group )
        {
            case ellipsoidal_harmonic_k_terms:
                return kTerms.at( degree );
            case ellipsoidal_harmonic_l_terms:
                return lTerms.at( degree );
            case ellipsoidal_harmonic_m_terms:
                return mTerms.at( degree );
            default:
                return nTerms.at( degree );
        }
    }

    EntryType& getGroup( const int degree, const EllipsoidalHarmonicTermGroup group )
    {
        switch( group )
        {
            case ellipsoidal_harmonic_k_terms:
                return kTerms.at( degree );
            case ellipsoidal_harmonic_l_terms:
                return lTerms.at( degree );
            case ellipsoidal_harmonic_m_terms:
                return mTerms.at( degree );
            default:
                return nTerms.at( degree );
        }
    }
};

//! One number per term, used for the field coefficients and for the radial function values.
typedef EllipsoidalHarmonicTermSet< Eigen::VectorXd > EllipsoidalHarmonicCoefficients;

//! One matrix per group per degree. Used for the Lame polynomial coefficients, where column j
//! holds the polynomial of term j in ascending powers, and for the square transform matrices.
typedef EllipsoidalHarmonicTermSet< Eigen::MatrixXd > EllipsoidalHarmonicMatrices;

//! Allocates a coefficient set of the given maximum degree, with all entries set to zero.
EllipsoidalHarmonicCoefficients createEllipsoidalHarmonicCoefficients( const int maximumDegree );

//! Converts a Cartesian position to the ellipsoidal coordinates used by this expansion.
/*!
 * The returned coordinates are the ellipsoidal radius rho, the angle phi measured away from the
 * plane normal to the longest axis, and the angle gamma about that axis. They relate to
 * Cartesian position as
 *
 *      x = rho * sin( phi )
 *      y = sqrt( rho^2 - h^2 ) * cos( phi ) * cos( gamma )
 *      z = sqrt( rho^2 - k^2 ) * cos( phi ) * sin( gamma )
 *
 * so surfaces of constant rho are ellipsoids with semi-axes rho, sqrt(rho^2-h^2) and
 * sqrt(rho^2-k^2). Note that the polar axis is x, the longest axis of the body, not z.
 *
 * \param bodyFixedPosition Position in the body-fixed frame.
 * \param h Focal parameter sqrt( a^2 - b^2 ) of the reference ellipsoid.
 * \param k Focal parameter sqrt( a^2 - c^2 ) of the reference ellipsoid, larger than h.
 * \return Vector holding rho, phi and gamma, with the angles in radians.
 */
Eigen::Vector3d convertCartesianToEllipsoidalCoordinates( const Eigen::Vector3d& bodyFixedPosition, const double h, const double k );

//! Nodes and weights of the Gauss-Legendre quadrature rule of the given order on [-1, 1].
/*!
 * \param order Number of nodes.
 * \param nodes Returned nodes, in ascending order.
 * \param weights Returned weights.
 */
void computeGaussLegendreNodesAndWeights( const int order, Eigen::VectorXd& nodes, Eigen::VectorXd& weights );

//! Lame function of the first kind, the radial part of a term inside the reference ellipsoid.
/*!
 * Only the term group matters here, not which term within the group. The group selects which
 * algebraic factors multiply the polynomial.
 *
 * \param polynomialCoefficients Polynomial of this term, in ascending powers of (rho/h)^2.
 * \param degree Degree of the term.
 * \param group Group the term belongs to.
 * \param rho Ellipsoidal radius at which to evaluate.
 * \param h First focal parameter of the reference ellipsoid.
 * \param k Second focal parameter of the reference ellipsoid.
 */
double evaluateLameFunctionFirstKind( const Eigen::VectorXd& polynomialCoefficients,
                                      const int degree,
                                      const EllipsoidalHarmonicTermGroup group,
                                      const double rho,
                                      const double h,
                                      const double k );

//! Radial part of every term at one ellipsoidal radius.
/*!
 * Fills radialFunctions with the second-kind Lame function of each term. These depend on rho
 * alone, which is why they are computed in one pass and shared by all terms of the expansion.
 *
 * \param lameFunctionCoefficients Lame polynomial coefficients of every term.
 * \param rho Ellipsoidal radius at which to evaluate. Must be larger than k.
 * \param h First focal parameter of the reference ellipsoid.
 * \param k Second focal parameter of the reference ellipsoid.
 * \param quadratureNodes Gauss-Legendre nodes on [-1, 1].
 * \param quadratureWeights Matching Gauss-Legendre weights.
 * \param radialFunctions Output, allocated by createEllipsoidalHarmonicCoefficients.
 */
void computeRadialFunctions( const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                             const double rho,
                             const double h,
                             const double k,
                             const Eigen::VectorXd& quadratureNodes,
                             const Eigen::VectorXd& quadratureWeights,
                             EllipsoidalHarmonicCoefficients& radialFunctions );

//! Legendre order that a term pairs with.
/*!
 * The four groups walk through the Legendre orders differently: K takes the even orders from
 * zero, L and M the odd orders, and N the even orders from two.
 *
 * \param group Group the term belongs to.
 * \param indexInGroup Zero-based position of the term within its group.
 */
int getTermLegendreOrder( const EllipsoidalHarmonicTermGroup group, const int indexInGroup );

//! Whether a term pairs with the sine of its angle rather than the cosine.
bool termUsesSineOfAngle( const EllipsoidalHarmonicTermGroup group );

//! Derivative of the first-kind Lame function with respect to the ellipsoidal radius.
/*!
 * Arguments are as for evaluateLameFunctionFirstKind.
 */
double evaluateLameFunctionFirstKindDerivative( const Eigen::VectorXd& polynomialCoefficients,
                                                const int degree,
                                                const EllipsoidalHarmonicTermGroup group,
                                                const double rho,
                                                const double h,
                                                const double k );

//! Radial part of every term and its derivative, at one ellipsoidal radius.
/*!
 * Arguments are as for computeRadialFunctions, with the derivatives written to a second output.
 *
 * \param radialFunctions Output, allocated by createEllipsoidalHarmonicCoefficients.
 * \param radialFunctionDerivatives Output for the derivatives, allocated the same way.
 */
void computeRadialFunctionsAndDerivatives( const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                           const double rho,
                                           const double h,
                                           const double k,
                                           const Eigen::VectorXd& quadratureNodes,
                                           const Eigen::VectorXd& quadratureWeights,
                                           EllipsoidalHarmonicCoefficients& radialFunctions,
                                           EllipsoidalHarmonicCoefficients& radialFunctionDerivatives );

//! Derivatives of Cartesian position with respect to the ellipsoidal coordinates.
/*!
 * Column j holds the derivative of position with respect to the j-th ellipsoidal coordinate.
 * The matrix is singular where the long axis of the body meets the surface of constant
 * ellipsoidal radius, which is the one place the gradient cannot be formed.
 *
 * \param ellipsoidalCoordinates Coordinates rho, phi and gamma of the point.
 * \param h First focal parameter of the reference ellipsoid.
 * \param k Second focal parameter of the reference ellipsoid.
 */
Eigen::Matrix3d computeEllipsoidalCoordinateJacobian( const Eigen::Vector3d& ellipsoidalCoordinates, const double h, const double k );

//! Gradient of the potential in the body-fixed frame, which is the gravitational acceleration.
/*!
 * \param coefficients Coefficients after normaliseCoefficientsToReferenceSurface.
 * \param radialFunctions Radial parts at this point.
 * \param radialFunctionDerivatives Their derivatives with respect to the ellipsoidal radius.
 * \param transformFactors Transform matrices that mix the terms within each group.
 * \param legendreCache Legendre cache updated with sin( phi ), with derivatives switched on.
 * \param ellipsoidalCoordinates Coordinates rho, phi and gamma of the point.
 * \param h First focal parameter of the reference ellipsoid.
 * \param k Second focal parameter of the reference ellipsoid.
 * \return Gradient of the potential, in the body-fixed frame.
 */
Eigen::Vector3d computeEllipsoidalHarmonicPotentialGradient( const EllipsoidalHarmonicCoefficients& coefficients,
                                                             const EllipsoidalHarmonicCoefficients& radialFunctions,
                                                             const EllipsoidalHarmonicCoefficients& radialFunctionDerivatives,
                                                             const EllipsoidalHarmonicMatrices& transformFactors,
                                                             const LegendreCache& legendreCache,
                                                             const Eigen::Vector3d& ellipsoidalCoordinates,
                                                             const double h,
                                                             const double k );

//! Radial parts with their first and second derivatives with respect to the ellipsoidal radius.
/*!
 * The second derivative is obtained by numerically differencing the analytic first derivative.
 *
 * Arguments are as for computeRadialFunctionsAndDerivatives, with a third output added.
 */
void computeRadialFunctionsAndSecondDerivatives( const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                                 const double rho,
                                                 const double h,
                                                 const double k,
                                                 const Eigen::VectorXd& quadratureNodes,
                                                 const Eigen::VectorXd& quadratureWeights,
                                                 EllipsoidalHarmonicCoefficients& radialFunctions,
                                                 EllipsoidalHarmonicCoefficients& radialFunctionDerivatives,
                                                 EllipsoidalHarmonicCoefficients& radialFunctionSecondDerivatives );

//! Second derivatives of Cartesian position with respect to the ellipsoidal coordinates.
/*!
 * Entry m of the output is the symmetric matrix of second derivatives of the m-th Cartesian
 * component. Together with the Jacobian these give the Hessian of anything expressed in
 * ellipsoidal coordinates, without needing the second derivatives of the inverse map.
 *
 * \param ellipsoidalCoordinates Coordinates rho, phi and gamma of the point.
 * \param h First focal parameter of the reference ellipsoid.
 * \param k Second focal parameter of the reference ellipsoid.
 * \param secondDerivatives Output, resized to three matrices.
 */
void computeEllipsoidalCoordinateSecondDerivatives( const Eigen::Vector3d& ellipsoidalCoordinates,
                                                    const double h,
                                                    const double k,
                                                    std::vector< Eigen::Matrix3d >& secondDerivatives );

//! Hessian of the potential in the body-fixed frame.
/*!
 * This is what the partial of the acceleration with respect to position needs.
 *
 * \param coefficients Coefficients after normaliseCoefficientsToReferenceSurface.
 * \param radialFunctions Radial parts at this point.
 * \param radialFunctionDerivatives Their first derivatives with respect to the ellipsoidal radius.
 * \param radialFunctionSecondDerivatives Their second derivatives.
 * \param transformFactors Transform matrices that mix the terms within each group.
 * \param legendreCache Legendre cache updated with sin( phi ), with both derivatives switched on.
 * \param ellipsoidalCoordinates Coordinates rho, phi and gamma of the point.
 * \param h First focal parameter of the reference ellipsoid.
 * \param k Second focal parameter of the reference ellipsoid.
 * \return Hessian of the potential, in the body-fixed frame.
 */
Eigen::Matrix3d computeEllipsoidalHarmonicPotentialHessian( const EllipsoidalHarmonicCoefficients& coefficients,
                                                            const EllipsoidalHarmonicCoefficients& radialFunctions,
                                                            const EllipsoidalHarmonicCoefficients& radialFunctionDerivatives,
                                                            const EllipsoidalHarmonicCoefficients& radialFunctionSecondDerivatives,
                                                            const EllipsoidalHarmonicMatrices& transformFactors,
                                                            const LegendreCache& legendreCache,
                                                            const Eigen::Vector3d& ellipsoidalCoordinates,
                                                            const double h,
                                                            const double k );

//! Total number of terms up to and including the given degree, which is ( degree + 1 )^2.
int getTotalNumberOfTerms( const int maximumDegree );

//! The radial parts tabulated against the ellipsoidal radius, for fast repeated evaluation.
/*!
 * The radial parts depend on the ellipsoidal radius only, so they are computed once on a grid of
 * radii and interpolated afterwards with cubic Hermite interpolation in the reciprocal radius.
 * This avoids one numerical integration per term at every evaluation. Outside the tabulated
 * range the caller should use the exact routines.
 */
class EllipsoidalHarmonicRadialTable
{
public:
    //! Constructor, which builds the table.
    /*!
     * \param lameFunctionCoefficients Lame polynomial coefficients of every term.
     * \param smallestRadius Smallest ellipsoidal radius to cover, normally the longest semi-axis
     *      of the reference ellipsoid, since the expansion is not valid inside it.
     * \param largestRadius Largest ellipsoidal radius to cover.
     * \param h First focal parameter of the reference ellipsoid.
     * \param k Second focal parameter of the reference ellipsoid.
     * \param quadratureOrder Number of Gauss-Legendre nodes used when building the table.
     * \param numberOfGridPoints Number of radii in the table.
     */
    EllipsoidalHarmonicRadialTable( const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                    const double smallestRadius,
                                    const double largestRadius,
                                    const double h,
                                    const double k,
                                    const int quadratureOrder,
                                    const int numberOfGridPoints );

    //! Whether the given ellipsoidal radius lies within the tabulated range.
    bool isRadiusInRange( const double rho ) const
    {
        return rho >= smallestRadius_ && rho <= largestRadius_;
    }

    //! Radial part of every term at the given radius, by interpolation.
    /*!
     * \param rho Ellipsoidal radius, which must lie within the tabulated range.
     * \param radialFunctions Output, allocated by createEllipsoidalHarmonicCoefficients.
     */
    void evaluate( const double rho, EllipsoidalHarmonicCoefficients& radialFunctions ) const;

    //! Radial part of every term and its derivative, by interpolation.
    /*!
     * \param rho Ellipsoidal radius, which must lie within the tabulated range.
     * \param radialFunctions Output for the values.
     * \param radialFunctionDerivatives Output for the derivatives with respect to the radius.
     */
    void evaluateWithDerivatives( const double rho,
                                  EllipsoidalHarmonicCoefficients& radialFunctions,
                                  EllipsoidalHarmonicCoefficients& radialFunctionDerivatives ) const;

private:
    //! Interpolates the scaled radial parts and their derivatives with respect to the interpolation variable, at the given radius.
    void interpolate( const double rho, Eigen::VectorXd& scaledValues, Eigen::VectorXd& scaledSlopes ) const;

    const int maximumDegree_;

    const double smallestRadius_;

    const double largestRadius_;

    const int numberOfGridPoints_;

    //! Value of the interpolation variable at the first grid point, and the spacing between them.
    double firstVariable_;

    double variableSpacing_;

    //! Degree of each term, in the flat order the table stores them in.
    Eigen::VectorXi degreeOfTerm_;

    //! Entry ( i, j ) is the scaled radial part of term j at grid point i.
    Eigen::MatrixXd scaledValues_;

    //! Entry ( i, j ) is its derivative with respect to the interpolation variable.
    Eigen::MatrixXd scaledSlopes_;
};

//! Divides each coefficient by the radial part of its own term at the reference ellipsoid.
/*!
 * The coefficient files are normalised to the surface of the reference ellipsoid, so this has
 * to be undone once before the coefficients can be used.
 *
 * \param fileCoefficients Coefficients as read from file.
 * \param radialFunctionsAtReferenceSurface Radial parts evaluated at rho equal to the longest
 *      semi-axis of the reference ellipsoid.
 */
EllipsoidalHarmonicCoefficients normaliseCoefficientsToReferenceSurface(
        const EllipsoidalHarmonicCoefficients& fileCoefficients,
        const EllipsoidalHarmonicCoefficients& radialFunctionsAtReferenceSurface );

//! Sums the expansion at one point, given its radial parts and angles.
/*!
 * \param coefficients Coefficients after normaliseCoefficientsToReferenceSurface.
 * \param radialFunctions Radial parts at this point, from computeRadialFunctions.
 * \param transformFactors Transform matrices that mix the terms within each group.
 * \param legendreCache Legendre cache already updated with sin( phi ). Only values are read,
 *      so the cache may have derivative computation switched off.
 * \param gamma Angle gamma of this point, in radians.
 * \return Gravitational potential.
 */
double sumEllipsoidalHarmonicTerms( const EllipsoidalHarmonicCoefficients& coefficients,
                                    const EllipsoidalHarmonicCoefficients& radialFunctions,
                                    const EllipsoidalHarmonicMatrices& transformFactors,
                                    const LegendreCache& legendreCache,
                                    const double gamma );

}  // namespace basic_mathematics

}  // namespace tudat

#endif  // TUDAT_ELLIPSOIDAL_HARMONICS_H
