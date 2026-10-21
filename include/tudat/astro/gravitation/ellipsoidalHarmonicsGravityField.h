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

#ifndef TUDAT_ELLIPSOIDAL_HARMONICS_GRAVITY_FIELD_H
#define TUDAT_ELLIPSOIDAL_HARMONICS_GRAVITY_FIELD_H

/*  Gravity field of a body described by ellipsoidal harmonics, valid outside the reference
 *  ellipsoid. See tudat/math/basic/ellipsoidalHarmonics.h for the expansion itself.
 */

#include <memory>
#include <string>

#include <Eigen/Core>

#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/math/basic/ellipsoidalHarmonics.h"
#include "tudat/math/basic/legendrePolynomials.h"

namespace tudat
{

namespace gravitation
{

//! Stores the parts of an ellipsoidal harmonic evaluation that depend only on position.
/*!
 * The potential, the gradient and the Hessian all need the same intermediate values, so they
 * are computed once per position and kept here.
 */
class EllipsoidalHarmonicsGravityCache
{
public:
    //! Constructor.
    /*!
     * \param lameFunctionCoefficients Lame polynomial coefficients of every term.
     * \param firstFocalParameter Focal parameter sqrt( a^2 - b^2 ) of the reference ellipsoid.
     * \param secondFocalParameter Focal parameter sqrt( a^2 - c^2 ) of the reference ellipsoid.
     * \param quadratureOrder Number of Gauss-Legendre nodes used for the radial parts.
     */
    EllipsoidalHarmonicsGravityCache( const basic_mathematics::EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                      const double firstFocalParameter,
                                      const double secondFocalParameter,
                                      const int quadratureOrder,
                                      const std::shared_ptr< basic_mathematics::EllipsoidalHarmonicRadialTable > radialTable = nullptr );

    //! Recomputes what the potential needs, unless the position is the one already cached.
    void update( const Eigen::Vector3d& currentBodyFixedPosition );

    //! Recomputes what the gradient needs as well.
    void updateIncludingDerivatives( const Eigen::Vector3d& currentBodyFixedPosition );

    //! Recomputes what the Hessian needs as well.
    void updateIncludingSecondDerivatives( const Eigen::Vector3d& currentBodyFixedPosition );

    //! Supplies the tabulated radial parts, replacing the exact integration where it applies.
    /*!
     * Anything already cached is discarded, so the next update recomputes.
     */
    void setRadialTable( const std::shared_ptr< basic_mathematics::EllipsoidalHarmonicRadialTable > radialTable );

    //! Ellipsoidal radius of the cached position.
    double getEllipsoidalRadius( ) const
    {
        return ellipsoidalCoordinates_( 0 );
    }

    //! Ellipsoidal coordinates rho, phi and gamma of the cached position.
    const Eigen::Vector3d& getEllipsoidalCoordinates( ) const
    {
        return ellipsoidalCoordinates_;
    }

    //! Angle of the cached position about the long axis.
    double getAngleAboutLongAxis( ) const
    {
        return ellipsoidalCoordinates_( 2 );
    }

    //! Radial part of every term at the cached position.
    const basic_mathematics::EllipsoidalHarmonicCoefficients& getRadialFunctions( ) const
    {
        return radialFunctions_;
    }

    //! Derivatives of the radial parts with respect to the ellipsoidal radius.
    const basic_mathematics::EllipsoidalHarmonicCoefficients& getRadialFunctionDerivatives( ) const
    {
        return radialFunctionDerivatives_;
    }

    //! Legendre functions at the cached position, values only.
    const basic_mathematics::LegendreCache& getLegendreCache( ) const
    {
        return legendreCache_;
    }

    //! Legendre functions at the cached position, with derivatives.
    const basic_mathematics::LegendreCache& getLegendreCacheWithDerivatives( ) const
    {
        return legendreCacheWithDerivatives_;
    }

    //! Second derivatives of the radial parts with respect to the ellipsoidal radius.
    const basic_mathematics::EllipsoidalHarmonicCoefficients& getRadialFunctionSecondDerivatives( ) const
    {
        return radialFunctionSecondDerivatives_;
    }

    //! Legendre functions at the cached position, with first and second derivatives.
    const basic_mathematics::LegendreCache& getLegendreCacheWithSecondDerivatives( ) const
    {
        return legendreCacheWithSecondDerivatives_;
    }

private:
    //! Radial parts and their first derivatives at one radius, from the table where it applies and by exact integration otherwise.
    void evaluateRadialFunctions( const double rho,
                                  basic_mathematics::EllipsoidalHarmonicCoefficients& values,
                                  basic_mathematics::EllipsoidalHarmonicCoefficients& derivatives ) const;

    //! Position the cached values belong to. Starts as NaN so the first update always runs.
    Eigen::Vector3d currentBodyFixedPosition_;

    //! Position the cached derivatives belong to, tracked separately from the values.
    Eigen::Vector3d positionOfCachedDerivatives_;

    //! Position the cached second derivatives belong to.
    Eigen::Vector3d positionOfCachedSecondDerivatives_;

    Eigen::Vector3d ellipsoidalCoordinates_;

    basic_mathematics::EllipsoidalHarmonicMatrices lameFunctionCoefficients_;

    basic_mathematics::EllipsoidalHarmonicCoefficients radialFunctions_;

    basic_mathematics::EllipsoidalHarmonicCoefficients radialFunctionDerivatives_;

    basic_mathematics::LegendreCache legendreCache_;

    basic_mathematics::EllipsoidalHarmonicCoefficients radialFunctionSecondDerivatives_;

    //! Scratch space for differencing the first derivatives to get the second ones.
    basic_mathematics::EllipsoidalHarmonicCoefficients forwardRadialFunctions_;

    basic_mathematics::EllipsoidalHarmonicCoefficients forwardRadialFunctionDerivatives_;

    basic_mathematics::EllipsoidalHarmonicCoefficients backwardRadialFunctions_;

    basic_mathematics::EllipsoidalHarmonicCoefficients backwardRadialFunctionDerivatives_;

    basic_mathematics::LegendreCache legendreCacheWithDerivatives_;

    basic_mathematics::LegendreCache legendreCacheWithSecondDerivatives_;

    //! Tabulated radial parts, used in place of the integration where it covers the radius.
    std::shared_ptr< basic_mathematics::EllipsoidalHarmonicRadialTable > radialTable_;

    const double firstFocalParameter_;

    const double secondFocalParameter_;

    Eigen::VectorXd quadratureNodes_;

    Eigen::VectorXd quadratureWeights_;
};

//! Gravity field of a body described by ellipsoidal harmonics.
class EllipsoidalHarmonicsGravityField : public GravityFieldModel
{
public:
    //! Constructor.
    /*!
     * The gravitational parameter is taken from the degree zero coefficient.
     *
     * \param referenceEllipsoidSemiAxes The three semi-axes of the reference ellipsoid, longest
     *      first and shortest last. The expansion is valid outside this ellipsoid.
     * \param coefficients Field coefficients as read from file, still normalised to the surface
     *      of the reference ellipsoid.
     * \param lameFunctionCoefficients Lame polynomial coefficients of every term.
     * \param transformFactors Factors that mix the terms within each group.
     * \param quadratureOrder Number of Gauss-Legendre nodes used for the radial parts.
     * \param radialTableGridPoints Number of radii at which to tabulate the radial parts. Zero
     *      integrates them exactly every time instead.
     * \param radialTableLargestRadius Largest ellipsoidal radius the table covers. Beyond it the
     *      radial parts are integrated exactly. Zero uses fifty times the longest semi-axis.
     * \param fixedReferenceFrame Identifier of the body-fixed frame the field is given in.
     * \param updateInertiaTensor Function called to update the inertia tensor, usually supplied
     *      by the Body class.
     */
    EllipsoidalHarmonicsGravityField( const Eigen::Vector3d& referenceEllipsoidSemiAxes,
                                      const basic_mathematics::EllipsoidalHarmonicCoefficients& coefficients,
                                      const basic_mathematics::EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                      const basic_mathematics::EllipsoidalHarmonicMatrices& transformFactors,
                                      const int quadratureOrder = 24,
                                      const int radialTableGridPoints = 3200,
                                      const double radialTableLargestRadius = 0.0,
                                      const std::string& fixedReferenceFrame = "",
                                      const std::function< void( ) > updateInertiaTensor = std::function< void( ) >( ) );

    //! Gravitational potential at a body-fixed position.
    /*!
     * Throws if the position is inside the reference ellipsoid, where the expansion does not
     * converge.
     */
    double getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition ) override;

    //! Gradient of the potential at a body-fixed position, which is the gravitational acceleration.
    /*!
     * Throws if the position is inside the reference ellipsoid, where the expansion does not
     * converge.
     *
     * Close to the long axis of the body, where the analytic gradient is singular, the gradient
     * is obtained by differencing the potential instead.
     */
    Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& bodyFixedPosition ) override;

    //! Second derivatives of the potential at a body-fixed position.
    /*!
     * Close to the long axis of the body the gradient is differenced instead.
     */
    Eigen::Matrix3d getHessianOfPotential( const Eigen::Vector3d& bodyFixedPosition );

    //! Laplacian of the potential, which is zero outside the body.
    /*!
     * Computed as the trace of the Hessian.
     */
    double getLaplacianOfPotential( const Eigen::Vector3d& bodyFixedPosition ) override;

    //! Identifier of the body-fixed reference frame the field is given in.
    std::string getFixedReferenceFrame( )
    {
        return fixedReferenceFrame_;
    }

private:
    //! Factor that brings the coefficients to the current gravitational parameter.
    /*!
     * Scaling by this factor makes a reset of the gravitational parameter change the field.
     */
    double getCoefficientScaling( ) const
    {
        return gravitationalParameter_ / gravitationalParameterOfCoefficients_;
    }

    //! Gradient obtained by differencing the potential, used close to the long axis.
    Eigen::Vector3d computeGradientByDifferencingPotential( const Eigen::Vector3d& bodyFixedPosition );

    //! Hessian obtained by differencing the gradient, used close to the long axis.
    Eigen::Matrix3d computeHessianByDifferencingGradient( const Eigen::Vector3d& bodyFixedPosition );

    //! Throws if the position lies inside the reference ellipsoid.
    void checkPositionIsOutsideReferenceEllipsoid( const Eigen::Vector3d& bodyFixedPosition, const double ellipsoidalRadius ) const;

    const Eigen::Vector3d referenceEllipsoidSemiAxes_;

    const double firstFocalParameter_;

    const double secondFocalParameter_;

    const int maximumDegree_;

    const int quadratureOrder_;

    //! Coefficients with the normalisation to the reference surface undone.
    basic_mathematics::EllipsoidalHarmonicCoefficients coefficients_;

    //! Gravitational parameter contained in the coefficients, as read from file.
    double gravitationalParameterOfCoefficients_;

    basic_mathematics::EllipsoidalHarmonicMatrices lameFunctionCoefficients_;

    basic_mathematics::EllipsoidalHarmonicMatrices transformFactors_;

    std::shared_ptr< EllipsoidalHarmonicsGravityCache > gravityCache_;

    //! Cosine of the angle from the long axis within which the potential is differenced instead.
    const double longAxisToleranceForGradient_ = 1.0E-4;

    std::string fixedReferenceFrame_;
};

}  // namespace gravitation

}  // namespace tudat

#endif  // TUDAT_ELLIPSOIDAL_HARMONICS_GRAVITY_FIELD_H
