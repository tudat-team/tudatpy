/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ELLIPSOIDAL_HARMONICS_GRAVITY_FIELD_H
#define TUDAT_ELLIPSOIDAL_HARMONICS_GRAVITY_FIELD_H

/*  Gravity field of a body described by ellipsoidal harmonics. Unlike spherical harmonics, which
 *  are only valid outside the sphere enclosing the body, this expansion is valid outside the
 *  enclosing ellipsoid, so it stays usable much closer to an elongated body. See
 *  tudat/math/basic/ellipsoidalHarmonics.h for the expansion itself.
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

/*! Stores the parts of an ellipsoidal harmonic evaluation that depend only on position.
 *
 * Converting to ellipsoidal coordinates and evaluating the radial part of every term is the
 * expensive step, and the potential, the acceleration and their derivatives all need the same
 * values. They are therefore computed once per position and kept here, as the polyhedron and
 * ring gravity fields do with their own intermediates.
 */
class EllipsoidalHarmonicsGravityCache
{
public:
    /*! Constructor.
     *
     * @param lameFunctionCoefficients Lame polynomial coefficients of every term.
     * @param firstFocalParameter Focal parameter sqrt( a^2 - b^2 ) of the reference ellipsoid.
     * @param secondFocalParameter Focal parameter sqrt( a^2 - c^2 ) of the reference ellipsoid.
     * @param quadratureOrder Number of Gauss-Legendre nodes used for the radial parts.
     * @param computeLegendreDerivatives Whether the Legendre derivatives are needed, which they
     *      are for the acceleration but not for the potential. They are singular on the long
     *      axis of the body, so asking for them restricts where the field can be evaluated.
     */
    EllipsoidalHarmonicsGravityCache( const basic_mathematics::EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                      const double firstFocalParameter,
                                      const double secondFocalParameter,
                                      const int quadratureOrder,
                                      const bool computeLegendreDerivatives );

    //! Recomputes the cached values, unless the position is the one already cached.
    void update( const Eigen::Vector3d& currentBodyFixedPosition );

    //! Ellipsoidal radius of the cached position.
    double getEllipsoidalRadius( ) const
    {
        return ellipsoidalCoordinates_( 0 );
    }

    //! Angle of the cached position measured away from the plane normal to the long axis.
    double getAngleFromLongAxis( ) const
    {
        return ellipsoidalCoordinates_( 1 );
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

    //! Legendre functions at the cached position.
    const basic_mathematics::LegendreCache& getLegendreCache( ) const
    {
        return legendreCache_;
    }

private:
    //! Position the cached values belong to; starts as NaN so the first update always runs.
    Eigen::Vector3d currentBodyFixedPosition_;

    Eigen::Vector3d ellipsoidalCoordinates_;

    basic_mathematics::EllipsoidalHarmonicMatrices lameFunctionCoefficients_;

    basic_mathematics::EllipsoidalHarmonicCoefficients radialFunctions_;

    basic_mathematics::LegendreCache legendreCache_;

    const double firstFocalParameter_;

    const double secondFocalParameter_;

    Eigen::VectorXd quadratureNodes_;

    Eigen::VectorXd quadratureWeights_;
};

//! Gravity field of a body described by ellipsoidal harmonics.
class EllipsoidalHarmonicsGravityField : public GravityFieldModel
{
public:
    /*! Constructor.
     *
     * @param gravitationalParameter Gravitational parameter of the body.
     * @param referenceEllipsoidSemiAxes The three semi-axes of the reference ellipsoid, longest
     *      first and shortest last. The expansion is valid outside this ellipsoid.
     * @param coefficients Field coefficients as read from file, still normalised to the surface
     *      of the reference ellipsoid; the constructor undoes that normalisation.
     * @param lameFunctionCoefficients Lame polynomial coefficients of every term.
     * @param transformFactors Factors that mix the terms within each group.
     * @param quadratureOrder Number of Gauss-Legendre nodes used for the radial parts. The
     *      default is accurate to well below the precision of typical input coefficients.
     * @param fixedReferenceFrame Identifier of the body-fixed frame the field is given in.
     * @param updateInertiaTensor Function called to update the inertia tensor, usually supplied
     *      by the Body class.
     */
    EllipsoidalHarmonicsGravityField( const double gravitationalParameter,
                                      const Eigen::Vector3d& referenceEllipsoidSemiAxes,
                                      const basic_mathematics::EllipsoidalHarmonicCoefficients& coefficients,
                                      const basic_mathematics::EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                      const basic_mathematics::EllipsoidalHarmonicMatrices& transformFactors,
                                      const int quadratureOrder = 24,
                                      const std::string& fixedReferenceFrame = "",
                                      const std::function< void( ) > updateInertiaTensor = std::function< void( ) >( ) );

    /*! Gravitational potential at a body-fixed position.
     *
     * Throws if the position is inside the reference ellipsoid, where the expansion does not
     * converge.
     */
    double getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition ) override;

    /*! Not implemented yet; throws.
     *
     * The base class would otherwise silently return the gradient of a point mass, which would
     * look plausible and be wrong.
     */
    Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& bodyFixedPosition ) override;

    //! Identifier of the body-fixed reference frame the field is given in.
    std::string getFixedReferenceFrame( )
    {
        return fixedReferenceFrame_;
    }

    //! The three semi-axes of the reference ellipsoid, longest first.
    Eigen::Vector3d getReferenceEllipsoidSemiAxes( )
    {
        return referenceEllipsoidSemiAxes_;
    }

    //! Focal parameter sqrt( a^2 - b^2 ) of the reference ellipsoid.
    double getFirstFocalParameter( )
    {
        return firstFocalParameter_;
    }

    //! Focal parameter sqrt( a^2 - c^2 ) of the reference ellipsoid.
    double getSecondFocalParameter( )
    {
        return secondFocalParameter_;
    }

    //! Highest degree of the expansion.
    int getMaximumDegree( )
    {
        return maximumDegree_;
    }

    //! Number of Gauss-Legendre nodes used for the radial parts.
    int getQuadratureOrder( )
    {
        return quadratureOrder_;
    }

    //! Field coefficients, with the normalisation to the reference surface already undone.
    basic_mathematics::EllipsoidalHarmonicCoefficients getCoefficients( )
    {
        return coefficients_;
    }

    //! Lame polynomial coefficients of every term.
    basic_mathematics::EllipsoidalHarmonicMatrices getLameFunctionCoefficients( )
    {
        return lameFunctionCoefficients_;
    }

    //! Factors that mix the terms within each group.
    basic_mathematics::EllipsoidalHarmonicMatrices getTransformFactors( )
    {
        return transformFactors_;
    }

private:
    //! Throws if the position lies inside the reference ellipsoid.
    void checkPositionIsOutsideReferenceEllipsoid( const Eigen::Vector3d& bodyFixedPosition, const double ellipsoidalRadius ) const;

    const Eigen::Vector3d referenceEllipsoidSemiAxes_;

    const double firstFocalParameter_;

    const double secondFocalParameter_;

    const int maximumDegree_;

    const int quadratureOrder_;

    //! Coefficients with the normalisation to the reference surface undone.
    basic_mathematics::EllipsoidalHarmonicCoefficients coefficients_;

    basic_mathematics::EllipsoidalHarmonicMatrices lameFunctionCoefficients_;

    basic_mathematics::EllipsoidalHarmonicMatrices transformFactors_;

    std::shared_ptr< EllipsoidalHarmonicsGravityCache > gravityCache_;

    std::string fixedReferenceFrame_;
};

}  // namespace gravitation

}  // namespace tudat

#endif  // TUDAT_ELLIPSOIDAL_HARMONICS_GRAVITY_FIELD_H
