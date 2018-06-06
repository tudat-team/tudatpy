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

#ifndef TUDAT_TORUS_H
#define TUDAT_TORUS_H

#include <memory>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Mathematics/GeometricShapes/singleSurfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

//! Torus class.
/*!
 * Class that defines the torus shape. The parameters are the majorRadius_,
 * denoting the distance from the center of the torus to the center of the tube
 * and the tubeRadius_, denoting the radius of the tube. Independent variables
 * are the major and minor circumferential angles, with the former denoting
 * the angle by which the tube is 'revolved' and the latter denoting the point
 * on the circular cross-section of the tube.
 */
class Torus : public SingleSurfaceGeometry
{
public:

    //! Default constructor.
    /*!
     *  Default constructor. Default values set to full torus.
     * \param majorRadius Major radius of torus (i.e. radius from center of torus to outside, of
     *          torus).
     * \param minorRadius Minor radius of torus (i.e. radius of torus 'tube')
     * \param minimumMajorCircumferentialAngle Minimum of circumferential angle 'subtended by
     *          major radius'.
     * \param maximumMajorCircumferentialAngle Maximum of circumferential angle 'subtended by
     *          major radius'.
     * \param minimumMinorCircumferentialAngle Minimum of circumferential angle 'subtended by
     *          minor radius'.
     * \param maximumMinorCircumferentialAngle Maximum of circumferential angle 'subtended by
     *          minor radius'.
     */
    Torus( const double majorRadius, const double minorRadius,
           const double minimumMajorCircumferentialAngle = 0.0,
           const double maximumMajorCircumferentialAngle = 2.0 * mathematical_constants::PI,
           const double minimumMinorCircumferentialAngle = 0.0,
           const double maximumMinorCircumferentialAngle = 2.0 * mathematical_constants::PI );

    //! Get surface point on torus.
    /*!
     * Retrieves a surface point in Cartesian coordinates on the torus from
     * values of the two independent variables.
     * Function uses cartesianPositionVector_ member variable. Values of
     * this vector set previous to function call are irrelevant.
     * \param majorCircumferentialAngle Major circumferential angle of the
     *         torus at which to retrieve the surface point.
     * \param minorCircumferentialAngle Minor circumferential angle of the
     *         torus at which to retrieve the surface point.
     * \return Point on torus in Cartesian coordinates.
     */
    Eigen::VectorXd getSurfacePoint( const double majorCircumferentialAngle,
                                     const double minorCircumferentialAngle );

    //! Get surface derivative on torus.
    /*!
     * Retrieves the derivatives of the surface point with respect to the two
     * independent variables. For powerOfMajorCircumferentialAngle = 1 and
     * powerOfMinorCircumferentialAngle = 2 the function returns:
     * \f[
     *      \frac{ d^{ 3 } ( x, y, z ) } { du * dv^{ 2 } }
     * \f]
     * (with u and v the major and minor angles).
     * \param majorCircumferentialAngle Major circumferential angle.
     * \param minorCircumferentialAngle Minor circumferential angle.
     * \param powerOfMajorCircumferentialAngleDerivative Power of derivative
     *          with respect to major circumferential angle.
     * \param powerOfMinorCircumferentialAngleDerivative Power of the derivative
     *          with respect to the minor circumferential angle.
     * \return Surface derivative on torus.
     */
    Eigen::VectorXd getSurfaceDerivative( const double majorCircumferentialAngle,
                                          const double minorCircumferentialAngle,
                                          const int powerOfMajorCircumferentialAngleDerivative,
                                          const int powerOfMinorCircumferentialAngleDerivative );

    //! Get parameter of torus.
    /*!
     * Retrieves a parameter of the torus.
     * Function uses parameter_ member variable to prevent multiple
     * declarations.
     * \param index Index of parameter to return ( index = 0: returns major
     *          radius; index = 1: returns minor radius ).
     * \return Selected parameter.
     */
    double getParameter( const int index );

    //! Get major radius.
    /*!
     * Returns the major radius.
     * \return Major radius.
     */
    double getMajorRadius( ) { return majorRadius_; }

    //! Get minor radius.
    /*!
     * Returns the minor radius.
     * \return Minor radius.
     */
    double getMinorRadius( ) { return minorRadius_; }

    //! Get maximum of major circumferential angle.
    /*!
     * Returns the maximum value of the major circumferential angle.
     * \return Maximum value of major circumferential angle.
     */
    double getMaximumMajorCircumferentialAngle( ) { return getMaximumIndependentVariable( 1 ); }

    //! Get maximum of minor circumferential angle.
    /*!
     * Returns the maximum of the minor circumferential angle.
     * \return Maximum minor circumferential angle.
     */
    double getMaximumMinorCircumferentialAngle( ) { return getMaximumIndependentVariable( 2 ); }

    //! Get minimum of major circumferential angle.
    /*!
     * Returns the minimum value of the major circumferential angle.
     * \return Minimum value of major circumferential angle.
     */
    double getMinimumMajorCircumferentialAngle( ) { return getMinimumIndependentVariable( 1 ); }

    //! Get minimum of minor circumferential angle.
    /*!
     * Returns the minimum value of the minor circumferential angle.
     * \return Minimum value of minor circumferential angle.
     */
    double getMinimumMinorCircumferentialAngle( ) { return getMinimumIndependentVariable( 2 ); }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information, prints the class type,
     * the ranges for the minor and major circumferential angles, and the
     * major and minor radii.
     * \param stream Stream to which info is to be printed.
     * \param torus Torus frustum of which info is to be printed.
     * \return Stream with printed info.
     */
    friend std::ostream &operator << ( std::ostream &stream, Torus& torus );

protected:

private:

    //! Major radius.
    /*!
     * Major radius, i.e. radius of center of torus to center of tube.
     */
    double majorRadius_;

    //! Minor radius.
    /*!
     * Minor radius, i.e. radius of tube.
     */
    double minorRadius_;
};

//! Typedef for shared-pointer to Torus object.
typedef std::shared_ptr< Torus > TorusPointer;

} // namespace geometric_shapes
} // namespace tudat

#endif // TUDAT_TORUS_H
