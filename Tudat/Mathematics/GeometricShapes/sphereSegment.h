/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      100910    D. Dirkx          First version of file.
 *      100915    D. Dirkx          Modified to correct comments, 80-lines rule, etc.
 *      100928    D. Dirkx          Modifications following first checking iteration.
 *      100929    D. Dirkx          Creation of separate file for class.
 *      101125    D. Dirkx          Update of class, better get and set functions.
 *      110208    K. Kumar          Updated file header; Doxygen comments corrected; minor changes
 *                                  to functions.
 *      110209    D. Dirkx          Minor changes.
 *      110209    K. Kumar          Minor changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120323    D. Dirkx          Removed set functions; moved functionality to constructor.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_SPHERE_SEGMENT_H
#define TUDAT_SPHERE_SEGMENT_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Mathematics/GeometricShapes/singleSurfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

//! Sphere segment class.
/*!
 * Class that defines the sphere ( segment ) shape. Parametrization is based
 * on spherical coordinates, with azimuth and zenith angles as 1st and 2nd
 * variables (see basic_mathematics::convertSphericalToCartesian from Tudat Core).
 * In addition to the minimum and maximum of these two variables, the sphere radius is
 * required for defining the sphere.
 * Note that by using the scalingMatrix_ member variable, this class can also
 * represent ellipsoids.
 */
class SphereSegment : public SingleSurfaceGeometry
{
public:

    //! Default constructor.
    /*!
     *  Default constructor. Default angle values are set to a full sphere.
     * \param radius Radius of sphere.
     * \param minimumAzimuthAngle Minimum of azimuth angle.
     * \param maximumAzimuthAngle Maximum of azimuth angle.
     * \param minimumZenithAngle Minimum of zenith angle.
     * \param maximumZenithAngle Maximum of zenith angle.
     */
    SphereSegment( const double radius,
                   const double minimumAzimuthAngle = 0.0,
                   const double maximumAzimuthAngle
                   = 2.0 * mathematical_constants::PI,
                   const double minimumZenithAngle = 0.0,
                   const double maximumZenithAngle
                   = mathematical_constants::PI );

    //! Get surface point on sphere segment.
    /*!
     * Retrieves a surface point on the sphere, by performing the mapping of
     * independent variables to cartesian coordinates for the sphere segment by
     * spherical coordinate transformation.
     * Function uses cartesianPositionVector_ member variable. Values of
     * this vector set previous to function call are irrelevant.
     * \param azimuthAngle Azimuth angle.
     * \param zenithAngle Zenith angle.
     * \return Surface point on sphere.
     */
    Eigen::VectorXd getSurfacePoint( const double azimuthAngle, const double zenithAngle );

    //! Get surface derivative on sphere segment.
    /*!
     * Retrieves the derivatives of the surface point with respect to the two
     * independent variables. For powerOfZenithDerivative = 1 and
     * powerOfAzimuthDerivative = 2 the function returns:
     * \f[
     *      \frac{ d^{ 3 } ( x, y, z ) } { du * dv^{ 2 } }
     * \f]
     * ( with u and v the zenith and azimuth angles ).
     * \param azimuthAngle Azimuth angle.
     * \param zenithAngle Zenith angle.
     * \param powerOfZenithAngleDerivative Power of derivative with respect to
     *          azimuth angle.
     * \param powerOfAzimuthAngleDerivative Power of derivative with respect to
     *          zenith angle.
     * \return Surface derivative on sphere.
     */
    Eigen::VectorXd getSurfaceDerivative( const double azimuthAngle, const double zenithAngle,
                                          const int powerOfZenithAngleDerivative,
                                          const int powerOfAzimuthAngleDerivative );

    //! Get parameter of sphere segment.
    /*!
     * Returns a parameter of the sphere segment.
     * Function uses parameter_ member variable to prevent multiple
     * declarations.
     * \param index Index of parameter to return ( index = 0: returns radius ).
     * \return Selected parameter.
     */
    double getParameter( const int index );

    //! Get radius.
    /*!
     * Returns the radius of the sphere segment.
     * \return Radius of the sphere segment.
     */
    double getRadius( ) { return radius_; }

    //! Get maximum value of azimuth angle.
    /*!
     * Returns the maximum values of the azimuth angle.
     * \return Maximum value of azimuth angle.
     */
    double getMaximumAzimuthAngle( ) { return maximumIndependentVariable1_; }

    //! Get minimum value of azimuth angle.
    /*!
     * Returns the minimum value of the azimuth angle.
     * \return Minimum value of azimuth angle.
     */
    double getMinimumAzimuthAngle( ) { return minimumIndependentVariable1_; }

    //! Get maximum value of zenith angle.
    /*!
     * Returns the maximum value of the zenith angle.
     * \return Maximum value of zenith angle.
     */
    double getMaximumZenithAngle( ) {  return maximumIndependentVariable2_; }

    //! Get minimum value of the zenith angle.
    /*!
     * Returns the minimum value of the zenith.
     * \return Minimum value of zenith angle.
     */
    double getMinimumZenithAngle( ) { return minimumIndependentVariable2_; }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information, prints the class type,
     * the ranges for the azimuth and zenith angles and the radius.
     * \param stream Stream object.
     * \param sphereSegment Sphere segment.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, SphereSegment& sphereSegment );

protected:

private:

    //! Sphere radius.
    /*!
     * Sphere radius.
     */
    double radius_;
};

//! Typedef for shared-pointer to SphereSegment object.
typedef boost::shared_ptr< SphereSegment > SphereSegmentPointer;

} // namespace geometric_shapes
} // namespace tudat

#endif // TUDAT_SPHERE_SEGMENT_H
