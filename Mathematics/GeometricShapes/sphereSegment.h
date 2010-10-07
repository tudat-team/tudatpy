/*! \file sphereSegment.h
 *  This file contains the definition of the Sphere Segment class
 *
 *  Path              : Mathematics/Geometry/
 *  Version           : 1
 *  Check status      : Checked
 *
 *  Author            : Dominic Dirkx
 *  Affiliation       : TU Delft
 *  E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *  Checker           : J. Melman
 *  Affiliation       : Delft University of Technology
 *  E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *  Date created      : 29 September, 2010
 *  Last modified     : 29 September, 2010
 *
 *  References
 *
 *  Notes
 *    The setParameter and getParameter functions could benefit from
 *    taking an enum as input instead of an int.
 *
 *    The original file was split into several parts, so that each class is
 *    defined separately in a file. Hence, the this file was created after the
 *    first entries in the changelog and the file version is lower than the
 *    number of entries in changelog.
 *
 *  Copyright (c) 2010 Delft University of Technology.
 *
 *  This software is protected by national and international copyright.
 *  Any unauthorized use, reproduction or modification is unlawful and
 *  will be prosecuted. Commercial and non-private application of the
 *  software in any form is strictly prohibited unless otherwise granted
 *  by the authors.
 *
 *  The code is provided without any warranty; without even the implied
 *  warranty of merchantibility or fitness for a particular purpose.
 *
 *  Changelog
 *    100910   D. Dirkx                    First version of file
 *    100915   D. Dirkx                    Modified to correct comments, 80
 *                                         lines rule, etc.
 *    100928   D. Dirkx                    Modifications following first
 *                                         checking iteration.
 *    100929   D. Dirkx                    Creation of separate file for class
 *
 */

#ifndef SPHERESEGMENT_H
#define SPHERESEGMENT_H

#include "surfaceGeometry.h"

//! Sphere segment class.
/*!
 * Class that defines the sphere (segment) shape. Parametrization is based
 * on spherical coordinates, with independentVariable1 representing the azimuth
 * angle and independentVariable2 the zenith angle (see
 * mathematics::convertSphericalToCartesian). In addition to the minimum and
 * maximum of these two variables the sphere radius is required for defining
 * the sphere.
 */
class SphereSegment : public SurfaceGeometry
{
public:

    /*!
     *  Default constructor.
     */
    SphereSegment( );

    //! Constructor to initialize parameters.
    /*!
     * Constructor for radius and independent variable bounds (1: azimuth,
     * 2: zenith).
     * \param minimumindependentVariable1 Lower bound for azimuth.
     * \param maximumindependentVariable1 Upper bound for azimuth.
     * \param minimumindependentVariable2 Lower bound for zenith.
     * \param maximumindependentVariable2 Upper bound for zenith.
     * \param radius Sphere radius.
     */
    SphereSegment( const double& minimumindependentVariable1 ,
                   const double& maximumindependentVariable1 ,
                   const double& minimumindependentVariable2 ,
                   const double& maximumindependentVariable2 ,
                   const double& radius );

    /*!
     * Default destructor.
     */
    virtual ~SphereSegment( );

    //! Function to retrieve a surface point on the sphere.
    /*!
     * Performs the mapping of parameterized variables to cartesian coordinates
     * for the sphere by spherical coordinate transformation and use of offset
     * and rotation settings.
     * \param independentVariable1 Azimuth angle.
     * \param independentVariable2 Zenith angle.
     */
    virtual Vector getSurfacePoint( const double& independentVariable1 ,
                                    const double& independentVariable2 );

    //! Calculates the surface derivative on the sphere w.r.t. the independent
    //! variables.
    /*!
     * Function to retrieve the derivatives of the surface point w.r.t.
     * the two independent variables. derivative1 and derivative2
     * set the power of the derivative w.r.t. the two variables. For instance,
     * for (derivative1, derivative2)=(1,2), the function returns
     * d^3(x,y,z)/(du*dv^2) (with u and v the first and second independent
     * variable).
     * \param independentVariable1 Azimuth angle.
     * \param independentVariable2 Zenith angle.
     * \param powerOfDerivative1 Power of the derivative w.r.t. the azimuth angle.
     * \param powerOfDerivative2 Power of the derivative w.r.t. the zenith angle.
     */
    virtual Vector getSurfaceDerivative( const double& independentVariable1,
                                         const double& independentVariable2,
                                         const int& powerOfDerivative1,
                                         const int& powerOfDerivative2 );

    //! Function to retrieve a parameter value (radius).
    /*!
     * Function to retrieve the parameters of the geometry. For a sphere
     * only the radius exists as a parameter, so only for parameterIndex = 0
     * a parameter value is returned. For other values of parameterIndex, 0 is
     * returned and a warning is displayed.
     * \param parameterIndex Index of parameter.
     */
    virtual double getParameter( const int& parameterIndex );

    //! Function to set a parameters value (radius).
    /*!
     * Function to set the parameters of the geometry. For a sphere only the
     * radius exists as a parameter, so only for parameterIndex = 0 is a
     * parameter value set. For other values of parameterIndex, a warning is displayed.
     * \param parameterIndex Index of parameter.
     * \param value Value of parameter to set.
     */
    virtual void setParameter( const int& parameterIndex, const double& value );

    //! Overloaded ostream to print class information.
    /*!
     *  Overloaded ostream to print class information, prints the class type,
     *  the ranges for the azimuth and zenith angles and the radius.
     */
    friend std::ostream& operator<<( std::ostream& stream, SphereSegment& geometry );

protected:

private:

    //! Sphere radius.
    /*!
     *  Sphere radius.
     */
    double radius_;
};

#endif // SPHERESEGMENT_H

// End of file
