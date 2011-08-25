/*! \file sphereSegment.h
 *    This file contains the definition of the Sphere Segment class.
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 29 September, 2010
 *    Last modified     : 9 February, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100910    D. Dirkx          First version of file.
 *      100915    D. Dirkx          Modified to correct comments, 80-lines
 *                                  rule, etc.
 *      100928    D. Dirkx          Modifications following first checking
 *                                  iteration.
 *      100929    D. Dirkx          Creation of separate file for class.
 *      101125    D. Dirkx          Update of class, better get and set
 *                                  functions.
 *      110208    K. Kumar          Updated file header; Doxygen comments
 *                                  corrected; minor changes to functions.
 *      110209    D. Dirkx          Minor changes.
 */

#ifndef SPHERESEGMENT_H
#define SPHERESEGMENT_H

// Include statements.
#include "singleSurfaceGeometry.h"
#include "conicalFrustum.h"

//! Sphere segment class.
/*!
 * Class that defines the sphere ( segment ) shape. Parametrization is based
 * on spherical coordinates, with azimuth and zenith angles as 1st and 2nd
 * variables ( see mathematics::convertSphericalToCartesian ). In addition to
 * the minimum and maximum of these two variables, the sphere radius is
 * required for defining the sphere.
 * Note that by using the scalingMatrix_ member variable, this class can also
 * represent ellipsoids.
 */
class SphereSegment : public SingleSurfaceGeometry
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     */
    SphereSegment( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SphereSegment( );

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
    VectorXd getSurfacePoint( const double& azimuthAngle,
                              const double& zenithAngle );

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
    VectorXd getSurfaceDerivative( const double& azimuthAngle,
                                   const double& zenithAngle,
                                   const int& powerOfZenithAngleDerivative,
                                   const int& powerOfAzimuthAngleDerivative );

    //! Get parameter of sphere segment.
    /*!
     * Retrieves a parameter of the sphere segment.
     * Function uses parameter_ member variable to prevent multiple
     * declarations.
     * \param index Index of parameter to return ( index = 0: returns radius ).
     * \return Selected parameter.
     */
    double getParameter( const int& index );

    //! Set parameter of sphere segment.
    /*!
     * Set a parameter of the sphere segment. Here only the radius can be set.
     * \param index Index of parameter which is to be retrieved. Only
     *          index = 0 is valid, which sets the radius.
     * \param parameter Value of parameter which is set.
     */
    void setParameter( const int& index, const double& parameter );

    //! Get radius.
    /*!
     * Retrieves the radius of the sphere segment.
     * \return Radius of the sphere segment.
     */
    double& getRadius( );

    //! Set radius.
    /*!
     * Sets the radius of the sphere segment.
     * \param radius Radius of sphere segment.
     */
    void setRadius( const double& radius );

    //! Set maximum value of azimuth angle.
    /*!
     * Sets the maximum value of the azimuth angle.
     * \param maximumAzimuthAngle Maximum value of azimuth angle.
     */
    void setMaximumAzimuthAngle( const double& maximumAzimuthAngle );

    //! Set minimum value of azimuth angle.
    /*!
     * Sets the minimum value of the azimuth angle.
     * \param minimumAzimuthAngle Minimum value of azimuth angle.
     */
    void setMinimumAzimuthAngle( const double& minimumAzimuthAngle );
    
    //! Set maximum value of zenith angle.
    /*!
     * Sets the maximum value of the zenith angle.
     * \param maximumZenithAngle Maximum value of zenith angle.
     */
    void setMaximumZenithAngle( const double& maximumZenithAngle );
    
    //! Set minimum value of zenith angle.
    /*!
     * Sets the minimum value of the zenith angle.
     * \param minimumZenithAngle Minimum value of zenith angle.
     */
    void setMinimumZenithAngle( const double& minimumZenithAngle );
    
    //! Get maximum value of azimuth angle.
    /*!
     * Returns the maximum values of the azimuth angle.
     * \return Maximum value of azimuth angle.
     */
    double getMaximumAzimuthAngle( );

    //! Get minimum value of azimuth angle.
    /*!
     * Returns the minimum value of the azimuth angle.
     * \return Minimum value of azimuth angle.
     */
    double getMinimumAzimuthAngle( );
    
    //! Get maximum value of zenith angle.
    /*!
     * Returns the maximum value of the zenith angle.
     * \return Maximum value of zenith angle.
     */
    double getMaximumZenithAngle( );
    
    //! Get minimum value of the zenith angle.
    /*!
     * Returns the minimum value of the zenith.
     * \return Minimum value of zenith angle.
     */
    double getMinimumZenithAngle( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information, prints the class type,
     * the ranges for the azimuth and zenith angles and the radius.
     * \param stream Stream object.
     * \param sphereSegment Sphere segment.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     SphereSegment& sphereSegment );

protected:

private:

    //! Sphere radius.
    /*!
     * Sphere radius.
     */
    double radius_;
};

#endif // SPHERESEGMENT_H

// End of file.
