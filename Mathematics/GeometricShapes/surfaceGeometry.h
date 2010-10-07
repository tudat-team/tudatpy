/*! \file surfaceGeometry.h
 *  This file contains the definition of the Surface Geometry base class
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

#ifndef SURFACEGEOMETRY_H
#define SURFACEGEOMETRY_H

#include "geometricShape.h"

//! Surface geometry base class.
/*!
 * Base class for surface geometry representations in terms of two
 * parameterizing variables 1 and 2. Class contains the minimum and
 * maximum values of these two variables as well as a function to retrieve
 * these. Also the offset vector and rotation matrix, which are used
 * to transform the shape from its 'standard' position and
 * orientation to any postion and orientation that is desired,
 * as well as a function to set and get these quantities, are included.
 */
class SurfaceGeometry : public GeometricShape
{
public:

    //! Default constructor.
    /*!
     * Default constructor; sets the initial offset equal to zero and sets the
     * rotation matrix equal to the identity matrix.
     */
    SurfaceGeometry( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~SurfaceGeometry( );

    //! Pure virtual function to get surface point.
    /*!
     * Function to retrieve a surface point. This function is defined
     * separately for each derived class. It contains no functionality
     * for this base class.
     */
    virtual Vector getSurfacePoint( const double& independentVariable1,
                                    const double& independentVariable2 ) = 0 ;

    //! Pure virtual function to get surface derivative.
    /*!
     * Function to retrieve a surface derivative. This function is
     * defined separately for each derived class. It contains no functionality
     * for this base class.
     */
    virtual Vector getSurfaceDerivative( const double& independentVariable1 ,
                                         const double& independentVariable2 ,
                                         const int& powerOfDerivative1 ,
                                         const int& powerOfDerivative2 ) = 0;

    //! Function to set the rotation and translation values.
    /*!
     * Function to set the offset vector and rotation matrix of the geometry.
     * When retrieving a point from the geometry, first the rotation matrix is
     * used on the surface point, followed by translation over the offset vector.
     * \param offset Translation vector for surface geometry.
     * \param rotationMatrix Rotation matrix for vehicle geometry.
     */
    void setTransFormValues( const Vector& offset, const MatrixXd& rotationMatrix );

    //! Pure virtual function to retrieve parameter.
    /*!
     * Base class function to retrieve parameter from the shape. It returns 0
     * by default for the base class.
     */
    virtual double getParameter( const int& parameterIndex ) = 0;

    //! Pure virtual function to set parameter.
    /*!
     * Base class function to set parameter of the shape (i.e. radius for
     * derived SphereSegment, cone angle for derived Cone).
     */
    virtual void setParameter( const int& parameterIndex, const double& value ) = 0;

    //! Function to retrieve the rotation matrix from the shape.
    /*!
     * Function to retrieve rotation matrix from the shape.
     * \return Rotation matrix of shape.
     */
    MatrixXd getRotationMatrix( );

    //! Function to retrieve offset from the shape.
    /*!
     * Function to retrieve offset from the shape.
     * \return Translation vector for geometry.
     */
    Vector getOffset( );

    //! Function to retrieve minimum values of independent variables
    /*!
     * Function to retrieve minimum values of independent variables.
     * \param parameterIndex Index of independent variable from which to
     * retrieve minimum,
     * \return Minimum of independent variable.
     */
    double getMinimumIndependentVariable( const int& parameterIndex );

    //! Function to retrieve maximum values of independent variables.
    /*!
     * Function to retrieve maximum values of independent variables.
     * \param parameterIndex Index of independent variable from which to
     * retrieve minimum.
     * \return Maximum of independent variable.
     */
    double getMaximumIndependentVariable( const int& parameterIndex );

    //! Function to set the maximum values of independent variables.
    /*!
     * Function to set the maximum values of independent variables,
     * \param parameterIndex Index of independent variable from which to
     * set maximum.
     * \param value Upper bound for independent variable.
     */
    void setMaximumIndependentVariable( const int& parameterIndex,
                                        const double& value);

    //! Function to set the minimum values of independent variables,
    /*!
     * Function to set the minimum values of independent variables,
     * \param parameterIndex Index of independent variable from which to
     * set minimum.
     * \param value Lower bound for independent variable 1.
     */
    void setMinimumIndependentVariable( const int& parameterIndex,
                                        const double& value);

protected:

    //! Bounds of independent variables.
    /*!
     * Variables that denote the minimum and maximum values of the two
     * independent variables.
     */
    double minimumIndependentVariable1_, maximumIndependentVariable1_,
           minimumIndependentVariable2_, maximumIndependentVariable2_;

    //! Vector by which to translate the center of the geometric shape.
    /*!
     * Vector by which to translate the center of the geometric shape.
     */
    Vector offset_;

    /*!
     * Rotation matrix to be applied to geometry to obtain correct orientation.
     */
    MatrixXd rotationMatrix_;

private:
};

#endif // SURFACEGEOMETRY_H

// End of file
