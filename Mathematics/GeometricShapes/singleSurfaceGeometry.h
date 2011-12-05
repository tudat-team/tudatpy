/*! \file surfaceGeometry.h
 *    This file contains the definition of the SingleSurfaceGeometry base
 *    class.
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 10
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
 *    Last modified     : 5 September, 2011
 *
 *    References
 *
 *    Notes
 *      Contents of this file used to be in singleGeometry.h, but as this class
 *      has been split into single and composite surface geometry, the contents
 *      have been moved, with most of the SurfaceGeometry class now belonging to
 *      the SingleSurfaceGeometry class.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      101125    D. Dirkx          Migration of contents to this file.
 *      110124    K. Kumar          Updated path; minor comment and layout
 *                                  changes; removed inherited variables.
 *      110204    K. Kumar          Minor comment and layout modifications;
 *                                  corrected Doxygen comments.
 *      110207    D. Dirkx          Removed overloaded ostream operator.
 *      110209    D. Dirkx          Minor changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

#ifndef SINGLESURFACEGEOMETRY_H
#define SINGLESURFACEGEOMETRY_H

// Include statements.
#include "Mathematics/GeometricShapes/surfaceGeometry.h"
#include "Mathematics/LinearAlgebra/linearAlgebra.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Surface geometry base class.
/*!
 * Base class for single surface geometry representations in terms of two
 * parameterizing variables. Class contains the minimum and maximum values of
 * these two variables as well as a function to retrieve them.
 *
 * Also, the offset vector, rotation matrix, and scaling matrix, which are used
 * to transform the shape from its 'standard' position and orientation to any
 * postion and orientation that is desired, as well as functions to set and
 * get these quantities, are included.
 */
class SingleSurfaceGeometry: public SurfaceGeometry
{
public:

    //! Independent variables.
    /*!
     * Independent variables.
     */
    enum IndependentVariables
    {
        firstIndependentVariable = 1,
        secondIndependentVariable = 2
    };

    //! Default constructor.
    /*!
     * Default constructor; sets the initial offset equal to zero and sets the
     * rotation matrix and scaling matrix equal to the identity matrix.
     */
    SingleSurfaceGeometry( ) : minimumIndependentVariable1_( -0.0 ),
        maximumIndependentVariable1_( -0.0 ), minimumIndependentVariable2_( -0.0 ),
        maximumIndependentVariable2_( -0.0 ), parameter_( -0.0 ),
        independentVariable_( firstIndependentVariable ),
        cartesianPositionVector_( VectorXd::Zero( 3 ) ), offset_( VectorXd::Zero( 3 ) ),
        rotationMatrix_( MatrixXd::Identity( 3, 3 ) ),
        scalingMatrix_( MatrixXd::Identity( 3, 3 ) ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~SingleSurfaceGeometry( ) { }

    //! Set offset of the shape.
    /*!
     * Sets the offset of the shape.
     * \param offset Offset.
     */
    void setOffset( const VectorXd& offset ) { offset_ = offset; }

    //! Set rotation matrix of the shape.
    /*!
     * Sets the rotation matrix of the shape.
     * \param rotationMatrix Rotation matrix.
     */
    void setRotationMatrix( const MatrixXd& rotationMatrix ) { rotationMatrix_ = rotationMatrix; }

    //! Set scaling matrix of the shape.
    /*!
     * Sets the scaling matrix of the shape.
     * \param scalingMatrix Scaling matrix.
     */
    void setScalingMatrix( const MatrixXd& scalingMatrix ) { scalingMatrix_ = scalingMatrix; }

    //! Set minimum value of independent variable.
    /*!
     * Sets the minimum value of a given independent variable.
     * \param parameterIndex Index of independent variable for which to set
                minimum value.
     * \param minimumValue Minimum value for given independent variable.
     */
    void setMinimumIndependentVariable( int parameterIndex, double minimumValue );

    //! Set maximum value of independent variable.
    /*!
     * Sets the maximum value of a given independent variable.
     * \param parameterIndex Index of independent variable for which to set
     *          maximum value.
     * \param maximumValue Maximum value for given independent variable.
     */
    void setMaximumIndependentVariable( int parameterIndex, double maximumValue );

    //! Set parameter.
    /*!
     * Sets parameter of the shape ( i.e., radius for derived SphereSegment,
     * cone angle for derived Cone ).
     * \param parameterIndex Index of parameter to set for a specific shape.
     * \param parameterValue Parameter value to set.
     */
    virtual void setParameter( int parameterIndex, double parameterValue ) = 0;

    //! Get parameter.
    /*!
     * Returns parameter from the shape. It returns 0 by default for the base
     * class.
     * \param parameterIndex Index of parameter to set for a specific shape.
     * \return Specific parameter from a shape.
     */
    virtual double getParameter( int parameterIndex ) = 0;

    //! Get offset from the shape.
    /*!
     * Returns the offset from the shape.
     * \return Offset of the surface.
     */
    VectorXd getOffset( ) { return offset_; }

    //! Get rotation matrix from the shape.
    /*!
     * Returns the rotation matrix from the shape.
     * \return Rotation matrix of the surface.
     */
    MatrixXd getRotationMatrix( ) { return rotationMatrix_; }

    //! Get scaling matrix from the shape.
    /*!
     * Returns the scaling matrix from the shape.
     * \return Scaling matrix of the surface.
     */
    MatrixXd getScalingMatrix( ) { return scalingMatrix_; }

    //! Get minimum value of independent variable.
    /*!
     * Gets the minimum value of a given independent variable.
     * \param parameterIndex Index of independent variable from which to
     *          retrieve minimum.
     * \return minimumValue Minimum value for given independent variable.
     */
    double getMinimumIndependentVariable( int parameterIndex );

    //! Get maximum value of independent variable.
    /*!
     * Gets the maximum value of a given independent variable.
     * \param parameterIndex Index of independent variable from which to
     *          retrieve maximum.
     * \return maximumValue Maximum value for given independent variable.
     */
    double getMaximumIndependentVariable( int parameterIndex );

    //! Get surface point.
    /*!
     * Returns a surface point. This function is defined separately for
     * each derived class. It contains no functionality for this base class.
     * \param independentVariable1 Independent variable 1.
     * \param independentVariable2 Independent variable 2.
     * \return Surface point.
     */
    virtual VectorXd getSurfacePoint( double independentVariable1,
                                      double independentVariable2 ) = 0;

    //! Get surface derivative.
    /*!
     * Returns a surface derivative. This function is defined separately for
     * each derived class. It contains no functionality for this base class.
     * \param independentVariable1 Independent variable 1.
     * \param independentVariable2 Independent variable 2.
     * \param powerOfDerivative1 Power of derivative wrt independent variable 1.
     * \param powerOfDerivative2 Power of derivative wrt independent variable 2.
     * \return Surface derivative.
     */
    virtual VectorXd getSurfaceDerivative( double independentVariable1,
                                           double independentVariable2,
                                           int powerOfDerivative1,
                                           int powerOfDerivative2 ) = 0;

    //! Apply transformation to vehicle part.
    /*!
     * Applies the surface transformation to the given point. The following
     * transformations are observed: local translation, scaling, rotation,
     * global translation.
     * \param point Point which is transformed.
     */
    void transformPoint( VectorXd& point );

protected:

    //! Minimum value of independent variable 1.
    /*!
     * Minimum value of independent variable 1.
     */
    double minimumIndependentVariable1_;

    //! Maximum value of independent variable 1.
    /*!
     * Maximum value of independent variable 1.
     */
    double maximumIndependentVariable1_;

    //! Minimum value of independent variable 2.
    /*!
     * Minimum value of independent variable 2.
     */
    double minimumIndependentVariable2_;

    //! Maximum value of independent variable 2.
    /*!
     * Maximum value of independent variable 2.
     */
    double maximumIndependentVariable2_;

    //! Parameter for use in get parameter function.
    /*!
     * Parameter for use in get parameter function.
     * Variable is declared here to prevent it being created multiple times
     * in local scope.
     */
    double parameter_;

    //! Indepedent variable.
    /*!
     * Indepedent variable.
     */
    IndependentVariables independentVariable_;

    //! Cartesian position vector for use in getSurfacePoint function.
    /*!
     * Cartesian position vector for use in getSurfacePoint function.
     * Variable is declared here to prevent it being created multiple times
     * in local scope.
     */
    VectorXd cartesianPositionVector_;

    //! Offset vector.
    /*!
     * Vector by which to translate the center of the geometric shape.
     */
    VectorXd offset_;

    //! Rotation matrix.
    /*!
     * Rotation matrix to be applied to geometry to obtain correct orientation.
     */
    MatrixXd rotationMatrix_;

    //! Scaling matrix.
    /*!
     * Scaling matrix of surface points.
     */
    MatrixXd scalingMatrix_;

private:
};

}

#endif // SINGLESURFACEGEOMETRY_H

// End of file.
