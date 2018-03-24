/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      Contents of this file used to be in singleGeometry.h, but as this class has been split
 *      into single and composite surface geometry, the contents have been moved, with most of the
 *      SurfaceGeometry class now belonging to the SingleSurfaceGeometry class.
 *
 */

#ifndef TUDAT_SINGLE_SURFACE_GEOMETRY_H
#define TUDAT_SINGLE_SURFACE_GEOMETRY_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Mathematics/GeometricShapes/surfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
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
    SingleSurfaceGeometry( )
        : minimumIndependentVariable1_( TUDAT_NAN ),
          maximumIndependentVariable1_( TUDAT_NAN ),
          minimumIndependentVariable2_( TUDAT_NAN ),
          maximumIndependentVariable2_( TUDAT_NAN ),
          parameter_( TUDAT_NAN ),
          independentVariable_( firstIndependentVariable ),
          cartesianPositionVector_( Eigen::VectorXd::Zero( 3 ) ),
          offset_( Eigen::VectorXd::Zero( 3 ) ),
          rotationMatrix_( Eigen::MatrixXd::Identity( 3, 3 ) ),
          scalingMatrix_( Eigen::MatrixXd::Identity( 3, 3 ) )
    { }

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
    void setOffset( const Eigen::VectorXd& offset ) { offset_ = offset; }

    //! Set rotation matrix of the shape.
    /*!
     * Sets the rotation matrix of the shape.
     * \param rotationMatrix Rotation matrix.
     */
    void setRotationMatrix( const Eigen::MatrixXd& rotationMatrix )
    {
        rotationMatrix_ = rotationMatrix;
    }

    //! Set scaling matrix of the shape.
    /*!
     * Sets the scaling matrix of the shape.
     * \param scalingMatrix Scaling matrix.
     */
    void setScalingMatrix( const Eigen::MatrixXd& scalingMatrix )
    {
        scalingMatrix_ = scalingMatrix;
    }

    //! Get parameter.
    /*!
     * Returns parameter from the shape. It returns 0 by default for the base
     * class.
     * \param parameterIndex Index of parameter to set for a specific shape.
     * \return Specific parameter from a shape.
     */
    virtual double getParameter( const int parameterIndex ) = 0;

    //! Get offset from the shape.
    /*!
     * Returns the offset from the shape.
     * \return Offset of the surface.
     */
    Eigen::VectorXd getOffset( ) { return offset_; }

    //! Get rotation matrix from the shape.
    /*!
     * Returns the rotation matrix from the shape.
     * \return Rotation matrix of the surface.
     */
    Eigen::MatrixXd getRotationMatrix( ) { return rotationMatrix_; }

    //! Get scaling matrix from the shape.
    /*!
     * Returns the scaling matrix from the shape.
     * \return Scaling matrix of the surface.
     */
    Eigen::MatrixXd getScalingMatrix( ) { return scalingMatrix_; }

    //! Get minimum value of independent variable.
    /*!
     * Gets the minimum value of a given independent variable.
     * \param parameterIndex Index of independent variable from which to
     *          retrieve minimum.
     * \return minimumValue Minimum value for given independent variable.
     */
    double getMinimumIndependentVariable( const int parameterIndex );

    //! Get maximum value of independent variable.
    /*!
     * Gets the maximum value of a given independent variable.
     * \param parameterIndex Index of independent variable from which to
     *          retrieve maximum.
     * \return maximumValue Maximum value for given independent variable.
     */
    double getMaximumIndependentVariable( const int parameterIndex );

    //! Get surface point.
    /*!
     * Returns a surface point. This function is defined separately for
     * each derived class. It contains no functionality for this base class.
     * \param independentVariable1 Independent variable 1.
     * \param independentVariable2 Independent variable 2.
     * \return Surface point.
     */
    virtual Eigen::VectorXd getSurfacePoint( const double independentVariable1,
                                             const double independentVariable2 ) = 0;

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
    virtual Eigen::VectorXd getSurfaceDerivative( const double independentVariable1,
                                                  const double independentVariable2,
                                                  const int powerOfDerivative1,
                                                  const int powerOfDerivative2 ) = 0;

    //! Apply transformation to vehicle part.
    /*!
     * Applies the surface transformation to the given point. The following
     * transformations are observed: local translation, scaling, rotation,
     * global translation.
     * \param point Point which is transformed.
     */
    void transformPoint( Eigen::VectorXd& point );

protected:

    //! Set minimum value of independent variable.
    /*!
     * Sets the minimum value of a given independent variable.
     * \param parameterIndex Index of independent variable for which to set
                minimum value.
     * \param minimumValue Minimum value for given independent variable.
     */
    void setMinimumIndependentVariable( const int parameterIndex, const double minimumValue );

    //! Set maximum value of independent variable.
    /*!
     * Sets the maximum value of a given independent variable.
     * \param parameterIndex Index of independent variable for which to set
     *          maximum value.
     * \param maximumValue Maximum value for given independent variable.
     */
    void setMaximumIndependentVariable( const int parameterIndex, const double maximumValue );

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
    Eigen::VectorXd cartesianPositionVector_;

    //! Offset vector.
    /*!
     * Vector by which to translate the center of the geometric shape.
     */
    Eigen::VectorXd offset_;

    //! Rotation matrix.
    /*!
     * Rotation matrix to be applied to geometry to obtain correct orientation.
     */
    Eigen::MatrixXd rotationMatrix_;

    //! Scaling matrix.
    /*!
     * Scaling matrix of surface points.
     */
    Eigen::MatrixXd scalingMatrix_;

private:
};

//! Typedef for shared-pointer to SingleSurfaceGeometry object.
typedef boost::shared_ptr< SingleSurfaceGeometry > SingleSurfaceGeometryPointer;

} // namespace geometric_shapes
} // namespace tudat

#endif // TUDAT_SINGLE_SURFACE_GEOMETRY_H
