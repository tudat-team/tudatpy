/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */

#ifndef TUDAT_BODYSHAPEMODEL_H
#define TUDAT_BODYSHAPEMODEL_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

namespace tudat
{

namespace basic_astrodynamics
{

//! Base class for body shape models.
/*!
 *  This is the base class for shape models for (celestial) bodies, such as spheres, ellipsoids,
 *  or tri-axial ellipsoids. It can be used to calculate the altitude above such a shape or the
 *  local radius.
 */
class BodyShapeModel
{
public:

    //! Default constructor.
    BodyShapeModel( ) { }

    //! Virtual destructor.
    virtual ~BodyShapeModel( ) { }

    //! Calculates the altitude above the shape
    /*!
     *  Function to calculate the altitude above the shape from the a body fixed position.
     *  \param bodyFixedPosition Cartesian, body-fixed position of the point at which the
     *  altitude is to be determined.
     *  \return Altitude above the body.
     */
    virtual double getAltitude( const Eigen::Vector3d& bodyFixedPosition ) = 0;

    //! Function to return the mean radius of the shape model
    /*!
     *  Function to return the mean radius of the shape model, to be used for functions
     *  requiring moderate to low accuracy for shape model, i.e. using a 'nearest equivalent'
     *  sphere.
     *  \return Average radius of shape model.
     */
    virtual double getAverageRadius( ) = 0;

protected:

};

//! Function to calculate the altitude of a point over a central body from positions of both the
//! point and the body (in any frame)
/*!
 *  Function to calculate the altitude of a point over a central body from positions of both the
 *  point and the body (in any frame). The rotation to the body-fixed frame is provided to calculate
 *  the input for the altitude function of the bodyShapeModel, which is to be in a body-fixed frame.
 *  \param bodyShapeModel Shape model of central body
 *  \param position Position of point of which altitude is to be calculated
 *  \param bodyPosition Position of central body above which altitude is to be calculated
 *  \param toBodyFixedFrame Rotation from frame in which input vectors are given to body-fixed frame.
 *  \return Altitude above body shape.
 */
double getAltitudeFromNonBodyFixedPosition(
        const boost::shared_ptr< BodyShapeModel > bodyShapeModel,
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& bodyPosition,
        const Eigen::Quaterniond& toBodyFixedFrame );

//! Function to calculate the altitude of a point over a central body from positions of both the
//! point and the body (in any frame)
/*!
 *  Function to calculate the altitude of a point over a central body from positions of both the
 *  point and the body (in any frame).  The rotation to the body-fixed frame is provided to
 *  calculate the input for the altitude function of the bodyShapeModel. The position of the central
 *  body and rotation are not provided directly, but as function pointers to allow easy binding with
 *  e.g. a Body class.
 *  \param bodyShapeModel Shape model of central body
 *  \param position Position of point of which altitude is to be calculated
 *  \param bodyPositionFunction Function returning position of central body above which altitude is
 *  to be calculated
 *  \param toBodyFixedFrameFunction Function returning rotation from frame in which input vectors are
 *  given to body-fixed frame.
 *  \return Altitude above body shape.
 */
double getAltitudeFromNonBodyFixedPositionFunctions(
        const boost::shared_ptr< BodyShapeModel > bodyShapeModel,
        const Eigen::Vector3d& position,
        const boost::function< Eigen::Vector3d( ) > bodyPositionFunction,
        const boost::function< Eigen::Quaterniond( ) > toBodyFixedFrameFunction );


} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_BODYSHAPEMODEL_H
