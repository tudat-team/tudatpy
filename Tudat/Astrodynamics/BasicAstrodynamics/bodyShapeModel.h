#ifndef BODYSHAPEMODEL_H
#define BODYSHAPEMODEL_H

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
 *      150409    D. Dirkx          Migrated from personal code.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 *    Notes
 *
 */


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
 *  \param position Position of central body above which altitude is to be calculated
 *  \param position Rotation from frame in which input vectors are given to body-fixed frame.
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
 *  \param position Function returning position of central body above which altitude is to be
 *  calculated
 *  \param position Function returning rotation from frame in which input vectors are given to
 *  body-fixed frame.
 *  \return Altitude above body shape.
 */
double getAltitudeFromNonBodyFixedPositionFunctions(
        const boost::shared_ptr< BodyShapeModel > bodyShapeModel,
        const Eigen::Vector3d& position,
        const boost::function< Eigen::Vector3d( ) > bodyPositionFunction,
        const boost::function< Eigen::Quaterniond( ) > toBodyFixedFrameFunction );


}

}

#endif // BODYSHAPEMODEL_H
