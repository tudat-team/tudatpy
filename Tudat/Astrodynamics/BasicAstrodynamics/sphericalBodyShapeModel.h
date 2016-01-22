#ifndef SPHERICALBODYSHAPEMODEL_H
#define SPHERICALBODYSHAPEMODEL_H

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


#include <iostream>

#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"


namespace tudat
{
namespace basic_astrodynamics
{

//! Body shape model defined as a sphere
/*!
 *  Body shape model defined as a sphere, defined by only its radius.
 */
class SphericalBodyShapeModel: public BodyShapeModel
{
public:

    //! Constructor.
    /*!
     *  Constructor, defines body shape.
     *  \param radius Radius of sphere
     */
    SphericalBodyShapeModel ( const double radius ):
        radius_( radius ){ }

    //! Destructor
    ~SphericalBodyShapeModel( ){ }

    //! Calculates the altitude above the sphere
    /*!
     *  Function to calculate the altitude above the sphere from a body fixed position.
     *  \param bodyFixedPosition Cartesian, body-fixed position of the point at which the altitude
     *  is to be determined.
     *  \return Altitude above the sphere.
     */
    double getAltitude( const Eigen::Vector3d& bodyFixedPosition )
    {
        return bodyFixedPosition.norm( ) - radius_;
    }

    //! Function to return the mean radius of the shape model
    /*!
     *  Function to return the mean radius of the sphere, equal to its radius.
     *  \return Average radius of shape model (equal to radius).
     */
    double getAverageRadius( )
    {
        return radius_;
    }

private:

    //! Radius of sphere
    double radius_;
};

}

}

#endif // SPHERICALBODYSHAPEMODEL_H
