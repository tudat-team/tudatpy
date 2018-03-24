#ifndef TUDAT_SPHERICALBODYSHAPEMODEL_H
#define TUDAT_SPHERICALBODYSHAPEMODEL_H

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

} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_SPHERICALBODYSHAPEMODEL_H
