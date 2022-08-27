/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDATBUNDLE_HYBRIDBODYSHAPEMODEL_H
#define TUDATBUNDLE_HYBRIDBODYSHAPEMODEL_H

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"
#include <iostream>

namespace tudat
{
namespace basic_astrodynamics
{

// Hybrid body shape model, consisting of a low-resolution shape model used for high altitudes and high-resolution model
// used for low altitudes. Useful when the evaluation of the altitude with high resolution is computationally expensive
// (e.g. polyhedron model).
class HybridBodyShapeModel: public BodyShapeModel
{
public:

    /*! Constructor.
     *
     * Constructor.
     * @param lowResolutionBodyShapeModel Model used to compute the altitude for altitudes higher than switchoverAltitude.
     * @param highResolutionBodyShapeModel Model used to compute the altitude for altitudes lower than switchoverAltitude.
     * @param switchoverAltitude Altitude at which the model used to compute the altitude is changed.
     */
    HybridBodyShapeModel(
            const std::shared_ptr< BodyShapeModel > lowResolutionBodyShapeModel,
            const std::shared_ptr< BodyShapeModel > highResolutionBodyShapeModel,
            double switchoverAltitude):
        lowResolutionBodyShapeModel_ ( lowResolutionBodyShapeModel ),
        highResolutionBodyShapeModel_ ( highResolutionBodyShapeModel ),
        switchoverAltitude_ ( switchoverAltitude )
    { }

    //! Destructor
    ~HybridBodyShapeModel( ){ }

    //! Calculates the altitude above the body.
    /*!
     *  Function to calculate the altitude above the body from a body fixed position.
     *  Function first evaluates the altitude with the low-resolution model. If the computed altitude is above the
     *  selected switchover altitude, then the function return the low-resolution altitude. If the computed altitude
     *  is below the switchover altitude, the function computes and returns the high-resolution altitude.
     *  \param bodyFixedPosition Cartesian, body-fixed position of the point at which the altitude
     *  is to be determined.
     *  \return Altitude above the polyhedron.
     */
    double getAltitude( const Eigen::Vector3d& bodyFixedPosition )
    {
        double altitude;
        double lowResolutionAltitude = lowResolutionBodyShapeModel_->getAltitude( bodyFixedPosition );

        if ( lowResolutionAltitude > switchoverAltitude_ )
        {
            altitude = lowResolutionAltitude;
        }
        else
        {
            altitude = highResolutionBodyShapeModel_->getAltitude( bodyFixedPosition );
        }

        return altitude;
    }

    //! Function to return the average radius of the polyhedron.
    /*!
     *  Function to return the average radius of the polyhedron.
     *  \return Average radius of the polyhedron.
     */
    double getAverageRadius( )
    {
        return highResolutionBodyShapeModel_->getAverageRadius( );
    }

    // Function to return the altitude at which the model used to compute the altitude is swithced between the high-
    // and low-resolution one (or vice-versa).
    double getSwitchoverAltitude( )
    {
        return switchoverAltitude_;
    }

    // Function to return the low-resolution model used to compute the altitude for high-altitudes.
    std::shared_ptr< BodyShapeModel > getLowResolutionBodyShapeModel ( )
    {
        return lowResolutionBodyShapeModel_;
    }

    // Function to return the high-resolution model used to compute the altitude for low-altitudes.
    std::shared_ptr< BodyShapeModel > getHighResolutionBodyShapeModel ( )
    {
        return highResolutionBodyShapeModel_;
    }

private:

    // Model used to compute the altitude for altitudes higher than switchoverAltitude.
    std::shared_ptr< BodyShapeModel > lowResolutionBodyShapeModel_;

    // Model used to compute the altitude for altitudes lower than switchoverAltitude.
    std::shared_ptr< BodyShapeModel > highResolutionBodyShapeModel_;

    // Altitude at which the model used to compute the altitude is changed.
    double switchoverAltitude_;
};


} // namespace basic_astrodynamics
} // namespace tudat

#endif //TUDATBUNDLE_HYBRIDBODYSHAPEMODEL_H
