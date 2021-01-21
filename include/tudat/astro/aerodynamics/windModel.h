/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_WIND_MODEL_H
#define TUDAT_WIND_MODEL_H

#include <Eigen/Core>

#include <memory>
#include <functional>

#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

namespace tudat
{
namespace aerodynamics
{

//! Base class for a wind model.
/*!
 *  Base class for a wind model.  The wind vector is defined as the local velocity of the atmosphere expressed in the
 *  frame corotating with the body.
 */
class WindModel
{
public:

    //! Constructor.
    WindModel( const reference_frames::AerodynamicsReferenceFrames associatedFrame =
            reference_frames::vertical_frame ):
        associatedFrame_( associatedFrame )
    {
        if( !( associatedFrame == reference_frames::inertial_frame ||
               associatedFrame == reference_frames::corotating_frame ||
               associatedFrame == reference_frames::vertical_frame ) )
        {
            throw std::runtime_error(
                        "Error when creating wind model, definition must be in inertial, corotating or vertical frame" );
        }
    }

    //! Destructor.
    virtual  ~WindModel( ){ }

    //! Function (pure virtual) to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
    /*!
     * Function (pure virtual) to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     * \param currentAltitude Altitude at which wind vector is to be retrieved.
     * \param currentLongitude Longitude at which wind vector is to be retrieved.
     * \param currentLatitude Latitude at which wind vector is to be retrieved.
     * \param currentTime Time at which wind vector is to be retrieved.
     * \return Wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     */
    virtual Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity(
            const double currentAltitude,
            const double currentLongitude,
            const double currentLatitude,
            const double currentTime ) = 0;

    reference_frames::AerodynamicsReferenceFrames getAssociatedFrame( )
    {
        return associatedFrame_;
    }

protected:

    reference_frames::AerodynamicsReferenceFrames associatedFrame_;
};

class ConstantWindModel: public WindModel
{
public:
    ConstantWindModel(
            const Eigen::Vector3d constantWindVelocity,
            const reference_frames::AerodynamicsReferenceFrames associatedFrame =
            reference_frames::vertical_frame ):
        WindModel( associatedFrame ), constantWindVelocity_( constantWindVelocity ){ }

    Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity(
            const double currentAltitude,
            const double currentLongitude,
            const double currentLatitude,
            const double currentTime )
    {
        return constantWindVelocity_;
    }

private:
    Eigen::Vector3d constantWindVelocity_;
};

//! Class for computing the wind velocity vector from a custom, user-defined function.
class CustomWindModel: public WindModel
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param windFunction Function that returns wind vector as a function of altitude, longitude, latitude and time (in that
     * order).
     */
    CustomWindModel(
            const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
            const reference_frames::AerodynamicsReferenceFrames associatedFrame =
            reference_frames::vertical_frame ):
        WindModel( associatedFrame ), windFunction_( windFunction ){ }

    //! Destructor
    ~CustomWindModel( ){ }

    //! Function to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
    /*!
     * Function to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     * \param currentAltitude Altitude at which wind vector is to be retrieved.
     * \param currentLongitude Longitude at which wind vector is to be retrieved.
     * \param currentLatitude Latitude at which wind vector is to be retrieved.
     * \param currentTime Time at which wind vector is to be retrieved.
     * \return Wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     */
    Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity(
            const double currentAltitude,
            const double currentLongitude,
            const double currentLatitude,
            const double currentTime )
    {
        return windFunction_( currentAltitude, currentLongitude, currentLatitude, currentTime );
    }

private:

    //! Function that returns wind vector as a function of altitude, longitude, latitude and time (in that order).
    std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction_;

};

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_WIND_MODEL_H
