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
    WindModel( const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
               const bool includeCorotation = true,
               const bool useRadius = false ):
        associatedFrame_( associatedFrame ), includeCorotation_( includeCorotation ), useRadius_( useRadius )
    {
        if( !( associatedFrame == reference_frames::inertial_frame || associatedFrame == reference_frames::corotating_frame ||
               associatedFrame == reference_frames::vertical_frame ) )
        {
            throw std::runtime_error( "Error when creating wind model, definition must be in inertial, corotating or vertical frame" );
        }
    }

    //! Destructor.
    virtual ~WindModel( ) {}

    //! Function (pure virtual) to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
    /*!
     * Function (pure virtual) to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     * \param currentAltitude Altitude at which wind vector is to be retrieved.
     * \param currentLongitude Longitude at which wind vector is to be retrieved.
     * \param currentLatitude Latitude at which wind vector is to be retrieved.
     * \param currentTime Time at which wind vector is to be retrieved.
     * \return Wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     */
    virtual Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity( const double currentAltitude,
                                                                      const double currentLongitude,
                                                                      const double currentLatitude,
                                                                      const double currentTime ) = 0;

    /*!
     * \brief Get the reference frame associated with this wind model.
     * \return Reference frame in which the wind model is defined
     */
    reference_frames::AerodynamicsReferenceFrames getAssociatedFrame( )
    {
        return associatedFrame_;
    }

    /*!
     * \brief Get whether atmospheric co-rotation is included.
     * \return Boolean indicating if atmospheric co-rotation is included in aerodynamic computations
     */
    bool getIncludeCorotation( ) const
    {
        return includeCorotation_;
    }

    /*!
     * \brief Get whether radius is used instead of altitude.
     * \return Boolean indicating if radius is used for wind computation instead of altitude
     */
    bool getUseRadius( ) const
    {
        return useRadius_;
    }

    /*!
     * \brief Setter for useRadius
     * \param useRadius Boolean indicating whether radius is used for wind computation instead of altitude
     */
    void setUseRadius( const bool useRadius )
    {
        useRadius_ = useRadius;
    }

protected:
    reference_frames::AerodynamicsReferenceFrames associatedFrame_;

    //! Boolean flag indicating whether atmospheric co-rotation should be included in aerodynamic computations
    bool includeCorotation_;

    //! Boolean flag indicating whether radius is used for wind computation instead of altitude
    bool useRadius_;
};

/*!
 * \brief Simple wind model with constant wind velocity.
 */
class ConstantWindModel : public WindModel
{
public:
    /*!
     * \brief Constructor for constant wind model.
     * \param constantWindVelocity Constant wind velocity vector [m/s]
     * \param associatedFrame Reference frame in which the wind is defined
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     * \param useRadius Boolean indicating whether radius is used for wind computation instead of altitude
     */
    ConstantWindModel( const Eigen::Vector3d constantWindVelocity,
                       const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
                       const bool includeCorotation = true,
                       const bool useRadius = false ):
        WindModel( associatedFrame, includeCorotation, useRadius ), constantWindVelocity_( constantWindVelocity )
    {}

    /*!
     * \brief Get the constant wind velocity vector.
     * \param currentAltitude Altitude (unused for constant wind)
     * \param currentLongitude Longitude (unused for constant wind)
     * \param currentLatitude Latitude (unused for constant wind)
     * \param currentTime Time (unused for constant wind)
     * \return Constant wind velocity vector [m/s]
     */
    Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity( const double currentAltitude,
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
class CustomWindModel : public WindModel
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param windFunction Function that returns wind vector as a function of altitude, longitude, latitude and time (in that
     * order).
     * \param associatedFrame Reference frame in which the wind is defined
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     * \param useRadius Boolean indicating whether radius is used for wind computation instead of altitude
     */
    CustomWindModel( const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
                     const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
                     const bool includeCorotation = true,
                     const bool useRadius = false ):
        WindModel( associatedFrame, includeCorotation, useRadius ), windFunction_( windFunction )
    {}

    //! Destructor
    ~CustomWindModel( ) {}

    //! Function to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
    /*!
     * Function to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     * \param currentAltitude Altitude at which wind vector is to be retrieved.
     * \param currentLongitude Longitude at which wind vector is to be retrieved.
     * \param currentLatitude Latitude at which wind vector is to be retrieved.
     * \param currentTime Time at which wind vector is to be retrieved.
     * \return Wind velocity vector in body-fixed, body-centered frame of body with atmosphere
     */
    Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity( const double currentAltitude,
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

//! Class for empty wind model (no wind, only controls co-rotation behavior)
/*!
 * Class for empty wind model that always returns zero wind velocity.
 * This is used when no physical wind model is required but co-rotation behavior needs to be specified.
 */
class EmptyWindModel : public WindModel
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     * \param useRadius Boolean indicating whether radius is used for wind computation instead of altitude
     */
    EmptyWindModel( const bool includeCorotation = true, const bool useRadius = false ):
        WindModel( reference_frames::vertical_frame, includeCorotation, useRadius )
    {}

    //! Destructor
    ~EmptyWindModel( ) {}

    //! Function to retrieve wind velocity vector (always zero)
    /*!
     * Function to retrieve wind velocity vector in body-fixed, body-centered frame of body with atmosphere.
     * Always returns zero vector as this is an empty wind model.
     * \param currentAltitude Altitude at which wind vector is to be retrieved (unused).
     * \param currentLongitude Longitude at which wind vector is to be retrieved (unused).
     * \param currentLatitude Latitude at which wind vector is to be retrieved (unused).
     * \param currentTime Time at which wind vector is to be retrieved (unused).
     * \return Zero wind velocity vector in body-fixed, body-centered frame
     */
    Eigen::Vector3d getCurrentBodyFixedCartesianWindVelocity( const double currentAltitude,
                                                              const double currentLongitude,
                                                              const double currentLatitude,
                                                              const double currentTime )
    {
        return Eigen::Vector3d::Zero( );
    }
};

}  // namespace aerodynamics

}  // namespace tudat

#endif  // TUDAT_WIND_MODEL_H
