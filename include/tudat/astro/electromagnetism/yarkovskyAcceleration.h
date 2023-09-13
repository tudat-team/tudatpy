/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *          Pérez-Hernández, J. A., & Benet, L. (2022). Non-zero Yarkovsky acceleration for near-Earth 
 *          asteroid (99942) Apophis. Communications Earth & Environment, 3(1), Article 1. 
 *          DOI: https://doi.org/10.1038/s43247-021-00337-x
 */

#ifndef TUDAT_YARKOVSKYACCELERATION_H
#define TUDAT_YARKOVSKYACCELERATION_H

// FIXME: Remove unncecessary includes.
#include <functional>
#include <boost/lambda/lambda.hpp>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

namespace tudat
{

namespace electromagnetism
{

//! Class for calculating an Yarkovsky acceleration, based on (Pérez-Hernández & Benet, 2022). 
/*!
 * FIXME: Add the documentation.
 */
class YarkovskyAcceleration: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param yarkovskyParameter Yarkovsky parameter
     * \param bodyGravitationalParameterFunction Function that returns the state of the body.
     * \param centralBodyStateFunction Functon that returns the state of central body.
     */
    YarkovskyAcceleration(
            const double yarkovskyParameter,
            const std::function< Eigen::Vector6d( ) > bodyStateFunction,
            const std::function< Eigen::Vector6d( ) > centralBodyStateFunction =
            [ ]( ){ return Eigen::Vector6d::Zero( ); } ):
        bodyStateFunction_( bodyStateFunction ), centralBodyStateFunction_( centralBodyStateFunction )
    { 
        updateMembers( 0.0 );
    }

    //! Destructor
    ~YarkovskyAcceleration( ){ }

    //! Function to update constituent elements of Yarkovsky acceleration to current time
    /*!
     * Function to update constituent elements of Yarkovsky acceleration to current time
     * \param currentTime Time to which Yarkovsky acceleration elements are to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            // Calulate current relative state of accelerated body
            currentState_ = bodyStateFunction_( ) - centralBodyStateFunction_( );

            // Calculate current body-fixed state of accelerated body.
            currentRotationFromRswToInertialFrame_ = Eigen::Quaterniond(
                        reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix( currentState_ ) ).inverse( );

            // FIXME: Implement Yarkovsky acceleration model.
            currentLocalAcclereration_ = Eigen::Vector3d::Zero( );

            // Perform sanity check.
            if( currentLocalAcclereration_ != currentLocalAcclereration_ )
            {
                throw std::runtime_error( "Error when computing Yarkovsky acceleration, result is NaN" );
            }

            this->currentTime_ = currentTime;
            this->currentAcceleration_ = currentRotationFromRswToInertialFrame_ * currentLocalAcclereration_ ;
        }
    }


    //! Function to retrieve current state of the body that is undergoing the Yarkovsky acceleration, relative to central body
    /*!
     *  Function to retrieve Current state of the body that is undergoing the Yarkovsky acceleration, relative to central body,
     *  in propagation frame.
     *  \return Current state of the body that is undergoing the Yarkovsky acceleration, relative to central body, in global frame.
     */
    Eigen::Vector6d getCurrentState( )
    {
        return currentState_;
    }

    //! Function to retrieve quaternion defining the rotation from NTW to inertial frame.
    /*!
     * Function to retrieve quaternion defining the rotation from NTW to inertial frame.
     * \return Quaternion defining the rotation from NTW to inertial frame.
     */
    Eigen::Quaterniond getCurrentToInertialFrame( )
    {
        return currentRotationFromRswToInertialFrame_;
    }

    //! Function to retrieve current Yarkovsky acceleration in NTW frame.
    /*!
     * Function to retrieve current Yarkovsky acceleration in NTW frame.
     * \return Current Yarkovsky acceleration in NTW frame.
     */
    Eigen::Vector3d getCurrentLocalAcceleration( )
    {
        return currentLocalAcclereration_;
    }

    // //! Function to retrieve current true anomaly of accelerated body in its orbit about the central body.
    // /*!
    //  * Function to retrieve current true anomaly of accelerated body in its orbit about the central body.
    //  * \return Current true anomaly of accelerated body in its orbit about the central body.
    //  */
    // double getCurrentTrueAnomaly( )
    // {
    //     return currentTrueAnomaly_;
    // }


private:

    //! Yarkovsky Parameter
    double yarkovskyParameter_;

    //! State function of the body that is undergoing the Yarkovsky acceleration.
    std::function< Eigen::Vector6d( ) > bodyStateFunction_;

    //! State function of the body that is being orbited
    std::function< Eigen::Vector6d( ) > centralBodyStateFunction_;

    //! Current state of the body that is undergoing the Yarkovsky acceleration, relative to central body, in global frame.
    Eigen::Vector6d currentState_;

    //! Quaternion defining the rotation from NTW to inertial frame.
    Eigen::Quaterniond currentRotationFromRswToInertialFrame_;

    // //! Current true anomaly of accelerated body in its orbit about the central body.
    // double currentTrueAnomaly_;

    //! Current Yarkovsky acceleration in NTW frame.
    Eigen::Vector3d currentLocalAcclereration_;

};

}

}

#endif // TUDAT_YARKOVSKYACCELERATION_H
