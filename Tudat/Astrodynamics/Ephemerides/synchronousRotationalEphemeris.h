/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_SYNCHRONOUSROTATIONALEPHEMERIS_H
#define TUDAT_SYNCHRONOUSROTATIONALEPHEMERIS_H

#include <functional>

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

namespace ephemerides
{

//! Class to define a fully synchronous rotation model for a body
/*!
 *  Class to define a fully synchronous rotation model for a body. The body-fixed x-axis always points directly to the
 *  central body, and itts z-axis is along the orbital angular momentum vector (r x v).
 *  NOTE: The time derivative of the rotation matrix is not yet implemented.
 */
class SynchronousRotationalEphemeris: public RotationalEphemeris
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param relativeStateFunction Function returning the current state of the body relative to the central body, in
     * the base frame
     * \param centralBodyName Name of central body
     * \param baseFrameOrientation Name of base frame
     * \param targetFrameOrientation Name of target frame
     */
    SynchronousRotationalEphemeris(
            const std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction,
            const std::string& centralBodyName,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation ):
        RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
        relativeStateFunction_( relativeStateFunction ),
        isBodyInPropagation_( 0 ),
        centralBodyName_( centralBodyName ),
        warningPrinted_( false )
    { }

    //! Destructor
    ~SynchronousRotationalEphemeris( ){ }

    //! Calculate rotation quaternion from target frame to base frame.
    /*!
     * Function that calculates the rotation quaternion from target frame to base
     * frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \return Rotation quaternion computed.
     */
    Eigen::Quaterniond getRotationToBaseFrame( const double currentTime );

    //! Calculate rotation quaternion from base frame to target frame.
    /*!
     * Returns the rotation quaternion from base frame to target frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \return Rotation quaternion computed.
     */
    Eigen::Quaterniond getRotationToTargetFrame(
            const double currentTime )
    {
        return getRotationToBaseFrame( currentTime ).inverse( );
    }

    //! Function to calculate the derivative of the rotation matrix from target frame to base frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame to base
     *  frame at specified time.
     *  \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     *  \return Derivative of rotation from target to base frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double currentTime );

    //! Function to calculate the derivative of the rotation matrix from base frame to target frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from base frame to target
     *  frame at specified time.
     *  \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     *  \return Derivative of rotation from base to target frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double currentTime )
    {
        return getDerivativeOfRotationToBaseFrame( currentTime ).transpose( );
    }

    //! Function to define whether the body is currently being propagated, or not
    /*!
     *  Function to define whether the body is currently being propagated, or not
     *  \param isBodyInPropagation Boolean defining whether the body is currently being propagated, or not
     */
    void setIsBodyInPropagation( const bool isBodyInPropagation )
    {
        isBodyInPropagation_ = isBodyInPropagation;
    }

private:

    //! Function returning the current state of the body relative to the central body, in the base frame
    const std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction_;

    //!  Boolean defining whether the body is currently being propagated, or not
    bool isBodyInPropagation_;

    //! Name of central body
    std::string centralBodyName_;

    //!  Boolean defining whether the warning for the time-derivative of rotation amtrix has been printed.
    bool warningPrinted_;
};

}

}


#endif // TUDAT_SYNCHRONOUSROTATIONALEPHEMERIS_H
