/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_CONSTANTROTATIONALEPHEMERIS_H
#define TUDAT_CONSTANTROTATIONALEPHEMERIS_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"


namespace tudat
{

namespace ephemerides
{

//! Rotational ephemeris class that gives a constant (i.e. time independent) state.
class ConstantRotationalEphemeris : public RotationalEphemeris
{
public:

    //! Constructor of a constant Ephemeris object
    /*!
     *  Constructor
     *  \param constantState Constant state value.
     *  \param baseFrameOrientation Origin of reference frame in which state is defined.
     *  \param targetFrameOrientation Orientation of reference frame in which state is defined.
     */
    ConstantRotationalEphemeris( const Eigen::Vector7d constantState,
                                 const std::string& baseFrameOrientation = "",
                                 const std::string& targetFrameOrientation = "" ):
        RotationalEphemeris( baseFrameOrientation, targetFrameOrientation )
    {
        updateConstantState( constantState );
    }


    //! Get rotation quaternion from target frame to base frame.
    /*!
     * Function to calculate and return the rotation quaternion from target frame to base frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which rotational ephemeris is to be evaluated.
     * \return Rotation quaternion computed from target frame to base frame
     */
    Eigen::Quaterniond getRotationToBaseFrame( const double secondsSinceEpoch )
    {
        return currentRotationToLocalFrame_.inverse( );
    }

    //! Get rotation quaternion from base frame to target frame.
    /*!
     * Function to calculate and return the rotation quaternion from base frame to target frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which rotational ephemeris is to be evaluated.
     * \return Rotation quaternion computed from base frame to target frame
     */
    Eigen::Quaterniond getRotationToTargetFrame( const double secondsSinceEpoch )
    {
        return currentRotationToLocalFrame_;
    }

    //! Function to calculate the derivative of the rotation matrix from base frame to target frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from base frame to target frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which rotational ephemeris is to be evaluated.
     *  \return Derivative of rotation from base to target (body-fixed) frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( const double secondsSinceEpoch )
    {
        return currentRotationToLocalFrameDerivative_;
    }

    //! Function to calculate the derivative of the rotation matrix from target frame to base frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame to base frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which rotational ephemeris is to be evaluated.
     *  \return Derivative of rotation from target (body-fixed) to base frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double secondsSinceEpoch )
    {
        return currentRotationToLocalFrameDerivative_.transpose( );
    }

    //! Modifies the constant state
    /*!
     * Changes the constant state value to a new value.
     * \param newState New value for constant state.
     */
    void updateConstantState( const Eigen::Vector7d& newState )
    {
        constantState_ = newState;

        Eigen::Quaterniond currentRotationToGlobalFrame =
                 Eigen::Quaterniond( constantState_( 0 ),
                                     constantState_( 1 ),
                                     constantState_( 2 ),
                                     constantState_( 3 ) );
         currentRotationToGlobalFrame.normalize( );

         currentRotationToLocalFrame_ = currentRotationToGlobalFrame.inverse( );

         Eigen::Matrix3d currentRotationMatrixToLocalFrame = ( currentRotationToLocalFrame_ ).toRotationMatrix( );
         currentRotationToLocalFrameDerivative_ = linear_algebra::getCrossProductMatrix(
                     constantState_.block( 4, 0, 3, 1 ) ) * currentRotationMatrixToLocalFrame;
    }

private:

    //! Constant rotational state vector
    Eigen::Vector7d constantState_;

    //! Rotation quaternion to local frame
    Eigen::Quaterniond currentRotationToLocalFrame_;

    //! Current derivative of rotation matrix to local frame
    Eigen::Matrix3d currentRotationToLocalFrameDerivative_;


};

} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_CONSTANTROTATIONALEPHEMERIS_H
