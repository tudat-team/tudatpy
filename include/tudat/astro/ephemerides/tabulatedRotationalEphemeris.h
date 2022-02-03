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

#ifndef TUDAT_TABULATEDROTATIONALEPHEMERIS_H
#define TUDAT_TABULATEDROTATIONALEPHEMERIS_H


#include <map>
#include <vector>

#include <memory>

#include <Eigen/Core>

#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/basics/timeType.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{


//! Class that computes the current rotational state from tabulated values of the rotational state with the use of an interpolator
/*!
 *  Class that computes the current rotational state from tabulated values of the rotational state with the use of an interpolator
 *  The interpolated data consists of the four entries (w,x,y,z) of the quaternion from the target frame to the base frame, and
 *  body's angular velocity vector, expressed in its body-fixed frame (target frame). The quaternion is normalized to 1 before
 *  being set as the current rotation. Other properties of the rotation are derived from these two quantities.
 */
template< typename StateScalarType = double, typename TimeType = double >
class TabulatedRotationalEphemeris : public RotationalEphemeris
{
public:

    typedef Eigen::Matrix< StateScalarType, 7, 1 > StateType;
    typedef Eigen::Matrix< StateScalarType, 4, 1 > OrientationType;

    //! Constructor, sets rotational state interpolator and frame data.
    /*!
     *  Constructor, sets rotational state interpolator and frame data.
     *  \param interpolator Interpolator that returns the interpolated rotational state as a function of time.
     *  The interpolated data must consist of the four entries (w,x,y,z) of the quaternion from the target frame to the base
     *  frame, and body's angular velocity vector, expressed in its body-fixed frame (target frame).The quaternion is normalized
     *  to 1 before being set as the current rotation. Other properties of the rotation are derived from these two quantities.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    TabulatedRotationalEphemeris(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateType > >
            interpolator,
            const std::string& baseFrameOrientation = "ECLIPJ2000",
            const std::string& targetFrameOrientation = "" ):
        RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ), interpolator_( interpolator ),
    currentTime_( TUDAT_NAN ){  }

    //! Destructor
    ~TabulatedRotationalEphemeris( ){ }

    //! Function to reset the rotational state interpolator.
    /*!
     *  Function to reset the rotational state interpolator, for instance following an update of the rotational state of the body
     *  after a new numerical propagation of rotational equations of motion.
     *  \param interpolator New interpolator that returns the interpolated rotational state as a function of time.
     *  The interpolated data must consist of the four entries (w,x,y,z) of the quaternion from the target frame to the base
     *  frame, and body's angular velocity vector, expressed in its body-fixed frame (target frame).The quaternion is normalized
     *  to 1 before being set as the current rotation. Other properties of the rotation are derived from these two quantities.
     */
    void reset( const std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateType > > interpolator )
    {
        interpolator_ = interpolator;
    }

    //! Function to retrieve the rotational state interpolator.
    /*!
     * Function to retrieve the rotational state interpolator.
     * \return Interpolator that returns the interpolated rotational state as a function of time.
     *  The interpolated data must consist of the four entries (w,x,y,z) of the quaternion from the target frame to the base
     *  frame, and body's angular velocity vector, expressed in its body-fixed frame (target frame).
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateType > > getInterpolator( )
    {
        return interpolator_;
    }

    //! Get rotation quaternion from target frame to base frame.
    /*!
     * Function to calculate and return the rotation quaternion from target frame to base frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which rotational ephemeris is to be evaluated.
     * \return Rotation quaternion computed from target frame to base frame
     */
    Eigen::Quaterniond getRotationToBaseFrame( const double secondsSinceEpoch )
    {
        updateInterpolator( secondsSinceEpoch );
        return currentRotationToBaseFrame_.template cast< double >( );
    }

    //! Get rotation quaternion from base frame to target frame.
    /*!
     * Function to calculate and return the rotation quaternion from base frame to target frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which rotational ephemeris is to be evaluated.
     * \return Rotation quaternion computed from base frame to target frame
     */
    Eigen::Quaterniond getRotationToTargetFrame( const double secondsSinceEpoch )
    {
        return ( getRotationToBaseFrame( secondsSinceEpoch ) ).inverse( );
    }

    //! Function to retrieve the angular velocity vector of the body, expressed in the base frame.
    /*!
     * Function to retrieve the angular velocity vector of the body, expressed in the base frame.
     * \param secondsSinceEpoch Seconds since epoch at which rotational ephemeris is to be evaluated.
     * \return Angular velocity vector of body, expressed in base frame.
     */
    Eigen::Vector3d getRotationalVelocityVectorInBaseFrame( const double secondsSinceEpoch )
    {
        updateInterpolator( secondsSinceEpoch );

        return ( currentRotationToBaseFrame_ * currentRotationalVelocityVectorInTargetFrame_ ).template cast< double >( );
    }

    //! Function to retrieve the angular velocity vector of the body, expressed in the target (body-fixed) frame.
    /*!
     * Function to retrieve the angular velocity vector of the body, expressed in the target (body-fixed) frame.
     * \param secondsSinceEpoch Seconds since epoch at which rotational ephemeris is to be evaluated.
     * \return Angular velocity vector of body, expressed in target (body-fixed) frame.
     */
    Eigen::Vector3d getRotationalVelocityVectorInTargetFrame( const double secondsSinceEpoch )
    {
        updateInterpolator( secondsSinceEpoch );

        return currentRotationalVelocityVectorInTargetFrame_.template cast< double >( );
    }

    //! Function to calculate the derivative of the rotation matrix from base frame to target frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from base frame to target frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which rotational ephemeris is to be evaluated.
     *  \return Derivative of rotation from base to target (body-fixed) frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( const double secondsSinceEpoch )
    {
        updateInterpolator( secondsSinceEpoch );

        return getDerivativeOfRotationMatrixToFrame(
                    ( currentRotationToBaseFrame_.inverse( ) ).toRotationMatrix( ).template cast< double >( ),
                    ( currentRotationToBaseFrame_ * currentRotationalVelocityVectorInTargetFrame_ ).template cast< double >( ) );
    }

    //! Function to calculate the derivative of the rotation matrix from target frame to base frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame to base frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which rotational ephemeris is to be evaluated.
     *  \return Derivative of rotation from target (body-fixed) to base frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double secondsSinceEpoch )
    {
        updateInterpolator( secondsSinceEpoch );

        return getDerivativeOfRotationToTargetFrame( secondsSinceEpoch ).transpose( );
    }


private:

    //! Function to retrieve the current rotational state from the interpolator
    /*!
     * Function to retrieve the current rotational state from the interpolator
     * \param time Time at which to evaluate the interpoaltor
     */
    void updateInterpolator( const double time )
    {
        if( !( time == currentTime_ ) )
        {
            // Retrieve data from interpolator
            currentRotationalState_ = interpolator_->interpolate( time );

            // Normalize quaternion and set current rotation quaternion
            double quaternionNorm = ( currentRotationalState_.block( 0, 0, 4, 1 ) ).norm( );
            currentRotationalState_.block( 0, 0, 4, 1 ) = currentRotationalState_.block( 0, 0, 4, 1 ) / quaternionNorm;
            currentRotationToBaseFrame_ = Eigen::Quaternion< StateScalarType >(
                        currentRotationalState_( 0 ), currentRotationalState_( 1 ),
                        currentRotationalState_( 2 ), currentRotationalState_( 3 ) );

            // Set current angular velocity vector.
            currentRotationalVelocityVectorInTargetFrame_ = currentRotationalState_.block( 4, 0, 3, 1 );

            // Set time at which interpolator was evaluated.
            currentTime_ = time;
        }
    }

    //! Interpolator that returns the interpolated rotational state as a function of time.
    /*!
     * Interpolator that returns the interpolated rotational state as a function of time.
     * The interpolated data must consist of the four entries (w,x,y,z) of the quaternion from the target frame to the base
     * frame, and body's angular velocity vector, expressed in its body-fixed frame (target frame).
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateType > > interpolator_;

    //! Last time at which the updateInterpolator function was called
    double currentTime_;

    //! Vector retrieved from interpolator_ on last call to updateInterpolator.
    StateType currentRotationalState_;

    //! Angular velocity vector of body in body-fixed frame obtained at last call to updateInterpolator.
    Eigen::Matrix< StateScalarType, 3, 1 > currentRotationalVelocityVectorInTargetFrame_;

    //! Rotation from body-fixed frame to base frame obtained at last call to updateInterpolator.
    Eigen::Quaternion< StateScalarType > currentRotationToBaseFrame_;

};

//! Create a tabulated rotation model from a given rotation model and interpolation settings
/*!
 * Create a tabulated rotation model from a given rotation model and interpolation settings
 * \param ephemerisToInterrogate Rotation model from which the tabulated model is tpo be synthesized
 * \param startTime Start time for tabulated model
 * \param endTime End time for tabulated model
 * \param timeStep Constant time step for tabulated model
 * \param interpolatorSettings Interpolation settings for tabulated model
 * \return Tabulated rotation model, as synthesized from a given rotation model and interpolation settings
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< RotationalEphemeris > getTabulatedRotationalEphemeris(
        const std::shared_ptr< RotationalEphemeris > ephemerisToInterrogate,
        const TimeType startTime,
        const TimeType endTime,
        const TimeType timeStep,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) )
{
    typedef Eigen::Matrix< StateScalarType, 7, 1 > StateType;

    // Create state map that is to be interpolated
    std::map< TimeType, StateType >  stateMap;
    TimeType currentTime = startTime;
    while( currentTime <= endTime )
    {
        stateMap[ currentTime ] = ephemerisToInterrogate->getRotationStateVector( currentTime );
        currentTime += timeStep;
    }

    // Create tabulated ephemeris model
    return std::make_shared< TabulatedRotationalEphemeris< StateScalarType, TimeType > >(
                interpolators::createOneDimensionalInterpolator( stateMap, interpolatorSettings ),
                ephemerisToInterrogate->getBaseFrameOrientation( ),
                ephemerisToInterrogate->getTargetFrameOrientation( ) );

}

extern template class TabulatedRotationalEphemeris< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class TabulatedRotationalEphemeris< long double, double >;
extern template class TabulatedRotationalEphemeris< double, Time >;
extern template class TabulatedRotationalEphemeris< long double, Time >;
#endif

bool isTabulatedRotationalEphemeris( const std::shared_ptr< RotationalEphemeris > rotationalEphemeris );

}

}

#endif // TUDAT_TABULATEDROTATIONALEPHEMERIS_H

