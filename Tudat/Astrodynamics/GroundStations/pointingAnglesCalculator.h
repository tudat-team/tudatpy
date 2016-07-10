#ifndef POINTINGANGLESCALCULATOR_H
#define POINTINGANGLESCALCULATOR_H

#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/GroundStations/nominalGroundStationState.h"

namespace tudat
{

namespace ground_stations
{

//! Class to calculate the pointing angles (azimuth elevation) from a topocentric frame
/*!
 *  Class to calculate the pointing (i.e. viewing) angles (azimuth elevation) from a topocentric frame. The rotations between an inertial and
 *  body-fixed and body-fixed and tropocentric are provided to the constructor and member functions can be used to calculated the pointing
 *  angles from a pointing vector in an inertial frame and a time.
 */
class PointingAnglesCalculator
{
public:

    //! Constructor, takes functions defininf the sub-rotations from inertial to topocentric frame.
    /*!
     *  Constructor, takes functions defininf the sub-rotations from inertial to topocentric frame.
     *  \param rotationFromInertialToBodyFixedFrame Function returning the rotation from the inertial to the body-fixed frame at a specified time.
     *  \param rotationFromBodyFixedToTopoCentricFrame Function returning the rotation from the body-fixed to the topocentric frame at a
     *  specified time (note that this rotation is typically time-independent).
     */
    PointingAnglesCalculator(
            const boost::function< Eigen::Quaterniond( const double ) > rotationFromInertialToBodyFixedFrame,
            const boost::function< Eigen::Quaterniond( const double ) > rotationFromBodyFixedToTopoCentricFrame ):
        rotationFromInertialToBodyFixedFrame_( rotationFromInertialToBodyFixedFrame ),
        rotationFromBodyFixedToTopoCentricFrame_( rotationFromBodyFixedToTopoCentricFrame ){ }

    //! Function to calculate the elevation angle from body-fixed point to given point.
    /*!
     *  Function to calculate the elevation angle from body-fixed reference point (typically ground station) to given point. The elevation angle
     *  is calculatwed wrt the shape of the body on which the reference point is located.
     *  \param inertialVectorAwayFromStation Vector from reference point to target point (i.e. to which elevation angle is to be calculated),
     *  expressed in inertial frame.
     *  \param time Time at which elevation angle is to be calculated.
     *  \return Elevation angle from reference point to input point (inertialVectorAwayFromStation)
     */
    double calculateElevationAngle( const Eigen::Vector3d inertialVectorAwayFromStation, const double time );

    //! Function to calculate the azimuth angle from body-fixed point to given point.
    /*!
     *  Function to calculate the azimuth angle from body-fixed reference point (typically ground station) to given point. The azimuth angle
     *  is calculatwed wrt the shape of the body on which the reference point is located.
     *  \param inertialVectorAwayFromStation Vector from reference point to target point (i.e. to which azimuth angle is to be calculated),
     *  expressed in inertial frame.
     *  \param time Time at which azimuth angle is to be calculated.
     *  \return Azimuth angle from reference point to input point (inertialVectorAwayFromStation)
     */
    double calculationAzimuthAngle( const Eigen::Vector3d inertialVectorAwayFromStation, const double time );

    //! Function to calculate the elevation and azimuth angles from body-fixed point to given point.
    /*!
     *  Function to calculate the azimuth and elevationangle from body-fixed reference point to given point.
     *  \param inertialVectorAwayFromStation Vector from reference point to target point expressed in inertial frame.
     *  \param time Time at which angles are to be calculated.
     *  \return Elevation and ezimuth angle (in that order) from reference point to input point (inertialVectorAwayFromStation)
     */
    std::pair< double, double > calculatePointingAngles( const Eigen::Vector3d inertialVectorAwayFromStation, const double time );

    //! Function to convert vector in inertial frame to topocentric frame.
    /*!
     *  Function to convert vector in inertial frame to topocentric frame, with unit vectors in ENU (Earth-North-Up) order.
     *  \param inertialVectorAwayFromStation Vector from reference point to target point expressed in inertial frame.
     *  \param time Time at which rotation is to be calculated.
     *  \return Vector expressed in topocentric frame.
     */
    Eigen::Vector3d convertVectorFromInertialToTopocentricFrame( const Eigen::Vector3d& inertialVector, const double time );

private:

    //! Function returning the rotation from the inertial to the body-fixed frame.
    /*!
     *  Function returning the rotation from the inertial to the body-fixed frame at a specified time.
     */
    const boost::function< Eigen::Quaterniond( const double ) > rotationFromInertialToBodyFixedFrame_;

    //! Function returning the rotation from the body-fixed to the topocentric frame
    /*!
     *  Function returning the rotation from the body-fixed to the topocentric frame at a
     *  specified time (note that this rotation is typically time-independent).
     */
    const boost::function< Eigen::Quaterniond( const double ) > rotationFromBodyFixedToTopoCentricFrame_;
};

}

}

#endif // POINTINGANGLESCALCULATOR_H
