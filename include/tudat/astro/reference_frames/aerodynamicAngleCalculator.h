/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_AERODYNAMICANGLECALCULATOR_H
#define TUDAT_AERODYNAMICANGLECALCULATOR_H

#include <vector>
#include <map>

#include <functional>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/reference_frames/dependentOrientationCalculator.h"
#include "tudat/astro/aerodynamics/windModel.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"

namespace tudat
{

namespace reference_frames
{

//! Enum to define ids for various angles needed for converting between inertial and body-fixed
//! frame, using transformation chain via aerodynamic frame.
enum AerodynamicsReferenceFrameAngles
{
    latitude_angle = 0,
    longitude_angle = 1,
    heading_angle = 2,
    flight_path_angle = 3,
    angle_of_attack = 4,
    angle_of_sideslip = 5,
    bank_angle = 6
};

//! Function to get a string representing a 'named identification' of an aerodynamic angle
/*!
 * Function to get a string representing a 'named identification' of an aerodynamic angle
 * \param angle Type of aerodynamic angle
 * \return String withaerodynamic angle id.
 */
std::string getAerodynamicAngleName( const AerodynamicsReferenceFrameAngles angle );

//! Object to calculate aerodynamic orientation angles from current vehicle state.
/*!
 *  Object to calculate aerodynamic orientation angles from current vehicle state.
 */
class AerodynamicAngleCalculator: public DependentOrientationCalculator
{
public:

    //! Constructor
    /*!
     *  Constructor.
     *  \param bodyFixedStateFunction Vehicle state in a frame fixed w.r.t. the central body.  Note
     *  that this state is w.r.t. the body itself, not w.r.t. the local atmosphere.
     *  \param rotationFromCorotatingToInertialFrame Function returning the quaternion that rotates from the corotating to
     *  the inertial frame.
     *  \param centralBodyName Name of central body w.r.t. which the angles are computed.
     *  \param calculateVerticalToAerodynamicFrame Boolean to determine whether to determine
     *  vertical <-> aerodynamic frame conversion when calling update function.
     *  \param angleOfAttackFunction Function to determine the angle of attack of the vehicle.
     *  \param angleOfSideslipFunction Function to determine the angle of sideslip of the vehicle.
     *  \param bankAngleFunction Function to determine the bank angle of the vehicle.
     * \param angleUpdateFunction Function to update the aerodynamic angles to the current time (default none).
     */
    AerodynamicAngleCalculator(
            const std::function< Eigen::Vector6d( ) > bodyFixedStateFunction,
            const std::function< Eigen::Quaterniond( ) > rotationFromCorotatingToInertialFrame,
            const std::string centralBodyName,
            const bool calculateVerticalToAerodynamicFrame = false,
            const std::function< double( ) > angleOfAttackFunction = std::function< double( ) >( ),
            const std::function< double( ) > angleOfSideslipFunction = std::function< double( ) >( ),
            const std::function< double( ) > bankAngleFunction = std::function< double( ) >( ),
            const std::function< void( const double ) > angleUpdateFunction = std::function< void( const double ) >( ) ):
        DependentOrientationCalculator( ),
        bodyFixedStateFunction_( bodyFixedStateFunction ),
        rotationFromCorotatingToInertialFrame_( rotationFromCorotatingToInertialFrame ),
        centralBodyName_( centralBodyName ),
        calculateVerticalToAerodynamicFrame_( calculateVerticalToAerodynamicFrame ),
        angleOfAttackFunction_( angleOfAttackFunction ),
        angleOfSideslipFunction_( angleOfSideslipFunction ),
        bankAngleFunction_( bankAngleFunction ),
        angleUpdateFunction_( angleUpdateFunction ),
        currentBodyAngleTime_( TUDAT_NAN )
    {
        currentAerodynamicAngles_.resize( 7 );
    }

    //! Function to set the atmospheric wind model
    /*!
     * Function to set the atmospheric wind model
     * \param windModel Model that computes the atmospheric wind as a function of position and time
     * \param shapeModel Shape model of central body, used in computation of altitude that is required for wind calculation
     */
    void setWindModel(
            const std::shared_ptr< aerodynamics::WindModel > windModel,
            const std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel )
    {
        windModel_ = windModel;
        shapeModel_ = shapeModel;
    }

    //! Function to get the current rotation from the global (propagation/inertial) to the local (body-fixed) frame.
    /*!
     * Function to get the current rotation from the global (propagation/inertial) to the local (body-fixed) frame.
     * Implementation of pure virtual base class function, calls getRotationQuaternionBetweenFrames function.
     * \return Current rotation from the global (propagation/inertial) to the local (body-fixed) frame.
     */
    Eigen::Quaterniond getRotationToLocalFrame( )
    {
        return getRotationQuaternionBetweenFrames( inertial_frame, body_frame );
    }

    //! Function to get the current rotation from the local (body-fixed) to the global (propagation/inertial) frame.
    /*!
     * Function to get the current rotation from the local (body-fixed) to the global (propagation/inertial) frame.
     * Implementation of pure virtual base class function, calls getRotationQuaternionBetweenFrames function.
     * \return Current rotation from the local (body-fixed) to the global (propagation/inertial) frame.
     */
    Eigen::Quaterniond getRotationToGlobalFrame( )
    {
        return getRotationQuaternionBetweenFrames( body_frame, inertial_frame );
    }

    //! Function to update the orientation angles to the current state.
    /*!
     *  Function to update the orientation angles to the current state, implements base class pure virtual function;
     *  calls update function. Note, when calling this function, the bank angle, angle of attack and sideslip angle are
     *  NOT updated (due to possible cross-dependency with aerodynamic coefficients).
     *  \sa AerodynamicAngleCalculator::update
     */
    void updateCalculator( const double currentTime )
    {
        update( currentTime, true );
    }

    //! Function to update the orientation angles to the current state.
    /*!
     *  This function updates all requires orientation angles to the current state of the vehicle.
     *  The current state is retrieved from the bodyFixedStateFunction_ member variable
     *  function pointer.
     *  \param currentTime Time to which angle calculator is to be updated.
     *  \param updateBodyOrientation Boolean denoting whether the trajectory<->body-fixed angles are to be updated.
     */
    void update( const double currentTime, const bool updateBodyOrientation );

    void getRotationQuaternionReferenceBetweenFrames(
            Eigen::Quaterniond& rotationToFrame,
            const AerodynamicsReferenceFrames originalFrame,
            const AerodynamicsReferenceFrames targetFrame );

    //! Function to get the rotation quaternion between two frames
    /*!
     * Function to get the rotation quaternion between two frames. This function uses the values
     * calculated by the previous call of the update( ) function.
     * \param originalFrame Id for 'current' frame
     * \param targetFrame Id for frame to which transformation object should transfrom a vector
     * (from originalFrame).
     * \return Rotation quaternion from originalFrame to targetFrame.
     */
    Eigen::Quaterniond getRotationQuaternionBetweenFrames(
            const AerodynamicsReferenceFrames originalFrame,
            const AerodynamicsReferenceFrames targetFrame );

    Eigen::Matrix3d getRotationMatrixBetweenFrames(
            const AerodynamicsReferenceFrames originalFrame,
            const AerodynamicsReferenceFrames targetFrame )
    {
        return getRotationQuaternionBetweenFrames( originalFrame, targetFrame ).toRotationMatrix( );
    }

    //! Function to get a single orientation angle.
    /*!
     * Function to get a single orientation angle, as calculated by previous call to update( )
     * function.
     * \param angleId Id of requested angle.
     * \return Value of requested angle.
     */
    double getAerodynamicAngle( const AerodynamicsReferenceFrameAngles angleId );

    //! Function to set the trajectory<->body-fixed orientation angles.
    /*!
     * Function to set the trajectory<->body-fixed orientation angles.
     * \param angleOfAttackFunction Function to return the angle of attack.
     * \param angleOfSideslipFunction Function to return the angle of sideslip.
     * \param bankAngleFunction Function to return the bank angle.
     * \param angleUpdateFunction Function to update the angles to the current time.
     */
    void setOrientationAngleFunctions(
            const std::function< double( ) > angleOfAttackFunction = std::function< double( ) >( ),
            const std::function< double( ) > angleOfSideslipFunction = std::function< double( ) >( ),
            const std::function< double( ) > bankAngleFunction =  std::function< double( ) >( ),
            const std::function< void( const double ) > angleUpdateFunction = std::function< void( const double ) >( ),
            const bool silenceWarnings = false );

    //! Function to set constant trajectory<->body-fixed orientation angles.
    /*!
     * Function to set constant trajectory<->body-fixed orientation angles.
     * \param angleOfAttack Constant angle of attack (default NaN, used if no angle is to be defined).
     * \param angleOfSideslip Constant angle of sideslip (default NaN, used if no angle is to be defined).
     * \param bankAngle Constant bank angle (default NaN, used if no angle is to be defined).
     */
    void setOrientationAngleFunctions(
            const double angleOfAttack = TUDAT_NAN,
            const double angleOfSideslip = TUDAT_NAN,
            const double bankAngle = TUDAT_NAN,
            const bool silenceWarnings = false );

    //! Function to get the function returning the quaternion that rotates from the corotating to the inertial frame.
    /*!
     * Function to get the function returning the quaternion that rotates from the corotating to the inertial frame.
     * \return Function returning the quaternion that rotates from the corotating to the inertial frame.
     */
    std::function< Eigen::Quaterniond( ) >  getRotationFromCorotatingToInertialFrame(  )
    {
        return rotationFromCorotatingToInertialFrame_;
    }

    //! Function to get the name of central body w.r.t. which the angles are computed.
    /*!
     * Function to get the name of central body w.r.t. which the angles are computed.
     * \return Name of central body w.r.t. which the angles are computed.
     */
    std::string getCentralBodyName( )
    {
        return centralBodyName_;
    }

    //! Function to get the current airspeed-based body-fixed state of vehicle, as set by previous call to update( ).
    /*!
     * Function to get the current airspeed-based body-fixed state of vehicle, as set by previous call to update( ).
     * \return Current airspeed-based body-fixed state of vehicle, as set by previous call to update( ).
     */
    Eigen::Vector6d getCurrentAirspeedBasedBodyFixedState( )
    {
        return currentBodyFixedAirspeedBasedState_;
    }

    //! Function to get the current groundspeed-based body-fixed state of vehicle, as set by previous call to update( ).
    /*!
     * Function to get the current groundspeed-based body-fixed state of vehicle, as set by previous call to update( ).
     * \return Current groundspeed-based body-fixed state of vehicle, as set by previous call to update( ).
     */
    Eigen::Vector6d getCurrentGroundspeedBasedBodyFixedState( )
    {
        return currentBodyFixedGroundSpeedBasedState_;
    }

    //! Function to get the current groundspeed-based body-fixed velocity of vehicle, as set by previous call to update( ).
    /*!
     * Function to get the current groundspeed-based body-fixed velocity of vehicle, as set by previous call to update( ).
     * \return Current groundspeed-based body-fixed velocity of vehicle, as set by previous call to update( ).
     */
    Eigen::Vector3d getCurrentGroundspeedBasedBodyFixedVelocity( )
    {
        return currentBodyFixedGroundSpeedBasedState_.segment( 3, 3 );
    }

    //! Function to reset the value of the currentBodyAngleTime_ variable
    /*!
     * Function to reset the value of the currentBodyAngleTime_ variable. Typically used to reset the time to NaN,
     * signalling the need to recompute all quantities upon the next relevant function call.
     * \param currentTime New current time.
     */
    void resetDerivedClassTime( const double currentTime = TUDAT_NAN )
    {
        currentBodyAngleTime_ = currentTime;
    }

private:

    //! Model that computes the atmospheric wind as a function of position and time
    std::shared_ptr< aerodynamics::WindModel > windModel_;

    //! Shape model of central body, used in computation of altitude that is required for wind calculation
    std::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel_;

    //! Map of current angles, as calculated by previous call to update( ) function.
    std::vector< double > currentAerodynamicAngles_;

    //! Map of current transformation quaternions, as calculated since previous call to update( ) function.
    std::map< std::pair< AerodynamicsReferenceFrames, AerodynamicsReferenceFrames >,
    Eigen::Quaterniond > currentRotationMatrices_;

    //! Current airspeed-based body-fixed state of vehicle, as set by previous call to update( ).
    Eigen::Vector6d currentBodyFixedAirspeedBasedState_;

    //! Current groundspeed-based body-fixed state of vehicle, as set by previous call to update( ).
    Eigen::Vector6d currentBodyFixedGroundSpeedBasedState_;

    //! Current rotation from central-body-corotating to inertial frame, as set by previous call to update( ).
    Eigen::Quaterniond currentRotationFromCorotatingToInertialFrame_;

    //! Vehicle state in a frame fixed w.r.t. the central body.
    /*!
     *  Vehicle state in a frame fixed w.r.t. the central body.
     *  Note that this state is w.r.t. the body itself, not w.r.t. the local atmosphere
     */
    std::function< Eigen::Vector6d( ) > bodyFixedStateFunction_;

    //! Function returning the quaternion that rotates from the corotating to the inertial frame.
    std::function< Eigen::Quaterniond( ) > rotationFromCorotatingToInertialFrame_;

    //! Name of central body w.r.t. which the angles are computed.
    std::string centralBodyName_;

    //! Boolean to determine whether to determine vertical <-> aerodynamic frame conversion
    //! when calling update function.
    bool calculateVerticalToAerodynamicFrame_;

    //! Function to determine the angle of attack of the vehicle.
    std::function< double( ) > angleOfAttackFunction_;

    //! Function to determine the angle of sideslip of the vehicle.
    std::function< double( ) > angleOfSideslipFunction_;

    //! Function to determine the bank angle of the vehicle.
    std::function< double( ) > bankAngleFunction_;

    //! Function to update the bank, attack and sideslip angles to current time.
    std::function< void( const double ) > angleUpdateFunction_;

    //! Current time to which the bank, attack and sideslip angles have been updated.
    double currentBodyAngleTime_;

};

//! Get a function to transform aerodynamic force from local to propagation frame.
/*!
 *  Get a function to transform aerodynamic force from local to propagation frame. The returned
 *  function takes the force in the accelerationFrame, and returns this force in the
 *  propagationFrame.
 *  \param aerodynamicAngleCalculator Object to calculate aerodynamic angles and rotation
 *  quaternions.
 *  \param accelerationFrame Id of frame in which aerodynamic force is calculated
 *  \param bodyFixedToInertialFrameFunction Function to get the central-body corotating frame to
 *  inertial frame quaternion.
 *  \param propagationFrame Id of frame in which orbit propagation is done.
 *  \return Aerodynamic force conversion function.
 */
std::function< Eigen::Vector3d( const Eigen::Vector3d& ) >
getAerodynamicForceTransformationFunction(
        const std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator,
        const AerodynamicsReferenceFrames accelerationFrame,
        const std::function< Eigen::Quaterniond( ) > bodyFixedToInertialFrameFunction =
        [ ]( ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
        const AerodynamicsReferenceFrames propagationFrame = inertial_frame );

// EFFICIENCY TODO: MAKE THIS AND PREVIOUS FUNCTION THE SAME
std::function< void( Eigen::Vector3d&, const Eigen::Vector3d& ) >
getAerodynamicForceTransformationReferenceFunction(
        const std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator,
        const AerodynamicsReferenceFrames accelerationFrame,
        const std::function< Eigen::Quaterniond&( ) > bodyFixedToInertialFrameFunction,
        const AerodynamicsReferenceFrames propagationFrame = inertial_frame );


//! Wrapper class to set closure between an imposed orientation of a body and its bank, sideslip and attack angles.
/*!
 * Wrapper class to set closure between an imposed orientation of a body and its bank, sideslip and attack angles.
 * Based on a given rotation matrix function and AerodynamicAngleCalculator, this class computes the required values
 * for the ank, sideslip and attack angles so that the output from AerodynamicAngleCalculator is consistent with the
 * given rotation matrix.
 */
class AerodynamicAnglesClosure
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param imposedRotationFromInertialToBodyFixedFrame Inertial to body-fixed frame rotation to which the
     * aerodynamicAngleCalculator object is to be made consistent
     * \param aerodynamicAngleCalculator Object from which the aerodynamic angles are computed.
     */
    AerodynamicAnglesClosure(
            const std::function< Eigen::Quaterniond( const double ) > imposedRotationFromInertialToBodyFixedFrame,
            const std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator ):
        imposedRotationFromInertialToBodyFixedFrame_( imposedRotationFromInertialToBodyFixedFrame ),
        aerodynamicAngleCalculator_( aerodynamicAngleCalculator )
    { }

    //! Function to update the aerodynamic angles to current time.
    /*!
     * Function to update the aerodynamic angles to current time, using closure between
     * imposedRotationFromInertialToBodyFixedFrame_ and aerodynamicAngleCalculator_.
     * \param currentTime Time to which angles are to be updated.
     */
    void updateAngles( const double currentTime );

    //! Function returning the current angle of attack, as computed by last call to updateAngles function.
    /*!
     *  Function returning the current angle of attack, as computed by last call to updateAngles function.
     * \return Current angle of attack, as computed by last call to updateAngles function.
     */
    double getCurrentAngleOfAttack( )
    {
        return currentAngleOfAttack_;
    }

    //! Function returning the current angle of sideslip, as computed by last call to updateAngles function.
    /*!
     *  Function returning the current angle of sideslip, as computed by last call to updateAngles function.
     * \return Current angle of sideslip, as computed by last call to updateAngles function.
     */
    double getCurrentAngleOfSideslip( )
    {
        return currentAngleOfSideslip_;
    }

    //! Function returning the current bank angle, as computed by last call to updateAngles function.
    /*!
     *  Function returning the current bank angle, as computed by last call to updateAngles function.
     * \return Current bank angle, as computed by last call to updateAngles function.
     */
    double getCurrentBankAngle( )
    {
        return currentBankAngle_;
    }

private:

    //! Inertial to body-fixed frame rotation to which the aerodynamicAngleCalculator_ object is to be made consistent.
    std::function< Eigen::Quaterniond( const double ) > imposedRotationFromInertialToBodyFixedFrame_;

    //! Object from which the aerodynamic angles are computed.
    std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator_;

    //! Current angle of attack, as computed by last call to updateAngles function.
    double currentAngleOfAttack_;

    //! Current angle of sideslip, as computed by last call to updateAngles function.
    double currentAngleOfSideslip_;

    //! Current bank angle, as computed by last call to updateAngles function.
    double currentBankAngle_;

    //! Current rotation matrix from body-fixed to trajectory, as computed by last call to updateAngles function.
     Eigen::Matrix3d currentRotationFromBodyToTrajectoryFrame_;

};

//! Function to make aerodynamic angle computation consistent with imposed body-fixed to inertial rotation.
/*!
 * Function to make aerodynamic angle computation consistent with imposed body-fixed to inertial rotation.
 * \param imposedRotationFromInertialToBodyFixedFrame Inertial to body-fixed frame rotation to which the
 * aerodynamicAngleCalculator object is to be made consistent
 * \param aerodynamicAngleCalculator Object from which the aerodynamic angles are computed.
 */
void setAerodynamicDependentOrientationCalculatorClosure(
        const std::function< Eigen::Quaterniond( const double ) > imposedRotationFromInertialToBodyFixedFrame,
        std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator );

//! Function to make aerodynamic angle computation consistent with existing DependentOrientationCalculator
/*!
 * Function to make aerodynamic angle computation consistent with existing DependentOrientationCalculator
 * \param dependentOrientationCalculator Object computing the current orientation based on the current state of the
 * environment. Aerodynamic angles are to be computed from output given by this class.
 * \param aerodynamicAngleCalculator Object from which the aerodynamic angles are computed.
 */
void setAerodynamicDependentOrientationCalculatorClosure(
            std::shared_ptr< DependentOrientationCalculator > dependentOrientationCalculator,
            std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator );

//! Function to make aerodynamic angle computation consistent with existing rotational ephemeris
/*!
 * Function to make aerodynamic angle computation consistent with existing  rotational ephemeris
 * \param rotationalEphemeris Object computing the current orientation of the body. Aerodynamic angles are to be computed
 * from output given by this class.
 * \param aerodynamicAngleCalculator Object from which the aerodynamic angles are computed.
 */
void setAerodynamicDependentOrientationCalculatorClosure(
            std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
            std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator );

} // namespace reference_frames

} // namespace tudat

#endif // TUDAT_AERODYNAMICANGLECALCULATOR_H
