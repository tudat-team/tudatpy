/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ROTATIONMATRIXPARTIAL_H
#define TUDAT_ROTATIONMATRIXPARTIAL_H

#include <vector>

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace observation_partials
{

//! Function to calculate partial of a rotation matrix from a body-fixed to inertial frame w.r.t. a constant rotation rate.
/*!
 * Function to calculate partial of  a rotation matrix from a body-fixed to inertial frame, as computed by the
 * SimpleRotationalEphemeris class,  w.r.t. the constant rotation rate.
 * \param inertialBodyFixedToIntegrationFrame Rotation matrix at reference epoch
 * \param rotationRate Nominal rotation rate (about local z-axis)
 * \param timeSinceEpoch Elapsed time (in seconds) since reference epoch
 * \return Partial derivative of rotation matrix w.r.t. rotation rate.
 */
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameWrtConstantRotationRate(
        const Eigen::Quaterniond inertialBodyFixedToIntegrationFrame,
        const double rotationRate, const double timeSinceEpoch );

//! Function to calculate partial of a rotation matrix derivative (body-fixed to inertial) w.r.t. a constant rotation rate.
/*!
 * Function to calculate partial of  a rotation matrix derivative from a body-fixed to inertial frame, as computed by the
 * SimpleRotationalEphemeris class,  w.r.t. the constant rotation rate.
 * \param currentRotationFromLocalToGlobalFrame Rotation matrix to inertial frame at reference epoch
 * \param rotationRate Nominal rotation rate (about local z-axis)
 * \param timeSinceEpoch Elapsed time (in seconds) since reference epoch
 * \return Partial derivative of rotation matrix derivative w.r.t. rotation rate.
 */
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtConstantRotationRate(
        const Eigen::Matrix3d currentRotationFromLocalToGlobalFrame,
        const double rotationRate, const double timeSinceEpoch );

//! Function to calculate a partial of rotation matrix from a body-fixed to inertial frame w.r.t. a constant pole right
//! ascension and declination.
/*!
 * Function to calculate a partial of rotation matrix from a body-fixed to inertial frame, as computed by the
 * SimpleRotationalEphemeris class,  w.r.t. a constant pole right ascension and declination.
 * \param initialOrientationAngles Rotation Euler angles at reference epoch, in order right ascension, declination,
 * prime meridian
 * \param rotationRate Rotation rate (about local z-axis)
 * \param timeSinceEpoch Elapsed time (in seconds) since reference epoch
 * \return Partial derivatives of rotation matrix w.r.t. right ascension (entry 0) and declination (entry 1) of pole.
 */
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation(
        const Eigen::Vector3d initialOrientationAngles,
        const double rotationRate, const double timeSinceEpoch );

//! Function to calculate a partial of rotation matrix derivative from a body-fixed to inertial frame w.r.t. a constant
//! pole right  ascension and declination.
/*!
 * Function to calculate a partial of rotation matrix derivative from a body-fixed to inertial frame, as computed by the
 * SimpleRotationalEphemeris class,  w.r.t. a constant pole right ascension and declination.
 * \param initialOrientationAngles Rotation Euler angles at reference epoch, in order right ascension, declination,
 * prime meridian
 * \param rotationRate Rotation rate (about local z-axis)
 * \param timeSinceEpoch Elapsed time (in seconds) since reference epoch
 * \return Partial derivatives of rotation matrix w.r.t. right ascension (entry 0) and declination (entry 1) of pole.
 */
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtPoleOrientation(
        const Eigen::Vector3d initialOrientationAngles,
        const double rotationRate, const double timeSinceEpoch );

//! Base class for partial derivatives of rotation matrix from body fixed to inertial frame w.r.t. a parameter.
/*!
 *  Base class for partial derivatives of rotation matrix from body fixed to inertial frame w.r.t. a parameter.
 *  A derived class is implemented for each estimatable parameter which represents a property of a rotation matrix from a
 *  body-fixed frame (as defined by the isParameterRotationMatrixProperty function). Pointers to this class are used for
 *  both PositionPartials and the SphericalHarmonicsGravityPartial partial. In this manner, only a single partial object
 *  needs to be implemented for both the state and partial and the sh acceleration
 *  partial wrt a rotational parameter when such a parameter is added.
 */
class RotationMatrixPartial
{
public:

    RotationMatrixPartial( const std::shared_ptr< ephemerides::RotationalEphemeris > rotationModel ):
        rotationModel_( rotationModel ){ }

    //! Virtual destructor
    virtual ~RotationMatrixPartial( ){ }

    //! Function to calculate partials of rotation matrix wrt a parameter (for either double or VectorXd parameters).
    /*!
     *  Function to calculate partials of rotation matrix from the body-fixed to the inertial base frame wrt a parameter.
     *  For a double parameter, the size of the return vector is 1, for a VectorXd parameter, the size is equal to the size
     *  of the parameter and the entries of the vector returned from here correspond to partials of the same index in the
     *  parameter.
     *  \param time Time at which the partial(s) is(are) to be evaluated.
     *  \return Vector of partials of rotation matrix from body-fixed to inertial frame wrt parameter(s)
     */
    virtual std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixToBaseFrameWrParameter(
            const double time ) = 0;

    //! Function to calculate partials of rotation matrix derivative wrt a parameter (for double or VectorXd parameters).
    /*!
     *  Function to calculate partials of rotation matrix derivative from the body-fixed to the inertial base frame.
     *  For a double parameter, the size of the return vector is 1, for a VectorXd parameter, the size is equal to the size
     *  of the parameter and the entries of the vector returned from here correspond to partials of the same index in the
     *  parameter.
     *  \param time Time at which the partial(s) is(are) to be evaluated.
     *  \return Vector of partials of rotation matrix from body-fixed to inertial frame wrt parameter(s)
     */
    virtual std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
            const double time ) = 0;

    //! Function to calculate the partial of the position of a vector, which is given in a body-fixed frame, in the inertial
    //! frame wrt a parameter.
    /*!
     *  Function to calculate the partial of the position of a vector, which is given in a body-fixed frame, in the inertial
     *  frame wrt a parameter  denoting a property of the rotation between a body-fixed and an inertial frame. The type of
     *  parameter is defined by the derived class and a a separate derived class is impleneted for each parameter.
     *  An example is the change in position of a ground station,  as expressed in the inertial frame, wrt a rotational
     *  parameter. Each column of the return Matrix denotes a single entry of the parameter
     *  (so it is a Vector3d for a double parameter).
     *  \param time Time at which the partial is to be evaluated.
     *  \param vectorInLocalFrame Vector, expressed in the body-fixed frame, of which the partial in an inertial frame wrt
     *  parameter is to be determined.
     *  \return Partial of the value of the vector in an inertial frame wrt the parameter(s).
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfInertialPositionWrtParameter(
            const double time,
            const Eigen::Vector3d vectorInLocalFrame );

    //! Function to calculate the partial of the velocity of a vector, which is given in a body-fixed frame, in the inertial
    //! frame wrt a parameter.
    /*!
     *  Function to calculate the partial of the velocity of a vector, which is given in a body-fixed frame, in the inertial
     *  frame wrt a parameter  denoting a property of the rotation between a body-fixed and an inertial frame. The type of
     *  parameter is defined by the derived class and a a separate derived class is impleneted for each parameter.
     *  An example is the change in velocity of a ground station,  as expressed in the inertial frame, wrt a rotational
     *  parameter. Each column of the return Matrix denotes a single entry of the parameter
     *  (so it is a Vector3d for a double parameter).
     *  NOTE: This function currently only works for points that are static in a local frame (e.g. 0 velocity in this local
     *  frame).
     *  \param time Time at which the partial is to be evaluated.
     *  \param vectorInLocalFrame Vector, expressed in the body-fixed frame, of which the partial in an inertial frame wrt
     *  parameter is to be determined.
     *  \return Partial of the value of the vector in an inertial frame wrt the parameter(s).
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartialOfInertialVelocityWrtParameter(
            const double time,
            const Eigen::Vector3d vectorInLocalFrame );

    //! Function to return the secondary identifier of the estimated parameter
    /*!
     * Function to return the secondary identifier of the estimated parameter. This function returns an empty string by
     * default, and can be overridden in the derived class if the associated parameter has a secondary identifier.
     * \return  Secondary identifier of the estimated parameter
     */
    virtual std::string getSecondaryIdentifier( )
    {
        return "";
    }
protected:
    std::shared_ptr< ephemerides::RotationalEphemeris > rotationModel_;

};

//! Class to calculate a rotation matrix from a body-fixed to inertial frame w.r.t. a constant rotation rate.
/*!
 * Class to calculate a rotation matrix from a body-fixed to inertial frame, as computed by the SimpleRotationalEphemeris
 * class,  w.r.t. the constant rotation rate.
 */
class RotationMatrixPartialWrtConstantRotationRate: public RotationMatrixPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param bodyRotationModel Rotation model for which the partial derivative w.r.t. the rotation rate is to be taken.
     */
    RotationMatrixPartialWrtConstantRotationRate(
            const std::shared_ptr< ephemerides::SimpleRotationalEphemeris > bodyRotationModel ):
        RotationMatrixPartial( bodyRotationModel ),
        bodyRotationModel_( bodyRotationModel ){ }

    //! Destructor.
    ~RotationMatrixPartialWrtConstantRotationRate( ){ }

    //! Function to calculate partial of rotation matrix from the body-fixed to the inertial frame wrt a rotation rate.
    /*!
     *  Function to calculate partial of rotation matrix from the body-fixed to the inertial frame wrt a rotation rate,
     *  using the properties of the bodyRotationModel_ member.
     *  \param time Time at which the partial is to be evaluated.
     *  \return Vector of size 1 containing partials of rotation matrix from body-fixed to inertial frame wrt
     *  rotation rate.
     */
    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixToBaseFrameWrParameter(
            const double time )
    {
        return { calculatePartialOfRotationMatrixFromLocalFrameWrtConstantRotationRate(
                        bodyRotationModel_->getInitialRotationToTargetFrame( ).inverse( ),
                        bodyRotationModel_->getRotationRate( ),
                        time - bodyRotationModel_->getInitialSecondsSinceEpoch( ) ) };
    }

    //! Function to calculate partial of rotation matrix derivative from the body-fixed to inertial frame wrt a rotation rate.
    /*!
     *  Function to calculate partial of rotation matrix derivative from the body-fixed to inertial frame wrt a rotation rate,
     *  using the properties of the bodyRotationModel_ member.
     *  \param time Time at which the partial is to be evaluated.
     *  \return Vector of size 1 containing partials of rotation derivative matrix from body-fixed to inertial frame wrt
     *  rotation rate.
     */
    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
            const double time )
    {
        return { calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtConstantRotationRate(
                        bodyRotationModel_->getRotationToBaseFrame( time ).toRotationMatrix( ),
                        bodyRotationModel_->getRotationRate( ),
                        time - bodyRotationModel_->getInitialSecondsSinceEpoch( ) ) };
    }

private:

    //! Rotation model for which the partial derivative w.r.t. the rotation rate is to be taken.
    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > bodyRotationModel_;
};

//! Class to calculate a rotation matrix from a body-fixed to inertial frame w.r.t. a constant pole right ascension
//! and declination
/*!
 * Class to calculate a rotation matrix from a body-fixed to inertial frame, as computed by the SimpleRotationalEphemeris
 * class,  w.r.t. the constant pole right ascension and declination
 */
class RotationMatrixPartialWrtPoleOrientation: public RotationMatrixPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param bodyRotationModel Rotation model for which the partial derivative w.r.t. the pole position is to be taken.
     */
    RotationMatrixPartialWrtPoleOrientation(
            const std::shared_ptr< ephemerides::SimpleRotationalEphemeris > bodyRotationModel ):
        RotationMatrixPartial( bodyRotationModel ),
        bodyRotationModel_( bodyRotationModel ){ }

    //! Destructor.
    ~RotationMatrixPartialWrtPoleOrientation( ){ }

    //! Function to calculate partial of rotation matrix from the body-fixed to the inertial frame wrt pole
    //! right ascension and declination.
    /*!
     *  Function to calculate partial of rotation matrix from the body-fixed to the inertial frame  wrt pole
     *  right ascension and declination, using the properties of the bodyRotationModel_ member.
     *  \param time Time at which the partial is to be evaluated.
     *  \return Vector of size 2 containing partials of rotation matrix from body-fixed to inertial frame right ascension
     *  (entry 0) and declination (entry 1) of pole.
     */
    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixToBaseFrameWrParameter(
            const double time )
    {
        return calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation(
                    bodyRotationModel_->getInitialEulerAngles( ),
                    bodyRotationModel_->getRotationRate( ), time - bodyRotationModel_->getInitialSecondsSinceEpoch( ) );
    }

    //! Function to calculate partial of rotation matrix derivative from the body-fixed to the inertial frame wrt pole
    //! right ascension and declination.
    /*!
     *  Function to calculate partial of rotation matrix derivative from the body-fixed to the inertial frame  wrt pole
     *  right ascension and declination, using the properties of the bodyRotationModel_ member.
     *  \param time Time at which the partial is to be evaluated.
     *  \return Vector of size 2 containing partials of rotation matrix derivative from body-fixed to inertial frame right ascension
     *  ascension (entry 0) and declination (entry 1) of pole.
     */
    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
            const double time )
    {
        return calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtPoleOrientation(
                    bodyRotationModel_->getInitialEulerAngles( ),
                    bodyRotationModel_->getRotationRate( ), time - bodyRotationModel_->getInitialSecondsSinceEpoch( ) );

    }

private:

    //! Rotation model for which the partial derivative w.r.t. the rotation rate is to be taken.
    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > bodyRotationModel_;
};

//! Class to calculate a rotation matrix from a body-fixed to inertial frame w.r.t. the associated quaternion elements
class RotationMatrixPartialWrtQuaternion: public RotationMatrixPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param currentRotationToInertialFrameFunction Function returning the current quaternion for body-fixed to inertial rotation
     */
    RotationMatrixPartialWrtQuaternion(
            const std::function< Eigen::Quaterniond( ) > currentRotationToInertialFrameFunction ):
        RotationMatrixPartial( nullptr ),
        currentRotationToInertialFrameFunction_( currentRotationToInertialFrameFunction )
    {
        currentQuaternionPartials_.resize( 4 );
    }

    //! Destructor
    ~RotationMatrixPartialWrtQuaternion( ){ }

    //! Function to compute the required partial derivative of rotation matrix.
    /*!
     * Function to compute the partial derivative of rotation matrix from a body-fixed to inertial frame w.r.t.
     * the associated quaternion elements
     * \param time Time at which partials are to be computed
     * \return Vector of size 4 containing partials of rotation matrix from body-fixed to inertial frame w.r.t. the four
     * associated quaternion elements
     */
    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixToBaseFrameWrParameter(
            const double time )
    {
        linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion(
        linear_algebra::convertQuaternionToVectorFormat( currentRotationToInertialFrameFunction_( ) ),
                   currentQuaternionPartials_ );

        return currentQuaternionPartials_;
    }

    //! Function to compute the required partial derivative of rotation matrix derivative.
    /*!
     * Function to compute the partial derivative of derivative of rotation matrix from a body-fixed to inertial frame w.r.t.
     * the associated quaternion elements. NOTE: function not yet implemented
     * \param time Time at which partials are to be computed
     * \return Vector of size 4 containing partials of rotation matrix derivative from body-fixed to inertial frame w.r.t. the
     * four associated quaternion elements
     */
    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
            const double time )
    {
        throw std::runtime_error( "Error when calling RotationMatrixPartialWrtQuaternion::calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter, function not yet implemented." );

    }

private:

    //! Function returning the current quaternion for body-fixed to inertial rotation
    std::function< Eigen::Quaterniond( ) > currentRotationToInertialFrameFunction_;

    //! List of rotation matrix partial derivatives, as last computed by calculatePartialOfRotationMatrixToBaseFrameWrParameter
    std::vector< Eigen::Matrix3d > currentQuaternionPartials_;
};

//! Class to calculate a rotation matrix from a body-fixed to inertial frame w.r.t. the rotational state
/*!
 *  Class to calculate a rotation matrix from a body-fixed to inertial frame w.r.t. the rotational state, consisting of the
 *  associated quaternion elements and body-fixed angular velocity vector.
 */
class RotationMatrixPartialWrtRotationalState: public RotationMatrixPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param currentRotationToInertialFrameFunction Function returning the current quaternion for body-fixed to inertial rotation
     */
    RotationMatrixPartialWrtRotationalState(
            const std::function< Eigen::Quaterniond( const double ) > currentRotationToInertialFrameFunction ):
        RotationMatrixPartial( nullptr ),
        currentRotationToInertialFrameFunction_( currentRotationToInertialFrameFunction )
    {
        currentQuaternionPartials_.resize( 7 );
    }

    //! Destructor
    ~RotationMatrixPartialWrtRotationalState( ){ }

    //! Function to compute the required partial derivative of rotation matrix.
    /*!
     * Function to compute the partial derivative of rotation matrix from a body-fixed to inertial frame w.r.t.
     * the rotational state. Only derivatives w.r.t. quaternion elements are non-zero
     * \param time Time at which partials are to be computed
     * \return Vector of size 7 containing partials of rotation matrix from body-fixed to inertial frame w.r.t. the four
     * associated rotational state
     */
    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixToBaseFrameWrParameter(
            const double time )
    {
        linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion(
        linear_algebra::convertQuaternionToVectorFormat( currentRotationToInertialFrameFunction_( time ) ),
                   currentQuaternionPartials_ );

        for( int i = 0; i < 3; i++ )
        {
            currentQuaternionPartials_[ i + 4 ] =  Eigen::Matrix3d::Zero( );
        }
        return currentQuaternionPartials_;
    }

    //! Function to compute the required partial derivative of rotation matrix derivative.
    /*!
     * Function to compute the partial derivative of derivative of rotation matrix from a body-fixed to inertial frame w.r.t.
     * the rotational state vector. NOTE: function not yet implemented
     * \param time Time at which partials are to be computed
     * \return Vector of size 7 containing partials of rotation matrix derivative from body-fixed to inertial frame w.r.t. the
     * rotational state vector
     */
    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
            const double time )
    {
        throw std::runtime_error( "Error when calling RotationMatrixPartialWrtQuaternion::calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter, function not yet implemented." );

    }

private:

    //! Function returning the current quaternion for body-fixed to inertial rotation
    std::function< Eigen::Quaterniond( const double ) > currentRotationToInertialFrameFunction_;

    //! List of rotation matrix partial derivatives w.r.t. rotational state vector, as last computed by
    //! calculatePartialOfRotationMatrixToBaseFrameWrParameter
    std::vector< Eigen::Matrix3d > currentQuaternionPartials_;
};

//! Typedef of list of RotationMatrixPartial objects, ordered by parameter.
typedef std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
std::shared_ptr< observation_partials::RotationMatrixPartial > > RotationMatrixPartialNamedList;

}

}

#endif // TUDAT_ROTATIONMATRIXPARTIAL_H
