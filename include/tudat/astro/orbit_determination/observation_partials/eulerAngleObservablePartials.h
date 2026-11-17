/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EULERANGLEOBSERVABLEPARTIALS_H
#define TUDAT_EULERANGLEOBSERVABLEPARTIALS_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <functional>

#include "tudat/basics/basicTypedefs.h"

#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/math/basic/rotationRepresentations.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute the partial derivatives of the 3-1-3 Euler angle observable w.r.t. the quaternion of the rotation from body-fixed
//! to base frame (the rotational state).
/*!
 * Function to compute the partial derivatives of the 3-1-3 Euler angle observable w.r.t. the quaternion of the rotation from body-fixed
 * to base frame (the rotational state). The observable is a function of the rotation to the body-fixed frame, the inverse of the
 * rotational state (with equal real part and opposite imaginary part). The partial changes sign with the quaternion, so it is computed
 * with the sign of the quaternion that is provided.
 * \param rotationToBaseFrame Rotation from body-fixed to base frame
 * \return Partial derivatives of the Euler angles w.r.t. the quaternion entries (real part first) of the rotation to the base frame
 */
inline Eigen::Matrix< double, 3, 4 > calculateEulerAngle313WrtRotationToBaseFrameQuaternionPartial(
        const Eigen::Quaterniond& rotationToBaseFrame )
{
    Eigen::Matrix< double, 3, 4 > partial = basic_mathematics::calculateEulerAngle313WrtQuaternionPartial( rotationToBaseFrame.inverse( ) );
    partial.block( 0, 1, 3, 3 ) *= -1.0;
    return partial;
}

//! Class to compute the partial derivatives of a 3-1-3 Euler angle observable.
class EulerAngleObervationPartialWrtCurrentRotationalState : public ObservationPartial< 3 >
{
public:
    //! Local typedef for return type (list of partial matrices and associated evaluation times).
    typedef std::vector< std::pair< Eigen::Matrix< double, 3, Eigen::Dynamic >, double > > EulerAngleObservationPartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param rotationToBaseFrameFunction Function returning the rotation from body-fixed to base frame as a function of time, with
     * the quaternion entries of the rotational state.
     */
    EulerAngleObervationPartialWrtCurrentRotationalState(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::function< Eigen::Quaterniond( const double ) > rotationToBaseFrameFunction ):
        ObservationPartial< 3 >( parameterIdentifier ), rotationToBaseFrameFunction_( rotationToBaseFrameFunction )
    {}

    //! Destructor
    ~EulerAngleObervationPartialWrtCurrentRotationalState( ) {}

    //! Function to calculate the observation partial(s) at required time.
    /*!
     *  Function to calculate the observation partial(s) at required time. Time is
     *  typically obtained from evaluation of observation model.
     *  \param states Link end state. (Not used; vector contains NaN entry)
     *  \param times Link end time (size 1 std::vector).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes)
     *  \return Vector of pairs containing partial values and associated times.
     */
    EulerAngleObservationPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Vector3d& currentObservation = Eigen::Vector3d::Constant( TUDAT_NAN ) )
    {
        EulerAngleObservationPartialReturnType returnPartial;

        currentTime_ = times.at( 0 );
        currentPartial_.setZero( );
        // The sign of the quaternion must be that of the state
        currentPartial_.block( 0, 0, 3, 4 ) =
                calculateEulerAngle313WrtRotationToBaseFrameQuaternionPartial( rotationToBaseFrameFunction_( currentTime_ ) );
        returnPartial.push_back( std::make_pair( currentPartial_, currentTime_ ) );
        return returnPartial;
    }

protected:
    //! Function returning the rotation from body-fixed to base frame as a function of time
    std::function< Eigen::Quaterniond( const double ) > rotationToBaseFrameFunction_;

    //! Partial derivative, as computed by last call to calculatePartial function
    Eigen::Matrix< double, 3, 7 > currentPartial_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;
};

//! Class to compute the partial derivatives of a 3-1-3 Euler angle observable w.r.t. a parameter of the rotation of the observed body
/*!
 * Class to compute the partial derivatives of a 3-1-3 Euler angle observable w.r.t. a parameter of the rotation of the observed body
 * (a rotation model parameter, or the translational state for a rotation model that depends on it), from the partial derivatives of
 * the rotation matrix R from body-fixed to base frame. A change dR = [ a x ] R of this rotation (with a the rotation vector of the
 * change) changes its quaternion q by dq = ( 0, a ) q / 2, which is mapped to the Euler angles by their partial w.r.t. q.
 */
class EulerAngleObservationPartialWrtRotationParameter : public ObservationPartial< 3 >
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param rotationMatrixPartial Object computing the partial derivatives of the rotation matrix from body-fixed to base frame
     * w.r.t. the parameter
     * \param rotationToBaseFrameFunction Function returning the rotation from body-fixed to base frame as a function of time
     * \param partialScalingFactor Factor by which the partials are multiplied (-1 for the state of the central body of a rotation
     * model that depends on the relative state of the observed body)
     */
    EulerAngleObservationPartialWrtRotationParameter( const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
                                                      const std::shared_ptr< RotationMatrixPartial > rotationMatrixPartial,
                                                      const std::function< Eigen::Quaterniond( const double ) > rotationToBaseFrameFunction,
                                                      const double partialScalingFactor = 1.0 ):
        ObservationPartial< 3 >( parameterIdentifier ), rotationMatrixPartial_( rotationMatrixPartial ),
        rotationToBaseFrameFunction_( rotationToBaseFrameFunction ), partialScalingFactor_( partialScalingFactor )
    {}

    //! Function to calculate the observation partial(s) at required time.
    /*!
     *  Function to calculate the observation partial(s) at required time.
     *  \param states Link end state. (Not used; vector contains NaN entry)
     *  \param times Link end time (size 1 std::vector).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param ancillarySettings Ancillary settings (not used)
     *  \param currentObservation Value of the observation (not used)
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, 3, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Vector3d& currentObservation = Eigen::Vector3d::Constant( TUDAT_NAN ) )
    {
        const double currentTime = times.at( 0 );
        const Eigen::Quaterniond rotationToBaseFrame = rotationToBaseFrameFunction_( currentTime );
        const Eigen::Matrix3d rotationMatrixToBaseFrame = rotationToBaseFrame.toRotationMatrix( );
        const Eigen::Matrix< double, 3, 4 > quaternionPartial =
                calculateEulerAngle313WrtRotationToBaseFrameQuaternionPartial( rotationToBaseFrame );
        const std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartial_->calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime );

        Eigen::Matrix< double, 3, Eigen::Dynamic > partial =
                Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, rotationMatrixPartials.size( ) );
        for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
        {
            // Rotation vector a of the change, from the skew-symmetric matrix dR R^T = [ a x ]
            const Eigen::Matrix3d rotationChange = rotationMatrixPartials.at( i ) * rotationMatrixToBaseFrame.transpose( );
            const Eigen::Vector3d rotationVector = 0.5 *
                    Eigen::Vector3d( rotationChange( 2, 1 ) - rotationChange( 1, 2 ),
                                     rotationChange( 0, 2 ) - rotationChange( 2, 0 ),
                                     rotationChange( 1, 0 ) - rotationChange( 0, 1 ) );
            const Eigen::Quaterniond quaternionChange =
                    Eigen::Quaterniond( 0.0, rotationVector.x( ), rotationVector.y( ), rotationVector.z( ) ) * rotationToBaseFrame;
            partial.col( i ) = 0.5 * partialScalingFactor_ * quaternionPartial *
                    Eigen::Vector4d( quaternionChange.w( ), quaternionChange.x( ), quaternionChange.y( ), quaternionChange.z( ) );
        }
        return { std::make_pair( partial, currentTime ) };
    }

protected:
    //! Object computing the partial derivatives of the rotation matrix from body-fixed to base frame w.r.t. the parameter
    std::shared_ptr< RotationMatrixPartial > rotationMatrixPartial_;

    //! Function returning the rotation from body-fixed to base frame as a function of time
    std::function< Eigen::Quaterniond( const double ) > rotationToBaseFrameFunction_;

    //! Factor by which the partials are multiplied
    double partialScalingFactor_;
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_EULERANGLEOBSERVABLEPARTIALS_H
