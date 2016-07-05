/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110519    F.M. Engelen      File created.
 *      110628    K. Kumar          Minor comment and layout changes; changed
 *                                  input arguments to pass-by-reference.
 *      110701    K. Kumar          Updated file path.
 *      110718    F.M. Engelen      Repaired incorrect sign for angle in R2I and I2R
 *                                  Added Quaternion transformations and ItoE tranformation.
 *      110726    K. Kumar          Minor modifications.
 *      110809    F.M. Engelen      Applied the minus one correction for angleAxisD,
 *                                  changed to local vertical frame.
 *      120530    E.A.G. Heeren     Namespace update.
 *      130121    K. Kumar          Updated functions to be const-correct.
 *      130219    D. Dirkx          Migrated from personal code.
 *      130312    A. Ronse          Added V-T, TA-AA and AA-B transformations.
 *
 *    References
 *      Mooij, E. The Motion of a Vehicle in a Planetary Atmosphere, TU Delft, 1997.
 *      Seidelmann, P. K. (Ed.). (2005). Explanatory supplement to the astronomical almanac.
 *              Univ Science Books.
 *
 *    Notes
 *      Because potential speed improvement it was chosen to use AngleAxisd and quaternions
 *      but to get things working, the rotation angle inputted in angleAxisd need to be inverted.
 *      In the future it might be better to change it to write out the complete transformation for
 *      clarity, or work with directional cosine matrices.
 *
 */

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{
namespace reference_frames
{

//! Wrapper function to transform a vector to a different frame from a single rotation function.
Eigen::Vector3d transformVector(
        const Eigen::Vector3d& originalVector,
        const boost::function< Eigen::Quaterniond( ) > rotation )
{
    return rotation( ) * originalVector;
}

//! Wrapper function to transform a vector to a different frame from a single transformation function.
Eigen::Vector3d transformVector(
        const boost::function< Eigen::Vector3d( ) > originalVector,
        const boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > transformationFunction )
{
    return transformationFunction( originalVector( ) );
}

//! Wrapper function to transform a vector to a different frame from a list of transformation function.
Eigen::Vector3d transformVector(
        const Eigen::Vector3d& originalVector,
        const std::vector< boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > >& rotationsList )
{
    Eigen::Vector3d currentVector = originalVector;
    Eigen::Vector3d newVector;

    // Apply each of the required tranformations.
    for( unsigned int i = 0; i < rotationsList.size( ); i++ )
    {
        newVector = rotationsList.at( i )( currentVector );
        currentVector = newVector;
    }
    return currentVector;
}

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation matrix.
Eigen::Matrix3d
getRotatingPlanetocentricToInertialFrameTransformationMatrix( const double angleFromXItoXR )
{
    // Declare local variables.
    // Declare local matrix.
    Eigen::Matrix3d localMatrix_;

    // Set local matrix.
    localMatrix_ = reference_frames::
            getInertialToPlanetocentricFrameTransformationMatrix( angleFromXItoXR );

    // Return transformation matrix.
    return localMatrix_.transpose( );
}

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation quaternion.
Eigen::Quaterniond getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
        const double angleFromXItoXR )
{
    // Compute transformation quaternion
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd eigenRotationObject = Eigen::AngleAxisd( -1.0 * -angleFromXItoXR,
                                                               Eigen::Vector3d::UnitZ( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond( eigenRotationObject );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get rotation from planet-fixed to inertial frame.
Eigen::Quaterniond getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd firstRotationAroundZaxis =
            Eigen::AngleAxisd( longitudeOfPrimeMeridian, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd rotationAroundXaxis =
            Eigen::AngleAxisd(
                ( mathematical_constants::PI / 2.0 - declinationOfPole ),
                Eigen::Vector3d::UnitX( ) );
    Eigen::AngleAxisd secondRotationAroundZaxis = Eigen::AngleAxisd(
                rightAscensionOfPole
                + mathematical_constants::PI / 2.0, Eigen::Vector3d::UnitZ( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( secondRotationAroundZaxis * rotationAroundXaxis *  firstRotationAroundZaxis ) );
    return frameTransformationQuaternion;
}

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion matrix.
Eigen::Matrix3d getInertialToPlanetocentricFrameTransformationMatrix(
        const double angleFromXItoXR )
{
    // Compute rotation about Z-Axis.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd eigenRotationObject = Eigen::AngleAxisd( -1.0 * angleFromXItoXR,
                                                               Eigen::Vector3d::UnitZ( ) );

    // Return transformation matrix.
    return eigenRotationObject.toRotationMatrix( );
}

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion quaternion.
Eigen::Quaterniond getInertialToPlanetocentricFrameTransformationQuaternion(
        const double angleFromXItoXR )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd eigenRotationObject = Eigen::AngleAxisd( -1.0 * angleFromXItoXR,
                                                               Eigen::Vector3d::UnitZ( ) );

    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond( eigenRotationObject );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get rotation from inertial to planet-fixed frame.
Eigen::Quaterniond getInertialToPlanetocentricFrameTransformationQuaternion(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd secondRotationAroundZaxis =
            Eigen::AngleAxisd( -longitudeOfPrimeMeridian, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd rotationAroundXaxis =
            Eigen::AngleAxisd(
                -( mathematical_constants::PI / 2.0 - declinationOfPole ),
                Eigen::Vector3d::UnitX( ) );
    Eigen::AngleAxisd firstRotationAroundZaxis = Eigen::AngleAxisd(
                - ( rightAscensionOfPole + mathematical_constants::PI / 2.0 ),
                Eigen::Vector3d::UnitZ( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( secondRotationAroundZaxis * rotationAroundXaxis *  firstRotationAroundZaxis ) );
    return frameTransformationQuaternion;
}

//! Create a Quaterniond rotation state object from four quaternion values in a Vector4d
Eigen::Quaterniond getQuaternionObjectFromQuaternionValues(
        const Eigen::Vector4d& vectorWithQuaternion )
{
    // Set transformation quaternion.
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                vectorWithQuaternion( 0 ), vectorWithQuaternion( 1 ),
                vectorWithQuaternion( 2 ), vectorWithQuaternion( 3 ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get transformation matrix from Planetocentric (R) to the Local vertical (V) frame.
Eigen::Matrix3d getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(
    const double longitude, const double latitude )
{
    return getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
            longitude, latitude ).toRotationMatrix( );
}

//! Get transformation quaternion from Planetocentric (R) to the Local vertical (V) frame.
Eigen::Quaterniond getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
        const double longitude, const double latitude )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd RotationAroundZaxis = Eigen::AngleAxisd(
                -1.0 * longitude, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd RotationAroundYaxis = Eigen::AngleAxisd(
                -1.0 * ( -latitude - mathematical_constants::PI / 2.0 ),
                Eigen::Vector3d::UnitY( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( RotationAroundYaxis * RotationAroundZaxis ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get transformation matrix from local vertical (V) to the Planetocentric frame (R).
Eigen::Matrix3d getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(
    const double longitude, const double latitude )
{
    return getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(
            longitude, latitude ).transpose( );
}

//! Get transformation quaternion from local vertical (V) to the Planetocentric frame (R).
Eigen::Quaterniond getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
        const double longitude, const double latitude )
{
    return getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
            longitude, latitude ).inverse( );
}

//! Get transformation matrix from the TA/TG to the V-frame.
Eigen::Matrix3d getTrajectoryToLocalVerticalFrameTransformationMatrix(
        const double flightPathAngle, const double headingAngle )
{
    return getTrajectoryToLocalVerticalFrameTransformationQuaternion(
            flightPathAngle, headingAngle ).toRotationMatrix( );
}

//! Get transformation quaternion from the TA/TG to the V-frame.
Eigen::Quaterniond getTrajectoryToLocalVerticalFrameTransformationQuaternion(
        const double flightPathAngle, const double headingAngle )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd rotationAroundZaxis = Eigen::AngleAxisd(
                -1.0 * -headingAngle, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd rotationAroundYaxis = Eigen::AngleAxisd(
                -1.0 * -flightPathAngle, Eigen::Vector3d::UnitY( ) );

    return Eigen::Quaterniond( rotationAroundZaxis * rotationAroundYaxis );
}

//! Get transformation matrix from the local V- to TA/TG-frame.
Eigen::Matrix3d getLocalVerticalFrameToTrajectoryTransformationMatrix(
        const double flightPathAngle, const double headingAngle )
{
    return getTrajectoryToLocalVerticalFrameTransformationMatrix(
            flightPathAngle, headingAngle ).transpose( );
}

//! Get transformation quaternion from V- to the TA/TG-frame.
Eigen::Quaterniond getLocalVerticalFrameToTrajectoryTransformationQuaternion(
        const double flightPathAngle, const double headingAngle )
{
    return getTrajectoryToLocalVerticalFrameTransformationQuaternion(
            flightPathAngle, headingAngle ).inverse( );
}

//! Get transformation matrix from the TA- to the AA-frame.
Eigen::Matrix3d getTrajectoryToAerodynamicFrameTransformationMatrix(
        const double bankAngle )
{
    return getTrajectoryToAerodynamicFrameTransformationQuaternion(
            bankAngle ).toRotationMatrix( );
}

//! Get transformation quaternion from the TA- to the AA-frame.
Eigen::Quaterniond getTrajectoryToAerodynamicFrameTransformationQuaternion(
        const double bankAngle )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd rotationAroundXaxis 
        = Eigen::AngleAxisd( bankAngle, Eigen::Vector3d::UnitX( ) );
    return Eigen::Quaterniond( rotationAroundXaxis );
}

//! Get transformation matrix from the AA- to the TA-frame.
Eigen::Matrix3d getAerodynamicToTrajectoryFrameTransformationMatrix(
        const double bankAngle )
{
    return getTrajectoryToAerodynamicFrameTransformationMatrix( bankAngle ).transpose( );
}

//! Get transformation quaternion from the AA- to the TA-frame.
Eigen::Quaterniond getAerodynamicToTrajectoryFrameTransformationQuaternion(
        const double bankAngle )
{
    return getTrajectoryToAerodynamicFrameTransformationQuaternion( bankAngle ).inverse( );
}

//! Get transformation matrix fom the B- to the AA-frame.
Eigen::Matrix3d getBodyToAirspeedBasedAerodynamicFrameTransformationMatrix(
        const double angleOfAttack, const double angleOfSideslip )
{
    return getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
            angleOfAttack, angleOfSideslip ).toRotationMatrix( );
}

//! Get transformation quaternion fom the B- to the AA-frame.
Eigen::Quaterniond getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
        const double angleOfAttack, const double angleOfSideslip )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd rotationAroundZaxis 
        = Eigen::AngleAxisd( -1.0 * angleOfSideslip, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd rotationAroundYaxis 
        = Eigen::AngleAxisd( -1.0 * -angleOfAttack, Eigen::Vector3d::UnitY( ) );

    return Eigen::Quaterniond( rotationAroundZaxis * rotationAroundYaxis );
}

//! Get transformation matrix fom the AA- to the B-frame.
Eigen::Matrix3d getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
        const double angleOfAttack, const double angleOfSideslip )
{
    return getBodyToAirspeedBasedAerodynamicFrameTransformationMatrix(
            angleOfAttack, angleOfSideslip ).transpose( );
}

//! Get transformation quaternion fom the AA- to the B-frame.
Eigen::Quaterniond getAirspeedBasedAerodynamicToBodyFrameTransformationQuaternion(
        const double angleOfAttack, const double angleOfSideslip )
{
    return getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
            angleOfAttack, angleOfSideslip ).inverse( );
}

//! Calculate current heading angle.
double calculateHeadingAngle( const Eigen::Vector3d& velocityInVerticalFrame )
{
    return std::atan2( velocityInVerticalFrame( 1 ), velocityInVerticalFrame( 0 ) );
}

//! Calculatre current flight path angle.
double calculateFlightPathAngle( const Eigen::Vector3d& velocityInVerticalFrame )
{
    return -std::asin( velocityInVerticalFrame( 2 ) / velocityInVerticalFrame.norm( ) );
}


} // namespace reference_frames
} // namespace tudat
