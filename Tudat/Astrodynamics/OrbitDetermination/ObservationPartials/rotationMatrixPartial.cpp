#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"

namespace tudat
{

namespace observation_partials
{


Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameWrtConstantRotationRate(
        const Eigen::Quaterniond intertialBodyFixedToIntegrationFrame,
        const double rotationRate, const double timeSinceEpoch )
{
    double currentRotationAngle = rotationRate * timeSinceEpoch;

    Eigen::Matrix3d rotationMatrixDerivative;
    double sineOfAngle = sin( currentRotationAngle );
    double cosineOfAngle = cos( currentRotationAngle );
    rotationMatrixDerivative << -sineOfAngle, -cosineOfAngle, 0.0, cosineOfAngle, - sineOfAngle, 0.0, 0.0, 0.0, 0.0;

    return timeSinceEpoch * ( intertialBodyFixedToIntegrationFrame.toRotationMatrix( ) ) * rotationMatrixDerivative;

}

std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation(
        const Eigen::Vector3d initialOrientationAngles,//right ascension, declination, longitude of prime meridian.
        const double rotationRate, const double timeSinceEpoch )
{
    Eigen::Matrix3d commonTerm = Eigen::AngleAxisd( -1.0 * ( -initialOrientationAngles.z( ) - rotationRate * timeSinceEpoch ),
                                                    Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );


    double rightAscension = initialOrientationAngles.x( );
    double declination = initialOrientationAngles.y( );

    Eigen::Matrix3d rightAscensionPartial;
    rightAscensionPartial<< -std::sin( rightAscension + mathematical_constants::PI / 2.0 ), -std::cos( rightAscension + mathematical_constants::PI / 2.0 ), 0.0,
            std::cos( rightAscension + mathematical_constants::PI / 2.0 ), -std::sin( rightAscension + mathematical_constants::PI / 2.0 ), 0.0,
            0.0, 0.0, 0.0;

    Eigen::Matrix3d declinationPartial;
    declinationPartial<< 0.0, 0.0, 0.0,
            0.0, std::sin( mathematical_constants::PI / 2.0 - declination ), std::cos( mathematical_constants::PI / 2.0 - declination ),
            0.0, -std::cos( mathematical_constants::PI / 2.0 - declination ), std::sin( mathematical_constants::PI / 2.0 - declination );

    std::vector< Eigen::Matrix3d > rotationMatrixPartials;
    rotationMatrixPartials.push_back( rightAscensionPartial * Eigen::Matrix3d(
                                          Eigen::AngleAxisd( -( declination - mathematical_constants::PI / 2.0 ), Eigen::Vector3d::UnitX( ) ) ) * commonTerm );
    rotationMatrixPartials.push_back( Eigen::AngleAxisd( ( rightAscension + mathematical_constants::PI / 2.0 ), Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( ) *
                                      declinationPartial * commonTerm );

    return rotationMatrixPartials;
}

Eigen::Matrix< double, 3, Eigen::Dynamic > RotationMatrixPartial::calculatePartialOfRotatedVector(
        const double time,
        const Eigen::Vector3d vectorInLocalFrame )
{
    std::vector< Eigen::Matrix3d > rotationMatrixPartials = calculatePartialOfRotationMatrixToBaseFrameWrParameter( time );
    Eigen::Matrix< double, 3, Eigen::Dynamic > rotatedVectorPartial = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero(
                3, rotationMatrixPartials.size( ) );

    for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
    {
        rotatedVectorPartial.block( 0, i, 3, 1 ) = rotationMatrixPartials[ i ] * vectorInLocalFrame;
    }
    return rotatedVectorPartial;
}

}

}
