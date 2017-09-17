#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"

namespace tudat
{

namespace earth_orientation
{

double getApproximateTioLocator( double secondsSinceJ2000 )
{
    return iauSp00( basic_astrodynamics::JULIAN_DAY_ON_J2000, secondsSinceJ2000 / physical_constants::JULIAN_DAY );
}

Eigen::Quaterniond calculateRotationFromCirsToGcrs( double X, double Y, double s )
{
    double X2 = X * X;
    double Y2 = Y * Y;
    double XY = X * Y;

    double a = 0.5 + 0.125 * ( X2 + Y2 );
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix << 1.0 - a * X2, -a * XY, X, - a * XY, 1.0 - a * Y2, Y, -X, -Y, 1.0 - a * ( X2 + Y2 );
    return Eigen::Quaterniond( rotationMatrix ) * Eigen::Quaterniond( Eigen::AngleAxisd( -s, Eigen::Vector3d::UnitZ( ) ) );
}

Eigen::Quaterniond calculateRotationFromTirsToCirs( double earthRotationAngle )
{
    return Eigen::Quaterniond( Eigen::AngleAxisd( earthRotationAngle, Eigen::Vector3d::UnitZ( ) ) );
}

Eigen::Quaterniond calculateRotationFromItrsToTirs( double xPole, double yPole, double tioLocator )
{
    return Eigen::Quaterniond( Eigen::AngleAxisd( tioLocator, Eigen::Vector3d::UnitZ( ) ) *
                               Eigen::AngleAxisd( -xPole, Eigen::Vector3d::UnitY( ) ) *
                               Eigen::AngleAxisd( -yPole, Eigen::Vector3d::UnitX( ) ) );
}

Eigen::Matrix3d calculateRotationRateFromItrsToGcrs( double X, double Y, double cioLocator,
                                                     double earthRotationAngle,
                                                     double xPole, double yPole,
                                                     double tioLocator )
{
    Eigen::Matrix3d auxiliaryMatrix = Eigen::Matrix3d::Zero( );
    auxiliaryMatrix( 1, 0 ) = 1.0;
    auxiliaryMatrix( 0, 1 ) = -1.0;
    auxiliaryMatrix *= 2.0 * mathematical_constants::PI / 86400.0 * 1.002737909350795;
    //std::cout<<X<<" "<<Y<<" "<<cioLocator<<" "<<earthRotationAngle<<" "<<xPole<<" "<<yPole<<" "<<tioLocator<<std::endl;
    return  ( calculateRotationFromCirsToGcrs( X, Y, cioLocator )  *
              calculateRotationFromTirsToCirs( earthRotationAngle ) ).toRotationMatrix( ) * auxiliaryMatrix *
            calculateRotationFromItrsToTirs( xPole, yPole, tioLocator ).toRotationMatrix( );

}

Eigen::Matrix3d calculateRotationRateFromItrsToGcrs( Eigen::Vector6d rotationAngles,
                                                     double secondsSinceJ2000 )
{
    //std::cout<<secondsSinceJ2000<<" "<<rotationAngles.transpose( )<<std::endl;
    return calculateRotationRateFromItrsToGcrs(
                rotationAngles[ 0 ], rotationAngles[ 1 ], rotationAngles[ 2 ],
                rotationAngles[ 3 ], rotationAngles[ 4 ], rotationAngles[ 5 ],
                getApproximateTioLocator( secondsSinceJ2000 ) );
}

Eigen::Quaterniond calculateRotationFromItrsToGcrs( double X, double Y, double cioLocator,
                                                    double earthRotationAngle,
                                                    double xPole, double yPole,
                                                    double tioLocator )
{
    //std::cout<<X<<" "<<Y<<" "<<cioLocator<<" "<<earthRotationAngle<<" "<<xPole<<" "<<yPole<<" "<<tioLocator<<std::endl;
    return  calculateRotationFromCirsToGcrs( X, Y, cioLocator ) *
            calculateRotationFromTirsToCirs( earthRotationAngle ) *
            calculateRotationFromItrsToTirs( xPole, yPole, tioLocator );

}

Eigen::Quaterniond calculateRotationFromItrsToGcrs( Eigen::Vector6d rotationAngles,
                                                    double secondsSinceJ2000 )
{
    //std::cout<<secondsSinceJ2000<<" "<<rotationAngles.transpose( )<<std::endl;
    return calculateRotationFromItrsToGcrs(
                rotationAngles[ 0 ], rotationAngles[ 1 ], rotationAngles[ 2 ],
                rotationAngles[ 3 ], rotationAngles[ 4 ], rotationAngles[ 5 ],
                getApproximateTioLocator( secondsSinceJ2000 ) );
}


Eigen::Vector6d EarthOrientationAnglesCalculator::getRotationAnglesFromItrsToGcrs(
        const double& timeValue, basic_astrodynamics::TimeScales timeScale, double julianDayReferenceEpoch )
{
    double terrestrialTime = terrestrialTimeScaleConverter_->getCurrentTime(
                timeScale, basic_astrodynamics::tt_scale, timeValue, Eigen::Vector3d::Zero( ) );

    double utc = terrestrialTimeScaleConverter_->getCurrentTime(
                timeScale, basic_astrodynamics::utc_scale, timeValue, Eigen::Vector3d::Zero( ) );

    double ut1 = terrestrialTimeScaleConverter_->getCurrentTime(
                timeScale, basic_astrodynamics::ut1_scale, timeValue, Eigen::Vector3d::Zero( ) );

    std::pair< Eigen::Vector2d, double > positionOfCipInGcrs =
            precessionNutationCalculator_->getPositionOfCipInGcrs(
                terrestrialTime, utc, julianDayReferenceEpoch );

    Eigen::Vector2d positionOfCipInItrs = polarMotionCalculator_->getPositionOfCipInItrs(
                terrestrialTime, utc, julianDayReferenceEpoch );

    double earthRotationAngle = calculateUnnormalizedEarthRotationAngle( ut1, julianDayReferenceEpoch );

    Eigen::Vector6d rotationAngles;
    rotationAngles[ 0 ] = positionOfCipInGcrs.first.x( );
    rotationAngles[ 1 ] = positionOfCipInGcrs.first.y( );
    rotationAngles[ 2 ] = positionOfCipInGcrs.second;
    rotationAngles[ 3 ] = earthRotationAngle;
    rotationAngles[ 4 ] = positionOfCipInItrs.x( );
    rotationAngles[ 5 ] = positionOfCipInItrs.y( );
    return rotationAngles;
}



boost::shared_ptr< EarthOrientationAnglesCalculator > createStandardEarthOrientationCalculator( )
{
    using namespace interpolators;
    boost::shared_ptr< EOPReader > eopReader = boost::make_shared< EOPReader >( );

    boost::shared_ptr< LinearInterpolator< double, Eigen::Vector2d > > cipInItrsInterpolator =
            boost::make_shared< LinearInterpolator< double, Eigen::Vector2d > >(
                eopReader->getCipInItrsMapInSecondsSinceJ2000( ) ); // xPole, yPole

    boost::shared_ptr< LinearInterpolator< double, Eigen::Vector2d > > cipInGcrsCorrectionInterpolator =
            boost::make_shared< LinearInterpolator< double, Eigen::Vector2d > >(
                eopReader->getCipInGcrsCorrectionMapInSecondsSinceJ2000( ) ); // dX, dY

    boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator =
            getDefaultPolarMotionCorrectionCalculator( );

    boost::shared_ptr< PolarMotionCalculator > polarMotionCalculator = boost::make_shared< PolarMotionCalculator >
            ( cipInItrsInterpolator, shortPeriodPolarMotionCalculator );

    boost::shared_ptr< PrecessionNutationCalculator > precessionNutationCalculator =
            boost::make_shared< PrecessionNutationCalculator >( iau_2006, cipInGcrsCorrectionInterpolator );

    boost::shared_ptr< TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter = createDefaultTimeConverter( eopReader );

    return boost::make_shared< EarthOrientationAnglesCalculator >(
                polarMotionCalculator, precessionNutationCalculator, terrestrialTimeScaleConverter );
}

double calculateUnnormalizedEarthRotationAngle( const double ut1SinceEpoch,
                                                const double referenceJulianDay )
{
    return 2.0 * mathematical_constants::PI * ( 0.7790572732640 + 1.00273781191135448 * ( ut1SinceEpoch / physical_constants::JULIAN_DAY +
                                                                               ( referenceJulianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) );
}

boost::shared_ptr< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6,1 > > >
createInterpolatorForItrsToGcrsTransformation( double intervalStart,
                                               double intervalEnd,
                                               double timeStep,
                                               basic_astrodynamics::TimeScales timeScale,
                                               double julianDayReferenceEpoch,
                                               boost::shared_ptr< EarthOrientationAnglesCalculator > earthOrientationCalculator )
{
    using namespace interpolators;



    std::map< double, Eigen::Matrix< double, 6,1 > > orientationMap;

    double currentTime = intervalStart;
    while( currentTime < intervalEnd )
    {
        orientationMap[ currentTime ] = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs(
                    currentTime, timeScale, julianDayReferenceEpoch );
        currentTime += timeStep;
    }

    boost::shared_ptr< LagrangeInterpolator< double, Eigen::Matrix< double, 6,1 > > > interpolator =
            boost::make_shared< LagrangeInterpolator< double, Eigen::Matrix< double, 6,1 > > >( orientationMap, 6 );
    return interpolator;
}


}

}
