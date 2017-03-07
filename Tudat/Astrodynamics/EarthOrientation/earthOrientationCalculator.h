#ifndef EARTHORIENTATIONCALCULATOR_H
#define EARTHORIENTATIONCALCULATOR_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"

#include "Tudat/Astrodynamics/EarthOrientation/earthSiderealTimeCalculator.h"
#include "Tudat/Astrodynamics/EarthOrientation/polarMotionCalculator.h"
#include "Tudat/Astrodynamics/EarthOrientation/precessionNutationCalculator.h"
#include "Tudat/Astrodynamics/EarthOrientation/eopReader.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"


namespace tudat
{

namespace earth_orientation
{

double getApproximateTioLocator( double secondsSinceJ2000 );

//! Calculate rotation from CIRS to GCRS, i.e. applying CIO-based rotations due to nutation and precession.
Eigen::Quaterniond calculateRotationFromCirsToGcrs( double X, double Y, double s );

//! Calculate rotation from TIRS to CIRS, i.e. rotating over earth rotation angle.
Eigen::Quaterniond calculateRotationFromTirsToCirs( double earthRotationAngle );

//! Calculate rotation from ITRS to TIRS, i.e. applying rotations due to polar motion.
Eigen::Quaterniond calculateRotationFromItrsToTirs( double xPole, double yPole, double tioLocator );

Eigen::Matrix3d calculateRotationRateFromItrsToGcrs( double X, double Y, double cioLocator,
                                                     double earthRotationAngle,
                                                     double xPole, double yPole,
                                                     double tioLocator );

Eigen::Matrix3d calculateRotationRateFromItrsToGcrs( Eigen::Vector6d rotationAngles,
                                                     double secondsSinceJ2000 );

//! Calculate rotation from ITRS to GCRS.
Eigen::Quaterniond calculateRotationFromItrsToGcrs( double X, double Y, double cioLocator,
                                                    double earthRotationAngle,
                                                    double xPole, double yPole,
                                                    double tioLocator );

//! Calculate rotation from ITRS to GCRS.
Eigen::Quaterniond calculateRotationFromItrsToGcrs( Eigen::Vector6d rotationAngles,
                                                    double secondsSinceJ2000 );

//! Class to calculate earth orientation angles, i.e. those used for transforming from ITRS to GCRS
class EarthOrientationAnglesCalculator
{
public:

    //! Constructor from objects calculating three sub-parts of earth orientation.
    /*!
     *  Constructor from objects calculating sub-parts of earth orientation and time conversion,
     *  i.e. polar motion, precession/nutation and conversion between TT,TDB,UTC and UT1.
     *  \param polarMotionCalculator Pointer to object for calcu2012lating position of pole in ITRS (polarm motion).
     *  \param precessionNutationCalculator Pointer to object for calculating position of pole in GCRS (precession/nutation).
     *  \param earthSiderealTimeCalculator Pointer to object to convert between different time scales.
     */
    EarthOrientationAnglesCalculator(
            const boost::shared_ptr< PolarMotionCalculator > polarMotionCalculator,
            const boost::shared_ptr< PrecessionNutationCalculator > precessionNutationCalculator,
            const boost::shared_ptr< EarthSiderealTimeCalculator > earthSiderealTimeCalculator ):
        polarMotionCalculator_( polarMotionCalculator ),
        precessionNutationCalculator_( precessionNutationCalculator ),
        earthSiderealTimeCalculator_( earthSiderealTimeCalculator ) { }

    //! Calculate rotation angles from ITRS to GCRS at given time value.
    /*!
     *  Calculate rotation angles from ITRS to GCRS at given time value. Any time scale combined with any time value
     *  can be used as input. TIO locator is not included in output as its value is minute and cal be easily evaluated, with
     *  no regard for time conversions in input value.
     *  \param timeValue Number of seconds since julianDayReferenceEpoch at which orientation is to be evaluated.
     *  \param timeScale Time scale in which the timeValue is given. To be taken from TimeScales enum.
     *  \param julianDayReferenceEpoch Number of julian days after JD=0 to be used as reference (i.e zero) epoch for timeValue.
     *  \return Rotation angles for ITRS<->GCRS transformation at given epoch. Order is: X, Y, s, ERA, x_p, y_p
     */
    Eigen::Vector6d getRotationAnglesFromItrsToGcrs(
            const double& timeValue,
            basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tt_scale,
            double julianDayReferenceEpoch =
            basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    //! Function to get object that calculates polar motion.
    /*!
     *  Function to get object that calculates polar motion.
     */
    boost::shared_ptr< PolarMotionCalculator > getPolarMotionCalculator( )
    { return polarMotionCalculator_; }

    //! Function to get object that calculates precession/nutation.
    /*!
     *  Function to get object that calculates precession/nutation.
     */
    boost::shared_ptr< PrecessionNutationCalculator > getPrecessionNutationCalculator( )
    { return precessionNutationCalculator_; }

    //! Function to get object that converts between time scales.
    /*!
     *  Function to get object that converts between time scales.
     */
    boost::shared_ptr< EarthSiderealTimeCalculator > getEarthSiderealTimeCalculator( )
    { return earthSiderealTimeCalculator_; }


private:
    //! Pointer to object for calculating position of pole in ITRS (polarm motion).
    /*!
     *  Pointer to object for calculating position of pole in ITRS (polarm motion).
     */
    boost::shared_ptr< PolarMotionCalculator > polarMotionCalculator_;

    //! Pointer to object for calculating position of pole in GCRS (precession/nutation).
    /*!
     *  Pointer to object for calculating position of pole in GCRS (precession/nutation).
     */
    boost::shared_ptr< PrecessionNutationCalculator > precessionNutationCalculator_;

    //! Pointer to object to convert between different time scales.
    /*!
     *  Pointer to object to convert between different time scales.
     */
    boost::shared_ptr< EarthSiderealTimeCalculator > earthSiderealTimeCalculator_;
};

boost::shared_ptr< EarthOrientationAnglesCalculator > createStandardEarthOrientationCalculator( );

double calculateUnnormalizedEarthRotationAngle( const double ut1SinceEpoch,
                                                const double referenceJulianDay );

boost::shared_ptr< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6,1 > > >
createInterpolatorForItrsToGcrsTransformation( double intervalStart,
                                               double intervalEnd,
                                               double timeStep,
                                               basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tt_scale,
                                               double julianDayReferenceEpoch =
        basic_astrodynamics::JULIAN_DAY_ON_J2000,
                                               boost::shared_ptr< EarthOrientationAnglesCalculator > earthOrientationCalculator =
        createStandardEarthOrientationCalculator( ) );

}

}

#endif // EARTHORIENTATIONCALCULATOR_H
