#ifndef POLARMOTIONCALCULATOR_H
#define POLARMOTIONCALCULATOR_H

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Astrodynamics/EarthOrientation/shortPeriodEarthOrientationCorrectionCalculator.h"

namespace tudat
{

namespace earth_orientation
{

//! Function to calculate the mean Celestial Intermediate Pole (CIP) position in the ITRS.
/*!
 *  Function to calculate the mean Celestial Intermediate Pole (CIP) position in the ITRS.
 *  The implementation is according to the IERS 2010 conventions, Eq. 7.25 and Table 7.7
 *  \param Number of Julian years since J2000. N.b. due to the very slow variation of the returned quantity,
 *  the time scale that is used for input is generally inconsequential for the output.
 *  \return Vector containing the pole positions, typically denoted as x_{p} and y_{p}
 */
Eigen::Vector2d calculateMeanCipPositionInItrs( const double julianYearsSinceJ2000 );

class PolarMotionCalculator
{
public:
    //! Constructor for polar motion (x_{p} and y_{p}; motion of pole in ITRS) calculator.
    /*!
     *  Constructor for polar motion (motion of pole in ITRS) calculator. This object combines the
     *  measured daily data as published by the IERS with the short-period variation models
     *  of the pole position due to libration and ocean tides as described in Sections 5.5.1 and 8.2.
     *  \param dailyIersValueInterpolator Interpolator, with time in UTC since J2000 as input and
     *  interpolated measured daily pole offset values as output.
     *  \param shortPeriodPolarMotionCalculator Object calculating short period polar motion variations.
     */
    PolarMotionCalculator( const boost::shared_ptr< interpolators::OneDimensionalInterpolator
                           < double, Eigen::Vector2d > > dailyIersValueInterpolator,
                           const boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >
                           shortPeriodPolarMotionCalculator ):
        dailyIersValueInterpolator_( dailyIersValueInterpolator ),
        shortPeriodPolarMotionCalculator_( shortPeriodPolarMotionCalculator )
    { }

    //! Calculate the position of the Celestial Intermediate Pole in the ITRS
    /*!
     *  Calculate the position of the Celestial Intermediate Pole in the ITRS, i.e. x_{p} and y_{p}.
     *  The doodson arguments are calculated internally in this function. If these are available to
     *  the user from a more computationally efficient solution (such as an interpolator), the
     *  overloaded version of this function is adviced for use.
     *  \param ttSinceEpoch Terrestrial time since the epochShift variable.
     *  \param utcSinceEpoch UTC since the epochShift variable.
     *  \param epochShift. Shift in epoch t=0 in days from julain day 0.
     */
    Eigen::Vector2d getPositionOfCipInItrs(
            const double ttSinceEpoch,
            const double utcSinceEpoch,
            const double epochShift = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    Eigen::Vector2d getPositionOfCipInItrs(
            Eigen::Vector6d doodsonArguments,
            const double utcSinceEpoch,
            const double epochShift = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    Eigen::Vector2d getPolarMotionWobbleVariables(
            Eigen::Vector6d doodsonArguments,
            const double utcSinceEpoch,
            const double epochShift = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    Eigen::Vector2d getPolarMotionWobbleVariables(
            const double ttSinceEpoch,
            const double utcSinceEpoch,
            const double epochShift );

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > >
        getDailyIersValueInterpolator( )
    {
        return dailyIersValueInterpolator_;
    }

    boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > getShortPeriodPolarMotionCalculator( )
    {
        return shortPeriodPolarMotionCalculator_;
    }

    void resetDailyValueInterpolator( boost::shared_ptr< interpolators::OneDimensionalInterpolator
                                      < double, Eigen::Vector2d > > dailyValueInterpolator )
    {
        dailyIersValueInterpolator_ = dailyValueInterpolator;
    }

private:
    //! Interpolator for daily IERS-measured pole offsets.
    /*!
     *  Interpolator, with time in UTC since J2000 as input and
     *  interpolated measured daily pole offset values as output.
     */
    boost::shared_ptr< interpolators::OneDimensionalInterpolator
        < double, Eigen::Vector2d > > dailyIersValueInterpolator_;

    //! Object calculating short period polar motion variations.
    /*!
     *  Object calculating short period polar motion variations.
     */
    boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator_;
};

}

}

#endif // POLARMOTIONCALCULATOR_H
