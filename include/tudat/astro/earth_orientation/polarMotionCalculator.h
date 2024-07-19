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

#ifndef TUDAT_POLARMOTIONCALCULATOR_H
#define TUDAT_POLARMOTIONCALCULATOR_H

#include <functional>

#include <Eigen/Core>

#include "tudat/math/interpolators/oneDimensionalInterpolator.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"

namespace tudat
{

namespace earth_orientation
{


Eigen::Vector2d getSecularPolePositionInMas( const double time );

//! Object to compute polar motion variables x_{p} and y_{p}
/*!
 * Object to compute polar motion variables x_{p} and y_{p}, describing the motion of the Earth pole in ITRS.
 * This class combines the  measured daily data as published by the IERS with the short-period variation models due to libration
 * and ocean tides as described in Sections 5.5.1 and 8.2 of IERS 2010 Conventions.
 */
class PolarMotionCalculator
{
public:
    //! Constructor
    /*!
     *  Constructor
     *  \param dailyIersValueInterpolator Interpolator, with time in UTC since J2000 as input and
     *  interpolated measured daily pole offset values as output.
     *  \param shortPeriodPolarMotionCalculator Object calculating short period polar motion variations.
     */
    PolarMotionCalculator(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator < double, Eigen::Vector2d > >
            dailyIersValueInterpolator,
            const std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >
            shortPeriodPolarMotionCalculator ):
        dailyIersValueInterpolator_( dailyIersValueInterpolator ),
        shortPeriodPolarMotionCalculator_( shortPeriodPolarMotionCalculator )
    { }

    //! Calculate the position of the Celestial Intermediate Pole in the ITRS
    /*!
     *  Calculate the position of the Celestial Intermediate Pole in the ITRS, i.e. x_{p} and y_{p}.
     *  The fundamental arguments are calculated internally in this function. If these are available to
     *  the user from a more computationally efficient solution (such as an interpolator), the
     *  overloaded version of this function is adviced for use.
     *  \param ttSinceEpoch Terrestrial time since J2000
     *  \param utcSinceEpoch UTC since the J2000
     */
    Eigen::Vector2d getPositionOfCipInItrs(
            const double ttSinceEpoch,
            const double utcSinceEpoch );

    //! Calculate the position of the Celestial Intermediate Pole in the ITRS
    /*!
     *  Calculate the position of the Celestial Intermediate Pole in the ITRS, i.e. x_{p} and y_{p}.
     *  The fundamental arguments are provided as input to this function for computation of short-period variations.
     *  \param fundamentalArguments Fundamental arguments used for short-period polar motion corrections
     *  \param utcSinceEpoch UTC since the J2000
     */
    Eigen::Vector2d getPositionOfCipInItrs(
            Eigen::Vector6d fundamentalArguments,
            const double utcSinceEpoch);

    //! Function to retrieve interpolator for daily IERS-measured pole offsets
    /*!
     * Function to retrieve interpolator for daily IERS-measured pole offsets
     * \return Interpolator for daily IERS-measured pole offsets
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > >
    getDailyIersValueInterpolator( )
    {
        return dailyIersValueInterpolator_;
    }

    //! Function to retrieve object calculating short period polar motion variations.
    /*!
     * Function to retrieve object calculating short period polar motion variations.
     * \return Object calculating short period polar motion variations.
     */
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > getShortPeriodPolarMotionCalculator( )
    {
        return shortPeriodPolarMotionCalculator_;
    }

private:
    //! Interpolator for daily IERS-measured pole offsets.
    /*!
     *  Interpolator, with time in UTC since J2000 as input and interpolated measured daily pole offset values as output.
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator
    < double, Eigen::Vector2d > > dailyIersValueInterpolator_;

    //! Object calculating short period polar motion variations.
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator_;
};

}

}

#endif // TUDAT_POLARMOTIONCALCULATOR_H
