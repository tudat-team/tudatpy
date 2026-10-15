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
#include "tudat/astro/earth_orientation/earthOrientationDataEvaluator.h"

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
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > > dailyIersValueInterpolator,
            const std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator,
            const std::shared_ptr< EarthOrientationDataEvaluator > evaluator = nullptr,
            const std::function< Time( const Time& ) > ut1FromTai = {} ):
        evaluator_( evaluator ), ut1FromTai_( ut1FromTai ), dailyIersValueInterpolator_( dailyIersValueInterpolator ),
        shortPeriodPolarMotionCalculator_( shortPeriodPolarMotionCalculator )
    {}

    //! Calculate the position of the Celestial Intermediate Pole in the ITRS at a TT epoch
    /*!
     *  Calculate the position of the Celestial Intermediate Pole in the ITRS, i.e. x_{p} and y_{p}, at a TT epoch.
     *  Models constructed from EOP data evaluate the daily values at TAI, which is computed from TT without using UTC.
     *  \param tt Terrestrial time since J2000
     */
    template< typename TimeType >
    Eigen::Vector2d getPositionOfCipInItrsAtTt( const TimeType& tt )
    {
        // SOFA evaluates its harmonic arguments using scalar TT seconds.
        const double ttDouble = static_cast< double >( tt );

        if( !evaluator_ )
        {
            return getPositionOfCipInItrs( ttDouble, sofa_interface::convertTTtoUTC( ttDouble ) );
        }
        const Time tai = basic_astrodynamics::convertTTtoTAI( Time( tt ) );

        // Use the selected UT1 model when forming the short-period arguments.
        const Time ut1 = ut1FromTai_ ? ut1FromTai_( tai ) : tai + evaluator_->evaluate( tai ).ut1MinusTai;

        // Delaunay arguments include GMST formed from UT1, with TT for precession.
        const Eigen::Vector6d arguments =
                sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst( ttDouble, ttDouble, static_cast< double >( ut1 ) );

        return getPositionOfCipInItrsFromTai( tai, arguments );
    }

    //! Calculate the position of the Celestial Intermediate Pole in the ITRS from TAI
    /*!
     *  Calculate the position of the Celestial Intermediate Pole in the ITRS, i.e. x_{p} and y_{p}: the daily values from the
     *  EOP data at the given TAI, plus the short-period corrections for the given fundamental arguments. Only available for
     *  models constructed from EOP data.
     *  \param tai TAI since J2000
     *  \param arguments Fundamental arguments used for short-period polar motion corrections
     */
    Eigen::Vector2d getPositionOfCipInItrsFromTai( const Time& tai, const Eigen::Vector6d& arguments );

    //! Function to retrieve the interpolated EOP data (nullptr for models constructed from a daily interpolator)
    std::shared_ptr< EarthOrientationDataEvaluator > getEarthOrientationDataEvaluator( ) const
    {
        return evaluator_;
    }

    //! Calculate the position of the Celestial Intermediate Pole in the ITRS
    /*!
     *  Calculate the position of the Celestial Intermediate Pole in the ITRS, i.e. x_{p} and y_{p}.
     *  The fundamental arguments are calculated internally in this function. If these are available to
     *  the user from a more computationally efficient solution (such as an interpolator), the
     *  overloaded version of this function is adviced for use.
     *  \param ttSinceEpoch Terrestrial time since J2000
     *  \param utcSinceEpoch UTC since the J2000
     */
    Eigen::Vector2d getPositionOfCipInItrs( const double ttSinceEpoch, const double utcSinceEpoch );

    //! Calculate the position of the Celestial Intermediate Pole in the ITRS
    /*!
     *  Calculate the position of the Celestial Intermediate Pole in the ITRS, i.e. x_{p} and y_{p}.
     *  The fundamental arguments are provided as input to this function for computation of short-period variations.
     *  \param fundamentalArguments Fundamental arguments used for short-period polar motion corrections
     *  \param utcSinceEpoch UTC since the J2000
     */
    Eigen::Vector2d getPositionOfCipInItrs( Eigen::Vector6d fundamentalArguments, const double utcSinceEpoch );

    //! Function to retrieve interpolator for daily IERS-measured pole offsets
    /*!
     * Function to retrieve interpolator for daily IERS-measured pole offsets
     * \return Interpolator for daily IERS-measured pole offsets
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > > getDailyIersValueInterpolator( )
    {
        if( evaluator_ )
        {
            throw std::runtime_error(
                    "Models constructed from EOP data do not have a daily UTC interpolator; use "
                    "getEarthOrientationDataEvaluator()->evaluate(tai).polarMotion." );
        }
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
    //! Interpolated EOP data; nullptr for models constructed from a user-supplied daily interpolator.
    std::shared_ptr< EarthOrientationDataEvaluator > evaluator_;

    //! Function returning UT1 at a given TAI, used for the fundamental arguments of the short-period corrections. If empty,
    //! UT1 is computed from the daily EOP data only.
    std::function< Time( const Time& ) > ut1FromTai_;

    //! Interpolator for daily IERS-measured pole offsets.
    /*!
     *  Interpolator, with time in UTC since J2000 as input and interpolated measured daily pole offset values as output.
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > > dailyIersValueInterpolator_;

    //! Object calculating short period polar motion variations.
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator_;
};

}  // namespace earth_orientation

}  // namespace tudat

#endif  // TUDAT_POLARMOTIONCALCULATOR_H
