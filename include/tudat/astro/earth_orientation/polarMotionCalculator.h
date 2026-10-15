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
#include <stdexcept>

#include <Eigen/Core>

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
 * This class combines the measured daily data as published by the IERS (interpolated by an EarthOrientationDataEvaluator)
 * with the short-period variation models due to libration and ocean tides as described in Sections 5.5.1 and 8.2 of IERS
 * 2010 Conventions.
 */
class PolarMotionCalculator
{
public:
    //! Constructor
    /*!
     *  Constructor
     *  \param eopData Interpolated daily Earth orientation data, from which the measured pole position is taken.
     *  \param shortPeriodPolarMotionCalculator Object calculating short period polar motion variations (none if nullptr).
     *  \param ut1FromTai Function returning UT1 at a given TAI, used for the fundamental arguments of the short-period
     *  variations. If empty, UT1 is computed from the daily data only.
     */
    PolarMotionCalculator(
            const std::shared_ptr< EarthOrientationDataEvaluator > eopData,
            const std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator,
            const std::function< Time( const Time& ) > ut1FromTai = {} ):
        eopData_( eopData ), ut1FromTai_( ut1FromTai ), shortPeriodPolarMotionCalculator_( shortPeriodPolarMotionCalculator )
    {
        if( eopData_ == nullptr )
        {
            throw std::invalid_argument( "Polar motion calculator requires Earth orientation data" );
        }
    }

    //! Calculate the position of the Celestial Intermediate Pole in the ITRS at a TT epoch
    /*!
     *  Calculate the position of the Celestial Intermediate Pole in the ITRS, i.e. x_{p} and y_{p}, at a TT epoch.
     *  \param tt Terrestrial time since J2000
     *  \return Position of the Celestial Intermediate Pole in the ITRS (radians)
     */
    template< typename TimeType >
    Eigen::Vector2d getPositionOfCipInItrsAtTt( const TimeType& tt )
    {
        // SOFA evaluates its harmonic arguments using scalar TT seconds.
        const double ttDouble = static_cast< double >( tt );

        const Time tai = basic_astrodynamics::convertTTtoTAI( Time( tt ) );

        // UT1 for the short-period arguments.
        const Time ut1 = ut1FromTai_ ? ut1FromTai_( tai ) : tai + eopData_->evaluate( tai ).ut1MinusTai;

        // Delaunay arguments include GMST formed from UT1, with TT for precession.
        const Eigen::Vector6d arguments =
                sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst( ttDouble, ttDouble, static_cast< double >( ut1 ) );

        return getPositionOfCipInItrsFromTai( tai, arguments );
    }

    //! Calculate the position of the Celestial Intermediate Pole in the ITRS from TAI
    /*!
     *  Calculate the position of the Celestial Intermediate Pole in the ITRS, i.e. x_{p} and y_{p}: the daily values from the
     *  Earth orientation data at the given TAI, plus the short-period variations for the given fundamental arguments.
     *  \param tai TAI since J2000
     *  \param arguments Fundamental arguments used for short-period polar motion variations
     *  \return Position of the Celestial Intermediate Pole in the ITRS (radians)
     */
    Eigen::Vector2d getPositionOfCipInItrsFromTai( const Time& tai, const Eigen::Vector6d& arguments );

    //! Function to retrieve the interpolated daily Earth orientation data
    std::shared_ptr< EarthOrientationDataEvaluator > getEarthOrientationDataEvaluator( ) const
    {
        return eopData_;
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
    //! Interpolated daily Earth orientation data.
    std::shared_ptr< EarthOrientationDataEvaluator > eopData_;

    //! Function returning UT1 at a given TAI, used for the fundamental arguments of the short-period variations. If empty,
    //! UT1 is computed from the daily Earth orientation data only.
    std::function< Time( const Time& ) > ut1FromTai_;

    //! Object calculating short period polar motion variations.
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator_;
};

}  // namespace earth_orientation

}  // namespace tudat

#endif  // TUDAT_POLARMOTIONCALCULATOR_H
