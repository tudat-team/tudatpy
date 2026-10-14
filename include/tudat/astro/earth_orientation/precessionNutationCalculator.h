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

#ifndef TUDAT_PRECESSIONNUTATIONCALCULATOR_H
#define TUDAT_PRECESSIONNUTATIONCALCULATOR_H

#include <memory>
#include <functional>
#include "tudat/astro/earth_orientation/earthOrientationDataEvaluator.h"

#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/interface/sofa/earthOrientation.h"

namespace tudat
{

namespace earth_orientation
{

//! Function to compute precession/nutation angles,
/*!
 *  Function to compute precession/nutation angles, defined by position of CIP in GCRS (X and Y parameters of IERS 2010
 *  Conventions, Section 5). Computed by combining IAU theory output with daily IERS corrections.
 */
class PrecessionNutationCalculator
{
public:
    //! Constructor for (CIO-based) precession-nutation calculation object.
    /*!
     * Constructor for CIO-based precession-nutation calculation object. This constructor binds the SOFA function wrapper
     * for the indicated precession-nutation model internally. If one wishes to use a different function (such as)
     * an interpolator, the alternative constructor should be used.
     * \param precessionNutationTheory IAU precession-nutation theory that is to be used.
     * \param dailyCorrectionInterpolator Interpolator taking UTC time since J2000 as input and returning
     * interpolated daily measured values of precession-nutation corrections.
     */
    PrecessionNutationCalculator(
            const basic_astrodynamics::IAUConventions precessionNutationTheory,
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > > dailyCorrectionInterpolator,
            const std::shared_ptr< interpolators::InterpolatorGenerationSettings< double > > angleInterpolatorSettings = nullptr,
            const std::shared_ptr< EarthOrientationDataEvaluator > evaluator = nullptr );

    //! Return nominal TT X/Y/s in radians with observed dX/dY evaluated at precise TAI.
    Eigen::Vector3d getPositionOfCipInGcrsFromTai( const double terrestrialTime, const Time& tai );

    //! Nominal X/Y/s and d(X,Y,s)/dTT, before observed CPO. The rate is that of the direct SOFA model.
    std::pair< Eigen::Vector3d, Eigen::Vector3d > getNominalPositionAndRate( const Time& tt );

    std::shared_ptr< EarthOrientationDataEvaluator > getEarthOrientationDataEvaluator( ) const
    {
        return evaluator_;
    }

    //! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
    /*!
     *  Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
     *  Models constructed from EOP data use continuous TAI. Custom-interpolator models calculate UTC internally;
     *  use the other overload when their UTC time is already known.
     *  \param terrestrialTime TT at which calculation is to be performed.
     *  \return Pair of CIP position in GCRS (X and Y) and CIO-locator (s).
     */
    Eigen::Vector3d getPositionOfCipInGcrs( const double terrestrialTime );

    //! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
    /*!
     *  Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
     *  Custom-interpolator models use the supplied TT and UTC. Models constructed from EOP data
     *  derive continuous TAI from TT and do not use UTC.
     *  \param terrestrialTime TT at which calculation is to be performed.
     *  \param utc UTC at which calculation is to be performed.
     *  \return Pair of CIP position in GCRS (X and Y) and CIO-locator (s).
     */
    Eigen::Vector3d getPositionOfCipInGcrs( const double terrestrialTime, const double utc );

    //! Function to retrieve the interpolator for daily measured values of precession-nutation corrections.
    /*!
     * Function to retrieve the interpolator for daily measured values of precession-nutation corrections.
     * \return Interpolator for daily measured values of precession-nutation corrections.
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > > getDailyCorrectionInterpolator( )
    {
        if( evaluator_ )
        {
            throw std::runtime_error(
                    "Models constructed from EOP data do not have a daily UTC interpolator; use "
                    "getEarthOrientationDataEvaluator()->getCelestialPoleOffsets(tai)." );
        }
        return dailyCorrectionInterpolator_;
    }

    //! Function to retrieve the IAU precession-nutation theory that is to be used.
    /*!
     * Function to retrieve the IAU precession-nutation theory that is to be used.
     * \return IAU precession-nutation theory that is to be used.
     */
    basic_astrodynamics::IAUConventions getPrecessionNutationTheory( )
    {
        return precessionNutationTheory_;
    }

private:
    std::shared_ptr< EarthOrientationDataEvaluator > evaluator_;
    //! IAU precession-nutation theory that is to be used.
    basic_astrodynamics::IAUConventions precessionNutationTheory_;

    //! Interpolator for daily measured values of precession-nutation corrections.
    /*!
     *  Interpolator taking UTC time since J2000 as input and returning
     *  interpolated daily measured values of precession-nutation corrections.
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector2d > > dailyCorrectionInterpolator_;

    //! Function pointer returning the nominal CIP position  in the GCRS  and CIO locator.
    /*!
     *  Function pointer returning the nominal CIP position in the GCRS (X and Y, see IERS Conventions 2010),
     *  first pair argument, and CIO locator (s, see IERS Conventions 2010), second pair argument.
     */
    std::function< Eigen::Vector3d( const double ) > nominalCipPositionFunction_;
};

}  // namespace earth_orientation

}  // namespace tudat

#endif  // TUDAT_PRECESSIONNUTATIONCALCULATOR_H
