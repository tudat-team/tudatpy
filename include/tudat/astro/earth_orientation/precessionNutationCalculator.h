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
     * Constructor for CIO-based precession-nutation calculation object.
     * \param precessionNutationTheory IAU precession-nutation theory that is to be used.
     * \param eopData Interpolated daily Earth orientation data, from which the measured celestial pole offsets are taken.
     * \param angleInterpolatorSettings Settings to interpolate the X, Y and s of the IAU model, instead of computing them at
     * each epoch (nullptr if they are computed at each epoch).
     */
    PrecessionNutationCalculator(
            const basic_astrodynamics::IAUConventions precessionNutationTheory,
            const std::shared_ptr< EarthOrientationDataEvaluator > eopData,
            const std::shared_ptr< interpolators::InterpolatorGenerationSettings< double > > angleInterpolatorSettings = nullptr );

    //! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
    /*!
     *  Function to calculate the position of CIP in GCRS (X and Y) and CIO locator (s), in radians: the IAU model at the
     *  given TT, plus the celestial pole offsets dX and dY from the Earth orientation data at the given TAI.
     *  \param terrestrialTime TT at which calculation is to be performed.
     *  \param tai TAI at the same epoch.
     *  \return CIP position in GCRS (X and Y) and CIO-locator (s).
     */
    Eigen::Vector3d getPositionOfCipInGcrsFromTai( const double terrestrialTime, const Time& tai );

    //! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
    /*!
     *  Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator, at a given TT.
     *  \param terrestrialTime TT at which calculation is to be performed.
     *  \return CIP position in GCRS (X and Y) and CIO-locator (s).
     */
    Eigen::Vector3d getPositionOfCipInGcrs( const double terrestrialTime );

    //! X, Y and s from the IAU model (without celestial pole offsets) at the given TT, and their derivatives per TT second.
    //! The derivatives are computed numerically from the SOFA model, also when the values themselves are interpolated.
    std::pair< Eigen::Vector3d, Eigen::Vector3d > getNominalPositionAndRate( const Time& tt );

    //! Function to retrieve the interpolated daily Earth orientation data
    std::shared_ptr< EarthOrientationDataEvaluator > getEarthOrientationDataEvaluator( ) const
    {
        return eopData_;
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
    //! X, Y and s (radians) from the SOFA implementation of the IAU model at the given TT (seconds since J2000).
    Eigen::Vector3d getNominalPositionFromSofa( const Time& tt );

    //! Interpolated daily Earth orientation data.
    std::shared_ptr< EarthOrientationDataEvaluator > eopData_;

    //! IAU precession-nutation theory that is to be used.
    basic_astrodynamics::IAUConventions precessionNutationTheory_;

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
