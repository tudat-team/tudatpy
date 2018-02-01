/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/External/SofaInterface/earthOrientation.h"

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
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator < double, Eigen::Vector2d > >
            dailyCorrectionInterpolator );

    //! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
    /*!
     *  Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
     *  Only TT is taken as input, UTC is calculated internally by using the SOFA wrappers. If UTC is
     *  known when calling this function, use the overloaded version, instead.
     *  \param terrestrialTime TT at which calculation is to be performed.
     *  \return Pair of CIP position in GCRS (X and Y) and CIO-locator (s).
     */
    std::pair< Eigen::Vector2d, double > getPositionOfCipInGcrs(
            const double terrestrialTime );

    //! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
    /*!
     *  Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
     *  TT and UTC is taken as input,If UTC is not known when calling this function, use the overloaded version, instead.
     *  \param terrestrialTime TT at which calculation is to be performed.
     *  \param utc UTC at which calculation is to be performed.
     *  \return Pair of CIP position in GCRS (X and Y) and CIO-locator (s).
     */
    std::pair< Eigen::Vector2d, double > getPositionOfCipInGcrs(
            const double terrestrialTime,
            const double utc );

    boost::shared_ptr< interpolators::OneDimensionalInterpolator
    < double, Eigen::Vector2d > > getDailyCorrectionInterpolator( )
    {
        return dailyCorrectionInterpolator_;
    }

private:

    //! Interpolator for daily measured values of precession-nutation corrections.
    /*!
     *  Interpolator taking UTC time since J2000 as input and returning
     *  interpolated daily measured values of precession-nutation corrections.
     */
    boost::shared_ptr< interpolators::OneDimensionalInterpolator
    < double, Eigen::Vector2d > > dailyCorrectionInterpolator_;

    //! Function pointer returning the nominal CIP position  in the GCRS  and CIO locator.
    /*!
     *  Function pointer returning the nominal CIP position in the GCRS (X and Y, see IERS Conventions 2010),
     *  first pair argument, and CIO locator (s, see IERS Conventions 2010), second pair argument.
     */
    boost::function< std::pair< Eigen::Vector2d, double > ( const double ) > nominalCipPositionFunction_;

};

}

}

#endif // TUDAT_PRECESSIONNUTATIONCALCULATOR_H
