#ifndef PRECESSIONNUTATIONCALCULATOR_H
#define PRECESSIONNUTATIONCALCULATOR_H

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"

#include "Tudat/External/SofaInterface/earthOrientation.h"

namespace tudat
{

namespace earth_orientation
{

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
            const IAUConventions precessionNutationTheory,
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, Eigen::Vector2d > > dailyCorrectionInterpolator );

    //! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
    /*!
     *  Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
     *  Only TT is taken as input, UTC is calculated internally by using the SOFA wrappers. If UTC is
     *  known when calling this function, use the overloaded version, instead.
     *  \param terrestrialTime TT at which calculation is to be performed.
     *  \param terrestrialTimeDaysEpochShift Julian days since julian day = 0 at which TT=0 is defined in input.
     *  \return Pair of CIP position in GCRS (X and Y) and CIO-locator (s).
     */
    std::pair< Eigen::Vector2d, double > getPositionOfCipInGcrs(
            const double terrestrialTime,
            const double terrestrialTimeDaysEpochShift = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    //! Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
    /*!
     *  Function to calculate the position of CIP in GCRS (CIO-based precession-nutation) and CIO-locator.
     *  TT and UTC is taken as input,If UTC is not known when calling this function, use the overloaded version, instead.
     *  \param terrestrialTime TT at which calculation is to be performed.
     *  \param utc UTC at which calculation is to be performed.
     *  \param terrestrialTimeDaysEpochShift Julian days since julian day = 0 at which TT=0 is defined in input.
     *  \return Pair of CIP position in GCRS (X and Y) and CIO-locator (s).
     */
    std::pair< Eigen::Vector2d, double > getPositionOfCipInGcrs(
            const double terrestrialTime,
            const double utc,
            const double terrestrialTimeDaysEpochShift = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

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

#endif // PRECESSIONNUTATIONCALCULATOR_H
