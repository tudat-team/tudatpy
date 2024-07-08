/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Regarding method in general:
 *          Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *          Wakker, K. F. astro I + II. Lecture Notes AE4-874, Delft University of
 *              Technology, Delft, Netherlands.
 *      Regarding the choice of initial guess:
 *          Musegaas, P., Optimization of Space Trajectories Including Multiple Gravity Assists and
 *              Deep Space Maneuvers, MSc thesis report, Delft University of Technology, 2012.
 *              [unpublished so far]. Section available on tudat website (tudat.tudelft.nl)
 *              under issue #539.
 *
 *    Notes
 *      There are known to be some issues on some systems with near-parabolic orbits that are very
 *      close to eccentricity=1.0. For these orbits, the unit tests establish that they are able
 *      to generate consistent results by converting to and from eccentric anomaly to a precision
 *      of 1.0e-9. If this quality of solution is not adequate for specific applications, the user
 *      should investigate the code further.
 *
 */

#ifndef TUDAT_CONVERT_MEAN_ANOMALY_TO_ECCENTRIC_ANOMALY_H
#define TUDAT_CONVERT_MEAN_ANOMALY_TO_ECCENTRIC_ANOMALY_H


#include <memory>

#include <boost/math/special_functions/asinh.hpp>

#include <cmath>

#include "tudat/math/root_finders/createRootFinder.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/math/basic/functionProxy.h"



namespace tudat
{

namespace orbital_element_conversions
{

//! Compute Kepler's function for elliptical orbits.
/*!
 * Computes Kepler's function, given as:
 * \f[
 *      f( E ) = E - e * sin( E ) - M
 * \f]
 * for elliptical orbits, where \f$ E \f$ is the eccentric anomaly, \f$ e \f$ is the
 * eccentricity, \f$ M \f$ is the mean anomaly. All eccentricities >= 0.0 and < 1.0 are valid.
 * \param eccentricAnomaly Eccentric anomaly.
 * \param eccentricity Eccentricity.
 * \param meanAnomaly Mean anomaly.
 * \return Value of Kepler's function for elliptical orbits.
 */
template< typename ScalarType = double >
ScalarType computeKeplersFunctionForEllipticalOrbits( const ScalarType eccentricAnomaly,
                                                      const ScalarType eccentricity,
                                                      const ScalarType meanAnomaly )
{
    return eccentricAnomaly - eccentricity * std::sin( eccentricAnomaly ) - meanAnomaly;
}

//! Compute first-derivative of Kepler's function for elliptical orbits.
/*!
 * Computes the first-derivative of Kepler's function, given as:
 * \f[
 *      \frac{ df( E ) } { dE } = 1 - e * cos( E )
 * \f]
 * for elliptical orbits, where \f$ E \f$ is the eccentric anomaly, and \f$ e \f$ is the
 * eccentricity. All eccentricities >= 0.0 and < 1.0 are valid.
 * \param eccentricAnomaly Eccentric anomaly.
 * \param eccentricity Eccentricity.
 * \return Value of first-derivative of Kepler's function for elliptical orbits.
 */
template< typename ScalarType = double >
ScalarType computeFirstDerivativeKeplersFunctionForEllipticalOrbits(
        const ScalarType eccentricAnomaly,
        const ScalarType eccentricity )
{
    return mathematical_constants::getFloatingInteger< ScalarType >( 1 )-
            eccentricity * std::cos( eccentricAnomaly );
}

//! Compute Kepler's function for hyperbolic orbits.
/*!
 * Computes Kepler's function, given as:
 * \f[
 *      f( F ) = e * sinh( F ) - F - M
 * \f]
 * for hyperbolic orbits, where \f$ F \f$ is the hyperbolic eccentric anomaly, \f$ e \f$ is the
 * eccentricity, \f$ M \f$ is the mean anomaly. All eccentricities > 1.0 are valid.
 * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.
 * \param eccentricity Eccentricity.
 * \param hyperbolicMeanAnomaly Hyperbolic mean anomaly..
 * \return Value of Kepler's function for hyperbolic orbits.
 */
template< typename ScalarType = double >
ScalarType computeKeplersFunctionForHyperbolicOrbits( const ScalarType hyperbolicEccentricAnomaly,
                                                      const ScalarType eccentricity,
                                                      const ScalarType hyperbolicMeanAnomaly )
{
    return eccentricity * std::sinh( hyperbolicEccentricAnomaly )
            - hyperbolicEccentricAnomaly - hyperbolicMeanAnomaly;
}

//! Compute first-derivative of Kepler's function for hyperbolic orbits.
/*!
 * Computes the first-derivative of Kepler's function, given as:
 * \f[
 *      \frac{ df( F ) } { dF } = e * cosh( F ) - 1
 * \f]
 * for hyperbolic orbits, where \f$ F \f$ is the hyperbolic eccentric anomaly, and \f$ e \f$ is
 * the eccentricity. All eccentricities > 1.0 are valid.
 * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly
 * \param eccentricity Eccentricity.
 * \return Value of first-derivative of Kepler's function for hyperbolic orbits.
 */
template< typename ScalarType = double >
ScalarType computeFirstDerivativeKeplersFunctionForHyperbolicOrbits(
        const ScalarType hyperbolicEccentricAnomaly,
        const ScalarType eccentricity)
{
    return eccentricity * std::cosh( hyperbolicEccentricAnomaly ) - 1.0;
}

//! Convert mean anomaly to eccentric anomaly.
/*!
 * Converts mean anomaly to eccentric anomaly for elliptical orbits for all eccentricities >=
 * 0.0 and < 1.0. If the conversion fails or the eccentricity falls outside the valid range,
 * then ScalarType::NaN is returned. Calculated with an accuracy of 1.0e-13 for all cases, except
 * for some near-parabolic cases in which macine precision problems occur. These are tested
 * against an accuracy of 1.0e-9. Near-parabolic in this sense means e > 1.0-1.0e-11. Also
 * note that your mean anomaly is automatically transformed to fit within the 0 to 2.0*PI
 * spectrum. Numerical tests performed using double ScalarType.
 * \param eccentricity Eccentricity of the orbit [-].
 * \param aMeanAnomaly Mean anomaly to convert to eccentric anomaly [rad].
 * \param useDefaultInitialGuess Boolean specifying whether to use default initial guess [-].
 * \param userSpecifiedInitialGuess Initial guess for rootfinder [rad].
 * \param rootFinder Shared-pointer to the rootfinder that is to be used. Default is
 *          Newton-Raphson using 1000 iterations as maximum and apprximately 1.0e-13 absolute
 *          X-tolerance (for doubles; 500 times ScalarType resolution ).
 *          Higher precision may invoke machine precision problems for some values.
 * \return Eccentric anomaly [rad].
 */
template< typename ScalarType = double >
ScalarType convertMeanAnomalyToEccentricAnomaly(
        const ScalarType eccentricity, const ScalarType aMeanAnomaly,
        const bool useDefaultInitialGuess = true,
        const ScalarType userSpecifiedInitialGuess = TUDAT_NAN,
        std::shared_ptr< root_finders::RootFinder< ScalarType > > rootFinder =
        std::shared_ptr< root_finders::RootFinder< ScalarType > >( ) )
{
    using namespace mathematical_constants;
    using namespace root_finders;


    // Set mean anomaly to region between 0 and 2 PI.
    ScalarType meanAnomaly = basic_mathematics::computeModulo< ScalarType >(
                aMeanAnomaly, getFloatingInteger< ScalarType >( 2 ) *
                getPi< ScalarType >( ) );

    // Required because the make_shared in the function definition gives problems for MSVC.
    if ( !rootFinder.get( ) )
    {
        ScalarType tolerance = 10.0 * std::numeric_limits< ScalarType >::epsilon( );

        // Loosen tolerance for near-parabolic orbits
        if( std::fabs( eccentricity - getFloatingInteger< ScalarType >( 1 ) ) <
                1.0E5 * std::numeric_limits< ScalarType >::epsilon( ) )
        {
            tolerance *= 2.5;
        }

        rootFinder = root_finders::createRootFinder< ScalarType >(
                    root_finders::newtonRaphsonRootFinderSettings(
                        TUDAT_NAN, tolerance, TUDAT_NAN, 20, root_finders::throw_exception ) );
    }

    // Declare eccentric anomaly.
    ScalarType eccentricAnomaly = TUDAT_NAN;

    // Check if orbit is elliptical, and not near-parabolic.
    if ( eccentricity < getFloatingInteger< ScalarType >( 1 ) &&
         eccentricity >= getFloatingInteger< ScalarType >( 0 ) )
    {
        // Create an object containing the function of which we whish to obtain the root from.
        std::shared_ptr< basic_mathematics::FunctionProxy< ScalarType, ScalarType > > rootFunction
                = std::make_shared< basic_mathematics::FunctionProxy< ScalarType, ScalarType >  >(
                    std::bind( &computeKeplersFunctionForEllipticalOrbits< ScalarType >, std::placeholders::_1,
                                 eccentricity, meanAnomaly ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding(
                    -1, std::bind(
                        &computeFirstDerivativeKeplersFunctionForEllipticalOrbits< ScalarType >, std::placeholders::_1,
                        eccentricity ) );

        // Declare initial guess.
        ScalarType initialGuess = TUDAT_NAN;

        // Set the initial guess. Check if the default scheme is to be used or a user specified
        // value should be used.
        // !!!!!!!!!!!!!     IMPORTANT     !!!!!!!!!!!!!
        // If this scheme is changed, please run a very extensive test suite. The root finder
        // function tends to be chaotic for some very specific combinations of mean anomaly and
        // eccentricity. Various random tests of 100.000.000 samples were done to verify the
        // functionality of this one, and of another option for the starter: PI. [Musegaas,2012]
        if ( useDefaultInitialGuess )
        {
            if ( meanAnomaly > getPi< ScalarType >( ) )
            {
                initialGuess =  meanAnomaly - eccentricity;
            }
            else
            {
                initialGuess = meanAnomaly + eccentricity;
            }
        }
        else
        {
            initialGuess = userSpecifiedInitialGuess;
        }

        try
        {
            // Set eccentric anomaly based on result of Newton-Raphson root-finding algorithm.
            eccentricAnomaly = rootFinder->execute( rootFunction, initialGuess );
        }
        // Use bisection algorithm if root finder fails
        catch( std::runtime_error const& )
        {
            // Set tolerance
            ScalarType tolerance = 10.0 * std::numeric_limits< ScalarType >::epsilon( );

            // Set upper/lower bounds (loosely)
            ScalarType lowerBound, upperBound;
            if ( meanAnomaly > getPi< ScalarType >( ) )
            {
                lowerBound = 0.5 * getPi< ScalarType >( ) + 10.0 * std::numeric_limits< ScalarType >::epsilon( );
                upperBound = mathematical_constants::getFloatingInteger< ScalarType >( 2 ) * getPi< ScalarType >( );
            }
            else
            {
                lowerBound = mathematical_constants::getFloatingInteger< ScalarType >( 0 );
                upperBound = 1.5 * getPi< ScalarType >( );
            }

            // Create root finder
            std::shared_ptr< RootFinder< ScalarType > > bisectionRootfinder =
                    rootFinder = root_finders::createRootFinder< ScalarType >(
                                root_finders::bisectionRootFinderSettings(
                                    TUDAT_NAN, tolerance, TUDAT_NAN, 100, root_finders::accept_result ),
                         lowerBound, upperBound );

            initialGuess = meanAnomaly;
            eccentricAnomaly = bisectionRootfinder->execute( rootFunction, initialGuess );
        }
    }
    //  Eccentricity is invalid: eccentricity < 0.0 or eccentricity >= 1.0.
    else
    {
        throw std::runtime_error( "Invalid eccentricity. Valid range is 0.0 <= e < 1.0. Eccentricity was: " +
                            std::to_string( eccentricity ) );
    }

    // Return eccentric anomaly.
    return eccentricAnomaly;
}


//! Convert mean anomaly to hyperbolic eccentric anomaly.
/*!
 * Converts mean anomaly to hyperbolic eccentric anomaly for hyperbolic orbits for all
 * eccentricities > 1.0. If the conversion fails or the eccentricity falls outside the valid
 * range, then double::NaN is returned. Calculated with an accuracy of 1.0e-14 for all
 * reasonable cases (eccentricities up to 1.0e15, mean anomalies -1.2e12 to 1.2e12, test cases
 * using doubles).
 * \param eccentricity Eccentricity of the orbit [-].
 * \param hyperbolicMeanAnomaly Hyperbolic mean anomaly to convert to eccentric anomaly [rad].
 * \param useDefaultInitialGuess Boolean specifying whether to use default initial guess [-].
 * \param userSpecifiedInitialGuess Initial guess for rootfinder [rad].
 * \param aRootFinder Shared-pointer to the rootfinder that is to be used. Default is
 *          Newton-Raphson using 1000 iterations as maximum and apprximately 5.0e-15 absolute
 *          X-tolerance (for doubles; 25 times ScalarType resolution ).
 *          Higher precision may invoke machine precision problems for some values.
 * \return Hyperbolic eccentric anomaly.
 */
template< typename ScalarType = double >
ScalarType convertMeanAnomalyToHyperbolicEccentricAnomaly(
        const ScalarType eccentricity, const ScalarType hyperbolicMeanAnomaly,
        const bool useDefaultInitialGuess = true,
        const ScalarType userSpecifiedInitialGuess = TUDAT_NAN,
        std::shared_ptr< root_finders::RootFinder< ScalarType > > aRootFinder =
        std::shared_ptr< root_finders::RootFinder< ScalarType > >( ) )
{
    using namespace mathematical_constants;
    using namespace root_finders;


    std::shared_ptr< RootFinder< ScalarType > > rootFinder = aRootFinder;

    // Required because the make_shared in the function definition gives problems for MSVC.
    if ( !rootFinder.get( ) )
    {
        double toleranceMultiplier = ( eccentricity - 1.0 ) < 1.0E-3 ? 10.0 : 1.0;
        rootFinder = createRootFinder< ScalarType >(
                    newtonRaphsonRootFinderSettings(
                        TUDAT_NAN, toleranceMultiplier * 25.0 * std::numeric_limits< ScalarType >::epsilon( ),
                        TUDAT_NAN, 1000, root_finders::throw_exception ) );
    }
    // Declare hyperbolic eccentric anomaly.
    ScalarType hyperbolicEccentricAnomaly = TUDAT_NAN;

    // Check if orbit is hyperbolic.
    if ( eccentricity > getFloatingInteger< ScalarType >( 1 ) )
    {
        // Create an object containing the function of which we whish to obtain the root from.
        std::shared_ptr< basic_mathematics::FunctionProxy< ScalarType, ScalarType > > rootFunction
                = std::make_shared< basic_mathematics::FunctionProxy< ScalarType, ScalarType > >(
                    std::bind( &computeKeplersFunctionForHyperbolicOrbits< ScalarType >, std::placeholders::_1,
                                 eccentricity, hyperbolicMeanAnomaly ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding(
                    -1, std::bind(
                        &computeFirstDerivativeKeplersFunctionForHyperbolicOrbits< ScalarType >,
                        std::placeholders::_1, eccentricity ) );

        // Declare initial guess.
        ScalarType initialGuess = TUDAT_NAN;

        // Set the initial guess. Check if the default scheme is to be used or a user specified
        // value should be used. See [Wakker, 2007] for derivations of the default values.
        // Note that an error was detected in these starter values, as is discussed in
        // [Musegaas,2012].
        // !!!!!!!!!!!!!     IMPORTANT     !!!!!!!!!!!!!
        // If this scheme is changed, please run a very extensive test suite. The root finder
        // function tends to be chaotic for some very specific combinations of mean anomaly and
        // eccentricity. Various random tests of 100.000.000 samples were done to verify the
        // functionality of this one. [Musegaas,2012]
        if ( useDefaultInitialGuess )
        {
            if ( std::abs( hyperbolicMeanAnomaly ) <
                 getFloatingInteger< ScalarType >( 6 ) * eccentricity )
            {
                initialGuess =
                        std::sqrt( getFloatingInteger< ScalarType >( 8 ) *
                                   ( eccentricity - getFloatingInteger< ScalarType >( 1 ) ) /
                                   eccentricity ) *
                        std::sinh( getFloatingFraction< ScalarType >( 1, 3 ) * boost::math::asinh(
                                       getFloatingInteger< ScalarType >( 3 ) *
                                       hyperbolicMeanAnomaly /
                                       ( std::sqrt( getFloatingInteger< ScalarType >( 8 ) *
                                                    ( eccentricity -
                                                      getFloatingInteger< ScalarType >( 1 ) ) /
                                                    eccentricity ) *
                                         ( eccentricity - getFloatingInteger< ScalarType >( 1 )
                                           ) ) ) );
            }
            else if ( hyperbolicMeanAnomaly > getFloatingInteger< ScalarType >( 6 ) * eccentricity )
            {
                initialGuess = ( std::log( getFloatingInteger< ScalarType >( 2 ) *
                                           hyperbolicMeanAnomaly / eccentricity ) );
            }
            else
            {
                initialGuess = ( - std::log( -getFloatingInteger< ScalarType >( 2 ) *
                                             hyperbolicMeanAnomaly / eccentricity ) );
            }
        }
        else
        {
            initialGuess = userSpecifiedInitialGuess;
        }

        // Set hyperbolic eccentric anomaly based on result of Newton-Raphson root-finding
        // algorithm.
        try
        {
            hyperbolicEccentricAnomaly = rootFinder->execute( rootFunction, initialGuess );

        }
        catch( std::runtime_error const& )
        {
            rootFinder = createRootFinder< ScalarType >(
                        bisectionRootFinderSettings(
                            TUDAT_NAN, 20.0 * std::numeric_limits< ScalarType >::epsilon( ),
                            TUDAT_NAN, 100 ), getFloatingInteger< ScalarType >( 0 ),
                        getFloatingInteger< ScalarType >( 2 ) * getPi< ScalarType >( ), root_finders::accept_result );

            hyperbolicEccentricAnomaly = rootFinder->execute( rootFunction, initialGuess );
        }

    }

    // In this case the orbit is not hyperbolic.
    else
    {
        throw std::runtime_error( "Invalid eccentricity. Valid range is e > 1.0. Eccentricity was: " +
                            std::to_string( eccentricity ) );
    }

    // Return hyperbolic eccentric anomaly.
    return hyperbolicEccentricAnomaly;
}

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_CONVERT_MEAN_ANOMALY_TO_ECCENTRIC_ANOMALY_H
