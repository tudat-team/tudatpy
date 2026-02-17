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

#ifndef TUDAT_TIME_EPHEMERIS_H
#define TUDAT_TIME_EPHEMERIS_H

#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/relativity/relativisticPotentials.h"
#include "tudat/astro/relativity/relativisticTimeConversion.h"
#include "tudat/math/interpolators/createInterpolator.h"

namespace tudat
{

//! Combine chained time-difference functions with configurable evaluation epochs.
/*!
 *  This utility evaluates a sequence of conversion functions and accumulates the resulting
 *  time differences. Each conversion function can be evaluated at either the original epoch
 *  or at one of the intermediate corrected epochs.
 *  \param timeDifferenceFunctions Ordered list of individual conversion functions.
 *  \param evaluationStepIndices Index list defining at which intermediate epoch each
 *  conversion function is evaluated.
 *  \param time Base input epoch.
 *  \return Total time difference obtained by summing all conversion steps.
 */
double combineTimeDifferenceFunction( const std::vector< std::function< double( const double ) > >& timeDifferenceFunctions,
                                      const std::vector< int >& evaluationStepIndices,
                                      const double time );

class TimeEphemeris
{
public:

    explicit TimeEphemeris( const std::string& centralBodyName ) : centralBodyName_( centralBodyName )
    { }

    virtual ~TimeEphemeris( ) { }

    //! Function to retrieve the time difference at a given time between two scales.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param inputTime Epoch at which the conversion is evaluated.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Time difference \f$t_{\mathrm{output}}-t_{\mathrm{input}}\f$.
     */
    template< typename TimeType >
    TimeType getTimeDifference( const basic_astrodynamics::TimeScales inputScale,
                                const basic_astrodynamics::TimeScales outputScale,
                                const TimeType inputTime,
                                const std::string& pointIdentifier = "" )
    {
        return static_cast< TimeType >( getTimeDifferenceFunction( inputScale, outputScale, pointIdentifier )(
                    static_cast< double >( inputTime ) ) );
    }

    //! Double-precision convenience overload.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param inputTime Epoch at which the conversion is evaluated.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Time difference \f$t_{\mathrm{output}}-t_{\mathrm{input}}\f$ in seconds.
     */
    double getTimeDifference( const basic_astrodynamics::TimeScales inputScale,
                              const basic_astrodynamics::TimeScales outputScale,
                              const double inputTime,
                              const std::string& pointIdentifier = "" )
    {
        return getTimeDifference< double >( inputScale, outputScale, inputTime, pointIdentifier );
    }

    //! Function to retrieve multiple time differences between two scales.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param inputTimes Epochs at which the conversion is evaluated.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Map of input epoch to time difference \f$t_{\mathrm{output}}-t_{\mathrm{input}}\f$.
     */
    template< typename TimeType >
    std::map< TimeType, TimeType > getTimeDifferences( const basic_astrodynamics::TimeScales inputScale,
                                                        const basic_astrodynamics::TimeScales outputScale,
                                                        const std::vector< TimeType >& inputTimes,
                                                        const std::string& pointIdentifier = "" )
    {
        const std::function< double( const double ) > timeDifferenceFunction =
                getTimeDifferenceFunction( inputScale, outputScale, pointIdentifier );

        std::map< TimeType, TimeType > timeDifferences;
        for( const TimeType& inputTime : inputTimes )
        {
            timeDifferences[ inputTime ] = static_cast< TimeType >( timeDifferenceFunction( static_cast< double >( inputTime ) ) );
        }
        return timeDifferences;
    }

    //! Double-precision convenience overload.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param inputTimes Epochs at which the conversion is evaluated.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Map of input epoch to time difference in seconds.
     */
    std::map< double, double > getTimeDifferences( const basic_astrodynamics::TimeScales inputScale,
                                                   const basic_astrodynamics::TimeScales outputScale,
                                                   const std::vector< double >& inputTimes,
                                                   const std::string& pointIdentifier = "" )
    {
        return getTimeDifferences< double >( inputScale, outputScale, inputTimes, pointIdentifier );
    }

    //! Retrieve the callable that computes the requested conversion at a given epoch.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Function object evaluating \f$t_{\mathrm{output}}-t_{\mathrm{input}}\f$.
     */
    virtual std::function< double( const double ) > getTimeDifferenceFunction(
            const basic_astrodynamics::TimeScales inputScale,
            const basic_astrodynamics::TimeScales outputScale,
            const std::string& pointIdentifier = "" ) = 0;

protected:

    std::string centralBodyName_;
};

//! Relativistic time ephemeris using post-Newtonian conversion terms and optional direct corrections.
/*!
 *  This class combines pre-integrated coordinate-time differences with direct correction terms to build
 *  conversions between barycentric, body-centric, and topocentric proper time scales.
 *
 *  The barycentric/body-centric component follows the Soffel et al. (2003), Eq. (58)-type
 *  post-Newtonian model, including first-order and optional second-order terms.
 *
 *  The body-centric/topocentric component follows the Turyshev et al. (2013), Eq. (22)-type
 *  local proper-time model, including kinetic, central-potential, tidal, and acceleration terms.
 *
 *  In Tudat, the propagated variables are time differences (e.g. \f$\tau-t_C\f$), which are then
 *  interpolated and combined to provide conversions between requested scales.
 */
class TimeEphemerisFromPostNewtonianExpansion : public TimeEphemeris
{
public:
    typedef std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > TimeDifferenceInterpolator;

    //! Constructor.
    /*!
     *  \param barycenterToPlanetCenterCoordinateTimeInterpolator Interpolator for
     *  barycentric-to-bodycentric coordinate-time difference.
     *  \param planetCenterToBarycenterCoordinateTimeInterpolator Interpolator for
     *  bodycentric-to-barycentric coordinate-time difference.
     *  \param centralBodyName Name of the central body associated with this time ephemeris.
     *  \param groundStationPositionFunctions Map from reference point identifier to position function.
     *  \param planetCoordinateToProperTimeInterpolators Map from reference point identifier to
     *  bodycentric-to-proper-time interpolator.
     *  \param properTimeToPlanetCoordinateInterpolators Map from reference point identifier to
     *  proper-time-to-bodycentric interpolator.
     */
    TimeEphemerisFromPostNewtonianExpansion(
            const TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
            const TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator,
            const std::string& centralBodyName,
            const std::map< std::string, std::function< Eigen::Vector3d( const double ) > >& groundStationPositionFunctions =
            ( std::map< std::string, std::function< Eigen::Vector3d( const double ) > >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ) ) :
        TimeEphemeris( centralBodyName ),
        barycenterToPlanetCenterCoordinateTimeInterpolator_( barycenterToPlanetCenterCoordinateTimeInterpolator ),
        planetCenterToBarycenterCoordinateTimeInterpolator_( planetCenterToBarycenterCoordinateTimeInterpolator ),
        groundStationPositionFunctions_( groundStationPositionFunctions ),
        planetCoordinateToProperTimeInterpolators_( planetCoordinateToProperTimeInterpolators ),
        properTimeToPlanetCoordinateInterpolators_( properTimeToPlanetCoordinateInterpolators )
    { }

    virtual ~TimeEphemerisFromPostNewtonianExpansion( ) { }

    //! Reset barycentric/bodycentric conversion interpolators.
    /*!
     *  \param barycenterToPlanetCenterCoordinateTimeInterpolator Interpolator for
     *  barycentric-to-bodycentric conversion.
     *  \param planetCenterToBarycenterCoordinateTimeInterpolator Interpolator for
     *  bodycentric-to-barycentric conversion.
     */
    void resetBarycentricToBodycentricInterpolators(
            const TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
            const TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator )
    {
        barycenterToPlanetCenterCoordinateTimeInterpolator_ = barycenterToPlanetCenterCoordinateTimeInterpolator;
        planetCenterToBarycenterCoordinateTimeInterpolator_ = planetCenterToBarycenterCoordinateTimeInterpolator;
    }

    //! Reset bodycentric/topocentric conversion interpolators for a reference point.
    /*!
     *  \param planetCoordinateToProperTimeInterpolator Interpolator for bodycentric-to-proper conversion.
     *  \param properTimeToPlanetCoordinateInterpolator Interpolator for proper-to-bodycentric conversion.
     *  \param referencePoint Reference point identifier.
     *  \param referencePointPositionFunction Optional position function for the reference point.
     */
    void resetBodycentricToTopocentricInterpolators(
            const TimeDifferenceInterpolator planetCoordinateToProperTimeInterpolator,
            const TimeDifferenceInterpolator properTimeToPlanetCoordinateInterpolator,
            const std::string& referencePoint,
            const std::function< Eigen::Vector3d( const double ) > referencePointPositionFunction = nullptr )
    {
        planetCoordinateToProperTimeInterpolators_[ referencePoint ] = planetCoordinateToProperTimeInterpolator;
        properTimeToPlanetCoordinateInterpolators_[ referencePoint ] = properTimeToPlanetCoordinateInterpolator;
        if( !doesReferencePointTopocentricConverterExist( referencePoint ) )
        {
            if( referencePointPositionFunction == nullptr )
            {
                std::cerr << "Error when resetting bodycentric to topocentric time converter for point " << referencePoint
                          << ", must also provide point position function; point not yet known in TimeEphemeris" << std::endl;
            }
            else
            {
                groundStationPositionFunctions_[ referencePoint ] = referencePointPositionFunction;
            }
        }
    }

    //! Retrieve conversion function for a requested pair of scales.
    std::function< double( const double ) > getTimeDifferenceFunction(
            const basic_astrodynamics::TimeScales inputScale,
            const basic_astrodynamics::TimeScales outputScale,
            const std::string& pointIdentifier = "" ) override;

    //! Compute direct correction term from a position function.
    /*!
     *  \param positionVectorFunctionFromReferencePoint Function returning point position at an epoch.
     *  \param currentTime Evaluation epoch.
     *  \param multiplier Optional multiplicative factor (used for inverse conversion chains).
     *  \return Direct correction contribution.
     */
    double calculateDirectTimeDifferenceTermFromFunction(
            const std::function< Eigen::Vector3d( const double ) > positionVectorFunctionFromReferencePoint,
            const double currentTime,
            const double multiplier = 1.0 )
    {
        return multiplier * calculateDirectTimeDifferenceTerm( positionVectorFunctionFromReferencePoint( currentTime ), currentTime );
    }

    //! Compute direct correction term for a point position.
    /*!
     *  \param positionVectorFromReferencePoint Point position relative to the body center.
     *  \param currentTime Evaluation epoch.
     *  \return Direct correction contribution.
     */
    virtual double calculateDirectTimeDifferenceTerm(
            const Eigen::Vector3d positionVectorFromReferencePoint,
            const double currentTime ) = 0;

    //! Check if all topocentric converter components are available for a reference point.
    /*!
     *  \param referencePointName Name/identifier of the reference point.
     *  \return True if position and both conversion interpolators are available.
     */
    bool doesReferencePointTopocentricConverterExist( const std::string& referencePointName )
    {
        return ( groundStationPositionFunctions_.count( referencePointName ) > 0 ) &&
               ( planetCoordinateToProperTimeInterpolators_.count( referencePointName ) > 0 ) &&
               ( properTimeToPlanetCoordinateInterpolators_.count( referencePointName ) > 0 );
    }

protected:

    TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator_;
    TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator_;

    std::map< std::string, std::function< Eigen::Vector3d( const double ) > > groundStationPositionFunctions_;
    std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators_;
    std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators_;
};

//! Post-Newtonian converter with first-order direct local term.
/*!
 *  Adds the first-order direct correction term for body-centric <-> topocentric conversion on top of the
 *  interpolated integral terms, using the central-body barycentric velocity (Soffel/IAU first-order formulation).
 */
class TimeEphemerisWithFirstOrderDirectConversion : public TimeEphemerisFromPostNewtonianExpansion
{
public:
    typedef std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > TimeDifferenceInterpolator;

    //! Constructor.
    /*!
     *  \param barycenterToPlanetCenterCoordinateTimeInterpolator Interpolator for
     *  barycentric-to-bodycentric conversion.
     *  \param planetCenterToBarycenterCoordinateTimeInterpolator Interpolator for
     *  bodycentric-to-barycentric conversion.
     *  \param centralBodyName Name of the central body.
     *  \param centralBodyStateFunction Function returning barycentric central-body state.
     *  \param groundStationPositionFunctions Map from reference point identifier to position function.
     *  \param planetCoordinateToProperTimeInterpolators Map of bodycentric-to-proper interpolators.
     *  \param properTimeToPlanetCoordinateInterpolators Map of proper-to-bodycentric interpolators.
     */
    TimeEphemerisWithFirstOrderDirectConversion(
            const TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
            const TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator,
            const std::string& centralBodyName,
            const std::function< Eigen::Vector6d( const double ) > centralBodyStateFunction,
            const std::map< std::string, std::function< Eigen::Vector3d( const double ) > >& groundStationPositionFunctions =
            ( std::map< std::string, std::function< Eigen::Vector3d( const double ) > >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ) ) :
        TimeEphemerisFromPostNewtonianExpansion(
            barycenterToPlanetCenterCoordinateTimeInterpolator,
            planetCenterToBarycenterCoordinateTimeInterpolator,
            centralBodyName,
            groundStationPositionFunctions,
            planetCoordinateToProperTimeInterpolators,
            properTimeToPlanetCoordinateInterpolators ),
        centralBodyStateFunction_( centralBodyStateFunction )
    { }

    //! Evaluate first-order direct term.
    double calculateDirectTimeDifferenceTerm( const Eigen::Vector3d positionVectorFromReferencePoint,
                                              const double currentTime ) override
    {
        return relativity::calculateFirstOrderTcbToTcgDirectCorrection(
                    positionVectorFromReferencePoint, centralBodyStateFunction_( currentTime ).segment( 3, 3 ) );
    }

private:

    std::function< Eigen::Vector6d( const double ) > centralBodyStateFunction_;
};

//! Post-Newtonian converter with placeholder for second-order direct local term.
/*!
 *  This class is intended for an extension with second-order direct corrections in the body-centric <-> topocentric
 *  conversion chain. The direct term is currently not implemented and returns zero.
 */
class TimeEphemerisWithSecondOrderDirectConversion : public TimeEphemerisFromPostNewtonianExpansion
{
public:
    //! Constructor.
    /*!
     *  \param barycenterToPlanetCenterCoordinateTimeInterpolator Interpolator for
     *  barycentric-to-bodycentric conversion.
     *  \param planetCenterToBarycenterCoordinateTimeInterpolator Interpolator for
     *  bodycentric-to-barycentric conversion.
     *  \param centralBodyName Name of the central body.
     *  \param centralBodyStateFunction Function returning barycentric central-body state.
     *  \param groundStationPositionFunctions Map from reference point identifier to position function.
     *  \param planetCoordinateToProperTimeInterpolators Map of bodycentric-to-proper interpolators.
     *  \param properTimeToPlanetCoordinateInterpolators Map of proper-to-bodycentric interpolators.
     */
    TimeEphemerisWithSecondOrderDirectConversion(
            const TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
            const TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator,
            const std::string& centralBodyName,
            const std::function< Eigen::Vector6d( const double ) > centralBodyStateFunction,
            const std::map< std::string, std::function< Eigen::Vector3d( const double ) > >& groundStationPositionFunctions =
            ( std::map< std::string, std::function< Eigen::Vector3d( const double ) > >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ) ) :
        TimeEphemerisFromPostNewtonianExpansion(
            barycenterToPlanetCenterCoordinateTimeInterpolator,
            planetCenterToBarycenterCoordinateTimeInterpolator,
            centralBodyName,
            groundStationPositionFunctions,
            planetCoordinateToProperTimeInterpolators,
            properTimeToPlanetCoordinateInterpolators ),
        centralBodyStateFunction_( centralBodyStateFunction )
    { }

    //! Evaluate direct term (currently placeholder for future second-order implementation).
    double calculateDirectTimeDifferenceTerm( const Eigen::Vector3d positionVectorFromReferencePoint,
                                              const double currentTime ) override
    {
        std::cerr << "Error, second order time ephemeris direct difference term not yet implemented" << std::endl;
        return 0.0;
    }

private:
    std::function< Eigen::Vector6d( const double ) > centralBodyStateFunction_;
};

//! Time ephemeris using direct metric-based coordinate-time/proper-time interpolation.
/*!
 *  Implements direct barycentric-coordinate <-> proper-time conversion by using interpolated values generated
 *  from metric-based propagation of \f$d(\tau-t)/dt\f$, without an intermediate body-centric scale.
 *
 *  The propagated rate is obtained from the metric perturbation and coordinate velocity
 *  (see `evaluateProperTimeEquation`), with configurable series-expansion order for the square-root term.
 */
class TimeEphemerisDirectFromMetric : public TimeEphemeris
{
public:
    typedef std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > TimeDifferenceInterpolator;

    //! Constructor.
    /*!
     *  \param centralBodyName Name of central body associated with this time ephemeris.
     *  \param globalCoordinateToProperTimeInterpolators Map from point identifier to
     *  global-coordinate-to-proper-time interpolator.
     *  \param properTimeToGlobalCoordinateInterpolators Map from point identifier to
     *  proper-time-to-global-coordinate interpolator.
     */
    TimeEphemerisDirectFromMetric(
            const std::string& centralBodyName,
            const std::map< std::string, TimeDifferenceInterpolator > globalCoordinateToProperTimeInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > properTimeToGlobalCoordinateInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ) ) :
        TimeEphemeris( centralBodyName ),
        globalCoordinateToProperTimeInterpolators_( globalCoordinateToProperTimeInterpolators ),
        properTimeToGlobalCoordinateInterpolators_( properTimeToGlobalCoordinateInterpolators )
    { }

    //! Reset direct conversion interpolators for a specific reference point.
    /*!
     *  \param globalCoordinateToProperTimeInterpolator Interpolator for global-to-proper conversion.
     *  \param properTimeToGlobalCoordinateInterpolator Interpolator for proper-to-global conversion.
     *  \param referencePoint Reference point identifier.
     */
    void resetGlobalToProperTimeInterpolators(
            const TimeDifferenceInterpolator globalCoordinateToProperTimeInterpolator,
            const TimeDifferenceInterpolator properTimeToGlobalCoordinateInterpolator,
            const std::string& referencePoint )
    {
        globalCoordinateToProperTimeInterpolators_[ referencePoint ] = globalCoordinateToProperTimeInterpolator;
        properTimeToGlobalCoordinateInterpolators_[ referencePoint ] = properTimeToGlobalCoordinateInterpolator;
    }

    //! Retrieve conversion function for a requested pair of scales.
    std::function< double( const double ) > getTimeDifferenceFunction(
            const basic_astrodynamics::TimeScales inputScale,
            const basic_astrodynamics::TimeScales outputScale,
            const std::string& pointIdentifier = "" ) override;

    //! Retrieve global-coordinate-to-proper-time difference from interpolator.
    /*!
     *  \param pointId Reference point identifier.
     *  \param evaluationTime Epoch at which interpolator is evaluated.
     *  \return Time difference \f$\tau - t\f$.
     */
    double getGlobalCoordinateToProperTimeDifference( const std::string& pointId, const double evaluationTime )
    {
        return globalCoordinateToProperTimeInterpolators_.at( pointId )->interpolate( evaluationTime );
    }

    //! Retrieve proper-time-to-global-coordinate difference from interpolator.
    /*!
     *  \param pointId Reference point identifier.
     *  \param evaluationTime Epoch at which interpolator is evaluated.
     *  \return Time difference \f$t - \tau\f$.
     */
    double getProperToGlobalCoordinateTimeDifference( const std::string& pointId, const double evaluationTime )
    {
        return properTimeToGlobalCoordinateInterpolators_.at( pointId )->interpolate( evaluationTime );
    }

protected:

    std::map< std::string, TimeDifferenceInterpolator > globalCoordinateToProperTimeInterpolators_;
    std::map< std::string, TimeDifferenceInterpolator > properTimeToGlobalCoordinateInterpolators_;
};

}  // namespace tudat

#endif  // TUDAT_TIME_EPHEMERIS_H
