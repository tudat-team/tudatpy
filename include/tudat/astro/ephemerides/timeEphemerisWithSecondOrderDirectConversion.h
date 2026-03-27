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

#ifndef TUDAT_TIME_EPHEMERIS_WITH_SECOND_ORDER_DIRECT_CONVERSION_H
#define TUDAT_TIME_EPHEMERIS_WITH_SECOND_ORDER_DIRECT_CONVERSION_H

#include <functional>
#include <map>
#include <stdexcept>
#include <string>

#include <Eigen/Core>

#include "tudat/astro/ephemerides/timeEphemerisFromPostNewtonianExpansion.h"

namespace tudat
{

//! Post-Newtonian converter with placeholder for second-order direct local term.
template< typename TimeType = double, typename StateScalarType = double >
class TimeEphemerisWithSecondOrderDirectConversion : public TimeEphemerisFromPostNewtonianExpansion< TimeType, StateScalarType >
{
public:
    using TimeDifferenceInterpolator =
            typename TimeEphemerisFromPostNewtonianExpansion< TimeType, StateScalarType >::TimeDifferenceInterpolator;

    TimeEphemerisWithSecondOrderDirectConversion(
            const TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
            const TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator,
            const std::string& centralBodyName,
            const std::function< Eigen::Vector6d( const TimeType ) > centralBodyStateFunction,
            const std::map< std::string, std::function< Eigen::Vector3d( const TimeType ) > >& groundStationPositionFunctions =
            ( std::map< std::string, std::function< Eigen::Vector3d( const TimeType ) > >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ) );

    double calculateDirectTimeDifferenceTerm( const Eigen::Vector3d positionVectorFromReferencePoint,
                                              const TimeType currentTime ) override;

private:
    std::function< Eigen::Vector6d( const TimeType ) > centralBodyStateFunction_;
};

template< typename TimeType, typename StateScalarType >
TimeEphemerisWithSecondOrderDirectConversion< TimeType, StateScalarType >::TimeEphemerisWithSecondOrderDirectConversion(
        const TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
        const TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator,
        const std::string& centralBodyName,
        const std::function< Eigen::Vector6d( const TimeType ) > centralBodyStateFunction,
        const std::map< std::string, std::function< Eigen::Vector3d( const TimeType ) > >& groundStationPositionFunctions,
        const std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators,
        const std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators ) :
    TimeEphemerisFromPostNewtonianExpansion< TimeType, StateScalarType >(
        barycenterToPlanetCenterCoordinateTimeInterpolator,
        planetCenterToBarycenterCoordinateTimeInterpolator,
        centralBodyName,
        groundStationPositionFunctions,
        planetCoordinateToProperTimeInterpolators,
        properTimeToPlanetCoordinateInterpolators ),
    centralBodyStateFunction_( centralBodyStateFunction )
{ }

template< typename TimeType, typename StateScalarType >
double TimeEphemerisWithSecondOrderDirectConversion< TimeType, StateScalarType >::calculateDirectTimeDifferenceTerm(
        const Eigen::Vector3d positionVectorFromReferencePoint,
        const TimeType currentTime )
{
    throw std::runtime_error( "Error, second order time ephemeris direct difference term not yet implemented" );
    return 0.0;
}

}  // namespace tudat

#endif  // TUDAT_TIME_EPHEMERIS_WITH_SECOND_ORDER_DIRECT_CONVERSION_H
