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

#ifndef TUDAT_TIME_EPHEMERIS_WITH_FIRST_ORDER_DIRECT_CONVERSION_H
#define TUDAT_TIME_EPHEMERIS_WITH_FIRST_ORDER_DIRECT_CONVERSION_H

#include <functional>
#include <map>
#include <string>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ephemerides/timeEphemerisFromPostNewtonianExpansion.h"

namespace tudat
{

//! Post-Newtonian converter with first-order direct local term.
template< typename TimeType = double, typename StateScalarType = double >
class TimeEphemerisWithFirstOrderDirectConversion : public TimeEphemerisFromPostNewtonianExpansion< TimeType, StateScalarType >
{
public:
    using TimeDifferenceInterpolator =
            typename TimeEphemerisFromPostNewtonianExpansion< TimeType, StateScalarType >::TimeDifferenceInterpolator;
    using TimeDifferenceDataMap = typename TimeEphemerisFromPostNewtonianExpansion< TimeType, StateScalarType >::TimeDifferenceDataMap;

    TimeEphemerisWithFirstOrderDirectConversion(
            const std::string& centralBodyName,
            const std::function< Eigen::Vector6d( const TimeType ) > centralBodyStateFunction,
            const TimeDifferenceDataMap& barycenterToPlanetCenterCoordinateTimeDifferences = TimeDifferenceDataMap( ),
            const std::map< std::string, std::function< Eigen::Vector3d( const TimeType ) > >& groundStationPositionFunctions =
                    ( std::map< std::string, std::function< Eigen::Vector3d( const TimeType ) > >( ) ),
            const std::map< std::string, TimeDifferenceDataMap > planetCoordinateToProperTimeDifferences =
                    ( std::map< std::string, TimeDifferenceDataMap >( ) ),
            const std::shared_ptr< interpolators::InterpolatorSettings >& timeDifferenceInterpolatorSettings =
                    TimeEphemerisFromPostNewtonianExpansion< TimeType, StateScalarType >::getDefaultTimeDifferenceInterpolatorSettings( ) );

    double calculateDirectTimeDifferenceTerm( const Eigen::Vector3d positionVectorFromReferencePoint, const TimeType currentTime ) override;

private:
    std::function< Eigen::Vector6d( const TimeType ) > centralBodyStateFunction_;
};

template< typename TimeType, typename StateScalarType >
TimeEphemerisWithFirstOrderDirectConversion< TimeType, StateScalarType >::TimeEphemerisWithFirstOrderDirectConversion(
        const std::string& centralBodyName,
        const std::function< Eigen::Vector6d( const TimeType ) > centralBodyStateFunction,
        const TimeDifferenceDataMap& barycenterToPlanetCenterCoordinateTimeDifferences,
        const std::map< std::string, std::function< Eigen::Vector3d( const TimeType ) > >& groundStationPositionFunctions,
        const std::map< std::string, TimeDifferenceDataMap > planetCoordinateToProperTimeDifferences,
        const std::shared_ptr< interpolators::InterpolatorSettings >& timeDifferenceInterpolatorSettings ):
    TimeEphemerisFromPostNewtonianExpansion< TimeType, StateScalarType >( centralBodyName,
                                                                          barycenterToPlanetCenterCoordinateTimeDifferences,
                                                                          groundStationPositionFunctions,
                                                                          planetCoordinateToProperTimeDifferences,
                                                                          timeDifferenceInterpolatorSettings ),
    centralBodyStateFunction_( centralBodyStateFunction )
{}

template< typename TimeType, typename StateScalarType >
double TimeEphemerisWithFirstOrderDirectConversion< TimeType, StateScalarType >::calculateDirectTimeDifferenceTerm(
        const Eigen::Vector3d positionVectorFromReferencePoint,
        const TimeType currentTime )
{
    return -centralBodyStateFunction_( currentTime ).segment( 3, 3 ).dot( positionVectorFromReferencePoint ) *
            physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
}

}  // namespace tudat

#endif  // TUDAT_TIME_EPHEMERIS_WITH_FIRST_ORDER_DIRECT_CONVERSION_H
