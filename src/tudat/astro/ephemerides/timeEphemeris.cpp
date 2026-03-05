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

#include <functional>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <Eigen/Core>

#include "tudat/astro/ephemerides/timeEphemeris.h"

namespace tudat
{

using namespace basic_astrodynamics;

std::function< Time( const Time ) > TimeEphemerisFromPostNewtonianExpansion::getTimeDifferenceFunctionFromExtendedTime(
        const TimeScales inputScale, const TimeScales outputScale, const std::string& pointIdentifier )
{
    std::function< Time( const Time ) > timeDifferenceFunction;

    using namespace std::placeholders;

    if( !isTimeScaleRelativistic( inputScale ) || !isTimeScaleRelativistic( outputScale ) )
    {
        std::cerr<<"Error when getting relatvistic time conversion, scales are not both relativistic"<<std::endl;
    }
    else
    {
        const bool requiresBarycentricInterpolators =
                ( inputScale == barycentric_coordinate_time_scale || outputScale == barycentric_coordinate_time_scale );
        const bool hasDoubleBarycentricInterpolators =
                ( barycenterToPlanetCenterCoordinateTimeInterpolator_ != nullptr &&
                  planetCenterToBarycenterCoordinateTimeInterpolator_ != nullptr );
        const bool hasExtendedBarycentricInterpolators =
                ( barycenterToPlanetCenterCoordinateTimeInterpolatorExtended_ != nullptr &&
                  planetCenterToBarycenterCoordinateTimeInterpolatorExtended_ != nullptr );
        if( requiresBarycentricInterpolators && !( hasDoubleBarycentricInterpolators || hasExtendedBarycentricInterpolators ) )
        {
            throw std::runtime_error(
                        "Error in TimeEphemerisFromPostNewtonianExpansion: barycentric/bodycentric interpolators not initialized for " +
                        centralBodyName_ );
        }

        auto getBarycentricToBodycentricDifference = [=]( const Time currentTime )
        {
            if( barycenterToPlanetCenterCoordinateTimeInterpolatorExtended_ != nullptr )
            {
                return barycenterToPlanetCenterCoordinateTimeInterpolatorExtended_->interpolate( currentTime );
            }
            else
            {
                return barycenterToPlanetCenterCoordinateTimeInterpolator_->interpolate( currentTime.getSeconds< double >( ) );
            }
        };

        auto getBodycentricToBarycentricDifference = [=]( const Time currentTime )
        {
            if( planetCenterToBarycenterCoordinateTimeInterpolatorExtended_ != nullptr )
            {
                return planetCenterToBarycenterCoordinateTimeInterpolatorExtended_->interpolate( currentTime );
            }
            else
            {
                return planetCenterToBarycenterCoordinateTimeInterpolator_->interpolate( currentTime.getSeconds< double >( ) );
            }
        };

        auto getPlanetCoordinateToProperDifference = [=]( const std::string& point, const Time currentTime )
        {
            if( planetCoordinateToProperTimeInterpolatorsExtended_.count( point ) > 0 )
            {
                return planetCoordinateToProperTimeInterpolatorsExtended_.at( point )->interpolate( currentTime );
            }
            else
            {
                return planetCoordinateToProperTimeInterpolators_.at( point )->interpolate( currentTime.getSeconds< double >( ) );
            }
        };

        auto getProperToPlanetCoordinateDifference = [=]( const std::string& point, const Time currentTime )
        {
            if( properTimeToPlanetCoordinateInterpolatorsExtended_.count( point ) > 0 )
            {
                return properTimeToPlanetCoordinateInterpolatorsExtended_.at( point )->interpolate( currentTime );
            }
            else
            {
                return properTimeToPlanetCoordinateInterpolators_.at( point )->interpolate( currentTime.getSeconds< double >( ) );
            }
        };

        switch( inputScale )
        {
        case barycentric_coordinate_time_scale:
        {
            if( outputScale == body_centered_coordinate_time_scale )
            {
                timeDifferenceFunction = [=]( const Time currentTime )
                {
                    return Time( getBarycentricToBodycentricDifference( currentTime ) );
                };
            }
            else if( outputScale == local_proper_time_scale )
            {
                if( planetCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 &&
                    planetCoordinateToProperTimeInterpolatorsExtended_.count( pointIdentifier ) == 0 )
                {
                    std::cerr << "Error, body-point " << pointIdentifier << " not found in time ephemeris for body "
                              << centralBodyName_ << std::endl;
                }
                else
                {
                    std::vector< std::function< double( const Time ) > > timeDifferenceFunctions;
                    timeDifferenceFunctions.push_back( getBarycentricToBodycentricDifference );

                    std::function< Eigen::Vector3d( const double ) > stationPositionFunction =
                            groundStationPositionFunctions_.at( pointIdentifier );

                    timeDifferenceFunctions.push_back(
                                std::bind( &TimeEphemerisFromPostNewtonianExpansion::calculateDirectTimeDifferenceTermFromFunction, this,
                                           stationPositionFunction, _1, 1.0 ) );

                    timeDifferenceFunctions.push_back( [=]( const Time currentTime )
                    {
                        return getPlanetCoordinateToProperDifference( pointIdentifier, currentTime );
                    } );

                    std::vector< int > timeEvaluationIndex = { 0, 1, 1 };

                    timeDifferenceFunction = [=]( const Time currentTime )
                    {
                        return Time( combineTimeDifferenceFunction( timeDifferenceFunctions, timeEvaluationIndex, currentTime ) );
                    };
                }
            }
            else
            {
                std::cerr<<"Error relativistic output scale not found"<<std::endl;
            }
            break;
        }

        case body_centered_coordinate_time_scale:
        {
            if( outputScale == barycentric_coordinate_time_scale )
            {
                timeDifferenceFunction = [=]( const Time currentTime )
                {
                    return Time( getBodycentricToBarycentricDifference( currentTime ) );
                };
            }
            else if( outputScale == local_proper_time_scale )
            {
                if( planetCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 &&
                    planetCoordinateToProperTimeInterpolatorsExtended_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found"<<std::endl;
                }
                else
                {
                    timeDifferenceFunction = [=]( const Time currentTime )
                    {
                        return Time( getPlanetCoordinateToProperDifference( pointIdentifier, currentTime ) );
                    };
                }
            }
            else
            {
                std::cerr<<"Error relativistic output scale not found"<<std::endl;
            }
            break;
        }

        case local_proper_time_scale:
        {
            if( outputScale == barycentric_coordinate_time_scale )
            {
                if( properTimeToPlanetCoordinateInterpolators_.count( pointIdentifier ) == 0 &&
                    properTimeToPlanetCoordinateInterpolatorsExtended_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found"<<std::endl;
                }
                else
                {
                    std::vector< std::function< double( const Time ) > > timeDifferenceFunctions;
                    timeDifferenceFunctions.push_back( [=]( const Time currentTime )
                    {
                        return getProperToPlanetCoordinateDifference( pointIdentifier, currentTime );
                    } );

                    std::function< Eigen::Vector3d( const double ) > stationPositionFunction =
                            groundStationPositionFunctions_.at( pointIdentifier );

                    timeDifferenceFunctions.push_back(
                                std::bind( &TimeEphemerisFromPostNewtonianExpansion::calculateDirectTimeDifferenceTermFromFunction,
                                           this, stationPositionFunction, _1, -1.0 ) );

                    timeDifferenceFunctions.push_back( getBodycentricToBarycentricDifference );

                    std::vector< int > timeEvaluationIndex = { 0, 1, 2 };

                    timeDifferenceFunction = [=]( const Time currentTime )
                    {
                        return Time( combineTimeDifferenceFunction( timeDifferenceFunctions, timeEvaluationIndex, currentTime ) );
                    };
                }
            }
            else if( outputScale == body_centered_coordinate_time_scale )
            {
                if( properTimeToPlanetCoordinateInterpolators_.count( pointIdentifier ) == 0 &&
                    properTimeToPlanetCoordinateInterpolatorsExtended_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found"<<std::endl;
                }
                else
                {
                    timeDifferenceFunction = [=]( const Time currentTime )
                    {
                        return Time( getProperToPlanetCoordinateDifference( pointIdentifier, currentTime ) );
                    };
                }
            }
            else
            {
                std::cerr<<"Error relativistic output scale not found"<<std::endl;
            }
            break;
        }
            default:
            throw std::runtime_error("Error, case not found in getTimeDifferenceFunctionFromExtendedTime.");
        }
    }

    return timeDifferenceFunction;
}


std::function< Time( const Time ) > TimeEphemerisDirectFromMetric::getTimeDifferenceFunctionFromExtendedTime(
        const TimeScales inputScale, const TimeScales outputScale, const std::string& pointIdentifier )
{
    std::function< Time( const Time ) > timeDifferenceFunction;

    if( !isTimeScaleRelativistic( inputScale ) || !isTimeScaleRelativistic( outputScale ) )
    {
        std::cerr<<"Error when getting relatvistic time conversion, scales are not both relativistic"<<std::endl;
    }
    else
    {
        if( inputScale == body_centered_coordinate_time_scale )
        {
            std::cerr<<"Cannot do body-centered input time scale for direct-from-metric time ephemeris"<<std::endl;
        }

        if( outputScale == body_centered_coordinate_time_scale )
        {
            std::cerr<<"Cannot do body-centered output time scale for direct-from-metric time ephemeris"<<std::endl;
        }

        if( inputScale == barycentric_coordinate_time_scale )
        {
            if( outputScale == local_proper_time_scale )
            {
                if( globalCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 &&
                    globalCoordinateToProperTimeInterpolatorsExtended_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found in direct-from-metric time ephemeris"<<std::endl;
                }
                else
                {
                    timeDifferenceFunction = [=]( const Time currentTime )
                    {
                        if( globalCoordinateToProperTimeInterpolatorsExtended_.count( pointIdentifier ) > 0 )
                        {
                            return Time( globalCoordinateToProperTimeInterpolatorsExtended_.at( pointIdentifier )->interpolate( currentTime ) );
                        }
                        else
                        {
                            return Time( getGlobalCoordinateToProperTimeDifference( pointIdentifier, currentTime.getSeconds< double >( ) ) );
                        }
                    };
                }
            }
            else
            {
                std::cerr<<"Error A when using direct from metric time ephemeris"<<std::endl;
            }
        }
        else if( inputScale == local_proper_time_scale )
        {
            if( outputScale == barycentric_coordinate_time_scale )
            {
                if( properTimeToGlobalCoordinateInterpolators_.count( pointIdentifier ) == 0 &&
                    properTimeToGlobalCoordinateInterpolatorsExtended_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found in direct-from-metric time ephemeris"<<std::endl;
                }
                else
                {
                    timeDifferenceFunction = [=]( const Time currentTime )
                    {
                        if( properTimeToGlobalCoordinateInterpolatorsExtended_.count( pointIdentifier ) > 0 )
                        {
                            return Time( properTimeToGlobalCoordinateInterpolatorsExtended_.at( pointIdentifier )->interpolate( currentTime ) );
                        }
                        else
                        {
                            return Time( getProperToGlobalCoordinateTimeDifference( pointIdentifier, currentTime.getSeconds< double >( ) ) );
                        }
                    };
                }
            }
            else
            {
                std::cerr<<"Error B when using direct from metric time ephemeris"<<std::endl;
            }
        }
    }

    return timeDifferenceFunction;
}

} // namespace tudat
