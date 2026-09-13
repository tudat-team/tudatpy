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

#ifndef TUDAT_TIME_EPHEMERIS_FROM_POST_NEWTONIAN_EXPANSION_H
#define TUDAT_TIME_EPHEMERIS_FROM_POST_NEWTONIAN_EXPANSION_H

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/ephemerides/timeEphemeris.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/math/interpolators/oneDimensionalInterpolator.h"

namespace tudat
{

//! Relativistic time ephemeris using post-Newtonian conversion terms and optional direct corrections.
/*!
 *  This class combines pre-integrated coordinate-time differences with direct correction terms to build
 *  conversions between barycentric, body-centric, and topocentric proper time scales.
 */
template< typename TimeType = double, typename StateScalarType = double >
class TimeEphemerisFromPostNewtonianExpansion : public TimeEphemeris
{
public:
    using TimeDifferenceInterpolator = std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateScalarType > >;
    using TimeDifferenceDataMap = std::map< TimeType, StateScalarType >;

    static std::shared_ptr< interpolators::InterpolatorSettings > getDefaultTimeDifferenceInterpolatorSettings( )
    {
        return std::make_shared< interpolators::InterpolatorSettings >(
                interpolators::cubic_spline_interpolator, interpolators::huntingAlgorithm, false, interpolators::extrapolate_at_boundary );
    }

    TimeEphemerisFromPostNewtonianExpansion(
            const std::string& centralBodyName,
            const TimeDifferenceDataMap& barycenterToPlanetCenterCoordinateTimeDifferences = TimeDifferenceDataMap( ),
            const std::map< std::string, std::function< Eigen::Vector3d( const TimeType ) > >& groundStationPositionFunctions =
                    ( std::map< std::string, std::function< Eigen::Vector3d( const TimeType ) > >( ) ),
            const std::map< std::string, TimeDifferenceDataMap > planetCoordinateToProperTimeDifferences =
                    ( std::map< std::string, TimeDifferenceDataMap >( ) ),
            const std::shared_ptr< interpolators::InterpolatorSettings >& timeDifferenceInterpolatorSettings =
                    getDefaultTimeDifferenceInterpolatorSettings( ) ):
        TimeEphemeris( centralBodyName ), groundStationPositionFunctions_( groundStationPositionFunctions ),
        timeDifferenceInterpolatorSettings_( timeDifferenceInterpolatorSettings )
    {
        if( timeDifferenceInterpolatorSettings_ == nullptr )
        {
            timeDifferenceInterpolatorSettings_ = getDefaultTimeDifferenceInterpolatorSettings( );
        }

        if( barycenterToPlanetCenterCoordinateTimeDifferences.size( ) > 0 )
        {
            resetBarycentricToBodycentricInterpolators( barycenterToPlanetCenterCoordinateTimeDifferences );
        }

        for( const auto& pointDataPair : planetCoordinateToProperTimeDifferences )
        {
            if( groundStationPositionFunctions_.count( pointDataPair.first ) == 0 )
            {
                throw std::runtime_error( "Error when initializing bodycentric to topocentric time converter for point " +
                                          pointDataPair.first +
                                          ", must also provide point position function; point not yet known in TimeEphemeris" );
            }

            resetBodycentricToTopocentricInterpolators(
                    pointDataPair.second, pointDataPair.first, groundStationPositionFunctions_.at( pointDataPair.first ) );
        }
    }

    ~TimeEphemerisFromPostNewtonianExpansion( ) override = default;

    void resetBarycentricToBodycentricInterpolators(
            const TimeDifferenceDataMap& barycenterToPlanetCenterCoordinateTimeDifferences,
            const std::shared_ptr< interpolators::InterpolatorSettings >& timeDifferenceInterpolatorSettings = nullptr )
    {
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettingsToUse = ( timeDifferenceInterpolatorSettings == nullptr )
                ? timeDifferenceInterpolatorSettings_
                : timeDifferenceInterpolatorSettings;

        if( interpolatorSettingsToUse == nullptr )
        {
            throw std::runtime_error( "Error when resetting barycentric/bodycentric time interpolators, interpolator settings are null." );
        }

        std::pair< TimeDifferenceInterpolator, TimeDifferenceInterpolator > timeInterpolators =
                createRelativisticTimeInterpolators( barycenterToPlanetCenterCoordinateTimeDifferences, interpolatorSettingsToUse );

        barycenterToPlanetCenterCoordinateTimeInterpolator_ = timeInterpolators.first;
        planetCenterToBarycenterCoordinateTimeInterpolator_ = timeInterpolators.second;
    }

    void resetBodycentricToTopocentricInterpolators(
            const TimeDifferenceDataMap& planetCoordinateToProperTimeDifferences,
            const std::string& referencePoint,
            const std::function< Eigen::Vector3d( const TimeType ) > referencePointPositionFunction = nullptr,
            const std::shared_ptr< interpolators::InterpolatorSettings >& timeDifferenceInterpolatorSettings = nullptr )
    {
        std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettingsToUse = ( timeDifferenceInterpolatorSettings == nullptr )
                ? timeDifferenceInterpolatorSettings_
                : timeDifferenceInterpolatorSettings;

        if( interpolatorSettingsToUse == nullptr )
        {
            throw std::runtime_error( "Error when resetting bodycentric/topocentric time interpolators for point " + referencePoint +
                                      ", interpolator settings are null." );
        }

        std::pair< TimeDifferenceInterpolator, TimeDifferenceInterpolator > timeInterpolators =
                createRelativisticTimeInterpolators( planetCoordinateToProperTimeDifferences, interpolatorSettingsToUse );

        planetCoordinateToProperTimeInterpolators_[ referencePoint ] = timeInterpolators.first;
        properTimeToPlanetCoordinateInterpolators_[ referencePoint ] = timeInterpolators.second;

        if( !doesReferencePointTopocentricConverterExist( referencePoint ) )
        {
            if( referencePointPositionFunction == nullptr )
            {
                throw std::runtime_error( "Error when resetting bodycentric to topocentric time converter for point " + referencePoint +
                                          ", must also provide point position function; point not yet known in TimeEphemeris" );
            }
            else
            {
                groundStationPositionFunctions_[ referencePoint ] = referencePointPositionFunction;
            }
        }
    }

    std::function< Time( const Time ) > getTimeDifferenceFunctionFromExtendedTime( const basic_astrodynamics::TimeScales inputScale,
                                                                                   const basic_astrodynamics::TimeScales outputScale,
                                                                                   const std::string& pointIdentifier = "" ) override;

    double calculateDirectTimeDifferenceTermFromFunction(
            const std::function< Eigen::Vector3d( const TimeType ) > positionVectorFunctionFromReferencePoint,
            const TimeType currentTime,
            const double multiplier = 1.0 )
    {
        return multiplier * calculateDirectTimeDifferenceTerm( positionVectorFunctionFromReferencePoint( currentTime ), currentTime );
    }

    virtual double calculateDirectTimeDifferenceTerm( const Eigen::Vector3d positionVectorFromReferencePoint,
                                                      const TimeType currentTime ) = 0;

    bool doesReferencePointTopocentricConverterExist( const std::string& referencePointName )
    {
        const bool hasPosition = ( groundStationPositionFunctions_.count( referencePointName ) > 0 );
        const bool hasInterpolators = ( planetCoordinateToProperTimeInterpolators_.count( referencePointName ) > 0 ) &&
                ( properTimeToPlanetCoordinateInterpolators_.count( referencePointName ) > 0 );
        return hasPosition && hasInterpolators;
    }

protected:
    static TimeType addTimeDifferenceToEpoch( const TimeType& currentTime, const StateScalarType& timeDifference )
    {
        if constexpr( std::is_same_v< TimeType, Time > )
        {
            return currentTime + static_cast< long double >( timeDifference );
        }
        else
        {
            return static_cast< TimeType >( static_cast< long double >( currentTime ) + static_cast< long double >( timeDifference ) );
        }
    }

    std::pair< TimeDifferenceInterpolator, TimeDifferenceInterpolator > createRelativisticTimeInterpolators(
            const TimeDifferenceDataMap& originalToTargetTimeMap,
            const std::shared_ptr< interpolators::InterpolatorSettings >& interpolatorSettings )
    {
        if( originalToTargetTimeMap.size( ) == 0 )
        {
            return std::make_pair( TimeDifferenceInterpolator( ), TimeDifferenceInterpolator( ) );
        }

        std::map< TimeType, StateScalarType > targetToOriginalTimeMap;
        for( const auto& mapEntry : originalToTargetTimeMap )
        {
            targetToOriginalTimeMap[ addTimeDifferenceToEpoch( mapEntry.first, mapEntry.second ) ] = -mapEntry.second;
        }

        return std::make_pair( interpolators::createOneDimensionalInterpolator< TimeType, StateScalarType >( originalToTargetTimeMap,
                                                                                                             interpolatorSettings ),
                               interpolators::createOneDimensionalInterpolator< TimeType, StateScalarType >( targetToOriginalTimeMap,
                                                                                                             interpolatorSettings ) );
    }

    TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator_;
    TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator_;

    std::map< std::string, std::function< Eigen::Vector3d( const TimeType ) > > groundStationPositionFunctions_;
    std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators_;
    std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators_;
    std::shared_ptr< interpolators::InterpolatorSettings > timeDifferenceInterpolatorSettings_;
};

template< typename TimeType, typename StateScalarType >
std::function< Time( const Time ) >
TimeEphemerisFromPostNewtonianExpansion< TimeType, StateScalarType >::getTimeDifferenceFunctionFromExtendedTime(
        const basic_astrodynamics::TimeScales inputScale,
        const basic_astrodynamics::TimeScales outputScale,
        const std::string& pointIdentifier )
{
    std::function< Time( const Time ) > timeDifferenceFunction;

    if( !basic_astrodynamics::isTimeScaleRelativistic( inputScale ) || !basic_astrodynamics::isTimeScaleRelativistic( outputScale ) )
    {
        throw std::runtime_error( "Error when getting relatvistic time conversion, scales are not both relativistic" );
    }
    else
    {
        const bool requiresBarycentricInterpolators = ( inputScale == basic_astrodynamics::barycentric_coordinate_time_scale ||
                                                        outputScale == basic_astrodynamics::barycentric_coordinate_time_scale );
        const bool hasBarycentricInterpolators = ( barycenterToPlanetCenterCoordinateTimeInterpolator_ != nullptr &&
                                                   planetCenterToBarycenterCoordinateTimeInterpolator_ != nullptr );
        if( requiresBarycentricInterpolators && !hasBarycentricInterpolators )
        {
            throw std::runtime_error(
                    "Error in TimeEphemerisFromPostNewtonianExpansion: barycentric/bodycentric interpolators not initialized for " +
                    centralBodyName_ );
        }

        auto getBarycentricToBodycentricDifference = [ = ]( const Time currentTime ) {
            return Time( barycenterToPlanetCenterCoordinateTimeInterpolator_->interpolate(
                    time_ephemeris_detail::convertInterpolatorTimeFromExtendedTime< TimeType >( currentTime ) ) );
        };

        auto getBodycentricToBarycentricDifference = [ = ]( const Time currentTime ) {
            return Time( planetCenterToBarycenterCoordinateTimeInterpolator_->interpolate(
                    time_ephemeris_detail::convertInterpolatorTimeFromExtendedTime< TimeType >( currentTime ) ) );
        };

        auto getPlanetCoordinateToProperDifference = [ = ]( const std::string& point, const Time currentTime ) {
            return Time( planetCoordinateToProperTimeInterpolators_.at( point )->interpolate(
                    time_ephemeris_detail::convertInterpolatorTimeFromExtendedTime< TimeType >( currentTime ) ) );
        };

        auto getProperToPlanetCoordinateDifference = [ = ]( const std::string& point, const Time currentTime ) {
            return Time( properTimeToPlanetCoordinateInterpolators_.at( point )->interpolate(
                    time_ephemeris_detail::convertInterpolatorTimeFromExtendedTime< TimeType >( currentTime ) ) );
        };

        switch( inputScale )
        {
            case basic_astrodynamics::barycentric_coordinate_time_scale: {
                if( outputScale == basic_astrodynamics::body_centered_coordinate_time_scale )
                {
                    timeDifferenceFunction = [ = ]( const Time currentTime ) {
                        return getBarycentricToBodycentricDifference( currentTime );
                    };
                }
                else if( outputScale == basic_astrodynamics::local_proper_time_scale )
                {
                    if( planetCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 )
                    {
                        throw std::runtime_error( "Error, body-point " + pointIdentifier + " not found in time ephemeris for body " +
                                                  centralBodyName_ );
                    }
                    else
                    {
                        std::vector< std::function< Time( const Time ) > > timeDifferenceFunctions;
                        timeDifferenceFunctions.push_back( getBarycentricToBodycentricDifference );

                        std::function< Eigen::Vector3d( const TimeType ) > stationPositionFunction =
                                groundStationPositionFunctions_.at( pointIdentifier );

                        timeDifferenceFunctions.push_back( [ = ]( const Time currentTime ) {
                            return Time( this->calculateDirectTimeDifferenceTermFromFunction(
                                    stationPositionFunction,
                                    time_ephemeris_detail::convertInterpolatorTimeFromExtendedTime< TimeType >( currentTime ),
                                    1.0 ) );
                        } );

                        timeDifferenceFunctions.push_back( [ = ]( const Time currentTime ) {
                            return getPlanetCoordinateToProperDifference( pointIdentifier, currentTime );
                        } );

                        std::vector< int > timeEvaluationIndex = { 0, 1, 1 };

                        timeDifferenceFunction = [ = ]( const Time currentTime ) {
                            return combineTimeDifferenceFunction( timeDifferenceFunctions, timeEvaluationIndex, currentTime );
                        };
                    }
                }
                else
                {
                    throw std::runtime_error( "Error relativistic output scale not found" );
                }
                break;
            }

            case basic_astrodynamics::body_centered_coordinate_time_scale: {
                if( outputScale == basic_astrodynamics::barycentric_coordinate_time_scale )
                {
                    timeDifferenceFunction = [ = ]( const Time currentTime ) {
                        return getBodycentricToBarycentricDifference( currentTime );
                    };
                }
                else if( outputScale == basic_astrodynamics::local_proper_time_scale )
                {
                    if( planetCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 )
                    {
                        throw std::runtime_error( "Error, body-point " + pointIdentifier + " not found" );
                    }
                    else
                    {
                        timeDifferenceFunction = [ = ]( const Time currentTime ) {
                            return getPlanetCoordinateToProperDifference( pointIdentifier, currentTime );
                        };
                    }
                }
                else
                {
                    throw std::runtime_error( "Error relativistic output scale not found" );
                }
                break;
            }

            case basic_astrodynamics::local_proper_time_scale: {
                if( outputScale == basic_astrodynamics::barycentric_coordinate_time_scale )
                {
                    if( properTimeToPlanetCoordinateInterpolators_.count( pointIdentifier ) == 0 )
                    {
                        throw std::runtime_error( "Error, body-point " + pointIdentifier + " not found" );
                    }
                    else
                    {
                        std::vector< std::function< Time( const Time ) > > timeDifferenceFunctions;
                        timeDifferenceFunctions.push_back( [ = ]( const Time currentTime ) {
                            return getProperToPlanetCoordinateDifference( pointIdentifier, currentTime );
                        } );

                        std::function< Eigen::Vector3d( const TimeType ) > stationPositionFunction =
                                groundStationPositionFunctions_.at( pointIdentifier );

                        timeDifferenceFunctions.push_back( [ = ]( const Time currentTime ) {
                            return Time( this->calculateDirectTimeDifferenceTermFromFunction(
                                    stationPositionFunction,
                                    time_ephemeris_detail::convertInterpolatorTimeFromExtendedTime< TimeType >( currentTime ),
                                    -1.0 ) );
                        } );

                        timeDifferenceFunctions.push_back( getBodycentricToBarycentricDifference );

                        std::vector< int > timeEvaluationIndex = { 0, 1, 2 };

                        timeDifferenceFunction = [ = ]( const Time currentTime ) {
                            return combineTimeDifferenceFunction( timeDifferenceFunctions, timeEvaluationIndex, currentTime );
                        };
                    }
                }
                else if( outputScale == basic_astrodynamics::body_centered_coordinate_time_scale )
                {
                    if( properTimeToPlanetCoordinateInterpolators_.count( pointIdentifier ) == 0 )
                    {
                        throw std::runtime_error( "Error, body-point " + pointIdentifier + " not found" );
                    }
                    else
                    {
                        timeDifferenceFunction = [ = ]( const Time currentTime ) {
                            return getProperToPlanetCoordinateDifference( pointIdentifier, currentTime );
                        };
                    }
                }
                else
                {
                    throw std::runtime_error( "Error relativistic output scale not found" );
                }
                break;
            }
            default:
                throw std::runtime_error( "Error, case not found in getTimeDifferenceFunctionFromExtendedTime." );
        }
    }

    return timeDifferenceFunction;
}

}  // namespace tudat

#endif  // TUDAT_TIME_EPHEMERIS_FROM_POST_NEWTONIAN_EXPANSION_H
