#include <functional>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <Eigen/Core>

#include "tudat/astro/ephemerides/timeEphemeris.h"

namespace tudat
{

double combineTimeDifferenceFunction( const std::vector< std::function< double( const double ) > > timeDifferenceFunctions,
                                      const std::vector< int > evaluationStepIndices,
                                      const double time )
{
    double timeDifference = 0.0;
    std::vector< double > intermediateTimes;
    intermediateTimes.push_back( time );

    for( unsigned int i = 0; i < timeDifferenceFunctions.size( ); i++ )
    {
        timeDifference += timeDifferenceFunctions.at( i )( intermediateTimes.at( evaluationStepIndices.at( i ) ) );
        intermediateTimes.push_back( time + timeDifference );
    }
    return timeDifference;
}


std::function< double( const double ) > TimeEphemerisFromPostNewtonianExpansion::getTimeDifferenceFunction(
        const TimeScales inputScale, const TimeScales outputScale, const std::string pointIdentifier )
{
    std::function< double( const double ) > timeDifferenceFunction;

    typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;
    using namespace std::placeholders;

    if( !isTimeScaleRelativistic( inputScale ) || !isTimeScaleRelativistic( outputScale ) )
    {
        std::cerr<<"Error when getting relatvistic time conversion, scales are not both relativistic"<<std::endl;
    }
    else
    {
        const bool requiresBarycentricInterpolators =
                ( inputScale == barycentric_coordinate_time_scale || outputScale == barycentric_coordinate_time_scale );
        if( requiresBarycentricInterpolators &&
            ( barycenterToPlanetCenterCoordinateTimeInterpolator_ == nullptr ||
              planetCenterToBarycenterCoordinateTimeInterpolator_ == nullptr ) )
        {
            throw std::runtime_error(
                        "Error in TimeEphemerisFromPostNewtonianExpansion: barycentric/bodycentric interpolators not initialized for " +
                        centralBodyName_ );
        }

        switch( inputScale )
        {
        case barycentric_coordinate_time_scale:
        {
            if( outputScale == body_centered_coordinate_time_scale )
            {
                timeDifferenceFunction = std::bind(
                            static_cast< double( LocalInterpolator::* )( const double ) >
                            ( &LocalInterpolator::interpolate ), barycenterToPlanetCenterCoordinateTimeInterpolator_, _1 );
            }
            else if( outputScale == local_proper_time_scale )
            {
                if( planetCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found "<<
                               planetCoordinateToProperTimeInterpolators_.size( )<<" "<<planetCoordinateToProperTimeInterpolators_.begin( )->first<<std::endl;
                }
                else
                {
                    std::vector< std::function< double( const double ) > > timeDifferenceFunctions;
                    timeDifferenceFunctions.push_back(
                                std::bind(
                                    static_cast< double( LocalInterpolator::* )( const double ) >
                                    ( &LocalInterpolator::interpolate ), barycenterToPlanetCenterCoordinateTimeInterpolator_, _1 ) );

                    std::function< Eigen::Vector3d( const double ) > stationPositionFunction =
                            groundStationPositionFunctions_.at( pointIdentifier );

                    timeDifferenceFunctions.push_back(
                                std::bind( &TimeEphemerisFromPostNewtonianExpansion::calculateDirectTimeDifferenceTermFromFunction, this,
                                           stationPositionFunction, _1, 1.0 ) );

                    timeDifferenceFunctions.push_back(
                                std::bind(
                                    static_cast< double( LocalInterpolator::* )( const double ) >
                                    ( &LocalInterpolator::interpolate ), planetCoordinateToProperTimeInterpolators_.at( pointIdentifier ), _1 ) );

                    std::vector< int > timeEvaluationIndex = { 0, 1, 1 };

                    timeDifferenceFunction = std::bind( combineTimeDifferenceFunction, timeDifferenceFunctions, timeEvaluationIndex, _1 );
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
                timeDifferenceFunction = std::bind(
                            static_cast< double( LocalInterpolator::* )( const double ) >
                            ( &LocalInterpolator::interpolate ), planetCenterToBarycenterCoordinateTimeInterpolator_, _1 );
            }
            else if( outputScale == local_proper_time_scale )
            {
                if( planetCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found"<<std::endl;
                }
                else
                {
                    timeDifferenceFunction = std::bind(
                                static_cast< double( LocalInterpolator::* )( const double ) >
                                ( &LocalInterpolator::interpolate ), planetCoordinateToProperTimeInterpolators_.at( pointIdentifier ) , _1 );
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
                if( planetCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found"<<std::endl;
                }
                else
                {
                    std::vector< std::function< double( const double ) > > timeDifferenceFunctions;
                    timeDifferenceFunctions.push_back(
                                std::bind(
                                    static_cast< double( LocalInterpolator::* )( const double ) >
                                    ( &LocalInterpolator::interpolate ), properTimeToPlanetCoordinateInterpolators_.at( pointIdentifier ), _1 ) );

                    std::function< Eigen::Vector3d( const double ) > stationPositionFunction =
                            groundStationPositionFunctions_.at( pointIdentifier );

                    timeDifferenceFunctions.push_back(
                                std::bind( &TimeEphemerisFromPostNewtonianExpansion::calculateDirectTimeDifferenceTermFromFunction,
                                           this, stationPositionFunction, _1, -1.0 ) );

                    timeDifferenceFunctions.push_back(
                                std::bind(
                                    static_cast< double( LocalInterpolator::* )( const double ) >
                                    ( &LocalInterpolator::interpolate ), planetCenterToBarycenterCoordinateTimeInterpolator_, _1 ) );

                    std::vector< int > timeEvaluationIndex = { 0, 1, 2 };

                    timeDifferenceFunction = std::bind( combineTimeDifferenceFunction, timeDifferenceFunctions, timeEvaluationIndex, _1 );
                }
            }
            else if( outputScale == body_centered_coordinate_time_scale )
            {
                if( planetCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found"<<std::endl;
                }
                else
                {
                    timeDifferenceFunction = std::bind(
                                static_cast< double( LocalInterpolator::* )( const double ) >
                                ( &LocalInterpolator::interpolate ), properTimeToPlanetCoordinateInterpolators_.at( pointIdentifier ), _1 );
                }
            }
            else
            {
                std::cerr<<"Error relativistic output scale not found"<<std::endl;
            }
            break;
        }
        }
    }

    return timeDifferenceFunction;
}


std::function< double( const double ) > TimeEphemerisDirectFromMetric::getTimeDifferenceFunction(
        const TimeScales inputScale, const TimeScales outputScale, const std::string pointIdentifier )
{
    std::function< double( const double ) > timeDifferenceFunction;
    using namespace std::placeholders;

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
                if( globalCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found in direct-from-metric time ephemeris"<<std::endl;
                }
                else
                {
                    timeDifferenceFunction = std::bind(
                                &TimeEphemerisDirectFromMetric::getGlobalCoordinateToProperTimeDifference,
                                this, pointIdentifier, _1 );
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
                if( properTimeToGlobalCoordinateInterpolators_.count( pointIdentifier ) == 0 )
                {
                    std::cerr<<"Error, body-point "<<pointIdentifier<<" not found in direct-from-metric time ephemeris"<<std::endl;
                }
                else
                {
                    timeDifferenceFunction = std::bind(
                                &TimeEphemerisDirectFromMetric::getProperToGlobalCoordinateTimeDifference,
                                this, pointIdentifier, _1 );
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
