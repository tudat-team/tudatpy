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

#include <Eigen/Core>

#include "tudat/astro/ephemerides/timeEphemeris.h"
#include "tudat/math/interpolators/oneDimensionalInterpolator.h"

namespace tudat
{

//! Relativistic time ephemeris using post-Newtonian conversion terms and optional direct corrections.
/*!
 *  This class combines pre-integrated coordinate-time differences with direct correction terms to build
 *  conversions between barycentric, body-centric, and topocentric proper time scales.
 */
class TimeEphemerisFromPostNewtonianExpansion : public TimeEphemeris
{
public:
    using TimeDifferenceInterpolator = std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > >;
    using ExtendedTimeDifferenceInterpolator = std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, double > >;

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

    ~TimeEphemerisFromPostNewtonianExpansion( ) override = default;

    void resetBarycentricToBodycentricInterpolators(
            const TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
            const TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator )
    {
        barycenterToPlanetCenterCoordinateTimeInterpolator_ = barycenterToPlanetCenterCoordinateTimeInterpolator;
        planetCenterToBarycenterCoordinateTimeInterpolator_ = planetCenterToBarycenterCoordinateTimeInterpolator;
    }

    void resetBarycentricToBodycentricInterpolators(
            const ExtendedTimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
            const ExtendedTimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator )
    {
        barycenterToPlanetCenterCoordinateTimeInterpolatorExtended_ = barycenterToPlanetCenterCoordinateTimeInterpolator;
        planetCenterToBarycenterCoordinateTimeInterpolatorExtended_ = planetCenterToBarycenterCoordinateTimeInterpolator;
    }

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
                throw std::runtime_error(
                            "Error when resetting bodycentric to topocentric time converter for point " + referencePoint +
                            ", must also provide point position function; point not yet known in TimeEphemeris" );
            }
            else
            {
                groundStationPositionFunctions_[ referencePoint ] = referencePointPositionFunction;
            }
        }
    }

    void resetBodycentricToTopocentricInterpolators(
            const ExtendedTimeDifferenceInterpolator planetCoordinateToProperTimeInterpolator,
            const ExtendedTimeDifferenceInterpolator properTimeToPlanetCoordinateInterpolator,
            const std::string& referencePoint,
            const std::function< Eigen::Vector3d( const double ) > referencePointPositionFunction = nullptr )
    {
        planetCoordinateToProperTimeInterpolatorsExtended_[ referencePoint ] = planetCoordinateToProperTimeInterpolator;
        properTimeToPlanetCoordinateInterpolatorsExtended_[ referencePoint ] = properTimeToPlanetCoordinateInterpolator;
        if( !doesReferencePointTopocentricConverterExist( referencePoint ) )
        {
            if( referencePointPositionFunction == nullptr )
            {
                throw std::runtime_error(
                            "Error when resetting bodycentric to topocentric time converter for point " + referencePoint +
                            ", must also provide point position function; point not yet known in TimeEphemeris" );
            }
            else
            {
                groundStationPositionFunctions_[ referencePoint ] = referencePointPositionFunction;
            }
        }
    }

    std::function< Time( const Time ) > getTimeDifferenceFunctionFromExtendedTime(
            const basic_astrodynamics::TimeScales inputScale,
            const basic_astrodynamics::TimeScales outputScale,
            const std::string& pointIdentifier = "" ) override;

    double calculateDirectTimeDifferenceTermFromFunction(
            const std::function< Eigen::Vector3d( const double ) > positionVectorFunctionFromReferencePoint,
            const double currentTime,
            const double multiplier = 1.0 )
    {
        return multiplier * calculateDirectTimeDifferenceTerm( positionVectorFunctionFromReferencePoint( currentTime ), currentTime );
    }

    virtual double calculateDirectTimeDifferenceTerm(
            const Eigen::Vector3d positionVectorFromReferencePoint,
            const double currentTime ) = 0;

    bool doesReferencePointTopocentricConverterExist( const std::string& referencePointName )
    {
        const bool hasPosition = ( groundStationPositionFunctions_.count( referencePointName ) > 0 );
        const bool hasDoubleInterpolators =
                ( planetCoordinateToProperTimeInterpolators_.count( referencePointName ) > 0 ) &&
                ( properTimeToPlanetCoordinateInterpolators_.count( referencePointName ) > 0 );
        const bool hasExtendedInterpolators =
                ( planetCoordinateToProperTimeInterpolatorsExtended_.count( referencePointName ) > 0 ) &&
                ( properTimeToPlanetCoordinateInterpolatorsExtended_.count( referencePointName ) > 0 );
        return hasPosition && ( hasDoubleInterpolators || hasExtendedInterpolators );
    }

protected:
    TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator_;
    TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator_;
    ExtendedTimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolatorExtended_;
    ExtendedTimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolatorExtended_;

    std::map< std::string, std::function< Eigen::Vector3d( const double ) > > groundStationPositionFunctions_;
    std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators_;
    std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators_;
    std::map< std::string, ExtendedTimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolatorsExtended_;
    std::map< std::string, ExtendedTimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolatorsExtended_;
};

}  // namespace tudat

#endif  // TUDAT_TIME_EPHEMERIS_FROM_POST_NEWTONIAN_EXPANSION_H
