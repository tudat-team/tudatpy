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

#ifndef TUDAT_TIME_EPHEMERIS_DIRECT_FROM_METRIC_H
#define TUDAT_TIME_EPHEMERIS_DIRECT_FROM_METRIC_H

#include <functional>
#include <map>
#include <memory>
#include <string>

#include "tudat/astro/ephemerides/timeEphemeris.h"
#include "tudat/math/interpolators/oneDimensionalInterpolator.h"

namespace tudat
{

//! Time ephemeris using direct metric-based coordinate-time/proper-time interpolation.
class TimeEphemerisDirectFromMetric : public TimeEphemeris
{
public:
    using TimeDifferenceInterpolator = std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > >;
    using ExtendedTimeDifferenceInterpolator = std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, double > >;

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

    void resetGlobalToProperTimeInterpolators(
            const TimeDifferenceInterpolator globalCoordinateToProperTimeInterpolator,
            const TimeDifferenceInterpolator properTimeToGlobalCoordinateInterpolator,
            const std::string& referencePoint )
    {
        globalCoordinateToProperTimeInterpolators_[ referencePoint ] = globalCoordinateToProperTimeInterpolator;
        properTimeToGlobalCoordinateInterpolators_[ referencePoint ] = properTimeToGlobalCoordinateInterpolator;
    }

    void resetGlobalToProperTimeInterpolators(
            const ExtendedTimeDifferenceInterpolator globalCoordinateToProperTimeInterpolator,
            const ExtendedTimeDifferenceInterpolator properTimeToGlobalCoordinateInterpolator,
            const std::string& referencePoint )
    {
        globalCoordinateToProperTimeInterpolatorsExtended_[ referencePoint ] = globalCoordinateToProperTimeInterpolator;
        properTimeToGlobalCoordinateInterpolatorsExtended_[ referencePoint ] = properTimeToGlobalCoordinateInterpolator;
    }

    std::function< Time( const Time ) > getTimeDifferenceFunctionFromExtendedTime(
            const basic_astrodynamics::TimeScales inputScale,
            const basic_astrodynamics::TimeScales outputScale,
            const std::string& pointIdentifier = "" ) override;

    double getGlobalCoordinateToProperTimeDifference( const std::string& pointId, const double evaluationTime )
    {
        return globalCoordinateToProperTimeInterpolators_.at( pointId )->interpolate( evaluationTime );
    }

    double getProperToGlobalCoordinateTimeDifference( const std::string& pointId, const double evaluationTime )
    {
        return properTimeToGlobalCoordinateInterpolators_.at( pointId )->interpolate( evaluationTime );
    }

protected:
    std::map< std::string, TimeDifferenceInterpolator > globalCoordinateToProperTimeInterpolators_;
    std::map< std::string, TimeDifferenceInterpolator > properTimeToGlobalCoordinateInterpolators_;
    std::map< std::string, ExtendedTimeDifferenceInterpolator > globalCoordinateToProperTimeInterpolatorsExtended_;
    std::map< std::string, ExtendedTimeDifferenceInterpolator > properTimeToGlobalCoordinateInterpolatorsExtended_;
};

}  // namespace tudat

#endif  // TUDAT_TIME_EPHEMERIS_DIRECT_FROM_METRIC_H
