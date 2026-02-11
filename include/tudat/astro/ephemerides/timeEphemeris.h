#ifndef TIMEEPHEMERIS_H
#define TIMEEPHEMERIS_H

#include <map>

#include <boost/function.hpp>

#include "tudat/math/interpolators/createInterpolator.h"


#include "tudat/astro/basic_astro/timeConversions.h"

#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/math/interpolators/createInterpolator.h"

#include "tudat/astro/relativity/relativisticTimeConversion.h"
#include "tudat/astro/relativity/relativisticPotentials.h"


namespace tudat
{

double combineTimeDifferenceFunction( const std::vector< std::function< double( const double ) > > timeDifferenceFunctions,
                                      const std::vector< int > evaluationStepIndices,
                                      const double time );

using namespace tudat::basic_astrodynamics;

class TimeEphemeris
{
public:

    TimeEphemeris( const std::string& centralBodyName ):centralBodyName_( centralBodyName )
    { }


    virtual ~TimeEphemeris( ){ }

    //! Function to retrieve the time difference at a given time between two scales
    double getTimeDifference( const TimeScales inputScale, const TimeScales outputScale,
                              const double inputTime, const std::string pointIdentifier = ""  )
    {
        return getTimeDifferenceFunction( inputScale, outputScale, pointIdentifier )( inputTime );
    }

    std::map< double, double > getTimeDifferences( const TimeScales inputScale, const TimeScales outputScale,
                                                   const std::vector< double >& inputTimes, const std::string pointIdentifier = ""  )
    {
        std::function< double( const double ) > timeDifferenceFunction =
                getTimeDifferenceFunction( inputScale, outputScale, pointIdentifier );

        std::map< double, double > timeDifferences;
        for( unsigned int i = 0; i < inputTimes.size( ); i++ )
        {
            timeDifferences[ inputTimes[ i ] ] = timeDifferenceFunction( inputTimes[ i ] );
        }
        return timeDifferences;
    }

    virtual std::function< double( const double ) > getTimeDifferenceFunction(
            const TimeScales inputScale, const TimeScales outputScale, const std::string pointIdentifier = "" ) = 0;
protected:

    std::string centralBodyName_;
};

//! Class that handles all relativistic time conversions for a single body (including all the ground station proper times)
class TimeEphemerisFromPostNewtonianExpansion: public TimeEphemeris
{
public:
    typedef std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > TimeDifferenceInterpolator;
    TimeEphemerisFromPostNewtonianExpansion(
            const TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
            const TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator,
            const std::string& centralBodyName,
            const std::map< std::string, std::function< Eigen::Vector3d( const double ) > >& groundStationPositionFunctions =
            ( std::map< std::string, std::function< Eigen::Vector3d( const double ) > >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ) ): TimeEphemeris( centralBodyName ),
        barycenterToPlanetCenterCoordinateTimeInterpolator_( barycenterToPlanetCenterCoordinateTimeInterpolator ),
        planetCenterToBarycenterCoordinateTimeInterpolator_( planetCenterToBarycenterCoordinateTimeInterpolator ),
        groundStationPositionFunctions_( groundStationPositionFunctions ),
        planetCoordinateToProperTimeInterpolators_( planetCoordinateToProperTimeInterpolators ),
        properTimeToPlanetCoordinateInterpolators_( properTimeToPlanetCoordinateInterpolators )
    {

    }

    virtual ~TimeEphemerisFromPostNewtonianExpansion( ){ }

    void resetBarycentricToBodycentricInterpolators(
            const TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
            const TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator )
    {
        barycenterToPlanetCenterCoordinateTimeInterpolator_ = barycenterToPlanetCenterCoordinateTimeInterpolator;
        planetCenterToBarycenterCoordinateTimeInterpolator_ = planetCenterToBarycenterCoordinateTimeInterpolator;
    }

    void resetBodycentricToTopocentricInterpolators(
            const TimeDifferenceInterpolator planetCoordinateToProperTimeInterpolator,
            const TimeDifferenceInterpolator properTimeToPlanetCoordinateInterpolator,
            const std::string& referencePoint,
            const std::function< Eigen::Vector3d( const double ) > referencePointPositionFunction = NULL )
    {
        planetCoordinateToProperTimeInterpolators_[ referencePoint ] = planetCoordinateToProperTimeInterpolator;
        properTimeToPlanetCoordinateInterpolators_[ referencePoint ] = properTimeToPlanetCoordinateInterpolator;
        if( !doesReferencePointTopocentricConverterExist( referencePoint ) )
        {
            if( referencePointPositionFunction == NULL )
            {
                std::cerr<<"Error when resetting bodycentric to topocentic time converter for point "<<referencePoint<<
                           ", must also provide point position function; point not yet known in TimeEphemeris"<<std::endl;
            }
            else
            {
                groundStationPositionFunctions_[ referencePoint ] = referencePointPositionFunction;
            }
        }
    }

    std::function< double( const double ) > getTimeDifferenceFunction(
            const TimeScales inputScale, const TimeScales outputScale, const std::string pointIdentifier = "" );

    double calculateDirectTimeDifferenceTermFromFunction(
            const std::function< Eigen::Vector3d( const double ) > positionVectorFunctionFromReferencePoint,
            const double currentTime,
            const double multiplier = 1.0 )
    {
        return multiplier * calculateDirectTimeDifferenceTerm( positionVectorFunctionFromReferencePoint( currentTime ), currentTime );
    }

    virtual double calculateDirectTimeDifferenceTerm(
            const Eigen::Vector3d positionVectorFromReferencePoint, const double currentTime ) = 0;

    bool doesReferencePointTopocentricConverterExist( const std::string& referencePointName )
    {
        bool isFound = 1;
        if( ( groundStationPositionFunctions_.count( referencePointName ) == 0 ) ||
                ( planetCoordinateToProperTimeInterpolators_.count( referencePointName ) == 0 ) ||
                ( properTimeToPlanetCoordinateInterpolators_.count( referencePointName ) == 0 ) )
        {
            isFound = 0;
        }
        return isFound;
    }

protected:

    TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator_;

    TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator_;

    std::map< std::string, std::function< Eigen::Vector3d( const double ) > > groundStationPositionFunctions_;

    std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators_;

    std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators_;
};


class TimeEphemerisWithFirstOrderDirectConversion: public TimeEphemerisFromPostNewtonianExpansion
{
public:
    typedef std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > TimeDifferenceInterpolator;

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
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ) ):
        TimeEphemerisFromPostNewtonianExpansion(
            barycenterToPlanetCenterCoordinateTimeInterpolator, planetCenterToBarycenterCoordinateTimeInterpolator, centralBodyName,
            groundStationPositionFunctions, planetCoordinateToProperTimeInterpolators, properTimeToPlanetCoordinateInterpolators ),
        centralBodyStateFunction_( centralBodyStateFunction )
    {

    }

    double calculateDirectTimeDifferenceTerm( const Eigen::Vector3d positionVectorFromReferencePoint,
                                              const double currentTime )
    {
        return relativity::calculateFirstOrderTcbToTcgDirectCorrection(
                    positionVectorFromReferencePoint, centralBodyStateFunction_( currentTime ).segment( 3, 3 ) );
    }

private:

    std::function< Eigen::Vector6d( const double ) > centralBodyStateFunction_;
};

class TimeEphemerisWithSecondOrderDirectConversion: public TimeEphemerisFromPostNewtonianExpansion
{
public:
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
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ) ):
        TimeEphemerisFromPostNewtonianExpansion(
            barycenterToPlanetCenterCoordinateTimeInterpolator, planetCenterToBarycenterCoordinateTimeInterpolator, centralBodyName,
            groundStationPositionFunctions, planetCoordinateToProperTimeInterpolators, properTimeToPlanetCoordinateInterpolators ),
        centralBodyStateFunction_( centralBodyStateFunction )
    {

    }

    double calculateDirectTimeDifferenceTerm( const Eigen::Vector3d positionVectorFromReferencePoint,
                                              const double currentTime )
    {
        std::cerr<<"Error, second order time ephemeris direct difference term not yet implemented"<<std::endl;
        return 0.0;
    }
private:
    std::function< Eigen::Vector6d( const double ) > centralBodyStateFunction_;
};

class TimeEphemerisDirectFromMetric: public TimeEphemeris
{
public:
    typedef std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > TimeDifferenceInterpolator;
    TimeEphemerisDirectFromMetric(
            const std::string& centralBodyName,
            const std::map< std::string, TimeDifferenceInterpolator > globalCoordinateToProperTimeInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > properTimeToGlobalCoordinateInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ) ): TimeEphemeris( centralBodyName ),
        globalCoordinateToProperTimeInterpolators_( globalCoordinateToProperTimeInterpolators ),
        properTimeToGlobalCoordinateInterpolators_( properTimeToGlobalCoordinateInterpolators )
    {

    }

    void resetGlobalToProperTimeInterpolators(
            const TimeDifferenceInterpolator globalCoordinateToProperTimeInterpolator,
            const TimeDifferenceInterpolator properTimeToGlobalCoordinateInterpolator,
            const std::string& referencePoint )
    {
        globalCoordinateToProperTimeInterpolators_[ referencePoint ] = globalCoordinateToProperTimeInterpolator;
        properTimeToGlobalCoordinateInterpolators_[ referencePoint ] = properTimeToGlobalCoordinateInterpolator;
    }

    std::function< double( const double ) > getTimeDifferenceFunction(
            const TimeScales inputScale, const TimeScales outputScale, const std::string pointIdentifier = "" );


    double getGlobalCoordinateToProperTimeDifference( const std::string& pointId, const double evaluationTime )
    {
        return globalCoordinateToProperTimeInterpolators_.at( pointId )->interpolate( evaluationTime );
    }

    double getProperToGlobalCoordinateTimeDifference( const std::string& pointId, const double evaluationTime  )
    {
        return properTimeToGlobalCoordinateInterpolators_.at( pointId )->interpolate( evaluationTime );
    }

protected:

    std::map< std::string, TimeDifferenceInterpolator > globalCoordinateToProperTimeInterpolators_;

    std::map< std::string, TimeDifferenceInterpolator > properTimeToGlobalCoordinateInterpolators_;
};

}

#endif // TIMEEPHEMERIS_H
