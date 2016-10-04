#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCalculator.h"
#include "Tudat/Astrodynamics/GroundStations/groundStation.h"
#include "Tudat/Astrodynamics/Ephemerides/compositeEphemeris.h"

namespace tudat
{

namespace observation_models
{

template< >
boost::shared_ptr< ephemerides::Ephemeris > createReferencePointEphemeris< double, double >(
        boost::shared_ptr< simulation_setup::Body > bodyWithReferencePoint,
        boost::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel,
        boost::function< basic_mathematics::Vector6d( const double& ) > referencePointRelativeStateFunction )
{
    typedef Eigen::Matrix< double, 6, 1 > StateType;

    std::map< int, boost::function< basic_mathematics::Vector6d( const double& ) > > stationEphemerisVector;
    stationEphemerisVector[ 2 ] = boost::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris, bodyWithReferencePoint, _1 );
    stationEphemerisVector[ 0 ] = referencePointRelativeStateFunction;

    std::map< int, boost::function< StateType( const double, const StateType& ) > > stationRotationVector;
    stationRotationVector[ 1 ] = boost::bind( &ephemerides::transformStateToGlobalFrame< double, double >, _2, _1, bodyRotationModel );

    boost::shared_ptr< ephemerides::Ephemeris > ephemeris = boost::make_shared< ephemerides::CompositeEphemeris< double, double > >(
                stationEphemerisVector, stationRotationVector,  "SSB", "ECLIPJ2000" );

    return ephemeris;
}

template< >
boost::shared_ptr< ephemerides::Ephemeris > createReferencePointEphemeris< double, long double >(
        boost::shared_ptr< simulation_setup::Body > bodyWithReferencePoint,
        boost::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel,
        boost::function< basic_mathematics::Vector6d( const double& ) > referencePointRelativeStateFunction )
{
    typedef Eigen::Matrix< long double, 6, 1 > StateType;

    std::map< int, boost::function< Eigen::Matrix< long double, 6, 1 >( const double& ) > > stationEphemerisVector;
    stationEphemerisVector[ 2 ] = boost::bind( &simulation_setup::Body::getLongStateInBaseFrameFromEphemeris, bodyWithReferencePoint, _1 );
    stationEphemerisVector[ 0 ] = boost::bind( &convertLongDoubleStateFromDoubleStateFunction< double >,
                                               _1, referencePointRelativeStateFunction );

    std::map< int, boost::function< StateType( const double, const StateType& ) > > stationRotationVector;
    stationRotationVector[ 1 ] =  boost::bind( &ephemerides::transformStateToGlobalFrame< long double, double >, _2, _1, bodyRotationModel );

    boost::shared_ptr< ephemerides::Ephemeris > ephemeris = boost::make_shared< ephemerides::CompositeEphemeris< double, long double > >(
                stationEphemerisVector, stationRotationVector, "SSB", "ECLIPJ2000" );

    return ephemeris;
}


}

}

