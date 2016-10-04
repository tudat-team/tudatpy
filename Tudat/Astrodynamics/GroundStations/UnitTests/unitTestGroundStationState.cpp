#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/GroundStations/nominalGroundStationState.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::coordinate_conversions;
using namespace tudat::basic_astrodynamics;
using namespace tudat::ground_stations;


BOOST_AUTO_TEST_SUITE( test_ground_station_state )

BOOST_AUTO_TEST_CASE( test_GroundStationState )
{
    {
        double bodyEquatorialRadius = 6378.0;
        double bodyFlattening = 1.0 / 298.257;
        double bodyPolarRadius = bodyEquatorialRadius * ( 1.0 - bodyFlattening );

        boost::shared_ptr< BodyShapeModel > oblateSpheroidShapeModel = boost::make_shared< OblateSpheroidBodyShapeModel >(
                    bodyEquatorialRadius, bodyFlattening );


        double stationAltitude = 0.0;
        double stationLongitude = unit_conversions::convertDegreesToRadians( 58.0 );
        double stationGeodeticLatitude = unit_conversions::convertDegreesToRadians( 32.0 );

        Eigen::Vector3d stationGeodeticCoordinates;
        stationGeodeticCoordinates<<stationAltitude, stationGeodeticLatitude, stationLongitude;

        Eigen::Vector3d stationCartesianPosition = convertCartesianToGeodeticCoordinates(
                    stationGeodeticCoordinates, bodyEquatorialRadius, bodyFlattening, 1.0E-3 );

        boost::shared_ptr< NominalGroundStationState > stationStateFromOblateSpheroid = boost::make_shared< NominalGroundStationState >(
                    stationCartesianPosition, oblateSpheroidShapeModel );
        //boost::shared_ptr< GeodeticNominalGroundStationState > geodeticStationStateFromOblateSpheroid = boost::make_shared< GeodeticNominalGroundStationState >(
        //            stationCartesianPosition, oblateSpheroidShapeModel, Eigen::Vector3d::Zero( ) );

        Eigen::Matrix3d toTopocentricFramFromOblateSpheroid =
                Eigen::Matrix3d( stationStateFromOblateSpheroid->getRotationFromBodyFixedToTopocentricFrame( 0.0 ) );
        //Eigen::Matrix3d toTopocentricFramFromEllipsoid =
        //        Eigen::Matrix3d( stationStateFromEllipsoid->getRotationFromBodyFixedToTopocentricFrame( 0.0 ) );
        //Eigen::Matrix3d toTopocentricFramFromGeodeticState =
        //        Eigen::Matrix3d( geodeticStationStateFromOblateSpheroid->getRotationFromBodyFixedToTopocentricFrame( 0.0 ) );

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
               // BOOST_CHECK_SMALL( toTopocentricFramFromOblateSpheroid( i, j ) - toTopocentricFramFromEllipsoid( i, j ), 1.0E-11 );
               // BOOST_CHECK_SMALL( toTopocentricFramFromOblateSpheroid( i, j ) - toTopocentricFramFromGeodeticState( i, j ), 1.0E-11);
            }
        }

        //BOOST_CHECK_CLOSE_FRACTION( stationStateFromOblateSpheroid->getLongitude( ), geodeticStationStateFromOblateSpheroid->getLongitude( ),
        //                            std::numeric_limits< double >::epsilon( ) );
        //BOOST_CHECK_CLOSE_FRACTION( stationStateFromOblateSpheroid->getLongitude( ), stationStateFromEllipsoid->getLongitude( ),
        //                            std::numeric_limits< double >::epsilon( ) );

        //BOOST_CHECK_CLOSE_FRACTION( stationGeodeticLatitude, geodeticStationStateFromOblateSpheroid->getGeodeticLatitude( ), 1.0E-11 );
        //BOOST_CHECK_CLOSE_FRACTION( stationAltitude, geodeticStationStateFromOblateSpheroid->getNominalElevation( ), 1.0E-11 );
    }

    {
        std::string baseFolder = input_output::getSpiceKernelPath( ) + "SpiceKernels/";
        spice_interface::loadSpiceKernelInTudat( baseFolder + "topocentricFramePck.pck" );
        spice_interface::loadSpiceKernelInTudat( baseFolder + "topocentricFrameSpk.spk" );
        spice_interface::loadSpiceKernelInTudat( baseFolder + "pck00009.tpc" );

        Eigen::Vector3d shanghaiStationPosition;
        shanghaiStationPosition<< -2847.6980296, +4659.8725718, +3283.9585330;
        shanghaiStationPosition = shanghaiStationPosition * 1000.0;

        // Call spice to retrieve the radii of Earth.
        double radii[ 3 ];
        SpiceInt numberOfReturnedParameters;
        std::string earthName = "Earth";
        bodvrd_c( earthName.c_str( ), "RADII", 3, &numberOfReturnedParameters, radii );

        double bodyEquatorialRadius = radii[ 0 ] * 1.0E3;
        double bodyFlattening = ( radii[ 0 ] - radii[ 2 ] ) / radii[ 0 ];

        boost::shared_ptr< BodyShapeModel > spiceCompatibleShapeModel = boost::make_shared< OblateSpheroidBodyShapeModel >(
                    bodyEquatorialRadius, bodyFlattening );

        boost::shared_ptr< NominalGroundStationState > shanghaiStationState = boost::make_shared< NominalGroundStationState >(
                    shanghaiStationPosition, spiceCompatibleShapeModel );

        Eigen::Matrix3d calculatedRotationMatrix = Eigen::Matrix3d( shanghaiStationState->getRotationFromBodyFixedToTopocentricFrame( 0.0 ) );
        Eigen::Matrix3d spiceRotationMatrix = Eigen::Matrix3d(
                    spice_interface::computeRotationQuaternionBetweenFrames( "IAU_Earth", "SHANGHAI_TOPO", 1.0 ) );
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( calculatedRotationMatrix( i, j ) - spiceRotationMatrix( i, j ), 2.0E-14 );
            }
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
