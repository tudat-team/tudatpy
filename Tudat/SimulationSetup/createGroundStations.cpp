#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

#include "Tudat/Astrodynamics/GroundStations/basicTidalBodyDeformation.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/SimulationSetup/createGroundStations.h"

namespace tudat
{

namespace simulation_setup
{

void createGroundStations( const std::map< std::string, boost::shared_ptr< Body > >& bodyMap,
                           const std::map< std::pair< std::string, std::string >, Eigen::Vector3d >& groundStationsWithPosition )
{
    using namespace tudat::ephemerides;
    using namespace tudat::coordinate_conversions;

    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > currentBodyShapeModel;

    for( std::map< std::pair< std::string, std::string >, Eigen::Vector3d >::const_iterator stationIterator = groundStationsWithPosition.begin( );
         stationIterator != groundStationsWithPosition.end( ); stationIterator++ )
    {
        if( bodyMap.count( stationIterator->first.first ) > 0 )
        {
            currentBodyShapeModel =  bodyMap.at( stationIterator->first.first )->getShapeModel( );

            if( currentBodyShapeModel == NULL )
            {
                std::cerr<<"Error when making ground station "<<stationIterator->first.second<<" of "<<stationIterator->first.first<<", body "<<
                           " has no shape model."<<std::endl;
            }
            else
            {

                boost::shared_ptr< ground_stations::NominalGroundStationState > stationState;

                    stationState = boost::make_shared< ground_stations::NominalGroundStationState >(
                                stationIterator->second, currentBodyShapeModel  );

                boost::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
                        boost::make_shared< ground_stations::PointingAnglesCalculator >(
                            boost::bind( &RotationalEphemeris::getRotationToTargetFrame,
                                         bodyMap.at( stationIterator->first.first )->getRotationalEphemeris( ),
                                         _1, basic_astrodynamics::JULIAN_DAY_ON_J2000 ),
                            boost::bind( &ground_stations::NominalGroundStationState::getRotationFromBodyFixedToTopocentricFrame, stationState, _1 ) );

                bodyMap.at( stationIterator->first.first )->addGroundStation(
                            stationIterator->first.second, boost::make_shared< ground_stations::GroundStation >(
                                stationState, pointingAnglesCalculator, stationIterator->first.second ) );

            }
        }
    }
}


void setSingleBodyGroundStationPositionVariationFunctions( boost::shared_ptr< Body > body,
                                                           gravitation::BodyDeformationTypes variationType,
                                                           const double initialTime,
                                                           const double finalTime,
                                                           const double timeStep )
{
    std::map< std::string, boost::shared_ptr< ground_stations::GroundStation > > groundStationMap =
            body->getGroundStationMap( );

    for( std::map< std::string, boost::shared_ptr< ground_stations::GroundStation > >::iterator groundStationIterator =
         groundStationMap.begin( ); groundStationIterator != groundStationMap.end( ); groundStationIterator++ )
    {
        setGroundStationPositionVariationFunction( groundStationIterator->second,
                                                   body, variationType, initialTime, finalTime, timeStep );
    }
}

void setSingleBodyGroundStationPositionVariationFunctions( boost::shared_ptr< Body > body,
                                                           std::vector< gravitation::BodyDeformationTypes > variationTypes,
                                                           const double initialTime,
                                                           const double finalTime,
                                                           const double timeStep )
{
    std::map< std::string, boost::shared_ptr< ground_stations::GroundStation > > groundStationMap =
            body->getGroundStationMap( );

    for( std::map< std::string, boost::shared_ptr< ground_stations::GroundStation > >::iterator groundStationIterator =
         groundStationMap.begin( ); groundStationIterator != groundStationMap.end( ); groundStationIterator++ )
    {
        setGroundStationPositionVariationFunction( groundStationIterator->second,
                                                   body, variationTypes, initialTime, finalTime, timeStep );
    }
}

boost::function< Eigen::Vector3d( const double ) > getGroundStationPositionVariationFunction(
        boost::shared_ptr< ground_stations::GroundStation > groundStation,
        boost::shared_ptr< Body > bodyWithStation,
        gravitation::BodyDeformationTypes variationType,
        const double initialTime,
        const double finalTime,
        const double timeStep )
{
    using namespace tudat::ground_stations;
    using namespace tudat::site_displacements;

    boost::function< Eigen::Vector3d( const double ) >  variationFunction;

    switch( variationType )
    {

//    case basic_solid_body:
//    {
//        boost::shared_ptr< BasicTidalBodyDeformation > deformationModel =
//                bodyWithStation->getBasicTidalBodyDeformation( );

//        boost::shared_ptr< ground_stations::NominalGroundStationState > nominalSiteState =
//                groundStation->getNominalStationState( );

//        boost::function< Eigen::Vector3d( const double ) > displacementFunction =
//                boost::bind( &BasicTidalBodyDeformation::calculateSiteDisplacement,
//                             deformationModel, _1, nominalSiteState,  );

//        variationFunction = displacementFunction;

//        break;
//    }

    default:
    {
        std::cerr<<"Error, ground station position variation type not yet implemented."<<std::endl;
    }
    }
    return variationFunction;
}


void setGroundStationPositionVariationFunction( boost::shared_ptr< ground_stations::GroundStation > groundStation,
                                                boost::shared_ptr< Body > bodyWithStation,
                                                std::vector< gravitation::BodyDeformationTypes > variationTypes,
                                                const double initialTime,
                                                const double finalTime,
                                                const double timeStep )

{
    std::vector< boost::function< Eigen::Vector3d( const double ) > > displacementFunctions;

    for( unsigned int i = 0; i < variationTypes.size( ); i++ )
    {
        displacementFunctions.push_back( getGroundStationPositionVariationFunction( groundStation, bodyWithStation, variationTypes[ i ],
                                                                                    initialTime, finalTime, timeStep ) );

    }
    groundStation->getNominalStationState( )->setPositionVariationsUpdateFunctions(
                boost::make_shared< ground_stations::DirectGroundStationPositionVariationSettings >(
                    initialTime, finalTime, timeStep, displacementFunctions ) );

}

void setGroundStationPositionVariationFunction( boost::shared_ptr< ground_stations::GroundStation > groundStation,
                                                boost::shared_ptr< Body > bodyWithStation,
                                                gravitation::BodyDeformationTypes variationType,
                                                const double initialTime,
                                                const double finalTime,
                                                const double timeStep )
{
    std::vector< gravitation::BodyDeformationTypes > variationTypes;
    variationTypes.push_back( variationType );
    setGroundStationPositionVariationFunction( groundStation, bodyWithStation, variationTypes, initialTime, finalTime, timeStep );
}

}

}
