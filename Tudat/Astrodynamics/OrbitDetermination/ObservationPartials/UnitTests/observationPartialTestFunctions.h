#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationRate.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/numericalObservationPartial.h"
#include "Tudat/SimulationSetup/createBodyDeformationModel.h"
#include "Tudat/SimulationSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/defaultBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::site_displacements;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

std::vector< boost::shared_ptr< simulation_setup::BodyDeformationSettings > > getSunOnlyEarthDeformationModel( );

NamedBodyMap setupEnvironment( const std::vector< LinkEndId > groundStations,
                               const double initialEphemerisTime = 1.0E7,
                               const double finalEphemerisTime = 1.2E7,
                               const double stateEvaluationTime = 0.0 );

boost::shared_ptr< EstimatableParameterSet< double > > createEstimatableParameters(
        const NamedBodyMap& bodyMap, const double initialTime );

Eigen::Matrix< double, 1, 3 > calculatePartialWrtConstantBodyState(
        const std::string& bodyName, const NamedBodyMap& bodyMap, const Eigen::Vector3d& bodyPositionVariation,
        const boost::function< double( const double ) > observationFunction, const double observationTime );

std::vector< double > calculateNumericalPartialsWrtDoubleParameters(
        const std::vector< boost::shared_ptr< EstimatableParameter< double > > >& doubleParameters,
        const std::vector< boost::function< void( ) > > updateFunctions,
        const std::vector< double >& parameterPerturbations,
        const boost::function< double( const double ) > observationFunction,
        const double observationTime );

std::vector< std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > > calculateAnalyticalPartials(
        const std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > >& partialObjectList,
        const std::vector< basic_mathematics::Vector6d >& states,
        const std::vector< double >& times,
        const LinkEndType linkEndOfFixedTime );

std::vector< int > getSingleAnalyticalPartialSize(
        const LinkEnds& linkEnds,
        const std::pair< std::vector< std::string >, boost::shared_ptr< EstimatableParameterSet< double > > >& estimatedParameters );

std::vector< std::vector< double > > getAnalyticalPartialEvaluationTimes(
        const LinkEnds& linkEnds,
        const ObservableType observableType,
        const std::vector< double >& linkEndTimes,
        const boost::shared_ptr< EstimatableParameterSet< double > >& estimatedParameters );

inline void testObservationPartials(
        const boost::shared_ptr< ObservationModel< 1, double, double, double > > observationModel, NamedBodyMap& bodyMap,
        const boost::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet,
        const LinkEnds& linkEnds, const ObservableType observableType )
{
    std::vector< boost::shared_ptr< EstimatableParameter< double > > > doubleParameterVector =
            fullEstimatableParameterSet->getEstimatedDoubleParameters( );

    // Create observation partials.
    boost::shared_ptr< ObservationPartialCreator< 1, double > > observationPartialCreator;
    std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > >,
            boost::shared_ptr< PositionPartialScaling > > fullAnalyticalPartialSet =
            observationPartialCreator->createObservationPartials(
                observableType, boost::assign::list_of( linkEnds ), bodyMap, fullEstimatableParameterSet ).begin( )->second;
    boost::shared_ptr< PositionPartialScaling > positionPartialScaler = fullAnalyticalPartialSet.second;

    std::vector< std::string > bodiesWithEstimatedState = estimatable_parameters::getListOfBodiesToEstimate(
                fullEstimatableParameterSet );

    for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( );
         linkEndIterator++ )
    {
        std::cout<<"Current link end "<<linkEndIterator->first<<" "<<linkEnds.size( )<<std::endl;
        // Evaluate nominal observation values
        std::vector< basic_mathematics::Vector6d > vectorOfStates;
        std::vector< double > vectorOfTimes;
        double observationTime = 1.1E7;
        observationModel->computeObservationsWithLinkEndData(
                    observationTime, linkEndIterator->first, vectorOfTimes, vectorOfStates );

        // Calculate analytical observation partials.
        positionPartialScaler->update( vectorOfStates, vectorOfTimes, static_cast< LinkEndType >( linkEndIterator->first ) );
        typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > ObservationPartialReturnType;
        std::vector< ObservationPartialReturnType > analyticalObservationPartials =
                calculateAnalyticalPartials( fullAnalyticalPartialSet.first, vectorOfStates, vectorOfTimes, linkEndIterator->first );

        // Set and test expected partial size and time
        std::vector< std::vector< double > > expectedPartialTimes = getAnalyticalPartialEvaluationTimes(
                    linkEnds, observableType, vectorOfTimes, fullEstimatableParameterSet );
        for( unsigned int i = 0; i < analyticalObservationPartials.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( analyticalObservationPartials[ i ].size( ), expectedPartialTimes[ i ].size( ) );

            for( unsigned int j = 0; j < expectedPartialTimes[ i ].size( ); j++ )
            {
                BOOST_CHECK_EQUAL( analyticalObservationPartials[ i ][ j ].second, expectedPartialTimes[ i ][ j ] );
            }

        }

        // Settings for body state partials
        boost::function< double( const double ) > observationFunction = boost::bind(
                    &ObservationModel< 1, double, double, double >::computeObservationEntry, observationModel, _1, linkEndIterator->first, 1 );
        Eigen::Vector3d bodyPositionVariation;
        bodyPositionVariation << 1000.0E3, 1000.0E3, 1000.0E3;

        // Calculate numerical partials w.r.t. Earth state.
        Eigen::Matrix< double, 1, 3 > bodyPositionPartial = Eigen::Matrix< double, 1, 3 >::Zero( );
        for( unsigned int i = 0; i < bodiesWithEstimatedState.size( ); i++ )
        {
            Eigen::Matrix< double, 1, 3 > numericalPartialWrtBodyPosition = calculatePartialWrtConstantBodyState(
                        bodiesWithEstimatedState[ i ], bodyMap, bodyPositionVariation, observationFunction, observationTime );
            bodyPositionPartial.setZero( );
            for( unsigned int j = 0; j < analyticalObservationPartials[ i ].size( ); j++ )
            {
                //std::cout<<"An part "<<std::setprecision( 16 )<< analyticalObservationPartials[ i ][ j ].first<<" "<<j<<std::endl;
               bodyPositionPartial +=  analyticalObservationPartials[ i ][ j ].first;
            }
            //std::cout<<"Num part "<<std::setprecision( 16 )<< numericalPartialWrtBodyPosition<<std::endl<<std::endl;
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( bodyPositionPartial, ( -1.0 * numericalPartialWrtBodyPosition ), 1.0E-5 );
        }

        // Settings for parameter partial functions.
        std::vector< double > parameterPerturbations = boost::assign::list_of( 1.0E-8 )( 1.0E-8 );
        std::vector< boost::function< void( ) > > updateFunctionList;
        updateFunctionList.push_back( emptyFunction2 );
        updateFunctionList.push_back( emptyFunction2 );

        // Calculate and test analytical against numerical partials.
        std::vector< double > numericalPartialsWrtDoubleParameters = calculateNumericalPartialsWrtDoubleParameters(
                    doubleParameterVector, updateFunctionList, parameterPerturbations, observationFunction, observationTime );

        double currentParameterPartial = 0.0;
        int numberOfEstimatedBodies = bodiesWithEstimatedState.size( );
        for( unsigned int i = 0; i < numericalPartialsWrtDoubleParameters.size( ); i++ )
        {
            std::cout<<"Double parameter "<<i<<std::endl;
            currentParameterPartial = 0.0;
            for( unsigned int j = 0; j < analyticalObservationPartials[ i + numberOfEstimatedBodies ].size( ); j++ )
            {
                currentParameterPartial += analyticalObservationPartials[ i + numberOfEstimatedBodies ][ j ].first.x( );

            }
            BOOST_CHECK_CLOSE_FRACTION( currentParameterPartial, -1.0 * numericalPartialsWrtDoubleParameters[ i ], 5.0E-3 );
        }
    }

}

}

}
