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
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationalOrientation.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/createObservationPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/numericalObservationPartial.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

//std::vector< boost::shared_ptr< simulation_setup::BodyDeformationSettings > > getSunOnlyEarthDeformationModel( );

NamedBodyMap setupEnvironment( const std::vector< LinkEndId > groundStations,
                               const double initialEphemerisTime = 1.0E7,
                               const double finalEphemerisTime = 1.2E7,
                               const double stateEvaluationTime = 0.0,
                               const bool useConstantEphemerides = 1 );

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

std::vector< Eigen::MatrixXd > calculateNumericalPartialsWrtVectorParameters(
        const std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >& vectorParameters,
        const std::vector< boost::function< void( ) > > updateFunctions,
        const std::vector< Eigen::VectorXd >& parameterPerturbations,
        const boost::function< Eigen::VectorXd( const double ) > observationFunction,
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
        const LinkEnds& linkEnds, const ObservableType observableType,
        const bool testPositionPartial = 1,
        const bool testParameterPartial = 1 )
{
    std::vector< boost::shared_ptr< EstimatableParameter< double > > > doubleParameterVector =
            fullEstimatableParameterSet->getEstimatedDoubleParameters( );
    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameterVector =
            fullEstimatableParameterSet->getEstimatedVectorParameters( );

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

        BOOST_CHECK_EQUAL( analyticalObservationPartials.size( ), expectedPartialTimes.size( ) );

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
                    &ObservationModel< 1, double, double, double >::computeObservationEntry, observationModel, _1, linkEndIterator->first, 0 );
        Eigen::Vector3d bodyPositionVariation;

        if( testPositionPartial )
        {
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
                    bodyPositionPartial +=  analyticalObservationPartials[ i ][ j ].first;
                }
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( bodyPositionPartial, ( -1.0 * numericalPartialWrtBodyPosition ), 1.0E-6 );
            }
        }

        if( testParameterPartial )
        {

            {
                // Settings for parameter partial functions.
                std::vector< double > parameterPerturbations = boost::assign::list_of( 1.0E-10 )( 1.0E-10 );
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
                    currentParameterPartial = 0.0;
                    for( unsigned int j = 0; j < analyticalObservationPartials[ i + numberOfEstimatedBodies ].size( ); j++ )
                    {
                        currentParameterPartial += analyticalObservationPartials[ i + numberOfEstimatedBodies ][ j ].first.x( );

                    }

                    BOOST_CHECK_CLOSE_FRACTION( currentParameterPartial, -1.0 * numericalPartialsWrtDoubleParameters[ i ], 1.0E-6 );
                }
            }

            {
                boost::function< Eigen::VectorXd( const double ) > vectorObservationFunction = boost::bind(
                            &ObservationModel< 1, double, double, double >::computeObservations, observationModel, _1, linkEndIterator->first );

                // Settings for parameter partial functions.
                std::vector< Eigen::VectorXd > parameterPerturbations;
                parameterPerturbations.push_back( Eigen::Vector2d::Constant( 1.0E-4 ) );
                parameterPerturbations.push_back( Eigen::Vector2d::Constant( 1.0E-4 ) );

                std::vector< boost::function< void( ) > > updateFunctionList;
                updateFunctionList.push_back( emptyFunction2 );
                updateFunctionList.push_back( emptyFunction2 );

                // Calculate and test analytical against numerical partials.
                std::vector< Eigen::MatrixXd > numericalPartialsWrtVectorParameters = calculateNumericalPartialsWrtVectorParameters(
                            vectorParameterVector, updateFunctionList, parameterPerturbations, vectorObservationFunction, observationTime );

                Eigen::VectorXd currentParameterPartial;
                int startIndex = bodiesWithEstimatedState.size( ) + doubleParameterVector.size( );
                for( unsigned int i = 0; i < numericalPartialsWrtVectorParameters.size( ); i++ )
                {
                    currentParameterPartial = Eigen::VectorXd::Zero( vectorParameterVector.at( i )->getParameterSize( ) );

                    for( unsigned int j = 0; j < analyticalObservationPartials[ i + startIndex ].size( ); j++ )
                    {
                        currentParameterPartial += analyticalObservationPartials[ i + startIndex ][ j ].first;

                    }

                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( currentParameterPartial.transpose( ) ), ( -1.0 * numericalPartialsWrtVectorParameters[ i ] ), 1.0E-6 );
                }
            }
        }
    }

}

}

}
