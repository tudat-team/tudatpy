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
#include "Tudat/SimulationSetup/EstimationSetup/createObservationPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/numericalObservationPartial.h"
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

Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyState(
        const std::string& bodyName, const NamedBodyMap& bodyMap, const Eigen::Vector3d& bodyPositionVariation,
        const boost::function< Eigen::VectorXd( const double ) > observationFunction, const double observationTime, const int observableSize );

std::vector< Eigen::VectorXd > calculateNumericalPartialsWrtDoubleParameters(
        const std::vector< boost::shared_ptr< EstimatableParameter< double > > >& doubleParameters,
        const std::vector< boost::function< void( ) > > updateFunctions,
        const std::vector< double >& parameterPerturbations,
        const boost::function< Eigen::VectorXd( const double ) > observationFunction,
        const double observationTime );

std::vector< Eigen::MatrixXd > calculateNumericalPartialsWrtVectorParameters(
        const std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >& vectorParameters,
        const std::vector< boost::function< void( ) > > updateFunctions,
        const std::vector< Eigen::VectorXd >& parameterPerturbations,
        const boost::function< Eigen::VectorXd( const double ) > observationFunction,
        const double observationTime );

template< int ObservableSize >
std::vector< std::vector< std::pair< Eigen::Matrix< double, ObservableSize, Eigen::Dynamic >, double > > > calculateAnalyticalPartials(
        const std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< ObservableSize > > >& partialObjectList,
        const std::vector< basic_mathematics::Vector6d >& states,
        const std::vector< double >& times,
        const LinkEndType linkEndOfFixedTime )
{
    std::vector< std::vector< std::pair< Eigen::Matrix< double, ObservableSize, Eigen::Dynamic >, double > > > partialList;

    for( typename std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< ObservableSize > > >::const_iterator partialIterator =
         partialObjectList.begin( ); partialIterator != partialObjectList.end( ); partialIterator++ )
    {
        partialList.push_back( partialIterator->second->calculatePartial( states, times, linkEndOfFixedTime ) );
    }
    return partialList;
}

std::vector< int > getSingleAnalyticalPartialSize(
        const LinkEnds& linkEnds,
        const std::pair< std::vector< std::string >, boost::shared_ptr< EstimatableParameterSet< double > > >& estimatedParameters );

std::vector< std::vector< double > > getAnalyticalPartialEvaluationTimes(
        const LinkEnds& linkEnds,
        const ObservableType observableType,
        const std::vector< double >& linkEndTimes,
        const boost::shared_ptr< EstimatableParameterSet< double > >& estimatedParameters );

template< int ObservableSize = 1 >
inline void testObservationPartials(
        const boost::shared_ptr< ObservationModel< ObservableSize, double, double, double > > observationModel, NamedBodyMap& bodyMap,
        const boost::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet,
        const LinkEnds& linkEnds, const ObservableType observableType,
        const double tolerance = 1.0E-6,
        const bool testPositionPartial = 1,
        const bool testParameterPartial = 1 )
{
    std::vector< boost::shared_ptr< EstimatableParameter< double > > > doubleParameterVector =
            fullEstimatableParameterSet->getEstimatedDoubleParameters( );
    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameterVector =
            fullEstimatableParameterSet->getEstimatedVectorParameters( );

    // Create observation partials.
    boost::shared_ptr< ObservationPartialCreator< ObservableSize, double > > observationPartialCreator;

    std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< ObservableSize > > >,
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
        typedef std::vector< std::pair< Eigen::Matrix< double, ObservableSize, Eigen::Dynamic >, double > > ObservationPartialReturnType;
        std::vector< ObservationPartialReturnType > analyticalObservationPartials =
                calculateAnalyticalPartials< ObservableSize >( fullAnalyticalPartialSet.first, vectorOfStates, vectorOfTimes, linkEndIterator->first );

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
        boost::function< Eigen::VectorXd( const double ) > observationFunction = boost::bind(
                    &ObservationModel< ObservableSize, double, double, double >::computeObservations, observationModel, _1, linkEndIterator->first );
        Eigen::Vector3d bodyPositionVariation;

        if( testPositionPartial )
        {
            bodyPositionVariation << 1000.0E3, 1000.0E3, 1000.0E3;

            // Calculate numerical partials w.r.t. Earth state.
            Eigen::Matrix< double, Eigen::Dynamic, 3 > bodyPositionPartial = Eigen::Matrix< double, ObservableSize, 3 >::Zero( );
            for( unsigned int i = 0; i < bodiesWithEstimatedState.size( ); i++ )
            {
                Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalPartialWrtBodyPosition = calculatePartialWrtConstantBodyState(
                            bodiesWithEstimatedState[ i ], bodyMap, bodyPositionVariation, observationFunction, observationTime, ObservableSize );
                bodyPositionPartial.setZero( );
                for( unsigned int j = 0; j < analyticalObservationPartials[ i ].size( ); j++ )
                {
                    bodyPositionPartial +=  analyticalObservationPartials[ i ][ j ].first;
                }

                if( observableType != angular_position )
                {
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( bodyPositionPartial, ( numericalPartialWrtBodyPosition ), tolerance );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( bodyPositionPartial( 0, 2 ) - numericalPartialWrtBodyPosition( 0, 2 ) ), 1.0E-20 );
                    bodyPositionPartial( 0, 2 ) = 0.0;
                    numericalPartialWrtBodyPosition( 0, 2 ) = 0.0;

                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( bodyPositionPartial, ( numericalPartialWrtBodyPosition ), tolerance );

                }
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
                std::vector< Eigen::VectorXd > numericalPartialsWrtDoubleParameters = calculateNumericalPartialsWrtDoubleParameters(
                            doubleParameterVector, updateFunctionList, parameterPerturbations, observationFunction, observationTime );

                Eigen::VectorXd currentParameterPartial;
                int numberOfEstimatedBodies = bodiesWithEstimatedState.size( );
                for( unsigned int i = 0; i < numericalPartialsWrtDoubleParameters.size( ); i++ )
                {
                    currentParameterPartial.setZero( ObservableSize );
                    for( unsigned int j = 0; j < analyticalObservationPartials[ i + numberOfEstimatedBodies ].size( ); j++ )
                    {
                        currentParameterPartial += analyticalObservationPartials[ i + numberOfEstimatedBodies ][ j ].first;

                    }

                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentParameterPartial, ( numericalPartialsWrtDoubleParameters[ i ] ), tolerance );
                }
            }

            {
                boost::function< Eigen::VectorXd( const double ) > vectorObservationFunction = boost::bind(
                            &ObservationModel< ObservableSize, double, double, double >::computeObservations, observationModel, _1, linkEndIterator->first );

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

                Eigen::MatrixXd currentParameterPartial;
                int startIndex = bodiesWithEstimatedState.size( ) + doubleParameterVector.size( );
                for( unsigned int i = 0; i < numericalPartialsWrtVectorParameters.size( ); i++ )
                {
                    currentParameterPartial = Eigen::MatrixXd::Zero( ObservableSize, vectorParameterVector.at( i )->getParameterSize( ) );

                    for( unsigned int j = 0; j < analyticalObservationPartials[ i + startIndex ].size( ); j++ )
                    {
                        currentParameterPartial += analyticalObservationPartials[ i + startIndex ][ j ].first;

                    }

                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( currentParameterPartial ), ( numericalPartialsWrtVectorParameters[ i ] ), tolerance );
                }
            }
        }
    }

}

}

}
