/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationalOrientation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/ppnParameters.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

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

//! Function to create environment for general observation partial tests.
SystemOfBodies setupEnvironment( const std::vector< std::pair< std::string, std::string > > groundStations,
                                 const double initialEphemerisTime = 1.0E7,
                                 const double finalEphemerisTime = 1.2E7,
                                 const double stateEvaluationTime = 0.0,
                                 const bool useConstantEphemerides = 1,
                                 const double gravitationalParameterScaling = 1.0,
                                 const bool useConstantRotationalEphemeris = false,
                                 const bool moveMarsToMoon = false,
                                 const bool useOblateEarthShape = false );

//! Function to create estimated parameters for general observation partial tests.
std::shared_ptr< EstimatableParameterSet< double > > createEstimatableParameters( const SystemOfBodies& bodies,
                                                                                  const double initialTime,
                                                                                  const bool useEquivalencePrincipleParameter = false,
                                                                                  const bool useRotationalStateAsParameter = false );

//! Function to compute numerical partials w.r.t. constant body states for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyState(
        const std::string& bodyName,
        const SystemOfBodies& bodies,
        const Eigen::Vector3d& bodyPositionVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction,
        const double observationTime,
        const int observableSize );

//! Function to compute numerical partials w.r.t. constant body orientation for general observation partial tests.
Eigen::MatrixXd calculateChangeDueToConstantBodyOrientation( const std::string& bodyName,
                                                             const SystemOfBodies& bodies,
                                                             const Eigen::Vector4d& bodyQuaternionVariation,
                                                             const std::function< Eigen::VectorXd( const double ) > observationFunction,
                                                             const double observationTime,
                                                             const int observableSize,
                                                             std::vector< Eigen::Vector4d >& appliedQuaternionPerturbation );

//! Function to compute numerical partials w.r.t. constant body angular velocity for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyAngularVelocityVector(
        const std::string& bodyName,
        const SystemOfBodies& bodies,
        const Eigen::Vector3d& bodyRotationVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction,
        const double observationTime,
        const int observableSize );

//! Function to compute numerical partials w.r.t. constant body states for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyVelocity(
        const std::string& bodyName,
        const SystemOfBodies& bodies,
        const Eigen::Vector3d& bodyVelocityVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction,
        const double observationTime,
        const int observableSize );

//! Function to compute numerical partials w.r.t. double parameters for general observation partial tests.
std::vector< Eigen::VectorXd > calculateNumericalPartialsWrtDoubleParameters(
        const std::vector< std::shared_ptr< EstimatableParameter< double > > >& doubleParameters,
        const std::vector< std::function< void( ) > > updateFunctions,
        const std::vector< double >& parameterPerturbations,
        const std::function< Eigen::VectorXd( const double ) > observationFunction,
        const double observationTime );

//! Function to compute numerical partials w.r.t. vector parameters for general observation partial tests.
std::vector< Eigen::MatrixXd > calculateNumericalPartialsWrtVectorParameters(
        const std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >& vectorParameters,
        const std::vector< std::function< void( ) > > updateFunctions,
        const std::vector< Eigen::VectorXd >& parameterPerturbations,
        const std::function< Eigen::VectorXd( const double ) > observationFunction,
        const double observationTime );

//! Function to compute analytical partials w.r.t. estimated parameters for general observation partial tests.
template< int ObservableSize >
std::vector< std::vector< std::pair< Eigen::Matrix< double, ObservableSize, Eigen::Dynamic >, double > > > calculateAnalyticalPartials(
        const std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservableSize > > >& partialObjectList,
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const LinkEndType linkEndOfFixedTime,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr,
        const Eigen::Matrix< double, ObservableSize, Eigen::Dynamic > currentObservation =
                Eigen::Matrix< double, ObservableSize, Eigen::Dynamic >::Constant( ObservableSize, TUDAT_NAN ) )
{
    std::vector< std::vector< std::pair< Eigen::Matrix< double, ObservableSize, Eigen::Dynamic >, double > > > partialList;

    for( typename std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservableSize > > >::const_iterator
                 partialIterator = partialObjectList.begin( );
         partialIterator != partialObjectList.end( );
         partialIterator++ )
    {
        partialList.push_back(
                partialIterator->second->calculatePartial( states, times, linkEndOfFixedTime, ancillarySettings, currentObservation ) );
    }
    return partialList;
}

//! Function to retrieve times associated with partial derivatives for general observation partial tests.
std::vector< std::vector< double > > getAnalyticalPartialEvaluationTimes(
        const LinkDefinition& linkEnds,
        const ObservableType observableType,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< EstimatableParameterSet< double > >& estimatedParameters );

//! Generalized test function for partial derivatives of observations
template< int ObservableSize = 1, typename TimeType = double >
void testObservationPartials(
        const std::shared_ptr< ObservationModel< ObservableSize, double, TimeType > > observationModel,
        SystemOfBodies& bodies,
        const std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet,
        const LinkDefinition& linkEnds,
        const ObservableType observableType,
        const double tolerance = 1.0E-6,
        const bool testPositionPartial = 1,
        const bool testParameterPartial = 1,
        const double positionPerturbationMultiplier = 1.0,
        const Eigen::VectorXd parameterPerturbationMultipliers = Eigen::VectorXd::Constant( 4, 1.0 ),
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr,
        double observationTime = 1.1E7,
        const double gammaToleranceWeakening = 1.0 )
{
    bool isNormalized = false;
    if( std::dynamic_pointer_cast< AngularPositionObservationModel< double, TimeType > >( observationModel ) != nullptr )
    {
        isNormalized = std::dynamic_pointer_cast< AngularPositionObservationModel< double, TimeType > >( observationModel )
                               ->getNormalizeRightAscension( );
    }
    printEstimatableParameterEntries( fullEstimatableParameterSet );

    // Retrieve double and vector parameters and estimate body states
    std::vector< std::string > bodiesWithEstimatedTranslationalState =
            estimatable_parameters::getListOfBodiesToEstimate( fullEstimatableParameterSet ).at( propagators::translational_state );
    int numberOfBodiesWithEstimatedTranslationalState = bodiesWithEstimatedTranslationalState.size( );

    std::vector< std::string > bodiesWithEstimatedRotationalState;
    if( estimatable_parameters::getListOfBodiesToEstimate( fullEstimatableParameterSet ).count( propagators::rotational_state ) > 0 )
    {
        bodiesWithEstimatedRotationalState =
                estimatable_parameters::getListOfBodiesToEstimate( fullEstimatableParameterSet ).at( propagators::rotational_state );
    }

    std::vector< std::shared_ptr< EstimatableParameter< double > > > doubleParameterVector =
            fullEstimatableParameterSet->getEstimatedDoubleParameters( );
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameterVector =
            fullEstimatableParameterSet->getEstimatedVectorParameters( );

    // Create observation partials.
    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservableSize > > >,
               std::shared_ptr< PositionPartialScaling > >
            fullAnalyticalPartialSet = ObservationPartialCreator< ObservableSize, double, TimeType >::createObservationPartials(
                    observationModel, bodies, fullEstimatableParameterSet );

    std::shared_ptr< PositionPartialScaling > positionPartialScaler = fullAnalyticalPartialSet.second;

    // Iterate over link ends, compute and test partials for observable referenced at each link end.
    LinkEnds currentLinkEnds = linkEnds.linkEnds_;
    for( LinkEnds::const_iterator linkEndIterator = currentLinkEnds.begin( ); linkEndIterator != currentLinkEnds.end( ); linkEndIterator++ )
    {
        bool runSimulation = true;
        if( ( observableType == observation_models::relative_angular_position ) && ( linkEndIterator->first != receiver ) )
        {
            runSimulation = false;
        }

        if( ( observableType == observation_models::position_angle || observableType == observation_models::separation_distance ||
              observableType == observation_models::position_angle_and_separation ) &&
            ( linkEndIterator->first != receiver ) )
        {
            runSimulation = false;
        }

        if( ( observableType == observation_models::n_way_differenced_range ) && ( linkEndIterator->first == retransmitter ) )
        {
            runSimulation = false;
        }

        if( observableType == observation_models::dsn_n_way_averaged_doppler && linkEndIterator->first != receiver )
        {
            runSimulation = false;
        }

        if( observableType == observation_models::dsn_n_way_range && linkEndIterator->first != receiver )
        {
            runSimulation = false;
        }

        if( observableType == observation_models::differenced_time_of_arrival && linkEndIterator->first != receiver )
        {
            runSimulation = false;
        }

        if( observableType == observation_models::doppler_measured_frequency && linkEndIterator->first != receiver )
        {
            runSimulation = false;
        }

        if( observableType == observation_models::one_way_doppler_measured_frequency && linkEndIterator->first != receiver )
        {
            runSimulation = false;
        }

        if( observableType == observation_models::differenced_frequency_of_arrival && linkEndIterator->first != receiver )
        {
            runSimulation = false;
        }

        if( observableType == observation_models::pixel_coordinates && linkEndIterator->first != receiver )
        {
            runSimulation = false;
        }

        // Remove retransmission delay from the retransmitting reference link end: computation of multi-leg light currently doesn't support
        // retransmission delays at the reference link end
        std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > modifiedAncillarySettings;
        if( observableType == n_way_range && ( linkEndIterator->first != receiver && linkEndIterator->first != transmitter ) )
        {
            std::vector< double > delays = ancillarySettings->getAncillaryDoubleVectorData( link_ends_delays, false );
            // If delays are specified including transmission, retransmission and reception
            if( delays.size( ) == currentLinkEnds.size( ) )
            {
                delays.at( static_cast< int >( linkEndIterator->first ) ) = 0.0;
            }
            // If only retransmission delays are specified
            else
            {
                delays.at( static_cast< int >( linkEndIterator->first ) - 1 ) = 0.0;
            }
            modifiedAncillarySettings = std::make_shared< ObservationAncillarySimulationSettings >( );
            modifiedAncillarySettings->setAncillaryDoubleVectorData( link_ends_delays, delays );
        }
        else
        {
            modifiedAncillarySettings = ancillarySettings;
        }

        if( runSimulation )
        {
            // Evaluate nominal observation values
            std::vector< Eigen::Vector6d > vectorOfStates;
            std::vector< double > vectorOfTimes;
            Eigen::VectorXd currentObservation = observationModel->computeObservationsWithLinkEndData(
                    observationTime, linkEndIterator->first, vectorOfTimes, vectorOfStates, modifiedAncillarySettings );

            // Calculate analytical observation partials.
            if( positionPartialScaler != nullptr )
            {
                positionPartialScaler->update( vectorOfStates,
                                               vectorOfTimes,
                                               static_cast< LinkEndType >( linkEndIterator->first ),
                                               currentObservation,
                                               modifiedAncillarySettings );
            }

            typedef std::vector< std::pair< Eigen::Matrix< double, ObservableSize, Eigen::Dynamic >, double > >
                    ObservationPartialReturnType;
            std::vector< ObservationPartialReturnType > analyticalObservationPartials =
                    calculateAnalyticalPartials< ObservableSize >( fullAnalyticalPartialSet.first,
                                                                   vectorOfStates,
                                                                   vectorOfTimes,
                                                                   linkEndIterator->first,
                                                                   modifiedAncillarySettings,
                                                                   currentObservation );

            // Set and test expected partial size and time
            if( observableType != euler_angle_313_observable )
            {
                std::vector< std::vector< double > > expectedPartialTimes =
                        getAnalyticalPartialEvaluationTimes( linkEnds, observableType, vectorOfTimes, fullEstimatableParameterSet );

                // Test analytical partial times.
                BOOST_CHECK_EQUAL( analyticalObservationPartials.size( ), expectedPartialTimes.size( ) );

                for( unsigned int i = 0; i < analyticalObservationPartials.size( ); i++ )
                {
                    if( observableType == two_way_doppler || observableType == doppler_measured_frequency )
                    {
                        std::vector< double > currentTimes = expectedPartialTimes.at( i );
                        if( currentTimes.size( ) == 2 )
                        {
                            if( linkEndIterator->first == transmitter )
                            {
                                currentTimes.insert( currentTimes.begin( ) + 1, currentTimes.at( 1 ) );
                            }
                            else if( linkEndIterator->first == receiver )
                            {
                                currentTimes.insert( currentTimes.begin( ) + 1, currentTimes.at( 0 ) );
                            }
                            expectedPartialTimes.at( i ) = currentTimes;
                        }
                    }

                    // Associated times for partial derivatives w.r.t. gamma not yet fully consistent (no impact on estimation)
                    if( i < 2 )
                    {
                        BOOST_CHECK_EQUAL( analyticalObservationPartials.at( i ).size( ), expectedPartialTimes.at( i ).size( ) );
                    }

                    for( unsigned int j = 0; j < expectedPartialTimes.at( i ).size( ); j++ )
                    {
                        BOOST_CHECK_EQUAL( analyticalObservationPartials.at( i ).at( j ).second, expectedPartialTimes.at( i ).at( j ) );
                    }
                }
            }

            // Define observation function for current observable/link end
            std::function< Eigen::VectorXd( const double ) > observationFunction =
                    std::bind( &ObservationModel< ObservableSize, double, TimeType >::computeObservations,
                               observationModel,
                               std::placeholders::_1,
                               linkEndIterator->first,
                               modifiedAncillarySettings );

            if( testPositionPartial )
            {
                // Set position perturbation for numerical partial
                Eigen::Vector3d bodyPositionVariation;
                bodyPositionVariation << 1000.0E3, 1000.0E3, 1000.0E3;
                bodyPositionVariation *= positionPerturbationMultiplier;

                // Calculate numerical partials w.r.t. estimate body state.
                Eigen::Matrix< double, Eigen::Dynamic, 3 > bodyPositionPartial = Eigen::Matrix< double, ObservableSize, 3 >::Zero( );
                for( unsigned int i = 0; i < bodiesWithEstimatedTranslationalState.size( ); i++ )
                {
                    // Compute numerical position partial
                    Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalPartialWrtBodyPosition =
                            calculatePartialWrtConstantBodyState( bodiesWithEstimatedTranslationalState.at( i ),
                                                                  bodies,
                                                                  bodyPositionVariation,
                                                                  observationFunction,
                                                                  observationTime,
                                                                  ObservableSize );

                    // Set total analytical partial
                    bodyPositionPartial.setZero( );
                    for( unsigned int j = 0; j < analyticalObservationPartials.at( i ).size( ); j++ )
                    {
                        bodyPositionPartial += analyticalObservationPartials.at( i ).at( j ).first.block( 0, 0, ObservableSize, 3 );
                    }

                    // Test position partial
                    if( ( ( observableType != angular_position ) || ( isNormalized == true ) ) &&
                        ( observableType != relative_angular_position ) )
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

            Eigen::Vector4d quaternionVariation;
            quaternionVariation << 1.0E-6, 1.0E-6, 1.0E-6, 1.0E-6;
            quaternionVariation *= positionPerturbationMultiplier;

            Eigen::Vector3d angularVelocityVariation;
            angularVelocityVariation << 1.0E-9, 1.0E-9, 1.0E-9;
            angularVelocityVariation *= positionPerturbationMultiplier;

            Eigen::Matrix< double, Eigen::Dynamic, 7 > bodyRotationalStatePartial = Eigen::Matrix< double, ObservableSize, 7 >::Zero( );
            for( unsigned int i = 0; i < bodiesWithEstimatedRotationalState.size( ); i++ )
            {
                // Compute numerical position partial
                Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalPartialWrtAngularVelocityVector =
                        calculatePartialWrtConstantBodyAngularVelocityVector( bodiesWithEstimatedRotationalState.at( i ),
                                                                              bodies,
                                                                              angularVelocityVariation,
                                                                              observationFunction,
                                                                              observationTime,
                                                                              ObservableSize );

                std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
                Eigen::MatrixXd changeDueToQuaternionChange =
                        calculateChangeDueToConstantBodyOrientation( bodiesWithEstimatedRotationalState.at( i ),
                                                                     bodies,
                                                                     quaternionVariation,
                                                                     observationFunction,
                                                                     observationTime,
                                                                     ObservableSize,
                                                                     appliedQuaternionPerturbation );

                int indexToUse = i + numberOfBodiesWithEstimatedTranslationalState;
                if( observableType == euler_angle_313_observable )
                {
                    indexToUse = i;
                }
                // Set total analytical partial
                bodyRotationalStatePartial.setZero( );
                for( unsigned int j = 0; j < analyticalObservationPartials[ indexToUse ].size( ); j++ )
                {
                    bodyRotationalStatePartial +=
                            analyticalObservationPartials[ indexToUse ].at( j ).first.block( 0, 0, ObservableSize, 7 );
                }

                for( int j = 0; j < 4; j++ )
                {
                    Eigen::MatrixXd testPartial = ( bodyRotationalStatePartial.block( 0, 0, ObservableSize, 4 ) *
                                                    appliedQuaternionPerturbation.at( j ).segment( 0, 4 ) );
                    BOOST_CHECK_SMALL( std::fabs( testPartial( 0 ) - changeDueToQuaternionChange( j ) ), 1.0E-4 );
                }

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( bodyRotationalStatePartial.block( 0, 4, ObservableSize, 3 ) ),
                                                   numericalPartialWrtAngularVelocityVector,
                                                   std::numeric_limits< double >::epsilon( ) );
            }

            if( testParameterPartial )
            {
                // Test double parameter partials
                {
                    // Settings for parameter partial functions.
                    std::vector< double > parameterPerturbations;
                    if( bodiesWithEstimatedRotationalState.size( ) == 0 )
                    {
                        parameterPerturbations.push_back( 1.0E-10 * parameterPerturbationMultipliers( 0 ) );
                        parameterPerturbations.push_back( 1.0E-10 * parameterPerturbationMultipliers( 1 ) );
                        parameterPerturbations.push_back( 1.0E8 * parameterPerturbationMultipliers( 2 ) );
                    }
                    else
                    {
                        parameterPerturbations.push_back( 1.0E8 * parameterPerturbationMultipliers( 2 ) );
                    }

                    std::vector< std::function< void( ) > > updateFunctionList;

                    updateFunctionList.push_back( emptyVoidFunction );
                    updateFunctionList.push_back( emptyVoidFunction );
                    updateFunctionList.push_back( emptyVoidFunction );

                    // Compute numerical partials.
                    std::vector< Eigen::VectorXd > numericalPartialsWrtDoubleParameters = calculateNumericalPartialsWrtDoubleParameters(
                            doubleParameterVector, updateFunctionList, parameterPerturbations, observationFunction, observationTime );

                    // Compute analytical partial and test against numerical partial
                    Eigen::VectorXd currentParameterPartial;
                    int numberOfEstimatedBodies =
                            bodiesWithEstimatedTranslationalState.size( ) + bodiesWithEstimatedRotationalState.size( );

                    for( unsigned int i = 0; i < numericalPartialsWrtDoubleParameters.size( ); i++ )
                    {
                        currentParameterPartial.setZero( ObservableSize );
                        for( unsigned int j = 0; j < analyticalObservationPartials[ i + numberOfEstimatedBodies ].size( ); j++ )
                        {
                            currentParameterPartial += analyticalObservationPartials[ i + numberOfEstimatedBodies ].at( j ).first;
                        }
                        std::cout << "Current double partial " << i << " " << std::setprecision( 16 )
                                  << analyticalObservationPartials[ i + numberOfEstimatedBodies ].size( ) << ", " << currentParameterPartial
                                  << ", " << numericalPartialsWrtDoubleParameters.at( i ) << ", "
                                  << ( currentParameterPartial( 0 ) - numericalPartialsWrtDoubleParameters.at( i )( 0 ) ) /
                                        numericalPartialsWrtDoubleParameters.at( i )( 0 )
                                  << std::endl;

                        double toleranceToUse = tolerance;
                        if( i == 2 )
                        {
                            toleranceToUse *= gammaToleranceWeakening;
                        }
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                currentParameterPartial, ( numericalPartialsWrtDoubleParameters.at( i ) ), toleranceToUse );
                    }
                }

                // Test vector parameter partials
                {
                    std::function< Eigen::Matrix< double, ObservableSize, 1 >( TimeType ) > vectorObservationFunction =
                            std::bind( &ObservationModel< ObservableSize, double, TimeType >::computeObservations,
                                       observationModel,
                                       std::placeholders::_1,
                                       linkEndIterator->first,
                                       modifiedAncillarySettings );

                    // Settings for parameter partial functions.
                    std::vector< Eigen::VectorXd > parameterPerturbations;
                    parameterPerturbations.push_back( Eigen::Vector2d::Constant( 1.0E-4 * parameterPerturbationMultipliers( 3 ) ) );
                    parameterPerturbations.push_back( Eigen::Vector2d::Constant( 1.0E-4 * parameterPerturbationMultipliers( 3 ) ) );

                    std::vector< std::function< void( ) > > updateFunctionList;
                    updateFunctionList.push_back( emptyVoidFunction );
                    updateFunctionList.push_back( emptyVoidFunction );

                    // Compute numerical partials.
                    std::vector< Eigen::MatrixXd > numericalPartialsWrtVectorParameters = calculateNumericalPartialsWrtVectorParameters(
                            vectorParameterVector, updateFunctionList, parameterPerturbations, vectorObservationFunction, observationTime );

                    // Compute analytical partial and test against numerical partial
                    Eigen::MatrixXd currentParameterPartial;
                    int startIndex = bodiesWithEstimatedTranslationalState.size( ) + bodiesWithEstimatedRotationalState.size( ) +
                            doubleParameterVector.size( );

                    for( unsigned int i = 0; i < numericalPartialsWrtVectorParameters.size( ); i++ )
                    {
                        currentParameterPartial =
                                Eigen::MatrixXd::Zero( ObservableSize, vectorParameterVector.at( i )->getParameterSize( ) );

                        for( unsigned int j = 0; j < analyticalObservationPartials[ i + startIndex ].size( ); j++ )
                        {
                            currentParameterPartial += analyticalObservationPartials[ i + startIndex ].at( j ).first;
                        }

                        std::cout << "Current vector partial " << i << " " << currentParameterPartial << " "
                                  << numericalPartialsWrtVectorParameters.at( i ) << std::endl;

                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                ( currentParameterPartial ), ( numericalPartialsWrtVectorParameters.at( i ) ), tolerance );
                    }
                }
            }
        }
    }
}

//! Create the gravitational parameters of the given bodies, and optionally PPN gamma, as the only estimated parameters.
inline std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > createShapiroParameters(
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::string >& gravitationalParameterBodies,
        const bool estimateGamma,
        std::vector< double >& parameterPerturbations,
        const double perturbationScaling = 1.0 )
{
    using namespace tudat::estimatable_parameters;

    // All parameters enter the (bending-free) Shapiro delay linearly, so large perturbations are used to suppress
    // round-off in the numerical partials. For planets, the Shapiro Doppler is so small that the perturbation must exceed
    // the nominal GM by orders of magnitude to resolve it. Observables with a large nominal value (e.g. frequencies) need
    // a further perturbationScaling.
    std::vector< std::shared_ptr< EstimatableParameter< double > > > doubleParameters;
    parameterPerturbations.clear( );
    for( const std::string& bodyName : gravitationalParameterBodies )
    {
        doubleParameters.push_back(
                std::make_shared< GravitationalParameter >( bodies.at( bodyName )->getGravityFieldModel( ), bodyName ) );
        parameterPerturbations.push_back( perturbationScaling * std::max( 0.5 * doubleParameters.back( )->getParameterValue( ), 1.0E19 ) );
    }
    if( estimateGamma )
    {
        doubleParameters.push_back( std::make_shared< PPNParameterGamma >( bodies.getSpaceTimeProperties( )->getPpnParameterSet( ) ) );
        parameterPerturbations.push_back( perturbationScaling );
    }
    return std::make_shared< EstimatableParameterSet< double > >(
            doubleParameters, std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( ) );
}

//! Compare analytical partials w.r.t. double parameters with central-difference partials of the observable.
/*!
 * Intended for parameters that enter the observation model directly (e.g. gravitational parameters of light-time
 * perturbing bodies, PPN gamma), so the parameter set must contain only double parameters. The observable is evaluated
 * for each requested reference link end, and the largest relative error per parameter is returned (for diagnostics).
 */
template< int ObservableSize = 1, typename TimeType = double >
std::vector< double > testDirectDoubleParameterPartials(
        const std::shared_ptr< ObservationModel< ObservableSize, double, TimeType > > observationModel,
        const SystemOfBodies& bodies,
        const std::shared_ptr< EstimatableParameterSet< double > > parameterSet,
        const std::vector< LinkEndType >& referenceLinkEnds,
        const std::vector< double >& parameterPerturbations,
        const double tolerance,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr,
        const double observationTime = 1.1E7 )
{
    const std::vector< std::shared_ptr< EstimatableParameter< double > > > doubleParameters = parameterSet->getEstimatedDoubleParameters( );
    if( parameterSet->getParameterSetSize( ) != static_cast< int >( doubleParameters.size( ) ) ||
        parameterPerturbations.size( ) != doubleParameters.size( ) )
    {
        throw std::runtime_error( "Error in testDirectDoubleParameterPartials: only double parameters are supported." );
    }

    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservableSize > > >,
               std::shared_ptr< PositionPartialScaling > >
            partialSet = ObservationPartialCreator< ObservableSize, double, TimeType >::createObservationPartials(
                    observationModel, bodies, parameterSet );

    std::vector< double > maximumRelativeErrors( doubleParameters.size( ), 0.0 );
    for( const LinkEndType referenceLinkEnd : referenceLinkEnds )
    {
        std::vector< Eigen::Vector6d > linkEndStates;
        std::vector< double > linkEndTimes;
        Eigen::Matrix< double, ObservableSize, 1 > nominalObservation = observationModel->computeObservationsWithLinkEndData(
                observationTime, referenceLinkEnd, linkEndTimes, linkEndStates, ancillarySettings );
        if( partialSet.second != nullptr )
        {
            partialSet.second->update( linkEndStates, linkEndTimes, referenceLinkEnd, nominalObservation );
        }

        // Compute all analytical partials before any perturbed evaluation: light-time correction partials use the
        // correction cached by the last light-time computation.
        std::vector< Eigen::Matrix< double, ObservableSize, 1 > > analyticalPartials( doubleParameters.size( ),
                                                                                      Eigen::Matrix< double, ObservableSize, 1 >::Zero( ) );
        for( const auto& partialEntry : partialSet.first )
        {
            const auto partialList = partialEntry.second->calculatePartial(
                    linkEndStates, linkEndTimes, referenceLinkEnd, ancillarySettings, nominalObservation );
            for( const auto& partial : partialList )
            {
                analyticalPartials.at( partialEntry.first.first ) += partial.first;
            }
        }

        for( unsigned int i = 0; i < doubleParameters.size( ); i++ )
        {
            const Eigen::Matrix< double, ObservableSize, 1 >& analyticalPartial = analyticalPartials.at( i );
            const double nominalValue = doubleParameters.at( i )->getParameterValue( );
            doubleParameters.at( i )->setParameterValue( nominalValue + parameterPerturbations.at( i ) );
            const Eigen::Matrix< double, ObservableSize, 1 > upObservation =
                    observationModel->computeObservations( observationTime, referenceLinkEnd, ancillarySettings );
            doubleParameters.at( i )->setParameterValue( nominalValue - parameterPerturbations.at( i ) );
            const Eigen::Matrix< double, ObservableSize, 1 > downObservation =
                    observationModel->computeObservations( observationTime, referenceLinkEnd, ancillarySettings );
            doubleParameters.at( i )->setParameterValue( nominalValue );

            const Eigen::Matrix< double, ObservableSize, 1 > numericalPartial =
                    ( upObservation - downObservation ) / ( 2.0 * parameterPerturbations.at( i ) );
            const double relativeError = ( analyticalPartial - numericalPartial ).norm( ) / numericalPartial.norm( );
            maximumRelativeErrors.at( i ) = std::max( maximumRelativeErrors.at( i ), relativeError );

            std::cout << "Direct parameter partial " << getParameterTypeString( doubleParameters.at( i )->getParameterName( ).first ) << " "
                      << doubleParameters.at( i )->getParameterName( ).second.first << ", reference link end "
                      << getLinkEndTypeString( referenceLinkEnd ) << ": analytical " << std::setprecision( 16 )
                      << analyticalPartial.transpose( ) << ", numerical " << numericalPartial.transpose( ) << ", relative error "
                      << relativeError << std::endl;

            BOOST_TEST_CONTEXT( "Parameter " << doubleParameters.at( i )->getParameterDescription( ) << ", reference link end "
                                             << getLinkEndTypeString( referenceLinkEnd ) )
            {
                BOOST_CHECK( numericalPartial.norm( ) > 0.0 );
                BOOST_CHECK_SMALL( relativeError, tolerance );
            }
        }
    }
    return maximumRelativeErrors;
}

//! Check that two matrices agree to within a tolerance relative to the norm of the expected matrix.
inline void checkMatrixCloseRelativeToNorm( const Eigen::MatrixXd& computedMatrix,
                                            const Eigen::MatrixXd& expectedMatrix,
                                            const double tolerance )
{
    BOOST_REQUIRE_EQUAL( computedMatrix.rows( ), expectedMatrix.rows( ) );
    BOOST_REQUIRE_EQUAL( computedMatrix.cols( ), expectedMatrix.cols( ) );
    BOOST_CHECK_SMALL( ( computedMatrix - expectedMatrix ).norm( ) / expectedMatrix.norm( ), tolerance );
}

//! Sum the partial-derivative contributions returned by an observation partial.
template< int ObservationSize >
Eigen::MatrixXd sumPartialContributions(
        const std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > >& partialContributions )
{
    Eigen::MatrixXd summedPartial = partialContributions.at( 0 ).first;
    for( unsigned int i = 1; i < partialContributions.size( ); i++ )
    {
        summedPartial += partialContributions.at( i ).first;
    }
    return summedPartial;
}

//! Compare analytical and numerical partials of a direction observable (angular position, position angle and/or separation) for
//! a fast-moving receiver (Earth/Graz), with transmitters Mars/MSL (and Moon).
/*!
 *  For the receiver as reference link end, the receiver (Earth) is given a constant velocity of order 0.01c, so that, if the ancillary
 *  settings select stellar aberration, the aberration factors in the partials differ from their astrometric values by order 1.0E-2,
 *  well above the tolerance of the numerical-partial comparison. Without stellar aberration, the partials w.r.t. the receiver velocity
 *  must be zero.
 *  \param observationModelSettings Settings of the observation model
 *  \param ancillarySettings Ancillary settings of the observation (celestial frame and stellar aberration)
 *  \param referenceLinkEnd Link end at which the observation time is kept fixed
 */
template< int ObservationSize >
void testAberratedObservationPartials( const std::shared_ptr< ObservationModelSettings > observationModelSettings,
                                       const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings,
                                       const LinkEndType referenceLinkEnd = receiver )
{
    std::vector< std::pair< std::string, std::string > > groundStations = { std::make_pair( "Earth", "Graz" ),
                                                                            std::make_pair( "Mars", "MSL" ),
                                                                            std::make_pair( "Moon", "" ) };
    SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

    // With the transmitter as reference link end, the light-time partials use the receiver velocity, which is not consistent with the
    // (fixed) position of a constant ephemeris, so that the receiver velocity is only set for the receiver as reference link end
    if( referenceLinkEnd == receiver )
    {
        std::shared_ptr< ConstantEphemeris > earthEphemeris =
                std::dynamic_pointer_cast< ConstantEphemeris >( bodies.at( "Earth" )->getEphemeris( ) );
        Eigen::Vector6d earthState = earthEphemeris->getCartesianState( 0.0 );
        earthState.segment< 3 >( 3 ) = Eigen::Vector3d( 1.2E6, -2.1E6, 0.8E6 );
        earthEphemeris->updateConstantState( earthState );
        bodies.at( "Earth" )->recomputeStateOnNextCall( );
    }

    std::shared_ptr< ObservationModel< ObservationSize > > observationModel =
            ObservationModelCreator< ObservationSize, double, double >::createObservationModel( observationModelSettings, bodies );
    std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createEstimatableParameters( bodies, 1.1E7 );

    const auto partialSet = ObservationPartialCreator< ObservationSize, double, double >::createObservationPartials(
            observationModel, bodies, parameterSet );

    const double observationTime = 1.1E7;
    std::vector< Eigen::Vector6d > linkEndStates;
    std::vector< double > linkEndTimes;
    const Eigen::VectorXd observation = observationModel->computeObservationsWithLinkEndData(
            observationTime, referenceLinkEnd, linkEndTimes, linkEndStates, ancillarySettings );
    partialSet.second->update( linkEndStates, linkEndTimes, referenceLinkEnd, observation, ancillarySettings );

    std::function< Eigen::VectorXd( const double ) > observationFunction =
            std::bind( &ObservationModel< ObservationSize >::computeObservations,
                       observationModel,
                       std::placeholders::_1,
                       referenceLinkEnd,
                       ancillarySettings );

    // Initial states of Earth (receiver) and Mars (transmitter) are the first two parameters
    const std::vector< std::string > bodiesWithEstimatedState = { "Earth", "Mars" };
    for( unsigned int i = 0; i < bodiesWithEstimatedState.size( ); i++ )
    {
        const std::pair< int, int > parameterIndices = std::make_pair( 6 * i, 6 );
        BOOST_REQUIRE( partialSet.first.count( parameterIndices ) == 1 );
        const Eigen::MatrixXd analyticalPartial = sumPartialContributions< ObservationSize >(
                partialSet.first.at( parameterIndices )
                        ->calculatePartial( linkEndStates, linkEndTimes, referenceLinkEnd, ancillarySettings, observation ) );

        const Eigen::MatrixXd numericalPositionPartial = calculatePartialWrtConstantBodyState( bodiesWithEstimatedState.at( i ),
                                                                                               bodies,
                                                                                               Eigen::Vector3d::Constant( 1.0E5 ),
                                                                                               observationFunction,
                                                                                               observationTime,
                                                                                               ObservationSize );
        checkMatrixCloseRelativeToNorm( analyticalPartial.block( 0, 0, ObservationSize, 3 ), numericalPositionPartial, 1.0E-6 );

        // Only the receiver velocity enters, through stellar aberration
        const Eigen::MatrixXd numericalVelocityPartial = calculatePartialWrtConstantBodyVelocity( bodiesWithEstimatedState.at( i ),
                                                                                                  bodies,
                                                                                                  Eigen::Vector3d::Constant( 1.0E2 ),
                                                                                                  observationFunction,
                                                                                                  observationTime,
                                                                                                  ObservationSize );
        if( bodiesWithEstimatedState.at( i ) == "Earth" && getPositionAngleUseStellarAberration( ancillarySettings ) )
        {
            BOOST_CHECK_GT( numericalVelocityPartial.norm( ), 1.0E-10 );
            checkMatrixCloseRelativeToNorm( analyticalPartial.block( 0, 3, ObservationSize, 3 ), numericalVelocityPartial, 1.0E-6 );
        }
        else
        {
            BOOST_CHECK_SMALL( analyticalPartial.block( 0, 3, ObservationSize, 3 ).norm( ), std::numeric_limits< double >::min( ) );
        }
    }

    // Earth rotation rate changes the receiver velocity as well as its position
    std::shared_ptr< EstimatableParameter< double > > earthRotationRate;
    std::pair< int, int > rotationRateIndices;
    for( const auto& doubleParameter : parameterSet->getDoubleParameters( ) )
    {
        if( doubleParameter.second->getParameterName( ).first == constant_rotation_rate &&
            doubleParameter.second->getParameterName( ).second.first == "Earth" )
        {
            earthRotationRate = doubleParameter.second;
            rotationRateIndices = std::make_pair( doubleParameter.first, 1 );
        }
    }
    BOOST_REQUIRE( earthRotationRate != nullptr );
    BOOST_REQUIRE( partialSet.first.count( rotationRateIndices ) == 1 );
    const Eigen::MatrixXd analyticalRotationRatePartial = sumPartialContributions< ObservationSize >(
            partialSet.first.at( rotationRateIndices )
                    ->calculatePartial( linkEndStates, linkEndTimes, referenceLinkEnd, ancillarySettings, observation ) );
    const std::vector< Eigen::VectorXd > numericalRotationRatePartial = calculateNumericalPartialsWrtDoubleParameters(
            { earthRotationRate }, { []( ) {} }, { 1.0E-10 }, observationFunction, observationTime );
    checkMatrixCloseRelativeToNorm( analyticalRotationRatePartial, numericalRotationRatePartial.at( 0 ), 1.0E-5 );
}

}  // namespace unit_tests

}  // namespace tudat
