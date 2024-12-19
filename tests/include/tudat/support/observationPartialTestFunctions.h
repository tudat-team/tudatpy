/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <boost/lambda/lambda.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationalOrientation.h"
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
                                 const bool useConstantRotationalEphemeris = false );

//! Function to create estimated parameters for general observation partial tests.
std::shared_ptr< EstimatableParameterSet< double > > createEstimatableParameters(
        const SystemOfBodies& bodies, const double initialTime,
        const bool useEquivalencePrincipleParameter = false,
        const bool useRotationalStateAsParameter = false );

//! Function to compute numerical partials w.r.t. constant body states for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyState(
        const std::string& bodyName, const SystemOfBodies& bodies, const Eigen::Vector3d& bodyPositionVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction,
        const double observationTime, const int observableSize );

//! Function to compute numerical partials w.r.t. constant body orientation for general observation partial tests.
Eigen::MatrixXd calculateChangeDueToConstantBodyOrientation(
        const std::string& bodyName, const SystemOfBodies& bodies, const Eigen::Vector4d& bodyQuaternionVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction, const double observationTime,
        const int observableSize,
        std::vector< Eigen::Vector4d >& appliedQuaternionPerturbation );

//! Function to compute numerical partials w.r.t. constant body angular velocity for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyAngularVelocityVector(
        const std::string& bodyName, const SystemOfBodies& bodies, const Eigen::Vector3d& bodyRotationVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction, const double observationTime,
        const int observableSize );

//! Function to compute numerical partials w.r.t. constant body states for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyVelocity(
        const std::string& bodyName, const SystemOfBodies& bodies, const Eigen::Vector3d& bodyVelocityVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction, const double observationTime,
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
std::vector< std::vector< std::pair< Eigen::Matrix< double, ObservableSize, Eigen::Dynamic >, double > > >
calculateAnalyticalPartials(
        const std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservableSize > > >& partialObjectList,
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const LinkEndType linkEndOfFixedTime,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr,
        const Eigen::Matrix< double, ObservableSize, Eigen::Dynamic > currentObservation =
            Eigen::Matrix< double, ObservableSize, Eigen::Dynamic >::Constant( ObservableSize, TUDAT_NAN ) )
{
    std::vector< std::vector< std::pair< Eigen::Matrix< double, ObservableSize, Eigen::Dynamic >, double > > > partialList;

    for( typename std::map< std::pair< int, int >,
         std::shared_ptr< ObservationPartial< ObservableSize > > >::const_iterator partialIterator =
         partialObjectList.begin( ); partialIterator != partialObjectList.end( ); partialIterator++ )
    {
        partialList.push_back( partialIterator->second->calculatePartial( states, times, linkEndOfFixedTime, ancilliarySettings, currentObservation ) );
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
        const LinkDefinition& linkEnds, const ObservableType observableType,
        const double tolerance = 1.0E-6,
        const bool testPositionPartial = 1,
        const bool testParameterPartial = 1,
        const double positionPerturbationMultiplier = 1.0,
        const Eigen::VectorXd parameterPerturbationMultipliers = Eigen::VectorXd::Constant( 4, 1.0 ),
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr,
        double observationTime = 1.1E7,
        const double gammaToleranceWeakening = 1.0 )
{
    printEstimatableParameterEntries( fullEstimatableParameterSet );

    // Retrieve double and vector parameters and estimate body states
    std::vector<std::string> bodiesWithEstimatedTranslationalState =
            estimatable_parameters::getListOfBodiesToEstimate(
                fullEstimatableParameterSet).at(propagators::translational_state);
    int numberOfBodiesWithEstimatedTranslationalState =
            bodiesWithEstimatedTranslationalState.size();

    std::vector<std::string> bodiesWithEstimatedRotationalState;
    if (estimatable_parameters::getListOfBodiesToEstimate(
                fullEstimatableParameterSet).count(propagators::rotational_state) > 0) {
        bodiesWithEstimatedRotationalState =
                estimatable_parameters::getListOfBodiesToEstimate(
                    fullEstimatableParameterSet).at(propagators::rotational_state);
    }

    std::vector<std::shared_ptr<EstimatableParameter<double> > > doubleParameterVector =
            fullEstimatableParameterSet->getEstimatedDoubleParameters();
    std::vector<std::shared_ptr<EstimatableParameter<Eigen::VectorXd> > > vectorParameterVector =
            fullEstimatableParameterSet->getEstimatedVectorParameters();

    // Create observation partials.
    std::pair<std::map<std::pair<int, int>, std::shared_ptr<ObservationPartial<ObservableSize> > >,
            std::shared_ptr<PositionPartialScaling> > fullAnalyticalPartialSet =
            ObservationPartialCreator< ObservableSize, double, TimeType >::createObservationPartials(
                observationModel, bodies, fullEstimatableParameterSet );

    std::shared_ptr<PositionPartialScaling> positionPartialScaler = fullAnalyticalPartialSet.second;


    // Iterate over link ends, compute and test partials for observable referenced at each link end.
    LinkEnds currentLinkEnds = linkEnds.linkEnds_;
    for (LinkEnds::const_iterator linkEndIterator = currentLinkEnds.begin(); linkEndIterator != currentLinkEnds.end();
         linkEndIterator++)
    {
        bool runSimulation = true;
        if( ( observableType == observation_models::relative_angular_position ) &&
                ( linkEndIterator->first != receiver ) )
        {
            runSimulation = false;
        }

        if( ( observableType == observation_models::n_way_differenced_range ) &&
               ( linkEndIterator->first == retransmitter ) )
        {
            runSimulation = false;
        }

        if( observableType == observation_models::dsn_n_way_averaged_doppler && linkEndIterator->first != receiver )
        {
            runSimulation = false;
        }

        if( observableType == observation_models::doppler_measured_frequency && linkEndIterator->first != receiver )
        {
            runSimulation = false;
        }

        // Remove retransmission delay from the retransmitting reference link end: computation of multi-leg light currently doesn't support
        // retransmission delays at the reference link end
        std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > modifiedAncilliarySettings;
        if ( observableType == n_way_range && ( linkEndIterator->first != receiver && linkEndIterator->first != transmitter ) )
        {
            std::vector< double > delays = ancilliarySettings->getAncilliaryDoubleVectorData( link_ends_delays, false );
            // If delays are specified including transmission, retransmission and reception
            if ( delays.size( ) == currentLinkEnds.size( ) )
            {
                delays.at( static_cast< int >( linkEndIterator->first ) ) = 0.0;
            }
            // If only retransmission delays are specified
            else
            {
                delays.at( static_cast< int >( linkEndIterator->first ) - 1 ) = 0.0;
            }
            modifiedAncilliarySettings = std::make_shared< ObservationAncilliarySimulationSettings >( );
            modifiedAncilliarySettings->setAncilliaryDoubleVectorData( link_ends_delays, delays );
        }
        else
        {
            modifiedAncilliarySettings = ancilliarySettings;
        }

        if ( runSimulation )
        {

            // Evaluate nominal observation values
            std::vector<Eigen::Vector6d> vectorOfStates;
            std::vector<double> vectorOfTimes;
            Eigen::VectorXd currentObservation = observationModel->computeObservationsWithLinkEndData(
                        observationTime, linkEndIterator->first, vectorOfTimes, vectorOfStates, modifiedAncilliarySettings );

            // Calculate analytical observation partials.
            if (positionPartialScaler != nullptr) {
                positionPartialScaler->update(vectorOfStates, vectorOfTimes, static_cast< LinkEndType >( linkEndIterator->first ),
                                              currentObservation);
            }

            typedef std::vector<std::pair<Eigen::Matrix<double, ObservableSize, Eigen::Dynamic>, double> >
                    ObservationPartialReturnType;
            std::vector<ObservationPartialReturnType> analyticalObservationPartials =
                    calculateAnalyticalPartials<ObservableSize>(
                        fullAnalyticalPartialSet.first, vectorOfStates, vectorOfTimes, linkEndIterator->first, modifiedAncilliarySettings, currentObservation);

            // Set and test expected partial size and time
            if (observableType != euler_angle_313_observable)
            {
                std::vector<std::vector<double> > expectedPartialTimes = getAnalyticalPartialEvaluationTimes(
                            linkEnds, observableType, vectorOfTimes, fullEstimatableParameterSet);

                // Test analytical partial times.
                BOOST_CHECK_EQUAL(analyticalObservationPartials.size(), expectedPartialTimes.size());

                for (unsigned int i = 0; i < analyticalObservationPartials.size(); i++)
                {
                    if (observableType == two_way_doppler || observableType == doppler_measured_frequency )
                    {
                        std::vector<double> currentTimes = expectedPartialTimes.at( i );
                        if (currentTimes.size() == 2)
                        {
                            if (linkEndIterator->first == transmitter)
                            {
                                currentTimes.insert(currentTimes.begin() + 1, currentTimes.at(1));
                            }
                            else if (linkEndIterator->first == receiver)
                            {
                                currentTimes.insert(currentTimes.begin() + 1, currentTimes.at(0));
                            }
                            expectedPartialTimes.at( i ) = currentTimes;
                        }
                    }

                    // Associated times for partial derivatives w.r.t. gamma not yet fully consistent (no impact on estimation)
                    if (i < 2)
                    {
                        BOOST_CHECK_EQUAL(analyticalObservationPartials.at( i ).size(), expectedPartialTimes.at( i ).size());
                    }

                    for (unsigned int j = 0; j < expectedPartialTimes.at( i ).size(); j++)
                    {
                        BOOST_CHECK_EQUAL(analyticalObservationPartials.at( i ).at( j ).second, expectedPartialTimes.at( i ).at( j ));
                    }
                }
            }

            // Define observation function for current observable/link end
            std::function<Eigen::VectorXd(const double)> observationFunction = std::bind(
                        &ObservationModel< ObservableSize, double, TimeType >::computeObservations,
                        observationModel, std::placeholders::_1, linkEndIterator->first, modifiedAncilliarySettings );

            if (testPositionPartial)
            {
                // Set position perturbation for numerical partial
                Eigen::Vector3d bodyPositionVariation;
                bodyPositionVariation << 1000.0E3, 1000.0E3, 1000.0E3;
                bodyPositionVariation *= positionPerturbationMultiplier;

                // Calculate numerical partials w.r.t. estimate body state.
                Eigen::Matrix<double, Eigen::Dynamic, 3> bodyPositionPartial =
                        Eigen::Matrix<double, ObservableSize, 3>::Zero();
                for (unsigned int i = 0; i < bodiesWithEstimatedTranslationalState.size(); i++)
                {
                    // Compute numerical position partial
                    Eigen::Matrix<double, Eigen::Dynamic, 3> numericalPartialWrtBodyPosition =
                            calculatePartialWrtConstantBodyState(
                                bodiesWithEstimatedTranslationalState.at( i ), bodies, bodyPositionVariation,
                                observationFunction, observationTime, ObservableSize);

                    // Set total analytical partial
                    bodyPositionPartial.setZero();
                    for (unsigned int j = 0; j < analyticalObservationPartials.at( i ).size(); j++)
                    {
                        bodyPositionPartial += analyticalObservationPartials.at( i ).at( j ).first.block(0, 0, ObservableSize, 3);
                    }

                    // Test position partial
                    if ( ( observableType != angular_position ) && ( observableType != relative_angular_position ) )
                    {
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(bodyPositionPartial, (numericalPartialWrtBodyPosition), tolerance);
                    }
                    else
                    {
                        BOOST_CHECK_SMALL(std::fabs(bodyPositionPartial(0, 2) - numericalPartialWrtBodyPosition(0, 2)),
                                          1.0E-20);
                        bodyPositionPartial(0, 2) = 0.0;
                        numericalPartialWrtBodyPosition(0, 2) = 0.0;

                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(bodyPositionPartial, (numericalPartialWrtBodyPosition), tolerance);

                    }
                }
            }

            Eigen::Vector4d quaternionVariation;
            quaternionVariation << 1.0E-6, 1.0E-6, 1.0E-6, 1.0E-6;
            quaternionVariation *= positionPerturbationMultiplier;

            Eigen::Vector3d angularVelocityVariation;
            angularVelocityVariation << 1.0E-9, 1.0E-9, 1.0E-9;
            angularVelocityVariation *= positionPerturbationMultiplier;

            Eigen::Matrix<double, Eigen::Dynamic, 7> bodyRotationalStatePartial =
                    Eigen::Matrix<double, ObservableSize, 7>::Zero();
            for (unsigned int i = 0; i < bodiesWithEstimatedRotationalState.size(); i++)
            {
                // Compute numerical position partial
                Eigen::Matrix<double, Eigen::Dynamic, 3> numericalPartialWrtAngularVelocityVector =
                        calculatePartialWrtConstantBodyAngularVelocityVector(
                            bodiesWithEstimatedRotationalState.at( i ), bodies, angularVelocityVariation,
                            observationFunction, observationTime, ObservableSize);

                std::vector<Eigen::Vector4d> appliedQuaternionPerturbation;
                Eigen::MatrixXd changeDueToQuaternionChange =
                        calculateChangeDueToConstantBodyOrientation(
                            bodiesWithEstimatedRotationalState.at( i ), bodies, quaternionVariation,
                            observationFunction, observationTime, ObservableSize, appliedQuaternionPerturbation);

                int indexToUse = i + numberOfBodiesWithEstimatedTranslationalState;
                if (observableType == euler_angle_313_observable) {
                    indexToUse = i;
                }
                // Set total analytical partial
                bodyRotationalStatePartial.setZero();
                for (unsigned int j = 0; j < analyticalObservationPartials[indexToUse].size(); j++) {
                    bodyRotationalStatePartial += analyticalObservationPartials[indexToUse].at( j ).first.block(
                                0, 0, ObservableSize, 7);
                }

                for (int j = 0; j < 4; j++) {
                    Eigen::MatrixXd testPartial = (bodyRotationalStatePartial.block(0, 0, ObservableSize, 4) * appliedQuaternionPerturbation.at(j).segment(0, 4));
                    BOOST_CHECK_SMALL(
                                std::fabs(testPartial(0) - changeDueToQuaternionChange(j)), 1.0E-4);
                }

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            (bodyRotationalStatePartial.block(0, 4, ObservableSize, 3)), numericalPartialWrtAngularVelocityVector,
                            std::numeric_limits<double>::epsilon());
            }


            if (testParameterPartial)
            {

                // Test double parameter partials
                {
                    // Settings for parameter partial functions.
                    std::vector<double> parameterPerturbations;
                    if (bodiesWithEstimatedRotationalState.size() == 0) {
                        parameterPerturbations.push_back(1.0E-10 * parameterPerturbationMultipliers(0));
                        parameterPerturbations.push_back(1.0E-10 * parameterPerturbationMultipliers(1));
                        parameterPerturbations.push_back(1.0E8 * parameterPerturbationMultipliers(2));
                    } else {
                        parameterPerturbations.push_back(1.0E8 * parameterPerturbationMultipliers(2));
                    }

                    std::vector<std::function<void()> > updateFunctionList;

                    updateFunctionList.push_back(emptyVoidFunction);
                    updateFunctionList.push_back(emptyVoidFunction);
                    updateFunctionList.push_back(emptyVoidFunction);

                    // Compute numerical partials.
                    std::vector<Eigen::VectorXd> numericalPartialsWrtDoubleParameters =
                            calculateNumericalPartialsWrtDoubleParameters(
                                doubleParameterVector, updateFunctionList, parameterPerturbations,
                                observationFunction, observationTime);

                    // Compute analytical partial and test against numerical partial
                    Eigen::VectorXd currentParameterPartial;
                    int numberOfEstimatedBodies =
                            bodiesWithEstimatedTranslationalState.size() + bodiesWithEstimatedRotationalState.size();

                    for (unsigned int i = 0; i < numericalPartialsWrtDoubleParameters.size(); i++)
                    {
                        currentParameterPartial.setZero(ObservableSize);
                        for (unsigned int j = 0; j < analyticalObservationPartials[i + numberOfEstimatedBodies].size(); j++)
                        {
                            currentParameterPartial += analyticalObservationPartials[i + numberOfEstimatedBodies].at( j ).first;

                        }
                        std::cout<<"Current double partial "<<i<<" "<<std::setprecision( 16 )<<analyticalObservationPartials[i + numberOfEstimatedBodies].size()<<" "<<currentParameterPartial<<" "<<numericalPartialsWrtDoubleParameters.at( i )<<" "
                        <<( currentParameterPartial( 0 ) - numericalPartialsWrtDoubleParameters.at( i )( 0 ) ) / numericalPartialsWrtDoubleParameters.at( i )( 0 )<<std::endl;

                        double toleranceToUse = tolerance;
                        if( i == 2 )
                        {
                            toleranceToUse *= gammaToleranceWeakening;
                        }
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                    currentParameterPartial, (numericalPartialsWrtDoubleParameters.at( i )), toleranceToUse);
                    }
                }

                // Test vector parameter partials
                {
                    std::function< Eigen::Matrix< double, ObservableSize, 1 >( TimeType ) > vectorObservationFunction =
                            std::bind( &ObservationModel< ObservableSize, double, TimeType >::computeObservations,
                                       observationModel, std::placeholders::_1, linkEndIterator->first, modifiedAncilliarySettings );

                    // Settings for parameter partial functions.
                    std::vector<Eigen::VectorXd> parameterPerturbations;
                    parameterPerturbations.push_back(Eigen::Vector2d::Constant(1.0E-4 * parameterPerturbationMultipliers(3)));
                    parameterPerturbations.push_back(Eigen::Vector2d::Constant(1.0E-4 * parameterPerturbationMultipliers(3)));

                    std::vector<std::function<void()> > updateFunctionList;
                    updateFunctionList.push_back(emptyVoidFunction);
                    updateFunctionList.push_back(emptyVoidFunction);

                    // Compute numerical partials.
                    std::vector<Eigen::MatrixXd> numericalPartialsWrtVectorParameters =
                            calculateNumericalPartialsWrtVectorParameters(
                                vectorParameterVector, updateFunctionList, parameterPerturbations,
                                vectorObservationFunction, observationTime);

                    // Compute analytical partial and test against numerical partial
                    Eigen::MatrixXd currentParameterPartial;
                    int startIndex = bodiesWithEstimatedTranslationalState.size() + bodiesWithEstimatedRotationalState.size() +
                            doubleParameterVector.size();

                    for (unsigned int i = 0; i < numericalPartialsWrtVectorParameters.size(); i++)
                    {
                        currentParameterPartial = Eigen::MatrixXd::Zero(
                                    ObservableSize, vectorParameterVector.at(i)->getParameterSize());

                        for (unsigned int j = 0; j < analyticalObservationPartials[i + startIndex].size(); j++)
                        {
                            currentParameterPartial += analyticalObservationPartials[i + startIndex].at( j ).first;

                        }

                        std::cout<<"Current vector partial "<<i<<" "<<currentParameterPartial<<" "<<numericalPartialsWrtVectorParameters.at( i )<<std::endl;

                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                    (currentParameterPartial), (numericalPartialsWrtVectorParameters.at( i )), tolerance);

                    }
                }
            }
        }
    }
}


}

}
