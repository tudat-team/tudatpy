/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Ephemerides/constantRotationalEphemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/observationPartialTestFunctions.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/ppnParameters.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/equivalencePrincipleViolationParameter.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"

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
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;

//! Function to create environment for general observation partial tests.
NamedBodyMap setupEnvironment( const std::vector< LinkEndId > groundStations,
                               const double initialEphemerisTime,
                               const double finalEphemerisTime,
                               const double stateEvaluationTime,
                               const bool useConstantEphemerides,
                               const double gravitationalParameterScaling,
                               const bool useConstantRotationalEphemeris )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies.
    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = std::make_shared< Body >( );
    bodyMap[ "Mars" ] = std::make_shared< Body >( );
    bodyMap[ "Moon" ] = std::make_shared< Body >( );
    bodyMap[ "Sun" ] = std::make_shared< Body >( );

    bodyMap[ "Earth" ]->setShapeModel( std::make_shared< basic_astrodynamics::SphericalBodyShapeModel >(
                                           spice_interface::getAverageRadius( "Earth" ) ) );
    bodyMap[ "Mars" ]->setShapeModel( std::make_shared< basic_astrodynamics::SphericalBodyShapeModel >(
                                          spice_interface::getAverageRadius( "Mars" ) ) );

    if( useConstantEphemerides )
    {
        Eigen::Vector6d bodyState = Eigen::Vector6d::Zero( );
        bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                    "Earth", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Earth" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
        bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                    "Mars", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Mars" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
        bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                    "Moon", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Moon" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
        bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                    "Sun", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Sun" ]->setEphemeris( std::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
    }
    else
    {
        bodyMap[ "Earth" ]->setEphemeris( std::make_shared< SpiceEphemeris >( "Earth", "SSB", false, false ) );
        bodyMap[ "Mars" ]->setEphemeris( std::make_shared< SpiceEphemeris >( "Mars", "SSB", false, false ) );
        bodyMap[ "Moon" ]->setEphemeris( std::make_shared< SpiceEphemeris >( "Moon", "SSB", false, false ) );
        bodyMap[ "Sun" ]->setEphemeris( std::make_shared< SpiceEphemeris >( "Sun", "SSB", false, false ) );
    }

    ( bodyMap[ "Sun" ] )->setGravityFieldModel(
                std::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Sun" ) ) );
    ( bodyMap[ "Moon" ] )->setGravityFieldModel(
                std::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Moon" ) ) );
    ( bodyMap[ "Mars" ] )->setGravityFieldModel(
                std::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Mars" ) *
                                                         gravitationalParameterScaling ) );
    ( bodyMap[ "Earth" ] )->setGravityFieldModel(
                std::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Earth" ) *
                                                         gravitationalParameterScaling ) );


    ( bodyMap[ "Earth" ] )->setRotationalEphemeris(
                createRotationModel(
                    std::make_shared< SimpleRotationModelSettings >(
                        "ECLIPJ2000", "IAU_Earth",
                        spice_interface::computeRotationQuaternionBetweenFrames(
                            "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                        initialEphemerisTime, 2.0 * mathematical_constants::PI /
                        physical_constants::JULIAN_DAY ), "Earth" ) );

    std::shared_ptr< ephemerides::RotationalEphemeris > marsRotationModel = createRotationModel(
                std::make_shared< SimpleRotationModelSettings >(
                    "ECLIPJ2000", "IAU_Mars",
                    spice_interface::computeRotationQuaternionBetweenFrames(
                        "ECLIPJ2000", "IAU_Mars", initialEphemerisTime ),
                    initialEphemerisTime, 2.0 * mathematical_constants::PI /
                    ( physical_constants::JULIAN_DAY + 40.0 * 60.0 ) ), "Mars" );

    if( !useConstantRotationalEphemeris )
    {
        bodyMap[ "Mars" ]->setRotationalEphemeris( marsRotationModel );
    }
    else
    {
        bodyMap[ "Mars" ]->setRotationalEphemeris(
                    std::make_shared< ephemerides::ConstantRotationalEphemeris >(
                        marsRotationModel->getRotationStateVector( 0.0 ), "ECLIPJ2000", "IAU_Mars" ) );
    }


    // Define and create ground stations.
    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStationsToCreate;
    groundStationsToCreate[ std::make_pair( "Earth", "Graz" ) ] =
            ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( );
    groundStationsToCreate[ std::make_pair( "Mars", "MSL" ) ] =
            ( Eigen::Vector3d( ) << -2.5E5, 3.2E6, -2.65E4 ).finished( );


    createGroundStations( bodyMap, groundStationsToCreate );


    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    return bodyMap;
}

//! Function to create estimated parameters for general observation partial tests.
std::shared_ptr< EstimatableParameterSet< double > > createEstimatableParameters(
        const NamedBodyMap& bodyMap, const double initialTime,
        const bool useEquivalencePrincipleParameter,
        const bool useRotationalStateAsParameter )
{
    std::vector< std::shared_ptr< EstimatableParameter< double > > > estimatableDoubleParameters;
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatableVectorParameters;

    if( !useRotationalStateAsParameter )
    {
        std::shared_ptr< RotationRate > earthRotationRate = std::make_shared< RotationRate >(
                    std::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                        bodyMap.at( "Earth" )->getRotationalEphemeris( ) ), "Earth" );
        std::shared_ptr< RotationRate > marsRotationRate = std::make_shared< RotationRate >(
                    std::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                        bodyMap.at( "Mars" )->getRotationalEphemeris( ) ), "Mars" );
        std::shared_ptr< ConstantRotationalOrientation > earthRotationOrientation =
                std::make_shared< ConstantRotationalOrientation >(
                    std::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                        bodyMap.at( "Earth" )->getRotationalEphemeris( ) ), "Earth" );
        std::shared_ptr< ConstantRotationalOrientation > marsRotationOrientation =
                std::make_shared< ConstantRotationalOrientation >(
                    std::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                        bodyMap.at( "Mars" )->getRotationalEphemeris( ) ), "Mars" );

        estimatableDoubleParameters.push_back( earthRotationRate );
        estimatableDoubleParameters.push_back( marsRotationRate );

        estimatableVectorParameters.push_back( earthRotationOrientation );
        estimatableVectorParameters.push_back( marsRotationOrientation );

        std::shared_ptr< EstimatableParameter< double > > relativisticParameter;
        if( useEquivalencePrincipleParameter )
        {
            relativisticParameter = std::make_shared< EquivalencePrincipleLpiViolationParameter >( );
        }
        else
        {
            relativisticParameter = std::make_shared< PPNParameterGamma >( );
        }
        estimatableDoubleParameters.push_back( relativisticParameter );
    }

    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatedInitialStateParameters;
    estimatedInitialStateParameters.push_back(
                std::make_shared< InitialTranslationalStateParameter< double > >(
                    "Earth", propagators::getInitialStateOfBody(
                        "Earth", "SSB", bodyMap, initialTime ) ) );
    estimatedInitialStateParameters.push_back(
                std::make_shared< InitialTranslationalStateParameter< double > >(
                    "Mars", propagators::getInitialStateOfBody(
                        "Mars", "SSB", bodyMap, initialTime ) ) );
    if( useRotationalStateAsParameter )
    {
        estimatedInitialStateParameters.push_back(
                    std::make_shared< InitialRotationalStateParameter< double > >(
                        "Mars", propagators::getInitialRotationalStateOfBody(
                            "Mars", "ECLIPJ2000", bodyMap, initialTime ),
                        std::bind( &simulation_setup::Body::getBodyInertiaTensor, bodyMap.at( "Mars" ) ) ) );
    }

    return std::make_shared< EstimatableParameterSet< double > >(
                estimatableDoubleParameters,
                estimatableVectorParameters,
                estimatedInitialStateParameters );
}

//! Function to compute numerical partials w.r.t. constant body states for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyState(
        const std::string& bodyName, const NamedBodyMap& bodyMap, const Eigen::Vector3d& bodyPositionVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction, const double observationTime,
        const int observableSize )
{
    // Calculate numerical partials w.r.t. body state.
    std::shared_ptr< ConstantEphemeris > bodyEphemeris = std::dynamic_pointer_cast< ConstantEphemeris >(
                bodyMap.at( bodyName )->getEphemeris( ) );
    Eigen::Vector6d bodyUnperturbedState = bodyEphemeris->getCartesianState( 0.0 );
    Eigen::Vector6d perturbedBodyState;

    Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalPartialWrtBodyPosition =
            Eigen::Matrix< double, Eigen::Dynamic, 3 >::Zero( observableSize, 3 );
    for( int i = 0; i < 3; i++ )
    {
        perturbedBodyState = bodyUnperturbedState;
        perturbedBodyState( i ) += bodyPositionVariation( i );
        bodyEphemeris->updateConstantState( perturbedBodyState );
        bodyMap.at( bodyName )->recomputeStateOnNextCall( );
        Eigen::VectorXd upPerturbedObservation = observationFunction( observationTime );

        perturbedBodyState = bodyUnperturbedState;
        perturbedBodyState( i ) -= bodyPositionVariation( i );
        bodyEphemeris->updateConstantState( perturbedBodyState );
        bodyMap.at( bodyName )->recomputeStateOnNextCall( );
        Eigen::VectorXd downPerturbedObservation = observationFunction( observationTime );

        numericalPartialWrtBodyPosition.block( 0, i, observableSize, 1  ) = ( upPerturbedObservation - downPerturbedObservation ) /
                ( 2.0 * bodyPositionVariation( i ) );
    }
    bodyEphemeris->updateConstantState( bodyUnperturbedState );
    bodyMap.at( bodyName )->recomputeStateOnNextCall( );

    return numericalPartialWrtBodyPosition;
}

//! Function to compute numerical partials w.r.t. constant body orientation for general observation partial tests.
Eigen::MatrixXd calculateChangeDueToConstantBodyOrientation(
        const std::string& bodyName, const NamedBodyMap& bodyMap, const Eigen::Vector4d& bodyQuaternionVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction, const double observationTime,
        const int observableSize,
        std::vector< Eigen::Vector4d >& appliedQuaternionPerturbation )
{

    Eigen::VectorXd nominalObservation = observationFunction( observationTime );
    appliedQuaternionPerturbation.resize( 4 );
    for( int i = 0; i < 4; i++ )
    {
        appliedQuaternionPerturbation[ i ].setZero( );
    }

    std::shared_ptr< ephemerides::ConstantRotationalEphemeris > rotationalEphemeris =
            std::dynamic_pointer_cast< ephemerides::ConstantRotationalEphemeris >(
                bodyMap.at( bodyName )->getRotationalEphemeris( ) );

    Eigen::Vector7d bodyUnperturbedState = rotationalEphemeris->getRotationStateVector( observationTime );
    Eigen::Vector7d perturbedBodyState;

    Eigen::Matrix< double, Eigen::Dynamic, 4 > numericalChangeDueToQuaternionChange =
            Eigen::Matrix< double, Eigen::Dynamic, 4 >::Zero( observableSize, 4 );
    for( int i = 0; i < 4; i++ )
    {
        perturbedBodyState = bodyUnperturbedState;
        perturbedBodyState( i ) += bodyQuaternionVariation( i );

        perturbedBodyState( 0 ) = ( bodyUnperturbedState( 0 ) > 0 ? 1.0 : -1.0 ) *
                std::sqrt( 1.0 - std::pow( perturbedBodyState.segment( 1, 3 ).norm( ), 2 ) );
        appliedQuaternionPerturbation[ i ] = perturbedBodyState.segment( 0, 4 ).normalized( ) -
                bodyUnperturbedState.segment( 0, 4 );

        rotationalEphemeris->updateConstantState( perturbedBodyState );
        bodyMap.at( bodyName )->recomputeStateOnNextCall( );
        Eigen::VectorXd upPerturbedObservation = observationFunction( observationTime );

        upPerturbedObservation -= nominalObservation;

        numericalChangeDueToQuaternionChange.block( 0, i, observableSize, 1  ) = upPerturbedObservation;

    }

    rotationalEphemeris->updateConstantState( bodyUnperturbedState );

    bodyMap.at( bodyName )->recomputeStateOnNextCall( );

    return numericalChangeDueToQuaternionChange;
}

//! Function to compute numerical partials w.r.t. constant body angular velocity for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyAngularVelocityVector(
        const std::string& bodyName, const NamedBodyMap& bodyMap, const Eigen::Vector3d& bodyRotationVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction, const double observationTime,
        const int observableSize )
{
    // Calculate numerical partials w.r.t. body state.
    Eigen::Vector7d bodyUnperturbedState = bodyMap.at( bodyName )->getRotationalStateVector( );
    Eigen::Vector7d perturbedBodyState;

    Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalPartialWrtBodyAngularVelocity =
            Eigen::Matrix< double, Eigen::Dynamic, 3 >::Zero( observableSize, 3 );
    for( int i = 0; i < 3; i++ )
    {
        perturbedBodyState = bodyUnperturbedState;
        perturbedBodyState( i + 4 ) += bodyRotationVariation( i );
        bodyMap.at( bodyName )->setCurrentRotationalStateToLocalFrame( perturbedBodyState );
        bodyMap.at( bodyName )->recomputeStateOnNextCall( );
        Eigen::VectorXd upPerturbedObservation = observationFunction( observationTime );

        perturbedBodyState = bodyUnperturbedState;
        perturbedBodyState( i + 4 ) -= bodyRotationVariation( i );
        bodyMap.at( bodyName )->setCurrentRotationalStateToLocalFrame( perturbedBodyState );
        bodyMap.at( bodyName )->recomputeStateOnNextCall( );
        Eigen::VectorXd downPerturbedObservation = observationFunction( observationTime );

        numericalPartialWrtBodyAngularVelocity.block( 0, i, observableSize, 1  ) = ( upPerturbedObservation - downPerturbedObservation ) /
                ( 2.0 * bodyRotationVariation( i ) );
    }
    bodyMap.at( bodyName )->setCurrentRotationalStateToLocalFrame( bodyUnperturbedState );
    bodyMap.at( bodyName )->recomputeStateOnNextCall( );

    return numericalPartialWrtBodyAngularVelocity;
}

//! Function to compute numerical partials w.r.t. constant body states for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyVelocity(
        const std::string& bodyName, const NamedBodyMap& bodyMap, const Eigen::Vector3d& bodyVelocityVariation,
        const std::function< Eigen::VectorXd( const double ) > observationFunction, const double observationTime,
        const int observableSize )
{
    // Calculate numerical partials w.r.t. body state.
    std::shared_ptr< ConstantEphemeris > bodyEphemeris = std::dynamic_pointer_cast< ConstantEphemeris >(
                bodyMap.at( bodyName )->getEphemeris( ) );
    Eigen::Vector6d bodyUnperturbedState = bodyEphemeris->getCartesianState( 0.0 );
    Eigen::Vector6d perturbedBodyState;

    Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalPartialWrtBodyPosition =
            Eigen::Matrix< double, Eigen::Dynamic, 3 >::Zero( observableSize, 3 );
    for( int i = 0; i < 3; i++ )
    {
        perturbedBodyState = bodyUnperturbedState;
        perturbedBodyState( i + 3 ) += bodyVelocityVariation( i );
        bodyEphemeris->updateConstantState( perturbedBodyState );
        bodyMap.at( bodyName )->recomputeStateOnNextCall( );
        Eigen::VectorXd upPerturbedObservation = observationFunction( observationTime );

        perturbedBodyState = bodyUnperturbedState;
        perturbedBodyState( i + 3 ) -= bodyVelocityVariation( i );
        bodyEphemeris->updateConstantState( perturbedBodyState );
        bodyMap.at( bodyName )->recomputeStateOnNextCall( );
        Eigen::VectorXd downPerturbedObservation = observationFunction( observationTime );

        numericalPartialWrtBodyPosition.block( 0, i, observableSize, 1  ) = ( upPerturbedObservation - downPerturbedObservation ) /
                ( 2.0 * bodyVelocityVariation( i ) );
    }
    bodyEphemeris->updateConstantState( bodyUnperturbedState );
    bodyMap.at( bodyName )->recomputeStateOnNextCall( );

    return numericalPartialWrtBodyPosition;
}


//! Function to compute numerical partials w.r.t. double parameters for general observation partial tests.
std::vector< Eigen::VectorXd > calculateNumericalPartialsWrtDoubleParameters(
        const std::vector< std::shared_ptr< EstimatableParameter< double > > >& doubleParameters,
        const std::vector< std::function< void( ) > > updateFunctions,
        const std::vector< double >& parameterPerturbations,
        const std::function< Eigen::VectorXd( const double ) > observationFunction,
        const double observationTime )
{
    std::vector< Eigen::VectorXd > partialVector;

    for( unsigned int i = 0; i < doubleParameters.size( ); i++ )
    {
        partialVector.push_back( calculateNumericalObservationParameterPartial(
                                     doubleParameters.at( i ), parameterPerturbations.at( i ), observationFunction,
                                     observationTime, updateFunctions.at( i ) ) );
    }

    return partialVector;
}

std::vector< Eigen::MatrixXd > calculateNumericalPartialsWrtVectorParameters(
        const std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >& vectorParameters,
        const std::vector< std::function< void( ) > > updateFunctions,
        const std::vector< Eigen::VectorXd >& parameterPerturbations,
        const std::function< Eigen::VectorXd( const double ) > observationFunction,
        const double observationTime )
{
    std::vector< Eigen::MatrixXd > partialVector;

    for( unsigned int i = 0; i < vectorParameters.size( ); i++ )
    {
        partialVector.push_back( calculateNumericalObservationParameterPartial(
                                     vectorParameters.at( i ), parameterPerturbations.at( i ), observationFunction,
                                     observationTime, updateFunctions.at( i ) ) );
    }

    return partialVector;
}

//! Function to retrieve times associated with analytical partial derivatives for general observation partial tests.
std::vector< std::vector< double > > getAnalyticalPartialEvaluationTimes(
        const LinkEnds& linkEnds,
        const ObservableType observableType,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< EstimatableParameterSet< double > >& estimatedParameters )
{
    std::vector< std::vector< double > > partialTimes;
    std::vector< double > currentPartialTimes;
    std::vector< int > currentPartialTimeIndices;

    std::vector< std::string > bodiesWithTranslationalState =
            estimatable_parameters::getListOfBodiesToEstimate(
                estimatedParameters ).at( propagators::translational_state );

    for( unsigned int i = 0; i < bodiesWithTranslationalState.size( ); i++ )
    {
        currentPartialTimes.clear( );
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( );
             linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second.first == bodiesWithTranslationalState.at( i ) )
            {
                currentPartialTimeIndices =
                        getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first, linkEnds.size( ) );

                for( unsigned int j = 0; j < currentPartialTimeIndices.size( ); j++ )
                {
                    currentPartialTimes.push_back( linkEndTimes.at( currentPartialTimeIndices.at( j ) ) );
                }
            }
        }
        partialTimes.push_back( currentPartialTimes );
    }

    std::vector< std::string > bodiesWithRotationalState;
    if( estimatable_parameters::getListOfBodiesToEstimate( estimatedParameters ).count( propagators::rotational_state ) > 0 )
    {
        bodiesWithRotationalState = estimatable_parameters::getListOfBodiesToEstimate(
                    estimatedParameters ).at( propagators::rotational_state );
    }

    for( unsigned int i = 0; i < bodiesWithRotationalState.size( ); i++ )
    {
        currentPartialTimes.clear( );
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( );
             linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second.first == bodiesWithRotationalState.at( i ) )
            {
                currentPartialTimeIndices =
                        getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first, linkEnds.size( ) );

                for( unsigned int j = 0; j < currentPartialTimeIndices.size( ); j++ )
                {
                    currentPartialTimes.push_back( linkEndTimes.at( currentPartialTimeIndices.at( j ) ) );
                }
            }
        }
        partialTimes.push_back( currentPartialTimes );
    }


    bool checkStationId;
    bool addContribution;
    for( unsigned int i = 0; i < estimatedParameters->getEstimatedDoubleParameters( ).size( ); i++ )
    {
        checkStationId = 0;
        std::pair< std::string, std::string > currentAssociatedLinkEndId =
                estimatedParameters->getEstimatedDoubleParameters( ).at( i )->getParameterName( ).second;
        if( currentAssociatedLinkEndId.second != "" )
        {
            checkStationId = 1;
        }
        currentPartialTimes.clear( );
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( );
             linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second.first == currentAssociatedLinkEndId.first )
            {
                addContribution = 0;
                if( checkStationId )
                {
                    if( linkEndIterator->second.second == currentAssociatedLinkEndId.second )
                    {
                        currentPartialTimeIndices =
                                getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first, linkEnds.size( ) );
                        addContribution = 1;
                    }
                }
                else
                {
                    currentPartialTimeIndices =
                            getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first, linkEnds.size( ) );
                    addContribution = 1;
                }
                if( addContribution )
                {
                    for( unsigned int j = 0; j < currentPartialTimeIndices.size( ); j++ )
                    {
                        currentPartialTimes.push_back( linkEndTimes.at( currentPartialTimeIndices.at( j ) ) );
                    }
                }
            }
        }
        partialTimes.push_back( currentPartialTimes );
    }

    for( unsigned int i = 0; i < estimatedParameters->getEstimatedVectorParameters( ).size( ); i++ )
    {
        checkStationId = 0;
        std::pair< std::string, std::string > currentAssociatedLinkEndId =
                estimatedParameters->getEstimatedVectorParameters( ).at( i )->getParameterName( ).second;
        if( currentAssociatedLinkEndId.second != "" )
        {
            checkStationId = 1;
        }
        currentPartialTimes.clear( );
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( );
             linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second.first == currentAssociatedLinkEndId.first )
            {
                addContribution = 0;
                if( checkStationId )
                {
                    if( linkEndIterator->second.second == currentAssociatedLinkEndId.second )
                    {
                        currentPartialTimeIndices =
                                getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first, linkEnds.size( ) );
                        addContribution = 1;
                    }
                }
                else
                {
                    currentPartialTimeIndices =
                            getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first, linkEnds.size( ) );
                    addContribution = 1;
                }
                if( addContribution )
                {
                    for( unsigned int j = 0; j < currentPartialTimeIndices.size( ); j++ )
                    {
                        currentPartialTimes.push_back( linkEndTimes.at( currentPartialTimeIndices.at( j ) ) );
                    }
                }
            }
        }
        partialTimes.push_back( currentPartialTimes );
    }

    return partialTimes;
}

template void testObservationPartials< 1 >(
        const std::shared_ptr< ObservationModel< 1, double, double > > observationModel,
        NamedBodyMap& bodyMap,
        const std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet,
        const LinkEnds& linkEnds, const ObservableType observableType,
        const double tolerance,
        const bool testPositionPartial,
        const bool testParameterPartial,
        const double positionPerturbationMultiplier,
        const Eigen::VectorXd parameterPerturbationMultipliers );

template void testObservationPartials< 2 >(
        const std::shared_ptr< ObservationModel< 2, double, double > > observationModel,
        NamedBodyMap& bodyMap,
        const std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet,
        const LinkEnds& linkEnds, const ObservableType observableType,
        const double tolerance,
        const bool testPositionPartial,
        const bool testParameterPartial,
        const double positionPerturbationMultiplier,
        const Eigen::VectorXd parameterPerturbationMultipliers );


}

}
