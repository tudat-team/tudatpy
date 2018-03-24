/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

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
                               const double gravitationalParameterScaling )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies.
    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = boost::make_shared< Body >( );
    bodyMap[ "Mars" ] = boost::make_shared< Body >( );
    bodyMap[ "Moon" ] = boost::make_shared< Body >( );
    bodyMap[ "Sun" ] = boost::make_shared< Body >( );

    bodyMap[ "Earth" ]->setShapeModel( boost::make_shared< basic_astrodynamics::SphericalBodyShapeModel >(
                                           spice_interface::getAverageRadius( "Earth" ) ) );
    bodyMap[ "Mars" ]->setShapeModel( boost::make_shared< basic_astrodynamics::SphericalBodyShapeModel >(
                                          spice_interface::getAverageRadius( "Mars" ) ) );

    if( useConstantEphemerides )
    {
        Eigen::Vector6d bodyState = Eigen::Vector6d::Zero( );
        bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                    "Earth", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Earth" ]->setEphemeris( boost::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
        bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                    "Mars", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Mars" ]->setEphemeris( boost::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
        bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                    "Moon", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Moon" ]->setEphemeris( boost::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
        bodyState.segment( 0, 3 ) = getBodyCartesianPositionAtEpoch(
                    "Sun", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Sun" ]->setEphemeris( boost::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
    }
    else
    {
        bodyMap[ "Earth" ]->setEphemeris( boost::make_shared< SpiceEphemeris >( "Earth", "SSB", false, false ) );
        bodyMap[ "Mars" ]->setEphemeris( boost::make_shared< SpiceEphemeris >( "Mars", "SSB", false, false ) );
        bodyMap[ "Moon" ]->setEphemeris( boost::make_shared< SpiceEphemeris >( "Moon", "SSB", false, false ) );
        bodyMap[ "Sun" ]->setEphemeris( boost::make_shared< SpiceEphemeris >( "Sun", "SSB", false, false ) );
    }

    ( bodyMap[ "Sun" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Sun" ) ) );
    ( bodyMap[ "Moon" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Moon" ) ) );
    ( bodyMap[ "Mars" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Mars" ) *
                                                         gravitationalParameterScaling ) );
    ( bodyMap[ "Earth" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Earth" ) *
                                                         gravitationalParameterScaling ) );


    ( bodyMap[ "Earth" ] )->setRotationalEphemeris(
                createRotationModel(
                    boost::make_shared< SimpleRotationModelSettings >(
                        "ECLIPJ2000", "IAU_Earth",
                        spice_interface::computeRotationQuaternionBetweenFrames(
                            "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                        initialEphemerisTime, 2.0 * mathematical_constants::PI /
                        physical_constants::JULIAN_DAY ), "Earth" ) );
    ( bodyMap[ "Mars" ] )->setRotationalEphemeris(
                createRotationModel(
                    boost::make_shared< SimpleRotationModelSettings >(
                        "ECLIPJ2000", "IAU_Mars",
                        spice_interface::computeRotationQuaternionBetweenFrames(
                            "ECLIPJ2000", "IAU_Mars", initialEphemerisTime ),
                        initialEphemerisTime, 2.0 * mathematical_constants::PI /
                        ( physical_constants::JULIAN_DAY + 40.0 * 60.0 ) ),
                    "Mars" ) );


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
boost::shared_ptr< EstimatableParameterSet< double > > createEstimatableParameters(
        const NamedBodyMap& bodyMap, const double initialTime,
        const bool useEquivalencePrincipleParameter )
{
    boost::shared_ptr< RotationRate > earthRotationRate = boost::make_shared< RotationRate >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap.at( "Earth" )->getRotationalEphemeris( ) ), "Earth" );
    boost::shared_ptr< RotationRate > marsRotationRate = boost::make_shared< RotationRate >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap.at( "Mars" )->getRotationalEphemeris( ) ), "Mars" );
    boost::shared_ptr< ConstantRotationalOrientation > earthRotationOrientation =
            boost::make_shared< ConstantRotationalOrientation >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap.at( "Earth" )->getRotationalEphemeris( ) ), "Earth" );
    boost::shared_ptr< ConstantRotationalOrientation > marsRotationOrientation =
            boost::make_shared< ConstantRotationalOrientation >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap.at( "Mars" )->getRotationalEphemeris( ) ), "Mars" );
    boost::shared_ptr< EstimatableParameter< double > > relativisticParameter;
    if( useEquivalencePrincipleParameter )
    {
        relativisticParameter = boost::make_shared< EquivalencePrincipleLpiViolationParameter >( );
    }
    else
    {
        relativisticParameter = boost::make_shared< PPNParameterGamma >( );
    }




    std::vector< boost::shared_ptr< EstimatableParameter< double > > > estimatableDoubleParameters;
    estimatableDoubleParameters.push_back( earthRotationRate );
    estimatableDoubleParameters.push_back( marsRotationRate );
    estimatableDoubleParameters.push_back( relativisticParameter );

    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatableVectorParameters;
    estimatableVectorParameters.push_back( earthRotationOrientation );
    estimatableVectorParameters.push_back( marsRotationOrientation );


    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatedInitialStateParameters;
    estimatedInitialStateParameters.push_back(
                boost::make_shared< InitialTranslationalStateParameter< double > >(
                    "Earth", propagators::getInitialStateOfBody(
                        "Earth", "SSB", bodyMap, initialTime ) ) );
    estimatedInitialStateParameters.push_back(
                boost::make_shared< InitialTranslationalStateParameter< double > >(
                    "Mars", propagators::getInitialStateOfBody(
                        "Mars", "SSB", bodyMap, initialTime ) ) );

    return boost::make_shared< EstimatableParameterSet< double > >(
                estimatableDoubleParameters,
                estimatableVectorParameters,
                estimatedInitialStateParameters );
}

//! Function to compute numerical partials w.r.t. constant body states for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyState(
        const std::string& bodyName, const NamedBodyMap& bodyMap, const Eigen::Vector3d& bodyPositionVariation,
        const boost::function< Eigen::VectorXd( const double ) > observationFunction, const double observationTime,
        const int observableSize )
{
    // Calculate numerical partials w.r.t. body state.
    boost::shared_ptr< ConstantEphemeris > bodyEphemeris = boost::dynamic_pointer_cast< ConstantEphemeris >(
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

//! Function to compute numerical partials w.r.t. constant body states for general observation partial tests.
Eigen::Matrix< double, Eigen::Dynamic, 3 > calculatePartialWrtConstantBodyVelocity(
        const std::string& bodyName, const NamedBodyMap& bodyMap, const Eigen::Vector3d& bodyVelocityVariation,
        const boost::function< Eigen::VectorXd( const double ) > observationFunction, const double observationTime,
        const int observableSize )
{
    // Calculate numerical partials w.r.t. body state.
    boost::shared_ptr< ConstantEphemeris > bodyEphemeris = boost::dynamic_pointer_cast< ConstantEphemeris >(
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
        const std::vector< boost::shared_ptr< EstimatableParameter< double > > >& doubleParameters,
        const std::vector< boost::function< void( ) > > updateFunctions,
        const std::vector< double >& parameterPerturbations,
        const boost::function< Eigen::VectorXd( const double ) > observationFunction,
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
        const std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >& vectorParameters,
        const std::vector< boost::function< void( ) > > updateFunctions,
        const std::vector< Eigen::VectorXd >& parameterPerturbations,
        const boost::function< Eigen::VectorXd( const double ) > observationFunction,
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
        const boost::shared_ptr< EstimatableParameterSet< double > >& estimatedParameters )
{
    std::vector< std::vector< double > > partialTimes;
    std::vector< double > currentPartialTimes;
    std::vector< int > currentPartialTimeIndices;

    std::vector< std::string > bodiesWithEstimatedState = estimatable_parameters::getListOfBodiesToEstimate(
                estimatedParameters );

    for( unsigned int i = 0; i < bodiesWithEstimatedState.size( ); i++ )
    {
        currentPartialTimes.clear( );
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( );
             linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second.first == bodiesWithEstimatedState.at( i ) )
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

}

}
