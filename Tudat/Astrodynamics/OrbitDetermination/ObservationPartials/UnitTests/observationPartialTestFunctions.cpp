#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/observationPartialTestFunctions.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
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


NamedBodyMap setupEnvironment( const std::vector< LinkEndId > groundStations,
                               const double initialEphemerisTime,
                               const double finalEphemerisTime,
                               const double stateEvaluationTime,
                               const bool useConstantEphemerides )
{
    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");

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
        basic_mathematics::Vector6d bodyState = basic_mathematics::Vector6d::Zero( );
        bodyState.segment( 0, 6 ) = getBodyCartesianPositionAtEpoch( "Earth", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Earth" ]->setEphemeris( boost::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
        bodyState.segment( 0, 6 ) = getBodyCartesianPositionAtEpoch( "Mars", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Mars" ]->setEphemeris( boost::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
        bodyState.segment( 0, 6 ) = getBodyCartesianPositionAtEpoch( "Moon", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
        bodyMap[ "Moon" ]->setEphemeris( boost::make_shared< ConstantEphemeris >( bodyState, "SSB", "ECLIPJ2000" ) );
        bodyState.segment( 0, 6 ) = getBodyCartesianPositionAtEpoch( "Sun", "SSB", "ECLIPJ2000", "NONE", stateEvaluationTime );
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
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Mars" ) ) );
    ( bodyMap[ "Earth" ] )->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Earth" ) ) );


    ( bodyMap[ "Earth" ] )->setRotationalEphemeris(
                createRotationModel(
                    boost::make_shared< SimpleRotationModelSettings >(
                        "ECLIPJ2000", "IAU_Earth",
                        spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                        initialEphemerisTime, 2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY ), "Earth" ) );
    ( bodyMap[ "Mars" ] )->setRotationalEphemeris(
                createRotationModel(
                    boost::make_shared< SimpleRotationModelSettings >(
                        "ECLIPJ2000", "IAU_Mars",
                        spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Mars", initialEphemerisTime ),
                        initialEphemerisTime, 2.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY + 40.0 * 60.0 ) ),
                    "Mars" ) );


    // Define and create ground stations.
    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStationsToCreate;
    groundStationsToCreate[ std::make_pair( "Earth", "Graz" ) ] =
            ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( );
    groundStationsToCreate[ std::make_pair( "Mars", "MSL" ) ] =
            ( Eigen::Vector3d( ) <<-2.5E5, 3.2E6, -2.65E4 ).finished( );


    createGroundStations( bodyMap, groundStationsToCreate );


    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    return bodyMap;
}

boost::shared_ptr< EstimatableParameterSet< double > > createEstimatableParameters(
        const NamedBodyMap& bodyMap, const double initialTime )
{
    boost::shared_ptr< RotationRate > earthRotationRate = boost::make_shared< RotationRate >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap.at( "Earth" )->getRotationalEphemeris( ) ), "Earth" );
    boost::shared_ptr< RotationRate > marsRotationRate = boost::make_shared< RotationRate >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap.at( "Mars" )->getRotationalEphemeris( ) ), "Mars" );
    boost::shared_ptr< ConstantRotationalOrientation > earthRotationOrientation = boost::make_shared< ConstantRotationalOrientation >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap.at( "Earth" )->getRotationalEphemeris( ) ), "Earth" );
    boost::shared_ptr< ConstantRotationalOrientation > marsRotationOrientation = boost::make_shared< ConstantRotationalOrientation >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap.at( "Mars" )->getRotationalEphemeris( ) ), "Mars" );


    std::vector< boost::shared_ptr< EstimatableParameter< double > > > estimatableDoubleParameters;
    estimatableDoubleParameters.push_back( earthRotationRate );
    estimatableDoubleParameters.push_back( marsRotationRate );

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

Eigen::Matrix< double, 1, 3 > calculatePartialWrtConstantBodyState(
        const std::string& bodyName, const NamedBodyMap& bodyMap, const Eigen::Vector3d& bodyPositionVariation,
        const boost::function< double( const double ) > observationFunction, const double observationTime )
{
    // Calculate numerical partials w.r.t. body state.
    boost::shared_ptr< ConstantEphemeris > bodyEphemeris = boost::dynamic_pointer_cast< ConstantEphemeris >(
                bodyMap.at( bodyName )->getEphemeris( ) );
    basic_mathematics::Vector6d bodyUnperturbedState = bodyEphemeris->getCartesianState( 0.0 );
    basic_mathematics::Vector6d perturbedBodyState;

    Eigen::Matrix< double, 1, 3 > numericalPartialWrtBodyPosition = Eigen::Matrix< double, 1, 3 >::Zero( );
    for( int i = 0; i < 3; i++ )
    {
        perturbedBodyState = bodyUnperturbedState;
        perturbedBodyState( i ) += bodyPositionVariation( i );
        bodyEphemeris->updateConstantState( perturbedBodyState );
        bodyMap.at( bodyName )->recomputeStateOnNextCall( );
        double upPerturbedObservation = observationFunction( observationTime );

        perturbedBodyState = bodyUnperturbedState;
        perturbedBodyState( i ) -= bodyPositionVariation( i );
        bodyEphemeris->updateConstantState( perturbedBodyState );
        bodyMap.at( bodyName )->recomputeStateOnNextCall( );
        double downPerturbedObservation = observationFunction( observationTime );

        numericalPartialWrtBodyPosition( 0, i ) = ( upPerturbedObservation - downPerturbedObservation ) /
                ( 2.0 * bodyPositionVariation( i ) );
    }
    bodyEphemeris->updateConstantState( bodyUnperturbedState );
    bodyMap.at( bodyName )->recomputeStateOnNextCall( );

    return numericalPartialWrtBodyPosition;
}

std::vector< double > calculateNumericalPartialsWrtDoubleParameters(
        const std::vector< boost::shared_ptr< EstimatableParameter< double > > >& doubleParameters,
        const std::vector< boost::function< void( ) > > updateFunctions,
        const std::vector< double >& parameterPerturbations,
        const boost::function< double( const double ) > observationFunction,
        const double observationTime )
{
    std::vector< double > partialVector;

    for( unsigned int i = 0; i < doubleParameters.size( ); i++ )
    {
        partialVector.push_back( calculateNumericalObservationParameterPartial(
                                     doubleParameters.at( i ), parameterPerturbations.at( i ), observationFunction,
                                     observationTime, updateFunctions.at( i ) )( 0 ) );
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

std::vector< std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > > calculateAnalyticalPartials(
        const std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > >& partialObjectList,
        const std::vector< basic_mathematics::Vector6d >& states,
        const std::vector< double >& times,
        const LinkEndType linkEndOfFixedTime )
{
    std::vector< std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > > partialList;

    for( std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > >::const_iterator partialIterator =
         partialObjectList.begin( ); partialIterator != partialObjectList.end( ); partialIterator++ )
    {
        partialList.push_back( partialIterator->second->calculatePartial( states, times, linkEndOfFixedTime ) );
    }
    return partialList;
}

std::vector< int > getSingleAnalyticalPartialSize(
        const LinkEnds& linkEnds,
        const std::pair< std::vector< std::string >, boost::shared_ptr< EstimatableParameterSet< double > > >& estimatedParameters )
{
    std::vector< int > partialSize;
    int currentSize;
    for( unsigned int i = 0; i < estimatedParameters.first.size( ); i++ )
    {
        currentSize = 0;
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second.first == estimatedParameters.first.at( i ) )
            {
                currentSize++;
            }
        }
        partialSize.push_back( currentSize );
    }

    bool checkStationId;
    for( unsigned int i = 0; i < estimatedParameters.second->getEstimatedDoubleParameters( ).size( ); i++ )
    {
        currentSize = 0;
        checkStationId = 0;
        std::pair< std::string, std::string > currentAssociatedLinkEndId =
                estimatedParameters.second->getEstimatedDoubleParameters( ).at( i )->getParameterName( ).second;
        if( currentAssociatedLinkEndId.second != "" )
        {
            checkStationId = 1;
        }

        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second.first == currentAssociatedLinkEndId.first )
            {
                if( checkStationId )
                {
                    if( linkEndIterator->second.second == currentAssociatedLinkEndId.second )
                    {
                        currentSize++;
                    }
                }
                else
                {
                    currentSize++;
                }
            }
        }
        partialSize.push_back( currentSize );
    }

    return partialSize;
}

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
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second.first == bodiesWithEstimatedState.at( i ) )
            {
                currentPartialTimeIndices = getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first );
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
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second.first == currentAssociatedLinkEndId.first )
            {
                addContribution = 0;
                if( checkStationId )
                {
                    if( linkEndIterator->second.second == currentAssociatedLinkEndId.second )
                    {
                        currentPartialTimeIndices = getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first );
                        addContribution = 1;
                    }
                }
                else
                {
                    currentPartialTimeIndices = getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first );
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
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second.first == currentAssociatedLinkEndId.first )
            {
                addContribution = 0;
                if( checkStationId )
                {
                    if( linkEndIterator->second.second == currentAssociatedLinkEndId.second )
                    {
                        currentPartialTimeIndices = getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first );
                        addContribution = 1;
                    }
                }
                else
                {
                    currentPartialTimeIndices = getLinkEndIndicesForLinkEndTypeAtObservable( observableType, linkEndIterator->first );
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
