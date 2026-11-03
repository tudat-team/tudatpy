/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/createInverseAprioriCovariance.h"
#include "tudat/simulation/estimation_setup/executePlanetaryParameterEstimationTestCase.h"
#include "tudat/simulation/estimation_setup/executeEarthOrbiterParameterEstimationTestCase.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"

namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_estimation_input_output )

//! Test utility functions for constructing and updating inverse a priori covariance matrices.
BOOST_AUTO_TEST_CASE( test_InverseAprioriCovarianceUtilities )
{
    using estimatable_parameters::EstimatebleParameterIdentifier;

    // Create two scalar gravitational parameters and one custom vector parameter.
    std::shared_ptr< gravitation::GravityFieldModel > earthGravityField =
            std::make_shared< gravitation::GravityFieldModel >( 3.986004418E14 );
    std::shared_ptr< gravitation::GravityFieldModel > marsGravityField = std::make_shared< gravitation::GravityFieldModel >( 4.282837E13 );

    std::shared_ptr< estimatable_parameters::GravitationalParameter > earthMuParameter =
            std::make_shared< estimatable_parameters::GravitationalParameter >( earthGravityField, "Earth" );
    std::shared_ptr< estimatable_parameters::GravitationalParameter > marsMuParameter =
            std::make_shared< estimatable_parameters::GravitationalParameter >( marsGravityField, "Mars" );

    Eigen::VectorXd customParameterValue = ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( );
    std::shared_ptr< estimatable_parameters::CustomEstimatableParameter > customVectorParameter =
            std::make_shared< estimatable_parameters::CustomEstimatableParameter >(
                    "CustomA",
                    3,
                    [ &customParameterValue ]( ) { return customParameterValue; },
                    [ &customParameterValue ]( const Eigen::VectorXd& newParameterValue ) { customParameterValue = newParameterValue; } );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > scalarParameters = { earthMuParameter,
                                                                                                                  marsMuParameter };
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParameters = {
        customVectorParameter
    };

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
            std::make_shared< estimatable_parameters::EstimatableParameterSet< double > >( scalarParameters, vectorParameters );

    BOOST_CHECK_EQUAL( parameterSet->getEstimatedParameterSetSize( ), 5 );

    // Test enum-only lookup for a unique enum.
    const EstimatebleParameterIdentifier customEnumOnlyIdentifier =
            std::make_pair( estimatable_parameters::custom_estimated_parameter, std::make_pair( "", "" ) );
    std::vector< std::pair< int, int > > customIndices = parameterSet->getIndicesForParameterType( customEnumOnlyIdentifier );
    BOOST_CHECK_EQUAL( customIndices.size( ), 1 );
    BOOST_CHECK_EQUAL( customIndices.at( 0 ).first, 2 );
    BOOST_CHECK_EQUAL( customIndices.at( 0 ).second, 3 );

    std::vector< std::pair< std::pair< int, int >, std::shared_ptr< estimatable_parameters::EstimatableParameterBase > > > customEntries =
            parameterSet->getParametersAndIndicesForParameterIdentifier( customEnumOnlyIdentifier );
    BOOST_CHECK_EQUAL( customEntries.size( ), 1 );
    BOOST_CHECK( customEntries.at( 0 ).second != nullptr );
    BOOST_CHECK_EQUAL( customEntries.at( 0 ).first.first, customIndices.at( 0 ).first );
    BOOST_CHECK_EQUAL( customEntries.at( 0 ).first.second, customIndices.at( 0 ).second );
    BOOST_CHECK_CLOSE_FRACTION( customEntries.at( 0 ).second->getParameterValueBase( ).norm( ), customParameterValue.norm( ), 1.0E-15 );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterBase > > customParameters =
            parameterSet->getParametersForParameterIdentifier( customEnumOnlyIdentifier );
    BOOST_CHECK_EQUAL( customParameters.size( ), 1 );
    BOOST_CHECK( customParameters.at( 0 ) != nullptr );
    customParameters.at( 0 )->setParameterValueBase( ( Eigen::Vector3d( ) << 4.0, 5.0, 6.0 ).finished( ) );
    BOOST_CHECK_SMALL( ( customParameterValue - ( Eigen::Vector3d( ) << 4.0, 5.0, 6.0 ).finished( ) ).norm( ), 1.0E-15 );

    const EstimatebleParameterIdentifier earthMuIdentifier =
            std::make_pair( estimatable_parameters::gravitational_parameter, std::make_pair( "Earth", "" ) );
    const EstimatebleParameterIdentifier marsMuIdentifier =
            std::make_pair( estimatable_parameters::gravitational_parameter, std::make_pair( "Mars", "" ) );

    // Test matrix creation from inverse a priori covariance entries.
    std::vector< std::pair< EstimatebleParameterIdentifier, Eigen::VectorXd > > creationEntries;
    creationEntries.push_back( std::make_pair( earthMuIdentifier, Eigen::VectorXd::Constant( 1, 1.0 / 4.0 ) ) );
    creationEntries.push_back(
            std::make_pair( customEnumOnlyIdentifier, ( Eigen::Vector3d( ) << 1.0 / 100.0, 1.0 / 400.0, 1.0 / 1600.0 ).finished( ) ) );

    Eigen::MatrixXd createdInverseAprioriCovariance =
            simulation_setup::createCovarianceFromDiagonalEntries< double >( parameterSet, creationEntries );

    BOOST_CHECK_EQUAL( createdInverseAprioriCovariance.rows( ), 5 );
    BOOST_CHECK_EQUAL( createdInverseAprioriCovariance.cols( ), 5 );
    BOOST_CHECK_CLOSE_FRACTION( createdInverseAprioriCovariance( 0, 0 ), 1.0 / 4.0, 1.0E-15 );
    BOOST_CHECK_EQUAL( createdInverseAprioriCovariance( 1, 1 ), 0.0 );
    BOOST_CHECK_CLOSE_FRACTION( createdInverseAprioriCovariance( 2, 2 ), 1.0 / 100.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( createdInverseAprioriCovariance( 3, 3 ), 1.0 / 400.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( createdInverseAprioriCovariance( 4, 4 ), 1.0 / 1600.0, 1.0E-15 );
    BOOST_CHECK_EQUAL( createdInverseAprioriCovariance( 0, 2 ), 0.0 );

    // Test update while preserving non-updated entries.
    Eigen::MatrixXd initialInverseAprioriCovariance = 3.0 * Eigen::MatrixXd::Identity( 5, 5 );
    initialInverseAprioriCovariance( 0, 1 ) = 9.0;
    initialInverseAprioriCovariance( 1, 0 ) = 9.0;

    std::vector< std::pair< EstimatebleParameterIdentifier, Eigen::VectorXd > > updateEntries;
    updateEntries.push_back( std::make_pair( marsMuIdentifier, Eigen::VectorXd::Constant( 1, 1.0 / 25.0 ) ) );

    Eigen::MatrixXd updatedInverseAprioriCovariance =
            simulation_setup::addCovarianceDiagonalEntries< double >( initialInverseAprioriCovariance, parameterSet, updateEntries );

    BOOST_CHECK_CLOSE_FRACTION( updatedInverseAprioriCovariance( 1, 1 ), 1.0 / 25.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( updatedInverseAprioriCovariance( 0, 0 ), 3.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( updatedInverseAprioriCovariance( 0, 1 ), 9.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( updatedInverseAprioriCovariance( 1, 0 ), 9.0, 1.0E-15 );

    // Enum-only lookup returns all matching parameters.
    const EstimatebleParameterIdentifier ambiguousGravitationalEnumIdentifier =
            std::make_pair( estimatable_parameters::gravitational_parameter, std::make_pair( "", "" ) );
    std::vector< std::pair< int, int > > ambiguousIndices =
            parameterSet->getIndicesForParameterType( ambiguousGravitationalEnumIdentifier );
    BOOST_CHECK_EQUAL( ambiguousIndices.size( ), 2 );
    BOOST_CHECK_EQUAL( ambiguousIndices.at( 0 ).second, 1 );
    BOOST_CHECK_EQUAL( ambiguousIndices.at( 1 ).second, 1 );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterBase > > ambiguousParameters =
            parameterSet->getParametersForParameterIdentifier( ambiguousGravitationalEnumIdentifier );
    BOOST_CHECK_EQUAL( ambiguousParameters.size( ), 2 );
    BOOST_CHECK( ambiguousParameters.at( 0 ) != nullptr );
    BOOST_CHECK( ambiguousParameters.at( 1 ) != nullptr );
    BOOST_CHECK_EQUAL( ambiguousParameters.at( 0 )->getParameterSize( ), 1 );
    BOOST_CHECK_EQUAL( ambiguousParameters.at( 1 )->getParameterSize( ), 1 );

    // Non-existent identifier returns no matches.
    const EstimatebleParameterIdentifier nonExistentIdentifier =
            std::make_pair( estimatable_parameters::gravitational_parameter, std::make_pair( "Jupiter", "" ) );
    BOOST_CHECK_EQUAL( parameterSet->getIndicesForParameterType( nonExistentIdentifier ).size( ), 0 );

    std::vector< std::pair< EstimatebleParameterIdentifier, Eigen::VectorXd > > nonMatchingEntries;
    nonMatchingEntries.push_back( std::make_pair( nonExistentIdentifier, Eigen::VectorXd::Constant( 1, 1.0 ) ) );

    // By default, unmatched parameter identifiers in covariance-entry input throw.
    BOOST_CHECK_THROW( simulation_setup::createCovarianceFromDiagonalEntries< double >( parameterSet, nonMatchingEntries ),
                       std::runtime_error );
    BOOST_CHECK_THROW(
            simulation_setup::addCovarianceDiagonalEntries< double >( Eigen::MatrixXd::Identity( 5, 5 ), parameterSet, nonMatchingEntries ),
            std::runtime_error );

    // If strict matching is disabled, unmatched identifiers are ignored.
    Eigen::MatrixXd covarianceWithIgnoredUnknown =
            simulation_setup::createCovarianceFromDiagonalEntries< double >( parameterSet, nonMatchingEntries, false );
    BOOST_CHECK_SMALL( covarianceWithIgnoredUnknown.norm( ), 1.0E-15 );

    Eigen::MatrixXd referenceCovariance = Eigen::MatrixXd::Identity( 5, 5 );
    Eigen::MatrixXd unchangedCovariance =
            simulation_setup::addCovarianceDiagonalEntries< double >( referenceCovariance, parameterSet, nonMatchingEntries, false );
    BOOST_CHECK_SMALL( ( unchangedCovariance - referenceCovariance ).norm( ), 1.0E-15 );

    // Enum-only inverse a priori covariance entry applies to all matching parameters.
    std::vector< std::pair< EstimatebleParameterIdentifier, Eigen::VectorXd > > ambiguousEntries;
    ambiguousEntries.push_back( std::make_pair( ambiguousGravitationalEnumIdentifier, Eigen::VectorXd::Constant( 1, 1.0 ) ) );

    Eigen::MatrixXd ambiguousInverseAprioriCovariance =
            simulation_setup::createCovarianceFromDiagonalEntries< double >( parameterSet, ambiguousEntries );
    BOOST_CHECK_CLOSE_FRACTION( ambiguousInverseAprioriCovariance( 0, 0 ), 1.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( ambiguousInverseAprioriCovariance( 1, 1 ), 1.0, 1.0E-15 );
    BOOST_CHECK_EQUAL( ambiguousInverseAprioriCovariance( 2, 2 ), 0.0 );
}

//! Test whether the covariance is correctly computed as a function of time
BOOST_AUTO_TEST_CASE( test_CovarianceAsFunctionOfTime )
{
    std::pair< std::shared_ptr< EstimationOutput< double > >, std::shared_ptr< EstimationInput< double, double > > > podData;

    // Simulate covariances directly by propagating to different final times
    std::map< int, Eigen::MatrixXd > manualCovariances;
    for( unsigned int i = 1; i < 5; i++ )
    {
        executeEarthOrbiterParameterEstimation< double, double >( podData, 1.0E7, i, 0, false );
        manualCovariances[ i ] = podData.first->getUnnormalizedCovarianceMatrix( );
    }

    // Use final calculations to compute covariance as a function of time
    std::map< double, Eigen::MatrixXd > automaticCovariances =
            simulation_setup::calculateCovarianceUsingDataUpToEpoch( podData.second, podData.first, 86400.0 - 1.0 );

    // Check consistency
    int counter = 1;
    for( std::map< double, Eigen::MatrixXd >::const_iterator covarianceIterator = automaticCovariances.begin( );
         covarianceIterator != automaticCovariances.end( );
         covarianceIterator++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceIterator->second, manualCovariances.at( counter ), 2.0E-8 );
        counter++;
    }
}

BOOST_AUTO_TEST_CASE( test_DesignMatrixSaving )
{
    // Simulate covariances directly by propagating to different final times
    for( unsigned int i = 0; i < 2; i++ )
    {
        std::vector< Eigen::MatrixXd > designMatrices;
        std::pair< std::shared_ptr< EstimationOutput< double > >, std::shared_ptr< EstimationInput< double, double > > > podData;
        executeEarthOrbiterParameterEstimation< double, double >( podData, 1.0E7, 1, 0, false, static_cast< bool >( i ) );
        designMatrices.push_back( podData.first->getNormalizedDesignMatrix( ) );
        designMatrices.push_back( podData.first->getUnnormalizedDesignMatrix( ) );
        designMatrices.push_back( podData.first->getNormalizedWeightedDesignMatrix( ) );
        designMatrices.push_back( podData.first->getUnnormalizedWeightedDesignMatrix( ) );

        if( i == 1 )
        {
            int numberOfParameters = podData.first->parameterEstimate_.rows( );
            int numberOfObservations = podData.second->getObservationCollection( )->getTotalObservableSize( );

            for( unsigned int j = 0; j < designMatrices.size( ); j++ )
            {
                BOOST_CHECK_EQUAL( designMatrices.at( j ).rows( ), numberOfObservations );
                BOOST_CHECK_EQUAL( designMatrices.at( j ).cols( ), numberOfParameters );
            }
        }

        if( i == 0 )
        {
            for( unsigned int j = 0; j < designMatrices.size( ); j++ )
            {
                BOOST_CHECK_EQUAL( designMatrices.at( j ).rows( ), 0 );
                BOOST_CHECK_EQUAL( designMatrices.at( j ).cols( ), 0 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
