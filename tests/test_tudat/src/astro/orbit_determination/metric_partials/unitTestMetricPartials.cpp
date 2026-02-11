/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <vector>

#include <Eigen/Core>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/ppnParameters.h"
#include "tudat/astro/orbit_determination/metric_partials/metricPartial.h"
#include "tudat/astro/orbit_determination/metric_partials/schwarzschildMetricPartial.h"
#include "tudat/astro/relativity/schwarzschildMetric.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::orbital_element_conversions;
using namespace tudat::orbit_determination;
using namespace tudat::orbit_determination::partial_derivatives;
using namespace tudat::relativity;

BOOST_AUTO_TEST_SUITE( test_metric_partials )

BOOST_AUTO_TEST_CASE( testFirstOrderSchwarzschildMetricPartials )
{
    const double nominalGravitationalParameter = 398600.4418E9;
    auto centralGravityField = std::make_shared< GravityFieldModel >( nominalGravitationalParameter );
    auto ppnParameterSet = std::make_shared< PPNParameterSet >( 1.0, 1.0 );
    auto schwarzschildMetric = std::make_shared< HarmonicSchwarzschildMetric >(
            std::bind( &GravityFieldModel::getGravitationalParameter, centralGravityField ),
            ppnParameterSet,
            "Earth",
            false );

    Eigen::Matrix< double, 6, 1 > initialKeplerElements;
    initialKeplerElements[ semiMajorAxisIndex ] = 8000.0E3;
    initialKeplerElements[ eccentricityIndex ] = 0.5;
    initialKeplerElements[ inclinationIndex ] = 30.0 * mathematical_constants::PI / 180.0;
    initialKeplerElements[ argumentOfPeriapsisIndex ] = 1.7;
    initialKeplerElements[ longitudeOfAscendingNodeIndex ] = 2.4;
    initialKeplerElements[ trueAnomalyIndex ] = 1.3 * mathematical_constants::PI;

    const Eigen::Matrix< double, 6, 1 > nominalState =
            convertKeplerianToCartesianElements( initialKeplerElements, centralGravityField->getGravitationalParameter( ) );

    schwarzschildMetric->update( nominalState, 0.0, true, false );

    auto metricPartial = std::make_shared< SchwarzschildMetricPartial >(
            schwarzschildMetric, std::make_pair( "Satellite", "" ) );
    metricPartial->update( );

    Eigen::Matrix< double, 6, 1 > statePerturbation;
    statePerturbation << 10.0, 10.0, 10.0, 0.1, 0.1, 0.1;

    std::vector< Eigen::Matrix< double, 4, 4 > > numericalMetricPartialWrtState;
    numericalMetricPartialWrtState.reserve( 6 );

    for( int i = 0; i < 6; ++i )
    {
        Eigen::Matrix< double, 6, 1 > perturbedState = nominalState;
        perturbedState( i ) += statePerturbation( i );
        schwarzschildMetric->update( perturbedState, 0.0, true, false );
        const Eigen::Matrix< double, 4, 4 > upPerturbedMetric =
                schwarzschildMetric->getCurrentCovariantMetricPeturbation( );

        perturbedState = nominalState;
        perturbedState( i ) -= statePerturbation( i );
        schwarzschildMetric->update( perturbedState, 0.0, true, false );
        const Eigen::Matrix< double, 4, 4 > downPerturbedMetric =
                schwarzschildMetric->getCurrentCovariantMetricPeturbation( );

        numericalMetricPartialWrtState.push_back(
                ( upPerturbedMetric - downPerturbedMetric ) / ( 2.0 * statePerturbation( i ) ) );
    }

    const auto partialFunctionPair = metricPartial->getDerivativeFunctionWrtStateOfIntegratedBody(
            std::make_pair( "Satellite", "" ), propagators::translational_state );

    BOOST_CHECK_EQUAL( partialFunctionPair.second, 6 );
    const auto analyticalMetricPartialWrtState = partialFunctionPair.first( );
    BOOST_CHECK_EQUAL( analyticalMetricPartialWrtState.size( ), 6 );

    for( int i = 0; i < 3; ++i )
    {
        for( int j = 0; j < 4; ++j )
        {
            for( int k = 0; k < 4; ++k )
            {
                if( j == k )
                {
                    BOOST_CHECK_CLOSE_FRACTION(
                            analyticalMetricPartialWrtState.at( i )( j, j ),
                            numericalMetricPartialWrtState.at( i )( j, j ),
                            1.0E-9 );
                }
                else
                {
                    BOOST_CHECK_EQUAL( analyticalMetricPartialWrtState.at( i )( j, k ), 0.0 );
                }
            }
        }
    }

    const double gravitationalParameterPerturbation = 1.0E12;

    centralGravityField->resetGravitationalParameter( nominalGravitationalParameter + gravitationalParameterPerturbation );
    schwarzschildMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix< double, 4, 4 > upMetricMu =
            schwarzschildMetric->getCurrentCovariantMetricPeturbation( );

    centralGravityField->resetGravitationalParameter( nominalGravitationalParameter - gravitationalParameterPerturbation );
    schwarzschildMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix< double, 4, 4 > downMetricMu =
            schwarzschildMetric->getCurrentCovariantMetricPeturbation( );

    const Eigen::Matrix< double, 4, 4 > numericalMetricPartialWrtMu =
            ( upMetricMu - downMetricMu ) / ( 2.0 * gravitationalParameterPerturbation );

    centralGravityField->resetGravitationalParameter( nominalGravitationalParameter );

    auto gravitationalParameterParameter = std::make_shared< estimatable_parameters::GravitationalParameter >(
            centralGravityField, "Earth" );

    const auto muPartialFunctionPair = metricPartial->getParameterPartialFunction( gravitationalParameterParameter );
    BOOST_CHECK_EQUAL( muPartialFunctionPair.second, 1 );

    const Eigen::Matrix< double, 4, 4 > analyticalMetricPartialWrtMu = muPartialFunctionPair.first( );

    for( int j = 0; j < 4; ++j )
    {
        for( int k = 0; k < 4; ++k )
        {
            if( j == k )
            {
                BOOST_CHECK_CLOSE_FRACTION(
                        analyticalMetricPartialWrtMu( j, j ),
                        numericalMetricPartialWrtMu( j, j ),
                        1.0E-9 );
            }
            else
            {
                BOOST_CHECK_EQUAL( analyticalMetricPartialWrtMu( j, k ), 0.0 );
            }
        }
    }

    const double ppnParameterPerturbation = 0.1;

    ppnParameterSet->setParameterGamma( 1.0 + ppnParameterPerturbation );
    schwarzschildMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix< double, 4, 4 > upMetricGamma =
            schwarzschildMetric->getCurrentCovariantMetricPeturbation( );

    ppnParameterSet->setParameterGamma( 1.0 - ppnParameterPerturbation );
    schwarzschildMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix< double, 4, 4 > downMetricGamma =
            schwarzschildMetric->getCurrentCovariantMetricPeturbation( );

    const Eigen::Matrix< double, 4, 4 > numericalMetricPartialWrtGamma =
            ( upMetricGamma - downMetricGamma ) / ( 2.0 * ppnParameterPerturbation );

    ppnParameterSet->setParameterGamma( 1.0 );

    auto gammaParameter = std::make_shared< estimatable_parameters::PPNParameterGamma >( ppnParameterSet );

    const auto gammaPartialFunctionPair = metricPartial->getParameterPartialFunction( gammaParameter );
    BOOST_CHECK_EQUAL( gammaPartialFunctionPair.second, 1 );
    const Eigen::Matrix< double, 4, 4 > analyticalMetricPartialWrtGamma = gammaPartialFunctionPair.first( );

    for( int j = 0; j < 4; ++j )
    {
        for( int k = 0; k < 4; ++k )
        {
            if( j == k )
            {
                BOOST_CHECK_CLOSE_FRACTION(
                        analyticalMetricPartialWrtGamma( j, j ),
                        numericalMetricPartialWrtGamma( j, j ),
                        1.0E-9 );
            }
            else
            {
                BOOST_CHECK_EQUAL( analyticalMetricPartialWrtGamma( j, k ), 0.0 );
            }
        }
    }

    ppnParameterSet->setParameterBeta( 1.0 + ppnParameterPerturbation );
    schwarzschildMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix< double, 4, 4 > upMetricBeta =
            schwarzschildMetric->getCurrentCovariantMetricPeturbation( );

    ppnParameterSet->setParameterBeta( 1.0 - ppnParameterPerturbation );
    schwarzschildMetric->update( nominalState, 0.0, true, false );
    const Eigen::Matrix< double, 4, 4 > downMetricBeta =
            schwarzschildMetric->getCurrentCovariantMetricPeturbation( );

    const Eigen::Matrix< double, 4, 4 > numericalMetricPartialWrtBeta =
            ( upMetricBeta - downMetricBeta ) / ( 2.0 * ppnParameterPerturbation );

    ppnParameterSet->setParameterBeta( 1.0 );

    auto betaParameter = std::make_shared< estimatable_parameters::PPNParameterBeta >( ppnParameterSet );
    const auto betaPartialFunctionPair = metricPartial->getParameterPartialFunction( betaParameter );
    BOOST_CHECK_EQUAL( betaPartialFunctionPair.second, 1 );
    const Eigen::Matrix< double, 4, 4 > analyticalMetricPartialWrtBeta = betaPartialFunctionPair.first( );

    for( int j = 0; j < 4; ++j )
    {
        for( int k = 0; k < 4; ++k )
        {
            if( j == k )
            {
                BOOST_CHECK_CLOSE_FRACTION(
                        analyticalMetricPartialWrtBeta( j, j ),
                        numericalMetricPartialWrtBeta( j, j ),
                        1.0E-5 );
            }
            else
            {
                BOOST_CHECK_EQUAL( analyticalMetricPartialWrtBeta( j, k ), 0.0 );
            }
        }
    }

    schwarzschildMetric->update( nominalState, 0.0, true, true );
    metricPartial->update( );

    const std::vector< Eigen::Matrix< double, 4, 4 > > directChristoffelSymbols =
            schwarzschildMetric->getCurrentChristoffelSymbols( );
    const std::vector< Eigen::Matrix< double, 4, 4 > > reconstructedChristoffelSymbols =
            calculateCurrentChristoffelSymbolsFromMetricPartials(
                    nominalState, 0.0, std::make_pair( "Satellite", "" ), schwarzschildMetric, metricPartial );

    for( unsigned int i = 0; i < directChristoffelSymbols.size( ); ++i )
    {
        for( int j = 0; j < 4; ++j )
        {
            for( int k = 0; k < 4; ++k )
            {
                // Allow slightly looser tolerance to accommodate numerical noise in double precision.
                const double tolerance = ( i == 0 || j == 0 || k == 0 ) ? 1.0E-20 : 1.0E-16;
                BOOST_CHECK_SMALL(
                        directChristoffelSymbols[ i ]( j, k ) - reconstructedChristoffelSymbols[ i ]( j, k ),
                        tolerance );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
