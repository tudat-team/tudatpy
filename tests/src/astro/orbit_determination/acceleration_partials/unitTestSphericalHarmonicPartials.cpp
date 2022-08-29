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

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/simulation/estimation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace gravitation;
using namespace aerodynamics;
using namespace ephemerides;
using namespace simulation_setup;
using namespace orbital_element_conversions;
using namespace unit_conversions;
using namespace orbit_determination;
using namespace acceleration_partials;
using namespace spice_interface;
using namespace orbit_determination;
using namespace estimatable_parameters;
using namespace electromagnetism;

BOOST_AUTO_TEST_SUITE( test_spherical_harmonic_partials )

BOOST_AUTO_TEST_CASE( testSphericalHarmonicPartials )
{
    // Short-cuts.
    using namespace gravitation;

    const double gravitationalParameter = 3.986004418e14;
    const double planetaryRadius = 6378137.0;

    const Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );
    const Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );

    Eigen::Vector3d position( 7.0e6, 8.0e6, 9.0e6 );
    Eigen::Vector3d nominalSphericalPosition = coordinate_conversions::
            convertCartesianToSpherical( position );
    nominalSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 -
            nominalSphericalPosition( 1 );

    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache
            = std::make_shared< basic_mathematics::SphericalHarmonicsCache >( 6, 6 );
    sphericalHarmonicsCache->update(
                nominalSphericalPosition( 0 ), std::sin( nominalSphericalPosition( 1 ) ), nominalSphericalPosition( 2 ),
                planetaryRadius );
    std::shared_ptr< basic_mathematics::LegendreCache > legendreCache = sphericalHarmonicsCache->getLegendreCache( );

    double currentLongitude = sphericalHarmonicsCache->getCurrentLongitude( );
    double currentPolynomialArgument = legendreCache->getCurrentPolynomialParameter( );

    Eigen::MatrixXd upPerturbedLegendrePolynomials = Eigen::MatrixXd( 6, 6 );
    Eigen::MatrixXd downPerturbedLegendrePolynomials = Eigen::MatrixXd( 6, 6 );

    Eigen::MatrixXd upPerturbedLegendrePolynomialPartials = Eigen::MatrixXd( 6, 6 );
    Eigen::MatrixXd downPerturbedLegendrePolynomialPartials = Eigen::MatrixXd( 6, 6 );

    Eigen::MatrixXd analyticalLegendrePolynomialPartials = Eigen::MatrixXd( 6, 6 );
    Eigen::MatrixXd analyticalLegendrePolynomialSecondPartials = Eigen::MatrixXd( 6, 6 );

    legendreCache->setComputeSecondDerivatives( 1 );
    legendreCache->update( currentPolynomialArgument + 0.1  );

    legendreCache->update( currentPolynomialArgument );

    for( unsigned int i = 0; i < 6; i++ )
    {
        for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
        {
            analyticalLegendrePolynomialPartials( i, j ) = legendreCache->getLegendrePolynomialDerivative( i, j );
            analyticalLegendrePolynomialSecondPartials( i, j ) = legendreCache->getLegendrePolynomialSecondDerivative( i, j );
        }

    }

    double polynomialArgumentPerturbation = 1.0E-6;
    {
        legendreCache->update( currentPolynomialArgument + polynomialArgumentPerturbation );
        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
            {
                upPerturbedLegendrePolynomials( i, j ) = legendreCache->getLegendrePolynomial( i, j );
                upPerturbedLegendrePolynomialPartials( i, j ) = legendreCache->getLegendrePolynomialDerivative( i, j );
            }
        }

        legendreCache->update( currentPolynomialArgument - polynomialArgumentPerturbation );
        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
            {
                downPerturbedLegendrePolynomials( i, j ) = legendreCache->getLegendrePolynomial( i, j );
                downPerturbedLegendrePolynomialPartials( i, j ) = legendreCache->getLegendrePolynomialDerivative( i, j );
            }
        }
    }

    Eigen::MatrixXd numericalLegendrePolynomialPartials =
            ( upPerturbedLegendrePolynomials - downPerturbedLegendrePolynomials ) /
            ( 2.0 * polynomialArgumentPerturbation );
    Eigen::MatrixXd numericalLegendrePolynomialSecondPartials =
            ( upPerturbedLegendrePolynomialPartials - downPerturbedLegendrePolynomialPartials ) /
            ( 2.0 * polynomialArgumentPerturbation );

    for( unsigned int i = 0; i < 6; i++ )
    {
        for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( numericalLegendrePolynomialPartials( i, j ) - analyticalLegendrePolynomialPartials( i, j ) ), 1.0E-8 );
            BOOST_CHECK_SMALL(
                        std::fabs( numericalLegendrePolynomialSecondPartials( i, j ) - analyticalLegendrePolynomialSecondPartials( i, j ) ), 1.0E-8 );
        }
    }


    std::vector< std::vector< Eigen::Vector3d > > sphericalPotentialGradients;

    std::vector< std::vector< Eigen::Matrix3d > > upPerturbedSphericalPotentialGradients;
    std::vector< std::vector< Eigen::Matrix3d > > downPerturbedSphericalPotentialGradients;

    std::vector< std::vector< Eigen::Matrix3d > > numericalSphericalPotentialHessian;
    std::vector< std::vector< Eigen::Matrix3d > > analyticalSphericalPotentialHessian;

    Eigen::Matrix3d normalization;
    normalization <<  nominalSphericalPosition( 0 ) * nominalSphericalPosition( 0 ), nominalSphericalPosition( 0 ), nominalSphericalPosition( 0 ),
            nominalSphericalPosition( 0 ), 1.0 , 1.0,
            nominalSphericalPosition( 0 ), 1.0, 1.0;
    for( unsigned int i = 0; i < 6; i++ )
    {
        std::vector< Eigen::Vector3d > singleTermPotentialGradients;
        singleTermPotentialGradients.resize( 6 );
        sphericalPotentialGradients.push_back( singleTermPotentialGradients );

        std::vector< Eigen::Matrix3d > singleTermPotentialPartials;
        singleTermPotentialPartials.resize( 6 );
        upPerturbedSphericalPotentialGradients.push_back( singleTermPotentialPartials );
        downPerturbedSphericalPotentialGradients.push_back( singleTermPotentialPartials );

        numericalSphericalPotentialHessian.push_back( singleTermPotentialPartials );
        analyticalSphericalPotentialHessian.push_back( singleTermPotentialPartials );


    }

    sphericalHarmonicsCache->update( position.norm( ), currentPolynomialArgument, currentLongitude, planetaryRadius );

    for( unsigned int i = 0; i < 6; i++ )
    {
        for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
        {
            sphericalPotentialGradients[ i ][ j ] = basic_mathematics::computePotentialGradient(
                        nominalSphericalPosition, gravitationalParameter / planetaryRadius, i, j,
                        cosineCoefficients( i, j ), sineCoefficients( i, j ), legendreCache->getLegendrePolynomial( i, j ),
                        legendreCache->getLegendrePolynomialDerivative( i, j ), sphericalHarmonicsCache );
        }
    }


    Eigen::Vector3d perturbedSphericalPosition;
    Eigen::Vector3d sphericalStatePerturbation;
    sphericalStatePerturbation << 10.0, 1.0E-7, 1.0E-8;

    for( unsigned parameter = 0; parameter < 3; parameter++ )
    {
        perturbedSphericalPosition = nominalSphericalPosition;
        perturbedSphericalPosition( parameter ) += sphericalStatePerturbation( parameter );

        sphericalHarmonicsCache->update(
                    perturbedSphericalPosition( 0 ), std::sin( perturbedSphericalPosition( 1 ) ),
                    perturbedSphericalPosition( 2 ), planetaryRadius );

        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
            {
                upPerturbedSphericalPotentialGradients[ i ][ j ].block( 0, parameter, 3, 1 ) = basic_mathematics::computePotentialGradient(
                            perturbedSphericalPosition, gravitationalParameter / planetaryRadius, i, j,
                            cosineCoefficients( i, j ), sineCoefficients( i, j ), legendreCache->getLegendrePolynomial( i, j ),
                            legendreCache->getLegendrePolynomialDerivative( i, j ), sphericalHarmonicsCache );
            }
        }

        perturbedSphericalPosition = nominalSphericalPosition;
        perturbedSphericalPosition( parameter ) -= sphericalStatePerturbation( parameter );

        sphericalHarmonicsCache->update(
                    perturbedSphericalPosition( 0 ), std::sin( perturbedSphericalPosition( 1 ) ),
                    perturbedSphericalPosition( 2 ), planetaryRadius );

        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
            {
                downPerturbedSphericalPotentialGradients[ i ][ j ].block( 0, parameter, 3, 1 ) = basic_mathematics::computePotentialGradient(
                            perturbedSphericalPosition, gravitationalParameter / planetaryRadius, i, j,
                            cosineCoefficients( i, j ), sineCoefficients( i, j ), legendreCache->getLegendrePolynomial( i, j ),
                            legendreCache->getLegendrePolynomialDerivative( i, j ), sphericalHarmonicsCache );
            }
        }

        for( unsigned int i = 0; i < 6; i++ )
        {
            for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
            {
                numericalSphericalPotentialHessian[ i ][ j ].block( 0, parameter, 3, 1 ) =
                        ( ( upPerturbedSphericalPotentialGradients[ i ][ j ].block( 0, parameter, 3, 1 ) -
                            downPerturbedSphericalPotentialGradients[ i ][ j ].block( 0, parameter, 3, 1 ) ) /
                          ( 2.0 * sphericalStatePerturbation( parameter ) ) );

            }

        }
    }

    for( unsigned int i = 1; i < 6; i++ )
    {
        for( unsigned int j = 0; ( j <= i && j < 6 ); j++ )
        {
            computePotentialSphericalHessian(
                        nominalSphericalPosition, gravitationalParameter / planetaryRadius, i, j,
                        cosineCoefficients( i, j ), sineCoefficients( i, j ), sphericalHarmonicsCache,
                        analyticalSphericalPotentialHessian[ i ][ j ] );

            analyticalSphericalPotentialHessian[ i ][ j ] = analyticalSphericalPotentialHessian[ i ][ j ].cwiseProduct(
                        normalization );
            numericalSphericalPotentialHessian[ i ][ j ] = numericalSphericalPotentialHessian[ i ][ j ].cwiseProduct(
                        normalization );
            for( unsigned int k = 0; k < 3; k++ )
            {
                for( unsigned int l = 0; l < 3; l++ )
                {
                    BOOST_CHECK_SMALL( analyticalSphericalPotentialHessian[ i ][ j ]( k, l ) -
                            numericalSphericalPotentialHessian[ i ][ j ]( k, l ), 2.5E-5 );
                }

            }

        }
    }

    Eigen::Matrix3d cumulativeSphericalHessian =  computeCumulativeSphericalHessian(
                nominalSphericalPosition, planetaryRadius, gravitationalParameter, cosineCoefficients, sineCoefficients,
                sphericalHarmonicsCache );

    Eigen::Matrix3d nominalGradientTransformationMatrix =
            coordinate_conversions::getSphericalToCartesianGradientMatrix( position );

    Eigen::Vector3d upPerturbedTotalGradient;
    Eigen::Vector3d downPerturbedTotalGradient;
    Eigen::Vector3d perturbedCartesianPosition;

    Eigen::Matrix3d numericalTotalSphericalGradient;

    sphericalStatePerturbation( 0 ) *= 10.0;
    sphericalStatePerturbation( 1 ) *= 100.0;
    sphericalStatePerturbation( 2 ) *= 1000.0;
    std::map< std::pair< int, int >, Eigen::Vector3d > dummyMap;

    for( unsigned parameter = 0; parameter < 3; parameter++ )
    {
        perturbedSphericalPosition = nominalSphericalPosition;
        perturbedSphericalPosition( parameter ) += sphericalStatePerturbation( parameter );

        perturbedSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - perturbedSphericalPosition( 1 );
        perturbedCartesianPosition = coordinate_conversions::convertSphericalToCartesian(
                    perturbedSphericalPosition );

        upPerturbedTotalGradient =
                coordinate_conversions::getSphericalToCartesianGradientMatrix( perturbedCartesianPosition ).inverse( ) *
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache, dummyMap );

        perturbedSphericalPosition = nominalSphericalPosition;
        perturbedSphericalPosition( parameter ) -= sphericalStatePerturbation( parameter );

        perturbedSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - perturbedSphericalPosition( 1 );
        perturbedCartesianPosition = coordinate_conversions::convertSphericalToCartesian(
                    perturbedSphericalPosition );

        downPerturbedTotalGradient =
                coordinate_conversions::getSphericalToCartesianGradientMatrix( perturbedCartesianPosition ).inverse( ) *
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache, dummyMap );

        numericalTotalSphericalGradient.block( 0, parameter, 3, 1 ) =
                ( upPerturbedTotalGradient - downPerturbedTotalGradient ) / ( 2.0 * sphericalStatePerturbation( parameter ) );

    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                cumulativeSphericalHessian, numericalTotalSphericalGradient, 1.0E-6 );

    Eigen::Vector3d cartesianStatePerturbation;
    cartesianStatePerturbation << 10.0, 10.0, 10.0;

    for( unsigned parameter = 0; parameter < 3; parameter++ )
    {
        perturbedCartesianPosition = position;
        perturbedCartesianPosition( parameter ) += cartesianStatePerturbation( parameter );

        upPerturbedTotalGradient =
                nominalGradientTransformationMatrix *
                coordinate_conversions::getSphericalToCartesianGradientMatrix( perturbedCartesianPosition ).inverse( ) *
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache, dummyMap );

        perturbedCartesianPosition = position;
        perturbedCartesianPosition( parameter ) -= cartesianStatePerturbation( parameter );

        downPerturbedTotalGradient =
                nominalGradientTransformationMatrix *
                coordinate_conversions::getSphericalToCartesianGradientMatrix( perturbedCartesianPosition ).inverse( ) *
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache, dummyMap );

        numericalTotalSphericalGradient.block( 0, parameter, 3, 1 ) =
                ( upPerturbedTotalGradient - downPerturbedTotalGradient ) / ( 2.0 * cartesianStatePerturbation( parameter ) );

    }

    Eigen::Matrix3d nominalGradientTransformationMatrixTranspose = nominalGradientTransformationMatrix.transpose( );
    Eigen::Matrix3d computedTotalSphericalGradient  =
            ( nominalGradientTransformationMatrix * cumulativeSphericalHessian  ) *
            nominalGradientTransformationMatrixTranspose;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                computedTotalSphericalGradient,
                numericalTotalSphericalGradient, 1.0E-6 );


    for( unsigned parameter = 0; parameter < 3; parameter++ )
    {
        perturbedCartesianPosition = position;
        perturbedCartesianPosition( parameter ) += cartesianStatePerturbation( parameter );

        upPerturbedTotalGradient =
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache, dummyMap );

        perturbedCartesianPosition = position;
        perturbedCartesianPosition( parameter ) -= cartesianStatePerturbation( parameter );

        downPerturbedTotalGradient =
                computeGeodesyNormalizedGravitationalAccelerationSum(
                    perturbedCartesianPosition, gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients,
                    sphericalHarmonicsCache, dummyMap );

        numericalTotalSphericalGradient.block( 0, parameter, 3, 1 ) =
                ( upPerturbedTotalGradient - downPerturbedTotalGradient ) / ( 2.0 * cartesianStatePerturbation( parameter ) );

    }


    Eigen::Matrix3d totalGradientCartesianPartial =
            computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
                position, planetaryRadius, gravitationalParameter, cosineCoefficients, sineCoefficients,
                sphericalHarmonicsCache );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                totalGradientCartesianPartial, numericalTotalSphericalGradient, 1.0E-6 );
}

//! Function to get tidal deformation model for Earth
std::vector< std::shared_ptr< GravityFieldVariationSettings > > getEarthGravityFieldVariationSettings( )
{
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariations;

    std::vector< std::string > deformingBodies;
    deformingBodies.push_back( "Moon" );

    std::map< int, std::vector< std::complex< double > > > loveNumbers;

    std::vector< std::complex< double > > degreeTwoLoveNumbers_;
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    std::vector< std::complex< double > > degreeThreeLoveNumbers_;
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    loveNumbers[ 2 ] = degreeTwoLoveNumbers_;
    loveNumbers[ 3 ] =degreeThreeLoveNumbers_;


    std::shared_ptr< GravityFieldVariationSettings > singleGravityFieldVariation =
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBodies, loveNumbers );
    gravityFieldVariations.push_back( singleGravityFieldVariation );
    return gravityFieldVariations;
}

BOOST_AUTO_TEST_CASE( testSphericalHarmonicAccelerationPartial )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create empty bodies, earth and vehicle.
    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    std::shared_ptr< Body > vehicle = std::make_shared< Body >( );

    const double gravitationalParameter = 3.986004418e14;
    const double planetaryRadius = 6378137.0;

    Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );


    Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );

    std::shared_ptr< GravityFieldSettings > earthGravityFieldSettings =
            std::make_shared< SphericalHarmonicsGravityFieldSettings >
            ( gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients, "IAU_Earth" );

    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings =
            getEarthGravityFieldVariationSettings( );


    SystemOfBodies bodies;
    bodies.addBody( earth, "Earth" );
    bodies.addBody( vehicle, "Vehicle" );;
    bodies.addBody( createSystemOfBodies( getDefaultBodySettings( { "Moon" } ) ).at( "Moon" ), "Moon" );

    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > simpleRotationalEphemeris =
            std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000" , "IAU_Earth", 0.0 ),
                2.0 * mathematical_constants::PI / 86400.0,
                1.0E7,
                "ECLIPJ2000" , "IAU_Earth" );
    earth->setRotationalEphemeris( simpleRotationalEphemeris );

    std::shared_ptr< tudat::gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
            std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField  >(
                tudat::simulation_setup::createGravityFieldModel(
                    earthGravityFieldSettings, "Earth", bodies, gravityFieldVariationSettings ) );
    earth->setGravityFieldModel( earthGravityField );
    bodies.at( "Earth" )->setGravityFieldVariationSet(
                createGravityFieldModelVariationsSet(
                    "Earth", bodies, gravityFieldVariationSettings ) );
    bodies.at( "Earth" )->updateConstantEphemerisDependentMemberQuantities( );

    


    // Set current state of vehicle and earth.
    double testTime = 1.0E6;
    earth->setState( Eigen::Vector6d::Zero( ) );
    earth->setCurrentRotationToLocalFrameFromEphemeris( testTime );
    bodies.at( "Moon" ) ->setState( tudat::spice_interface::getBodyCartesianStateAtEpoch( "Moon", "Earth", "ECLIPJ2000" ,"None", testTime ) );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    Eigen::Vector6d asterixInitialState = orbital_element_conversions::convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, gravitationalParameter );


    vehicle->setState( asterixInitialState );

    // Create acceleration due to vehicle on earth.
    std::shared_ptr< SphericalHarmonicAccelerationSettings > accelerationSettings =
            std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 );
    std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > gravitationalAcceleration =
            std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                createAccelerationModel( vehicle, earth, accelerationSettings, "Vehicle", "Earth" ) );

    gravitationalAcceleration->updateMembers( 0.0 );
    gravitationalAcceleration->getAcceleration( );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 10.0, 10.0, 10.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, earth, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > vehicleStateSetFunction =
            std::bind( &Body::setState, vehicle, std::placeholders::_1  );
    std::function< Eigen::Vector6d ( ) > earthStateGetFunction =
            std::bind( &Body::getState, earth );
    std::function< Eigen::Vector6d ( ) > vehicleStateGetFunction =
            std::bind( &Body::getState, vehicle );



    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialRotationalStateEstimatableParameterSettings< double > >(
                                  "Earth", 0.0, "ECLIPJ2000" ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", gravitational_parameter) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", constant_rotation_rate ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", rotation_pole_position ) );


    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 5, 4, "Earth", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 1, 5, 4, "Earth", spherical_harmonics_sine_coefficient_block ) );

    parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                                  "Earth", 2, std::vector< int >{ 2, 0, 1 }, "", false ) );
    parameterNames.push_back( std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                                  "Earth", 2, "Moon", true ) );
    parameterNames.push_back( std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                                  "Earth", 3, "", false ) );
    parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                                  "Earth", 3, std::vector< int >{ 0, 3 }, "", true ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodies );

    // Check if incompatible tidal parameters correctly throw an error
    {
        bool isExceptionCaught = false;
        std::vector< std::shared_ptr< EstimatableParameterSettings > > wrongParameterNames;
        wrongParameterNames.resize( 1 );
        wrongParameterNames[ 0 ] = std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                    "Earth", 2, std::vector< int >{ 2, 0, 1 }, "Sun", false );
        try
        {
            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
                    createParametersToEstimate( wrongParameterNames, bodies );
        }
        catch( std::runtime_error const& )
        {
            isExceptionCaught = true;
        }
        BOOST_CHECK_EQUAL( isExceptionCaught, true );

        std::vector< std::string > deformingBodyNames;
        deformingBodyNames.push_back( "Moon" );
        deformingBodyNames.push_back( "Sun" );

        wrongParameterNames[ 0 ] = std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                    "Earth", 2, std::vector< int >{ 2, 0, 1 }, deformingBodyNames, false );
        isExceptionCaught = false;
        try
        {
            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
                    createParametersToEstimate( wrongParameterNames, bodies );
        }
        catch( std::runtime_error const& )
        {
            isExceptionCaught = true;
        }
        BOOST_CHECK_EQUAL( isExceptionCaught, true );

        wrongParameterNames[ 0 ] = std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                    "Earth", 2, deformingBodyNames, true );
        isExceptionCaught = false;
        try
        {
            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
                    createParametersToEstimate( wrongParameterNames, bodies );
        }
        catch( std::runtime_error const& )
        {
            isExceptionCaught = true;
        }
        BOOST_CHECK_EQUAL( isExceptionCaught, true );

        wrongParameterNames[ 0 ] = std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                    "Earth", 3, "Sun", false );
        isExceptionCaught = false;
        try
        {
            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
                    createParametersToEstimate( wrongParameterNames, bodies );
        }
        catch( std::runtime_error const& )
        {
            isExceptionCaught = true;
        }
        BOOST_CHECK_EQUAL( isExceptionCaught, true );
    }



    // Create acceleration partial object.
    std::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartial =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityPartial > (
                createAnalyticalAccelerationPartial(
                    gravitationalAcceleration, std::make_pair( "Vehicle", vehicle ), std::make_pair( "Earth", earth ),
                    bodies, parameterSet ) );

    accelerationPartial->update( testTime );

    Eigen::MatrixXd partialWrtVehiclePosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtVehiclePosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtVehicleVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtVehicleVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );
    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtEarthOrientation = Eigen::MatrixXd::Zero( 3, 7 );
    accelerationPartial->wrtNonTranslationalStateOfAdditionalBody(
                partialWrtEarthOrientation.block( 0, 0, 3, 7 ), std::make_pair( "Earth", "" ), propagators::rotational_state );

    // Calculate numerical partials.
    testPartialWrtVehiclePosition = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, gravitationalAcceleration, vehicle->getState( ), positionPerturbation, 0 );
    testPartialWrtVehicleVelocity = calculateAccelerationWrtStatePartials(
                vehicleStateSetFunction, gravitationalAcceleration, vehicle->getState( ), velocityPerturbation, 3 );

    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), positionPerturbation, 0 );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), velocityPerturbation, 3 );

    Eigen::Vector4d orientationPerturbation;
    orientationPerturbation << 1.0E-8, 1.0E-8, 1.0E-8, 1.0E-8;

    std::function< void( Eigen::Vector7d ) > earthRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, earth, std::placeholders::_1 );
    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd accelerationDeviations = calculateAccelerationDeviationDueToOrientationChange(
                earthRotationalStateSetFunction, gravitationalAcceleration, earth->getRotationalStateVector( ), orientationPerturbation,
                appliedQuaternionPerturbation );

    // Calculate numerical partials.
    std::map< int, std::shared_ptr< EstimatableParameter< double > > > doubleParameters =
            parameterSet->getDoubleParameters( );
    std::map< int, std::shared_ptr< EstimatableParameter< double > > >::iterator doubleParametersIterator =
            doubleParameters.begin( );

    Eigen::Vector3d testPartialWrtEarthGravitationalParameter = calculateAccelerationWrtParameterPartials(
                doubleParametersIterator->second, gravitationalAcceleration, 1.0E12 );
    Eigen::Vector3d partialWrtEarthGravitationalParameter = accelerationPartial->wrtParameter(
                doubleParametersIterator->second );
    doubleParametersIterator++;

    Eigen::Vector3d partialWrtEarthRotationRate = accelerationPartial->wrtParameter(
                doubleParametersIterator->second );
    Eigen::Vector3d testPartialWrtEarthRotationRate = calculateAccelerationWrtParameterPartials(
                doubleParametersIterator->second, gravitationalAcceleration, 1.0E-12, &emptyFunction, testTime, std::bind(
                    &Body::setCurrentRotationToLocalFrameFromEphemeris, earth, std::placeholders::_1 ) );


    std::map< int, std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameters =
            parameterSet->getVectorParameters( );
    std::map< int, std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >::iterator vectorParametersIterator =
            vectorParameters.begin( );
    Eigen::MatrixXd partialWrtPolePosition = accelerationPartial->wrtParameter(
                vectorParametersIterator->second );
    Eigen::MatrixXd testPartialWrtPosition = calculateAccelerationWrtParameterPartials(
                vectorParametersIterator->second, gravitationalAcceleration, Eigen::Vector2d::Constant( 1.0E-6 ),
                &emptyFunction, testTime, std::bind(
                    &Body::setCurrentRotationToLocalFrameFromEphemeris, earth, std::placeholders::_1 ) );
    vectorParametersIterator++;

    std::function< void( ) > sphericalHarmonicFieldUpdate =
            std::bind( &tudat::gravitation::TimeDependentSphericalHarmonicsGravityField::update,
                       std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( earthGravityField ), testTime );

    Eigen::MatrixXd partialWrtCosineCoefficients = accelerationPartial->wrtParameter(
                vectorParametersIterator->second );
    Eigen::MatrixXd testPartialWrtCosineCoefficients = calculateAccelerationWrtParameterPartials(
                vectorParametersIterator->second, gravitationalAcceleration, vectorParametersIterator->second->getParameterValue( ) * 1.0E-2,
                sphericalHarmonicFieldUpdate );
    vectorParametersIterator++;

    accelerationPartial->getParameterPartialFunction(  vectorParametersIterator->second );
    Eigen::MatrixXd partialWrtSineCoefficients = accelerationPartial->wrtParameter(
                vectorParametersIterator->second );
    Eigen::MatrixXd testPartialWrtSineCoefficients = calculateAccelerationWrtParameterPartials(
                vectorParametersIterator->second, gravitationalAcceleration, vectorParametersIterator->second->getParameterValue( ) * 1.0E-2,
                sphericalHarmonicFieldUpdate );
    vectorParametersIterator++;

    Eigen::MatrixXd partialWrtDegreeTwoLoveNumberAtSeparateOrders = accelerationPartial->wrtParameter(
                vectorParametersIterator->second );
    Eigen::MatrixXd testPartialWrtDegreeTwoOrderTwoLoveNumberAtSeparateOrders = calculateAccelerationWrtParameterPartials(
                vectorParametersIterator->second, gravitationalAcceleration, Eigen::VectorXd::Constant( 3, 1.0 ), sphericalHarmonicFieldUpdate );
    vectorParametersIterator++;

    Eigen::MatrixXd partialWrtComplexDegreeTwoLoveNumber = accelerationPartial->wrtParameter(
                vectorParametersIterator->second );
    Eigen::MatrixXd testPartialWrtComplexDegreeTwoLoveNumber = calculateAccelerationWrtParameterPartials(
                vectorParametersIterator->second, gravitationalAcceleration, Eigen::VectorXd::Constant( 2, 1.0 ), sphericalHarmonicFieldUpdate );
    vectorParametersIterator++;

    Eigen::MatrixXd partialWrtDegreeThreeLoveNumber = accelerationPartial->wrtParameter(
                vectorParametersIterator->second );
    Eigen::MatrixXd testPartialWrtDegreeThreeLoveNumber = calculateAccelerationWrtParameterPartials(
                vectorParametersIterator->second, gravitationalAcceleration, Eigen::VectorXd::Constant( 1, 10.0 ), sphericalHarmonicFieldUpdate );
    vectorParametersIterator++;

    Eigen::MatrixXd partialWrtComplexDegreeThreeLoveNumberAtSeparateOrder = accelerationPartial->wrtParameter(
                vectorParametersIterator->second );
    Eigen::MatrixXd testPartialWrtComplexDegreeThreeLoveNumberAtSeparateOrder = calculateAccelerationWrtParameterPartials(
                vectorParametersIterator->second, gravitationalAcceleration, Eigen::VectorXd::Constant( 4, 10.0 ), sphericalHarmonicFieldUpdate );




    Eigen::VectorXd nominalTidalParameter = vectorParametersIterator->second->getParameterValue( );

    vectorParametersIterator->second->setParameterValue( nominalTidalParameter + Eigen::VectorXd::Constant( nominalTidalParameter.rows( ), 1.0 ) );
    earthGravityField->update( testTime );
    Eigen::MatrixXd upperturbedCosineCoefficients =
            earthGravityField->getCosineCoefficients( ).block( 0, 0, 3, 3 );
    Eigen::MatrixXd upperturbedSineCoefficients =
            earthGravityField->getSineCoefficients( ).block( 0, 0, 3, 3 );

    vectorParametersIterator->second->setParameterValue( nominalTidalParameter - Eigen::VectorXd::Constant( nominalTidalParameter.rows( ), 1.0 ) );
    earthGravityField->update( testTime );
    Eigen::MatrixXd downperturbedCosineCoefficients =
            earthGravityField->getCosineCoefficients( ).block( 0, 0, 3, 3 );
    Eigen::MatrixXd downperturbedSineCoefficients =
            earthGravityField->getSineCoefficients( ).block( 0, 0, 3, 3 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehiclePosition, partialWrtVehiclePosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtVehicleVelocity, partialWrtVehicleVelocity, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition, partialWrtEarthPosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity, partialWrtEarthVelocity, 1.0E-6 );

    //    // Compare numerical and analytical results.
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector3d numericalChangeInAcceleration = accelerationDeviations.block( 0, index - 1, 3, 1 );
        Eigen::Vector3d analyticalChangeInAcceleration =
                partialWrtEarthOrientation.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                partialWrtEarthOrientation.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInAcceleration,
                                           analyticalChangeInAcceleration, 1.0E-6 );
    }


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthGravitationalParameter, partialWrtEarthGravitationalParameter, 1.0E-12 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthRotationRate, partialWrtEarthRotationRate, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPosition, partialWrtPolePosition, 1.0E-6 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtCosineCoefficients, partialWrtCosineCoefficients, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtSineCoefficients, partialWrtSineCoefficients, 1.0E-6 );

    BOOST_CHECK_EQUAL( testPartialWrtCosineCoefficients.cols( ), 17 );
    BOOST_CHECK_EQUAL( testPartialWrtSineCoefficients.cols( ), 13 );


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtDegreeTwoLoveNumberAtSeparateOrders, testPartialWrtDegreeTwoOrderTwoLoveNumberAtSeparateOrders, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtComplexDegreeTwoLoveNumber, testPartialWrtComplexDegreeTwoLoveNumber, 1.0E-6 );


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtDegreeThreeLoveNumber, testPartialWrtDegreeThreeLoveNumber, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtComplexDegreeThreeLoveNumberAtSeparateOrder, testPartialWrtComplexDegreeThreeLoveNumberAtSeparateOrder, 1.0E-6 );

}

//! Unit test to check working onf spherical harmonic state partial for synchronously rotating body (and rotation depending on state)
BOOST_AUTO_TEST_CASE( testSphericalHarmonicAccelerationPartialWithSynchronousRotation )
{
    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings =
            std::make_shared< SynchronousRotationModelSettings >(
                "Moon", "ECLIPJ2000", "IAU_Earth" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    std::shared_ptr< tudat::simulation_setup::Body > earth = bodies.at( "Earth" );
    std::shared_ptr< tudat::simulation_setup::Body > moon = bodies.at( "Moon" );
    std::dynamic_pointer_cast< tudat::ephemerides::SynchronousRotationalEphemeris >(
                earth->getRotationalEphemeris( ) )->setIsBodyInPropagation( 1 );

    // Set translational and rotational state of bodies
    double testTime = 1.0E6;
    earth->setStateFromEphemeris( testTime );
    Eigen::Vector6d moonState = moon->getStateInBaseFrameFromEphemeris( testTime );
    moon->setState( moonState * 0.1 );

    earth->setCurrentRotationToLocalFrameFromEphemeris( testTime );
    moon->setCurrentRotationToLocalFrameFromEphemeris( testTime );

    // Create acceleration model
    std::shared_ptr< SphericalHarmonicAccelerationSettings > accelerationSettings =
            std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 );
    std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > gravitationalAcceleration =
            std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                createAccelerationModel( moon, earth, accelerationSettings, "Moon", "Earth" ) );
    gravitationalAcceleration->updateMembers( 0.0 );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMoonPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMoonVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in position for numerical partial/
    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1000.0, 1000.0, 1000.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0, 1.0E-1, 1.0;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector6d ) > earthStateSetFunction =
            std::bind( &Body::setState, earth, std::placeholders::_1  );
    std::function< void( Eigen::Vector6d ) > moonStateSetFunction =
            std::bind( &Body::setState, moon, std::placeholders::_1  );
    std::function< Eigen::Vector6d ( ) > earthStateGetFunction =
            std::bind( &Body::getState, earth );
    std::function< Eigen::Vector6d ( ) > moonStateGetFunction =
            std::bind( &Body::getState, moon );


    // Define estimated parameters
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialRotationalStateEstimatableParameterSettings< double > >(
                                  "Earth", 0.0, "ECLIPJ2000" ) );
    parameterNames.push_back( std::make_shared< InitialRotationalStateEstimatableParameterSettings< double > >(
                                  "Moon", 0.0, "ECLIPJ2000" ) );
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodies );


    // Create acceleration partial object.
    std::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartial =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityPartial > (
                createAnalyticalAccelerationPartial(
                    gravitationalAcceleration, std::make_pair( "Moon", moon ), std::make_pair( "Earth", earth ),
                    bodies, parameterSet ) );
    accelerationPartial->update( testTime );

    // Calculate analytical partials.
    Eigen::MatrixXd partialWrtMoonPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratedBody( partialWrtMoonPosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtMoonVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratedBody( partialWrtMoonVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtEarthPosition = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtPositionOfAcceleratingBody( partialWrtEarthPosition.block( 0, 0, 3, 3 ) );

    Eigen::MatrixXd partialWrtEarthVelocity = Eigen::Matrix3d::Zero( );
    accelerationPartial->wrtVelocityOfAcceleratingBody( partialWrtEarthVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    // Calculate numerical partials.
    std::function< void( ) > updateFunction =
            std::bind( &Body::setCurrentRotationToLocalFrameFromEphemeris, bodies.at( "Earth" ), testTime );

    testPartialWrtMoonPosition = calculateAccelerationWrtStatePartials(
                moonStateSetFunction, gravitationalAcceleration, moon->getState( ), positionPerturbation, 0,
                updateFunction );
    testPartialWrtMoonVelocity = calculateAccelerationWrtStatePartials(
                moonStateSetFunction, gravitationalAcceleration, moon->getState( ), velocityPerturbation, 3,
                updateFunction );

    testPartialWrtEarthPosition = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), positionPerturbation, 0,
                updateFunction );
    testPartialWrtEarthVelocity = calculateAccelerationWrtStatePartials(
                earthStateSetFunction, gravitationalAcceleration, earth->getState( ), velocityPerturbation, 3,
                updateFunction );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonPosition, partialWrtMoonPosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMoonVelocity, partialWrtMoonVelocity, 1.0E-3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthPosition, partialWrtEarthPosition, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtEarthVelocity, partialWrtEarthVelocity, 1.0E-3 );

}
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat


