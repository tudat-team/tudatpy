/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120209    K. Kumar          File created.
 *
 *    References
 *      Easy calculation. Newton's Law of Gravity Tutorial,
 *          http://easycalculation.com/physics/classical-physics/learn-newtons-law.php, last
 *          accessed: 12th February, 2012.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_gravity_field_variations )

using namespace tudat::gravitation;

boost::shared_ptr< GravityFieldVariationsSet > getTestGravityFieldVariations( )
{
    std::vector< std::string > deformingBodies;
    deformingBodies.push_back( "Io" );
    deformingBodies.push_back( "Europa" );

    std::vector< boost::function< basic_mathematics::Vector6d( const double ) > >
            deformingBodyStateFunctions;
    std::vector< boost::function< double( ) > > deformingBodyMasses;

    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
    {
        deformingBodyStateFunctions.push_back(
                    boost::bind( &spice_interface::getBodyCartesianStateAtEpoch,
                                 deformingBodies.at( i ), "SSB", "J2000",
                                 "None", _1 ) );
        deformingBodyMasses.push_back(
                    boost::bind( &spice_interface::getBodyGravitationalParameter, deformingBodies.at( i ) ) );
    }

    std::vector< std::vector< std::complex< double > > > loveNumbers;
    std::complex< double > constantLoveNumber = std::complex< double >( 0.5, 0.5E-3 );

    std::vector< std::complex< double > > constantSingleDegreeLoveNumber =
    { constantLoveNumber, constantLoveNumber, constantLoveNumber };
    loveNumbers.push_back( constantSingleDegreeLoveNumber );

    // Set up gravity field variation of Jupiter due to Galilean moons.
    boost::shared_ptr< GravityFieldVariations > solidBodyGravityFieldVariations =
            boost::make_shared< BasicSolidBodyTideGravityFieldVariations >(
                boost::bind( &spice_interface::getBodyCartesianStateAtEpoch,
                             "Jupiter", "SSB", "J2000", "None", _1 ),
                boost::bind( &spice_interface::computeRotationQuaternionBetweenFrames,
                             "J2000", "IAU_Jupiter", _1 ),
                deformingBodyStateFunctions,
                spice_interface::getAverageRadius( "Jupiter" ),
                boost::bind( &spice_interface::getBodyGravitationalParameter, "Jupiter" ),
                deformingBodyMasses, loveNumbers, deformingBodies );

    boost::shared_ptr< GravityFieldVariations > tabulatedGravityFieldVariations;

    boost::shared_ptr< GravityFieldVariationsSet > gravityFieldVariationsSet =
            boost::make_shared< GravityFieldVariationsSet >(
                boost::assign::list_of( solidBodyGravityFieldVariations )( tabulatedGravityFieldVariations ),
                boost::assign::list_of( basic_solid_body )( tabulated_variation ),
                boost::assign::list_of( "BasicTidal" )( "Tabulated" ) );

    return gravityFieldVariationsSet;
}
//! Test if gravitational force is computed correctly.
BOOST_AUTO_TEST_CASE( testGravityFieldVariations )
{
    double testTime = 1.0E7;

    double gravitationalParameter;
    double referenceRadius;
    Eigen::MatrixXd nominalCosineCoefficients;
    Eigen::MatrixXd nominalSineCoefficients;

    std::vector< boost::shared_ptr< GravityFieldVariations > > gravityFieldVariationsList =
            getTestGravityFieldVariations( )->getVariationObjects( );

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > directGravityFieldVariations;
    Eigen::Matrix< double, 5, 5 > expectedCosineCoefficientsCorrections;
    Eigen::Matrix< double, 5, 5 > expectedSineCoefficientsCorrections;

    int numberOfDegrees;
    int numberOfOrders;

    for( unsigned int i = 0; i < gravityFieldVariationsList.size( ); i++ )
    {
        numberOfDegrees = gravityFieldVariationsList.at( i )->getMaximumDegree( );
        numberOfOrders = gravityFieldVariationsList.at( i )->getMaximumOrder( );

        directGravityFieldVariations =
                    gravityFieldVariationsList.at( i )->
                    calculateSphericalHarmonicsCorrections( testTime );
        expectedCosineCoefficientsCorrections.block(
                    0, 0, numberOfDegrees, numberOfOrders ) += directGravityFieldVariations.first;
        expectedSineCoefficientsCorrections.block(
                    0, 0, numberOfDegrees, numberOfOrders ) += directGravityFieldVariations.second;
    }

    boost::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
            boost::make_shared< TimeDependentSphericalHarmonicsGravityField >(
                gravitationalParameter, referenceRadius, nominalCosineCoefficients,
                nominalSineCoefficients, getTestGravityFieldVariations( ) );

    timeDependentGravityField->update( testTime );

    Eigen::MatrixXd perturbedCosineCoefficients =
            timeDependentGravityField->getCosineCoefficients( );
    Eigen::MatrixXd perturbedSineCoefficients =
            timeDependentGravityField->getSineCoefficients( );

    Eigen::MatrixXd calculatedCosineCoefficientCorrections =
            ( perturbedCosineCoefficients - nominalCosineCoefficients ).block( 0, 0, 5, 5 );
    Eigen::MatrixXd calculatedSineCoefficientCorrections =
            ( perturbedSineCoefficients - nominalSineCoefficients ).block( 0, 0, 5, 5 );

    for( unsigned int i = 0; i < 5; i++ )
    {
        for( unsigned int j = 0; j < 5; j++ )
        {
            BOOST_CHECK_SMALL( calculatedCosineCoefficientCorrections( i, j ) -
                               expectedCosineCoefficientsCorrections( i, j ), 1.0E-20 );
            BOOST_CHECK_SMALL( calculatedSineCoefficientCorrections( i, j ) -
                               expectedSineCoefficientsCorrections( i, j ), 1.0E-20 );
        }
    }



}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

