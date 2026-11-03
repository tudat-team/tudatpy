/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/orbit_determination/observation_partials/eulerAngleObservablePartials.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

#include "tudat/support/observationPartialTestFunctions.h"

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

BOOST_AUTO_TEST_SUITE( test_euler_angle_partials )

//! Test partial derivatives of Euler angle observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testEulerAnglePartials )
{
    std::cout << " **************************************************************************************** " << std::endl;

    // Test partials with constant rotational ephemerides (with test of rotation state partial)
    {
        // Create environment
        SystemOfBodies bodies =
                setupEnvironment( std::vector< std::pair< std::string, std::string > >( ), 1.0E7, 1.2E7, 1.1E7, false, 1.0, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ observed_body ] = std::make_pair< std::string, std::string >( "Mars", "" );

        // Generate one-way range model
        std::shared_ptr< ObservationModel< 3 > > eulerAngleModel =
                observation_models::ObservationModelCreator< 3, double, double >::createObservationModel(
                        std::make_shared< observation_models::ObservationModelSettings >( euler_angle_313_observable, linkEnds ), bodies );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7, false, true );

        testObservationPartials< 3 >(
                eulerAngleModel, bodies, fullEstimatableParameterSet, linkEnds, euler_angle_313_observable, 1.0E-6, false, false );
    }
}

//! Test all components of the partial w.r.t. the rotational state (rotation to the base frame), for both signs of the quaternion, against
//! central differences of the observable function.
BOOST_AUTO_TEST_CASE( testEulerAnglePartialWrtRotationalState )
{
    const Eigen::Quaterniond rotationToBaseFrame =
            Eigen::Quaterniond( Eigen::AngleAxisd( 0.7, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitX( ) ) *
                                Eigen::AngleAxisd( -1.9, Eigen::Vector3d::UnitZ( ) ) );
    const double perturbation = 1.0E-6;

    for( double sign : { 1.0, -1.0 } )
    {
        const Eigen::Vector4d state = sign * linear_algebra::convertQuaternionToVectorFormat( rotationToBaseFrame );
        const Eigen::Quaterniond signedRotation( state( 0 ), state( 1 ), state( 2 ), state( 3 ) );

        EulerAngleObervationPartialWrtCurrentRotationalState partialObject(
                std::make_pair( initial_rotational_body_state, std::make_pair( "Mars", "" ) ),
                [ signedRotation ]( const double ) { return signedRotation; } );
        const Eigen::Matrix< double, 3, Eigen::Dynamic > partial =
                partialObject
                        .calculatePartial( std::vector< Eigen::Vector6d >( { Eigen::Vector6d::Zero( ) } ),
                                           std::vector< double >( { 0.0 } ),
                                           observed_body )
                        .at( 0 )
                        .first;

        // The state is a unit quaternion: compare along each state entry, projected on the tangent space
        for( int i = 0; i < 4; i++ )
        {
            Eigen::Vector4d direction = Eigen::Vector4d::Unit( i );
            direction -= direction.dot( state ) * state;

            const Eigen::Vector4d upperState = ( state + perturbation * direction ).normalized( );
            const Eigen::Vector4d lowerState = ( state - perturbation * direction ).normalized( );
            const Eigen::Vector3d upperAngles = basic_mathematics::get313EulerAnglesFromQuaternion(
                    Eigen::Quaterniond( upperState( 0 ), upperState( 1 ), upperState( 2 ), upperState( 3 ) ).inverse( ) );
            const Eigen::Vector3d lowerAngles = basic_mathematics::get313EulerAnglesFromQuaternion(
                    Eigen::Quaterniond( lowerState( 0 ), lowerState( 1 ), lowerState( 2 ), lowerState( 3 ) ).inverse( ) );
            const Eigen::Vector3d numericalChange = ( upperAngles - lowerAngles ) / ( 2.0 * perturbation );
            const Eigen::Vector3d analyticalChange = partial.block( 0, 0, 3, 4 ) * direction;
            for( int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( numericalChange( j ) - analyticalChange( j ), 1.0E-8 );
            }
        }

        // The observable does not depend on the angular velocity
        BOOST_CHECK_EQUAL( partial.block( 0, 4, 3, 3 ).norm( ), 0.0 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
