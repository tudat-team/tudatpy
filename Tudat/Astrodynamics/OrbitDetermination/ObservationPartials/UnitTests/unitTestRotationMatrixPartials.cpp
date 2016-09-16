
#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/createPositionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ephemerides;
using namespace tudat::estimatable_parameters;
using namespace tudat::observation_partials;

BOOST_AUTO_TEST_SUITE( test_rotation_matrix_partaisl )

BOOST_AUTO_TEST_CASE( testSimpleRotationalEphemerisPartials )
{
    // Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");

    double nominalRotationRate = 2.0 * mathematical_constants::PI / 86400.0;
    boost::shared_ptr< SimpleRotationalEphemeris > rotationalEphemeris =
            boost::make_shared< SimpleRotationalEphemeris >(
                spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", 1.0E7 ),
                nominalRotationRate, 1.0E7, basic_astrodynamics::JULIAN_DAY_ON_J2000, "ECLIPJ2000", "IAU_Earth" );

    {
        boost::shared_ptr< RotationMatrixPartialWrtConstantRotationRate > rotationMatrixPartialObject =
                boost::make_shared< RotationMatrixPartialWrtConstantRotationRate >( rotationalEphemeris );

        double testTime = 1.0E6;
        Eigen::Matrix3d rotationMatrixPartial = rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter(
                    testTime ).at( 0 );


        double perturbation = 1.0E-12;
        rotationalEphemeris->resetRotationRate( nominalRotationRate + perturbation );
        Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                    testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 ).toRotationMatrix( );

        rotationalEphemeris->resetRotationRate( nominalRotationRate - perturbation );
        Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                    testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 ).toRotationMatrix( );

        Eigen::Matrix3d numericalRotationMatrixPartial = ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
        Eigen::Matrix3d matrixDifference = rotationMatrixPartial - numericalRotationMatrixPartial;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 0.1 );
            }
        }
    }

    {
        boost::shared_ptr< RotationMatrixPartialWrtPoleOrientation > rotationMatrixPartialObject =
                boost::make_shared< RotationMatrixPartialWrtPoleOrientation >( rotationalEphemeris );

        double testTime = 1.0E6;
        std::vector< Eigen::Matrix3d > rotationMatrixPartials = rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter(
                    testTime );

        Eigen::Vector3d nominalEulerAngles = rotationalEphemeris->getInitialEulerAngles( );
        double perturbedAngle;

        double perturbation = 1.0E-6;
        {
            {
                perturbedAngle = nominalEulerAngles( 0 ) + perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( perturbedAngle, nominalEulerAngles( 1 ) );
                Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 ).toRotationMatrix( );

                perturbedAngle = nominalEulerAngles( 0 ) - perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( perturbedAngle, nominalEulerAngles( 1 ) );
                Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 ).toRotationMatrix( );

                Eigen::Matrix3d numericalRotationMatrixPartial = ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 0 ) - numericalRotationMatrixPartial;

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-8 );
                    }
                }
            }

            {
                perturbedAngle = nominalEulerAngles( 1 ) + perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( nominalEulerAngles( 0 ), perturbedAngle );
                Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 ).toRotationMatrix( );

                perturbedAngle = nominalEulerAngles( 1 ) - perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( nominalEulerAngles( 0 ), perturbedAngle );
                Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 ).toRotationMatrix( );

                Eigen::Matrix3d numericalRotationMatrixPartial = ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 1 ) - numericalRotationMatrixPartial;

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-8 );
                    }
                }
            }
        }
    }



}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





