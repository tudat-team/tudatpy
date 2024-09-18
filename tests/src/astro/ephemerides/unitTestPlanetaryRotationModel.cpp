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

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"



namespace tudat
{
    namespace unit_tests
    {
        
        using namespace tudat::ephemerides;
        using namespace tudat::simulation_setup;
        
        BOOST_AUTO_TEST_SUITE( test_planetary_rotation_model )
        
        BOOST_AUTO_TEST_CASE( testPlanetaryRotationModel )
        {
            double initialTime = 0.0;
            double finalTime = 1.0E8;
            
            spice_interface::loadStandardSpiceKernels( );
            
            SystemOfBodies bodies;
            bodies.createEmptyBody( "Mars" );
            
            std::shared_ptr< RotationModelSettings > defaultMarsRotationSettings =
                    getHighAccuracyMarsRotationModel( );
            
            std::shared_ptr< RotationalEphemeris > marsRotationModel =
                    createRotationModel( defaultMarsRotationSettings, "Mars" );
            
            std::shared_ptr< PlanetaryRotationModelSettings > modifiedMarsRotationModel =
                    std::dynamic_pointer_cast< PlanetaryRotationModelSettings >( defaultMarsRotationSettings );

            modifiedMarsRotationModel->setPeriodTermsToZero( );
            
            std::shared_ptr< RotationalEphemeris > marsSimplifiedRotationModel =
                    createRotationModel( modifiedMarsRotationModel, "Mars" );

            
            std::shared_ptr< RotationalEphemeris > spiceMarsRotationModel =
                    std::make_shared< SpiceRotationalEphemeris >( "ECLIPJ2000", "IAU_Mars" );

            
            double timeStep = 3600.0;
            double currentTime = initialTime + timeStep;
            double dt = 0.9;
            
            std::map< double, Eigen::MatrixXd > rotationDifferenceMap;
            std::map< double, Eigen::MatrixXd > DerivativeDifferenceMap;
            
            while( currentTime < finalTime - timeStep )
            {

                rotationDifferenceMap[ currentTime ] =  ( marsRotationModel->getRotationToTargetFrame( currentTime ) ).toRotationMatrix( ) -
                        ( marsSimplifiedRotationModel->getRotationToTargetFrame( currentTime ) ).toRotationMatrix( );

//                std::cout<<rotationDifferenceMap[ currentTime ]<<std::endl<<std::endl;

                rotationDifferenceMap[ currentTime ] =  ( marsRotationModel->getRotationToTargetFrame( currentTime ) ).toRotationMatrix( ) -
                        ( spiceMarsRotationModel->getRotationToTargetFrame( currentTime ) ).toRotationMatrix( );
                
//                std::cout<<rotationDifferenceMap[ currentTime ]<<std::endl<<std::endl;

//                rotationDifferenceMap[ currentTime ] =
//                        ( marsSimplifiedRotationModel->getRotationToTargetFrame( currentTime ) ).toRotationMatrix( ) -
//                        ( spiceMarsRotationModel->getRotationToTargetFrame( currentTime ) ).toRotationMatrix( );

//                std::cout<<rotationDifferenceMap[ currentTime ]<<std::endl<<std::endl;

                Eigen::MatrixXd numericalRotationMatrixderivative =
                        ( ( marsRotationModel->getRotationToBaseFrame( currentTime + dt ) ).toRotationMatrix( ) -
                          (marsRotationModel->getRotationToBaseFrame( currentTime - dt ) ).toRotationMatrix( ) ) / ( 2*dt );

                DerivativeDifferenceMap [ currentTime ] =
                        marsRotationModel->getDerivativeOfRotationToBaseFrame( currentTime ) - numericalRotationMatrixderivative;


                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( rotationDifferenceMap[ currentTime ]( i, j ), 1.0E-4 );
                        BOOST_CHECK_SMALL( DerivativeDifferenceMap[ currentTime ]( i, j ), 1.0E-10 );
                    }
                }

                currentTime += timeStep;
            }

//            std::cout << "analitical Rotation Matrix derivative" << std::endl;
//            std::cout <<  marsRotationModel->getDerivativeOfRotationToBaseFrame( (finalTime - timeStep) ) << std::endl;

//            Eigen::MatrixXd numericalRotationMatrixderivative = (
//                        ( marsRotationModel->getRotationToBaseFrame( (finalTime - timeStep) + dt ) ).toRotationMatrix( ) -
//                        ( marsRotationModel->getRotationToBaseFrame( (finalTime - timeStep) - dt ) ).toRotationMatrix( ) ) / ( 2*dt );
//            std::cout << "numerical Rotation Matrix derivative" << std::endl;
//            std::cout << numericalRotationMatrixderivative << std::endl;

//            Eigen::MatrixXd DerivativeDifference =
//                    marsRotationModel->getDerivativeOfRotationToBaseFrame( (finalTime - timeStep) ) - numericalRotationMatrixderivative;
//            std::cout << "rotation Matrix derivative difference" << std::endl;

//            std::cout << DerivativeDifference << std::endl;

        }

        BOOST_AUTO_TEST_SUITE_END( )

    }

}
