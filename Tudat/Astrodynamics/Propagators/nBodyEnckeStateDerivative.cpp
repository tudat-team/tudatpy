/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Propagators/nBodyEnckeStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to remove the central gravity acceleration from an AccelerationMap
std::vector< boost::function< double( ) > > removeCentralGravityAccelerations(
        const std::vector< std::string >& centralBodies, const std::vector< std::string >& bodiesToIntegrate,
        basic_astrodynamics::AccelerationMap& accelerationModelsPerBody )
{
    using namespace basic_astrodynamics;
    using namespace gravitation;

    std::vector< boost::function< double( ) > > centralBodyGravitationalParameters;
    centralBodyGravitationalParameters.resize( bodiesToIntegrate.size( ) );

    std::vector< boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > > listOfAccelerations;

    // Iterate over all central bodies.
    for( unsigned int i = 0; i < centralBodies.size( ); i++ )
    {
        // Check if current central body is exerting any accelerations on current body.
        if( accelerationModelsPerBody[ bodiesToIntegrate.at( i ) ].count( centralBodies.at( i ) ) == 0 )
        {
            std::string errorMessage =
                    "Error, cannot remove central point gravity of body " + bodiesToIntegrate.at( i ) +
                    " with central body " + centralBodies.at( i ) + " no accelerations due to requested central body found.";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            // Find central acceleration to central body, i.e. iterate over list of accelerations and find
            // central acceleration candiates (should be 1)
            int lastCandidate = -1;
            int numberOfCandidates = 0;
            bool isLastCandidateSphericalHarmonic = 0;
            listOfAccelerations =
                    accelerationModelsPerBody[ bodiesToIntegrate.at( i ) ][ centralBodies.at( i ) ];

            for( unsigned int j = 0; j < listOfAccelerations.size( ); j++ )
            {
                // Get type of current acceleration.
                AvailableAcceleration currentAccelerationType = getAccelerationModelType( listOfAccelerations[ j ] );

                // If central gravity, set as central acceleration candidate.
                if( currentAccelerationType == central_gravity )
                {
                    numberOfCandidates++;
                    lastCandidate = j;
                }
                else if( currentAccelerationType == spherical_harmonic_gravity )
                {
                    isLastCandidateSphericalHarmonic = 1;
                    numberOfCandidates++;
                    lastCandidate = j;
                }
                else if( ( currentAccelerationType == third_body_central_gravity ) ||
                         ( currentAccelerationType == third_body_spherical_harmonic_gravity ) )
                {
                    std::string errorMessage =
                            "Error when removing central body point gravity term, removal of 3rd body accelerations (of " +
                            centralBodies.at( i ) +
                            " on " + bodiesToIntegrate.at( i ) + ",) not yet supported";
                    throw std::runtime_error( errorMessage );
                }
            }

            // If no or multiple central acceleration candidates were found, give error.
            if( numberOfCandidates == 0 )
            {
                std::string errorMessage =
                        "Error when removing central body point gravity term on body " +
                        bodiesToIntegrate.at( i ) +
                        " with central body " + centralBodies.at( i ) + ", no central gravity found.";
                throw std::runtime_error( errorMessage );

            }
            else if( numberOfCandidates != 1 )
            {
                std::string errorMessage =
                        "Error when removing central body point gravity term on body " +
                        bodiesToIntegrate.at( i ) +
                        " with central body " + centralBodies.at( i ) + ", multiple central gravities found.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                if( !isLastCandidateSphericalHarmonic )
                {
                    // Set central body gravitational parameter (used for Kepler orbit propagation)
                    centralBodyGravitationalParameters.at( i ) =
                            boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                                listOfAccelerations.at( lastCandidate ) )->getGravitationalParameterFunction( );

                    // Remove central acceleration from list of accelerations that are evaluated at each time step.
                    listOfAccelerations.erase( listOfAccelerations.begin( ) + lastCandidate,
                                               listOfAccelerations.begin( ) + lastCandidate + 1 );
                    accelerationModelsPerBody[ bodiesToIntegrate.at( i ) ][ centralBodies.at( i ) ] =
                            listOfAccelerations;
                }
                else
                {
                    boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > originalAcceleration =
                            boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                                listOfAccelerations.at( lastCandidate ) );
                    centralBodyGravitationalParameters.at( i ) = originalAcceleration->getGravitationalParameterFunction( );

                    // Create 'additional' acceleration model which subtracts central gravity term.
                    accelerationModelsPerBody[ bodiesToIntegrate.at( i ) ][ centralBodies.at( i ) ].push_back(
                                boost::make_shared< CentralGravitationalAccelerationModel3d >
                                ( originalAcceleration->getStateFunctionOfBodyExertingAcceleration( ),
                                  originalAcceleration->getGravitationalParameterFunction( ),
                                  originalAcceleration->getStateFunctionOfBodyUndergoingAcceleration( ) ) );
                }
            }
        }
    }
    return centralBodyGravitationalParameters;
}

}

}
