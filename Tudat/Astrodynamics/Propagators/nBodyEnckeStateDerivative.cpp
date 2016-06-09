#include "Tudat/Astrodynamics/Propagators/nBodyEnckeStateDerivative.h"

namespace tudat
{

namespace propagators
{


std::vector< boost::function< double( ) > > removeCentralGravityAccelerations(
        const std::vector< std::string >& centralBodies, const std::vector< std::string >& bodiesToBeIntegratedNumerically,
        basic_astrodynamics::AccelerationMap& accelerationModelsPerBody )
{
    using namespace basic_astrodynamics;
    using namespace gravitation;

    std::vector< boost::function< double( ) > > centralBodyGravitationalParameters;
    centralBodyGravitationalParameters.resize( bodiesToBeIntegratedNumerically.size( ) );

    std::vector< boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > > listOfAccelerations;

    // Iterate over all central bodies.
    for( unsigned int i = 0; i < centralBodies.size( ); i++ )
    {
        // Check if current central body is exerting any accelerations on current body.
        if( accelerationModelsPerBody[ bodiesToBeIntegratedNumerically.at( i ) ].count( centralBodies.at( i ) ) == 0 )
        {
            std::cerr<<"Warning cannot remove central point gravity of body "<<bodiesToBeIntegratedNumerically.at( i )
                    <<" with central body "<<centralBodies.at( i )<<" no accelerations due to requested central body found. "<<std::endl;
        }
        else
        {
            // Find central acceleration to central body, i.e. iterate over list of accelerations and find
            // central acceleration candiates (should be 1)
            int lastCandidate = -1;
            int numberOfCandidates = 0;
            bool isLastCandidateSphericalHarmonic = 0;
            listOfAccelerations = accelerationModelsPerBody[ bodiesToBeIntegratedNumerically.at( i ) ][ centralBodies.at( i ) ];

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
                    std::cerr<<"Error when removing central body point gravity term, removal of 3rd body accelerations "
                            <<" (of "<<centralBodies.at( i )<<" on "<<bodiesToBeIntegratedNumerically.at( i )<<") not yet supported"<<std::endl;
                }
            }

            // If no or multiple central acceleration candidates were found, give error.
            if( numberOfCandidates == 0 )
            {
                std::cerr<<"Error when removing central body point gravity term on body "<<bodiesToBeIntegratedNumerically.at( i )
                        <<" with central body "<<centralBodies.at( i )<<", no central gravity found."<<std::endl;
            }
            else if( numberOfCandidates != 1 )
            {
                std::cerr<<"Error when removing central body point gravity term on body "<<bodiesToBeIntegratedNumerically.at( i )
                        <<" with central body "<<centralBodies.at( i )<<", "<<numberOfCandidates
                       <<"central gravities found, not 1."<<std::endl;
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
                    accelerationModelsPerBody[ bodiesToBeIntegratedNumerically.at( i ) ][ centralBodies.at( i ) ] =
                            listOfAccelerations;
                }
                else
                {
                    boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModelXd > originalAcceleration =
                            boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModelXd >(
                                listOfAccelerations.at( lastCandidate ) );
                    centralBodyGravitationalParameters.at( i ) = originalAcceleration->getGravitationalParameterFunction( );

                    // Create 'additional' acceleration model which subtracts central gravity term.
                    accelerationModelsPerBody[ bodiesToBeIntegratedNumerically.at( i ) ][ centralBodies.at( i ) ].push_back(
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
