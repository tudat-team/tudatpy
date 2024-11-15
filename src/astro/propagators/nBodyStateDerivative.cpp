/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/propagators/nBodyStateDerivative.h"

namespace tudat
{

namespace propagators
{

std::string getTranslationalPropagatorName( const TranslationalPropagatorType propagatorType )
{
    std::string propagatorName;
    switch( propagatorType )
    {
        case cowell:
            propagatorName = "Cowell";
            break;
        case encke:
            propagatorName = "Encke";
            break;
        case gauss_keplerian:
            propagatorName = "Kepler elements";
            break;
        case gauss_modified_equinoctial:
            propagatorName = "Modified equinoctial elements";
            break;
        case unified_state_model_quaternions:
            propagatorName = "USM/quaternions";
            break;
        case unified_state_model_modified_rodrigues_parameters:
            propagatorName = "USM/modified Rodrigues parameters";
            break;
        case unified_state_model_exponential_map:
            propagatorName = "USM/exponential map";
            break;
        default:
            throw std::runtime_error( "Error when getting translational propagator name, " + std::to_string( static_cast< int >( propagatorType ) ) + " not found" );
            break;
    }
    return propagatorName;
}

int getTranslationalStateSize( const TranslationalPropagatorType propagatorType )
{
    int stateSize;
    switch( propagatorType )
    {
        case cowell:
            stateSize = 6;
            break;
        case encke:
            stateSize = 6;
            break;
        case gauss_keplerian:
            stateSize = 6;
            break;
        case gauss_modified_equinoctial:
            stateSize = 6;
            break;
        case unified_state_model_quaternions:
            stateSize = 7;
            break;
        case unified_state_model_modified_rodrigues_parameters:
            stateSize = 7;
            break;
        case unified_state_model_exponential_map:
            stateSize = 7;
            break;
        default:
            throw std::runtime_error( "Error when getting translational propagator size, " + std::to_string( static_cast< int >( propagatorType ) ) + " not found" );
            break;
    }
    return stateSize;
}


//! Function to remove the central gravity acceleration from an AccelerationMap
std::vector< std::function< double( ) > > removeCentralGravityAccelerations(
        const std::vector< std::string >& centralBodies, const std::vector< std::string >& bodiesToIntegrate,
        basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
        std::map< std::string, std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > >& removedAcceleration )
{
    using namespace basic_astrodynamics;
    using namespace gravitation;

    std::vector< std::function< double( ) > > centralBodyGravitationalParameters;
    centralBodyGravitationalParameters.resize( bodiesToIntegrate.size( ) );

    std::vector< std::shared_ptr< AccelerationModel< Eigen::Vector3d > > > listOfAccelerations;

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
            bool isLastCandidatePolyhedron = 0;
            listOfAccelerations =
                    accelerationModelsPerBody[ bodiesToIntegrate.at( i ) ][ centralBodies.at( i ) ];

            for( unsigned int j = 0; j < listOfAccelerations.size( ); j++ )
            {
                // Get type of current acceleration.
                AvailableAcceleration currentAccelerationType = getAccelerationModelType( listOfAccelerations[ j ] );
                // If central gravity, set as central acceleration candidate.
                if( currentAccelerationType == point_mass_gravity )
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
                else if ( currentAccelerationType == polyhedron_gravity )
                {
                    isLastCandidatePolyhedron = 1;
                    numberOfCandidates++;
                    lastCandidate = j;
                }
                else if( ( currentAccelerationType == third_body_point_mass_gravity ) ||
                         ( currentAccelerationType == third_body_spherical_harmonic_gravity ) ||
                        ( currentAccelerationType == third_body_polyhedron_gravity ) )
                {
                    std::string errorMessage =
                            "Error when removing central body point gravity term, removal of 3rd body accelerations (of " +
                            centralBodies.at( i ) +
                            " on " + bodiesToIntegrate.at( i ) + "), not yet supported";
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
                if( !isLastCandidateSphericalHarmonic && !isLastCandidatePolyhedron )
                {
                    // Set central body gravitational parameter (used for Kepler orbit propagation)
                    removedAcceleration[ bodiesToIntegrate.at( i ) ] = std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                                listOfAccelerations.at( lastCandidate ) );
                    centralBodyGravitationalParameters.at( i ) =
                            removedAcceleration[ bodiesToIntegrate.at( i ) ]->getGravitationalParameterFunction( );

                    // Remove central acceleration from list of accelerations that are evaluated at each time step.
                    listOfAccelerations.erase( listOfAccelerations.begin( ) + lastCandidate,
                                               listOfAccelerations.begin( ) + lastCandidate + 1 );
                    accelerationModelsPerBody[ bodiesToIntegrate.at( i ) ][ centralBodies.at( i ) ] =
                            listOfAccelerations;
                }
                else if ( isLastCandidateSphericalHarmonic )
                {
                    std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > originalAcceleration =
                            std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                                listOfAccelerations.at( lastCandidate ) );
                    centralBodyGravitationalParameters.at( i ) = originalAcceleration->getGravitationalParameterFunction( );

                    // Create 'additional' acceleration model which subtracts central gravity term.
                    accelerationModelsPerBody[ bodiesToIntegrate.at( i ) ][ centralBodies.at( i ) ].push_back(
                                std::make_shared< CentralGravitationalAccelerationModel3d >
                                ( originalAcceleration->getStateFunctionOfBodyExertingAcceleration( ),
                                  originalAcceleration->getGravitationalParameterFunction( ),
                                  originalAcceleration->getStateFunctionOfBodyUndergoingAcceleration( ) ) );
                }
                else
                {
                    std::shared_ptr< PolyhedronGravitationalAccelerationModel > originalAcceleration =
                            std::dynamic_pointer_cast< PolyhedronGravitationalAccelerationModel >(
                                listOfAccelerations.at( lastCandidate ) );
                    centralBodyGravitationalParameters.at( i ) = originalAcceleration->getGravitationalParameterFunction( );

                    // Create 'additional' acceleration model which subtracts central gravity term.
                    accelerationModelsPerBody[ bodiesToIntegrate.at( i ) ][ centralBodies.at( i ) ].push_back(
                                std::make_shared< CentralGravitationalAccelerationModel3d >
                                ( originalAcceleration->getStateFunctionOfBodyExertingAcceleration( ),
                                  originalAcceleration->getGravitationalParameterFunction( ),
                                  originalAcceleration->getStateFunctionOfBodyUndergoingAcceleration( ) ) );
                }
            }
        }
    }
    return centralBodyGravitationalParameters;
}

//! Function to determine in which order the ephemerides are to be updated
std::vector< std::string > determineEphemerisUpdateorder( std::vector< std::string > integratedBodies,
                                                          std::vector< std::string > centralBodies,
                                                          std::vector< std::string > ephemerisOrigins )
{

    std::vector< std::string > updateOrder;

    // Declare variables
    bool isFinished = 0;
    int currentIndex = 0;
    int counter = 0;
    std::vector< std::string >::const_iterator centralBodyIterator;
    std::vector< std::string >::const_iterator ephemerisOriginIterator;

    // Continue iterating until all integratedBodies have been handled.
    while( !isFinished )
    {
        // Check if current central body or ephemeris origin is integratedBodies
        centralBodyIterator = std::find(
                    integratedBodies.begin( ), integratedBodies.end( ),
                    centralBodies.at( currentIndex ) );
        ephemerisOriginIterator = std::find(
                    integratedBodies.begin( ), integratedBodies.end( ),
                    ephemerisOrigins.at( currentIndex ) );

        // If neither is in list, there is no dependency, and current body can be added to update
        // list.
        if( centralBodyIterator == integratedBodies.end( )
            && ephemerisOriginIterator == integratedBodies.end( ) )
        {
            // Add to list.
            updateOrder.push_back( integratedBodies.at( currentIndex ) );

            // Remove from list of handled bodies.
            integratedBodies.erase( integratedBodies.begin( ) + currentIndex );
            centralBodies.erase( centralBodies.begin( ) + currentIndex );
            ephemerisOrigins.erase( ephemerisOrigins.begin( ) + currentIndex );

            // Handle first entry at next iteration.
            currentIndex = 0;

            // Check if any bodies left.
            if( integratedBodies.size( ) == 0 )
            {
                isFinished = 1;
            }
        }
        else
        {
            // If both are found in list, start next iteration at whichever is first in list.
            if( centralBodyIterator != integratedBodies.end( )
                && ephemerisOriginIterator != integratedBodies.end( ) )
            {

                currentIndex = std::min(
                            std::distance(
                                integratedBodies.begin( ),
                                std::find( integratedBodies.begin( ), integratedBodies.end( ),
                                           centralBodies.at( currentIndex ) ) ),
                            std::distance(
                                integratedBodies.begin( ),
                                std::find( integratedBodies.begin( ), integratedBodies.end( ),
                                           ephemerisOrigins.at( currentIndex ) ) ));

            }
            // If only central body is found, start with this central body in next iteration
            else if( centralBodyIterator != integratedBodies.end( ) )
            {
                currentIndex = std::distance(
                            integratedBodies.begin( ),
                            std::find( integratedBodies.begin( ), integratedBodies.end( ),
                                       centralBodies.at( currentIndex ) ) );
            }
            // If only ephemeris origin is found, start with this ephemeris origin in next iteration
            else if( ephemerisOriginIterator != integratedBodies.end( ) )
            {
                currentIndex = std::distance(
                            integratedBodies.begin( ),
                            std::find( integratedBodies.begin( ), integratedBodies.end( ),
                                       ephemerisOrigins.at( currentIndex ) ) );
            }
        }

        // Check to break circular dependency that occurs for inadmissible input data.
        counter++;
        if( counter > 10000 )
        {
            throw std::runtime_error( "Warning, ephemeris update order determination now at iteration " +
                                      std::to_string( counter ) );
        }
    }

    return updateOrder;
}

//template class NBodyStateDerivative< double, double >;

}

}
