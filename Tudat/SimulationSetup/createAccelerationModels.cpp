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
 *      150501    D. Dirkx          Ported from personal code
 *
 *    References
 *
 *    Notes
 *
 */

#include <boost/make_shared.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/SimulationSetup/accelerationModelTypes.h"
#include "Tudat/SimulationSetup/createAccelerationModels.h"

namespace tudat
{

namespace simulation_setup
{

using namespace aerodynamics;
using namespace gravitation;
using namespace basic_astrodynamics;
using namespace electro_magnetism;

//! Function to determine if a given frame is an inertial frame.
bool isFrameInertial( const std::string& frame )
{
    bool isFrameInertial_;
    if( frame == "SSB" || frame == "" || frame == "Inertial" )
    {
        isFrameInertial_ = true;
    }
    else
    {
        isFrameInertial_ = false;
    }
    return isFrameInertial_;
}


double evaluateDoubleFunctions(
        const boost::function< double( ) >& function1,
        const boost::function< double( ) >& function2 )
{
    return function1( ) + function2( );
}


//! Function to create central gravity acceleration model.
boost::shared_ptr< CentralGravitationalAccelerationModel3d > createCentralGravityAcceleratioModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const bool useCentralBodyFixedFrame )
{
    // Declare pointer to return object.
    boost::shared_ptr< CentralGravitationalAccelerationModel3d > accelerationModelPointer;

    // Check if body is endowed with a gravity field model (i.e. is capable of exerting
    // gravitation acceleration).
    if( bodyExertingAcceleration->getGravityFieldModel( ) == NULL )
    {
        throw std::runtime_error(
                    std::string( "Error, gravity field model not set when making central ") +
                    " gravitational acceleration of " + nameOfBodyExertingAcceleration + " on " +
                    nameOfBodyUndergoingAcceleration );
    }
    else
    {
        boost::function< double( ) > gravitationalParameterFunction;

        // Set correct value for gravitational parameter.
        if( useCentralBodyFixedFrame == 0  ||
                bodyUndergoingAcceleration->getGravityFieldModel( ) == NULL )
        {
            gravitationalParameterFunction =
                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                 bodyExertingAcceleration->getGravityFieldModel( ) );
        }
        else
        {
            boost::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                 bodyExertingAcceleration->getGravityFieldModel( ) );
            boost::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                 bodyUndergoingAcceleration->getGravityFieldModel( ) );
            gravitationalParameterFunction =
                    boost::bind( &evaluateDoubleFunctions,
                                 gravitationalParameterOfBodyExertingAcceleration,
                                 gravitationalParameterOfBodyUndergoingAcceleration );
        }

        // Create acceleration object.
        accelerationModelPointer =
                boost::make_shared< CentralGravitationalAccelerationModel3d >(
                    boost::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                    gravitationalParameterFunction,
                    boost::bind( &Body::getPosition, bodyExertingAcceleration ) );
    }


    return accelerationModelPointer;
}

//! Function to create spherical harmonic gravity acceleration model.
boost::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModelXd >
createSphericalHarmonicsGravityAcceleration(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const bool useCentralBodyFixedFrame )
{
    // Declare pointer to return object
    boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModelXd > accelerationModel;

    // Dynamic cast acceleration settings to required type and check consistency.
    boost::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicsSettings =
            boost::dynamic_pointer_cast< SphericalHarmonicAccelerationSettings >(
                accelerationSettings );
    if( sphericalHarmonicsSettings == NULL )
    {
        throw std::runtime_error(
                    std::string( "Error, acceleration settings inconsistent ") +
                    " making sh gravitational acceleration of " + nameOfBodyExertingAcceleration +
                    " on " + nameOfBodyUndergoingAcceleration );
    }
    else
    {
        // Get pointer to gravity field of central body and cast to required type.
        boost::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) );
        if( sphericalHarmonicsGravityField == NULL )
        {
            throw std::runtime_error(
                        std::string( "Error, spherical harmonic gravity field model not set when ")
                        + " making sh gravitational acceleration of " +
                        nameOfBodyExertingAcceleration +
                        " on " + nameOfBodyUndergoingAcceleration );
        }
        else
        {
            boost::function< double( ) > gravitationalParameterFunction;

            // Check if mutual acceleration is to be used.
            if( useCentralBodyFixedFrame == false ||
                    bodyUndergoingAcceleration->getGravityFieldModel( ) == NULL )
            {
                gravitationalParameterFunction =
                        boost::bind( &SphericalHarmonicsGravityField::getGravitationalParameter,
                                     sphericalHarmonicsGravityField );
            }
            else
            {
                // Create function returning summed gravitational parameter of the two bodies.
                boost::function< double( ) > gravitationalParameterOfBodyExertingAcceleration =
                        boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                     sphericalHarmonicsGravityField );
                boost::function< double( ) > gravitationalParameterOfBodyUndergoingAcceleration =
                        boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                     bodyUndergoingAcceleration->getGravityFieldModel( ) );
                gravitationalParameterFunction =
                        boost::bind( &evaluateDoubleFunctions,
                                     gravitationalParameterOfBodyExertingAcceleration,
                                     gravitationalParameterOfBodyUndergoingAcceleration );
            }

            // Create acceleration object.
            accelerationModel =
                    boost::make_shared< SphericalHarmonicsGravitationalAccelerationModelXd >
                    ( boost::bind( &Body::getPosition, bodyUndergoingAcceleration ),
                      gravitationalParameterFunction,
                      sphericalHarmonicsGravityField->getReferenceRadius( ),
                      boost::bind( &SphericalHarmonicsGravityField::getCosineCoefficients,
                                   sphericalHarmonicsGravityField,
                                   sphericalHarmonicsSettings->maximumDegree_,
                                   sphericalHarmonicsSettings->maximumOrder_ ),
                      boost::bind( &SphericalHarmonicsGravityField::getSineCoefficients,
                                   sphericalHarmonicsGravityField,
                                   sphericalHarmonicsSettings->maximumDegree_,
                                   sphericalHarmonicsSettings->maximumOrder_ ),
                      boost::bind( &Body::getPosition, bodyExertingAcceleration ) );
        }
    }
    return accelerationModel;
}


//! Function to create a third body central gravity acceleration model.
boost::shared_ptr< gravitation::ThirdBodyCentralGravityAcceleration >
createThirdBodyCentralGravityAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const std::string& nameOfCentralBody )
{
    // Declare pointer to return object.
    boost::shared_ptr< ThirdBodyCentralGravityAcceleration > accelerationModelPointer;

    // Create acceleration object.
    accelerationModelPointer =  boost::make_shared< ThirdBodyCentralGravityAcceleration >(
                boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                    createCentralGravityAcceleratioModel( bodyUndergoingAcceleration,
                                                          bodyExertingAcceleration,
                                                          nameOfBodyUndergoingAcceleration,
                                                          nameOfBodyExertingAcceleration, 0 ) ),
                boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                    createCentralGravityAcceleratioModel( centralBody, bodyExertingAcceleration,
                                                          nameOfCentralBody,
                                                          nameOfBodyExertingAcceleration, 0 ) ) );

    return accelerationModelPointer;
}

//! Function to create acceleration model object.
boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > createAccelerationModel(
        const boost::shared_ptr< Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< Body > bodyExertingAcceleration,
        const boost::shared_ptr< AccelerationSettings > accelerationSettings,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfCentralBody )
{
    // Declare pointer to return object.
    boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelPointer;

    // Switch to call correct acceleration model type factory function.
    switch( accelerationSettings->accelerationType_ )
    {
    case central_gravity:
        if( nameOfCentralBody == nameOfBodyExertingAcceleration ||
                isFrameInertial( nameOfCentralBody ) )
        {
            bool useCentralBodyFixedFrame = 0;
            if( nameOfCentralBody == nameOfBodyExertingAcceleration )
            {
                useCentralBodyFixedFrame = 1;
            }

            accelerationModelPointer = createCentralGravityAcceleratioModel(
                        bodyUndergoingAcceleration,
                        bodyExertingAcceleration,
                        nameOfBodyUndergoingAcceleration,
                        nameOfBodyExertingAcceleration, useCentralBodyFixedFrame );
        }
        else
        {

            accelerationModelPointer = createThirdBodyCentralGravityAccelerationModel(
                        bodyUndergoingAcceleration,
                        bodyExertingAcceleration,
                        centralBody,
                        nameOfBodyUndergoingAcceleration,
                        nameOfBodyExertingAcceleration,
                        nameOfCentralBody );
        }
        break;
    case spherical_harmonic_gravity:
        if( nameOfCentralBody == nameOfBodyExertingAcceleration ||
                isFrameInertial( nameOfCentralBody ) )
        {
            bool useCentralBodyFixedFrame = 0;
            if( nameOfCentralBody == nameOfBodyExertingAcceleration )
            {
                useCentralBodyFixedFrame = 1;
            }

            accelerationModelPointer = createSphericalHarmonicsGravityAcceleration(
                        bodyUndergoingAcceleration,
                        bodyExertingAcceleration,
                        nameOfBodyUndergoingAcceleration,
                        nameOfBodyExertingAcceleration,
                        accelerationSettings, useCentralBodyFixedFrame );

        }
        else
        {
            throw std::runtime_error(
                        "Error, cannot yet make third body spherical harmonic acceleration." );

        }
        break;
    default:
        throw std::runtime_error(
                    std::string( "Error, acceleration model ") +
                    boost::lexical_cast< std::string >( accelerationSettings->accelerationType_ ) +
                    " not recognized when making acceleration model of" +
                    nameOfBodyExertingAcceleration + " on " +
                    nameOfBodyUndergoingAcceleration );
        break;
    }
    return accelerationModelPointer;
}

//! Function to create a set of acceleration models from a map of bodies and acceleration model
//! types.
AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::map< std::string, std::string >& centralBodies )
{
    // Declare return map.
    AccelerationMap accelerationModelMap;

    // Iterate over all bodies which are undergoing acceleration
    for( SelectedAccelerationMap::const_iterator bodyIterator =
         selectedAccelerationPerBody.begin( ); bodyIterator != selectedAccelerationPerBody.end( );
         bodyIterator++ )
    {
        boost::shared_ptr< Body > currentCentralBody;

        // Retrieve name of body undergoing acceleration.
        std::string bodyUndergoingAcceleration = bodyIterator->first;

        // Retrieve name of current central body.
        std::string currentCentralBodyName = centralBodies.at( bodyUndergoingAcceleration );

        if( !isFrameInertial( currentCentralBodyName ) )
        {
            if( bodyMap.count( currentCentralBodyName ) == 0 )
            {
                throw std::runtime_error(
                            std::string( "Error, could not find non-inertial central body ") +
                            currentCentralBodyName + " of " + bodyUndergoingAcceleration +
                            " when making acceleration model." );
            }
            else
            {
                currentCentralBody = bodyMap.at( currentCentralBodyName );
            }
        }

        // Check if body undergoing acceleration is included in bodyMap
        if( bodyMap.count( bodyUndergoingAcceleration ) ==  0 )
        {
            throw std::runtime_error(
                        std::string( "Error when making acceleration models, requested forces" ) +
                                     "acting on body " + bodyUndergoingAcceleration  +
                                     ", but no such body found in map of bodies" );
        }

        // Declare map of acceleration models acting on current body.
        SingleBodyAccelerationMap mapOfAccelerationsForBody;

        // Retrieve list of required acceleration model types and bodies exerting accelerationd on
        // current body.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > >
                accelerationsForBody = bodyIterator->second;

        // Iterate over all bodies exerting an acceleration
        for( std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > >::
             iterator body2Iterator = accelerationsForBody.begin( );
             body2Iterator != accelerationsForBody.end( ); body2Iterator++ )
        {
            // Retrieve name of body exerting acceleration.
            std::string bodyExertingAcceleration = body2Iterator->first;

            // Check if body exerting acceleration is included in bodyMap
            if( bodyMap.count( bodyExertingAcceleration ) ==  0 )
            {
                throw std::runtime_error(
                            std::string( "Error when making acceleration models, requested forces ")
                            + "acting on body " + bodyUndergoingAcceleration  + " due to body " +
                            bodyExertingAcceleration +
                            ", but no such body found in map of bodies" );
            }

            // Retrieve list of accelerations due to current body.
            std::vector< boost::shared_ptr< AccelerationSettings > > accelerationList =
                    body2Iterator->second;

            for( unsigned int i = 0; i < accelerationList.size( ); i++ )
            {
                // Create acceleration model.
                mapOfAccelerationsForBody[ bodyExertingAcceleration ].push_back(
                            createAccelerationModel( bodyMap.at( bodyUndergoingAcceleration ),
                                                     bodyMap.at( bodyExertingAcceleration ),
                                                     accelerationList.at( i ),
                                                     bodyUndergoingAcceleration,
                                                     bodyExertingAcceleration,
                                                     currentCentralBody,
                                                     currentCentralBodyName ) );
            }
        }

        // Put acceleration models on current body in return map.
        accelerationModelMap[ bodyUndergoingAcceleration ] = mapOfAccelerationsForBody;
    }

    return accelerationModelMap;
}

}

}
