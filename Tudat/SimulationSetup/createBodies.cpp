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
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/SimulationSetup/createBodies.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;
using namespace gravitation;

std::vector< std::pair< std::string, boost::shared_ptr< BodySettings > > > determineBodyCreationOrder(
        const std::map< std::string, boost::shared_ptr< BodySettings > >& bodySettings )
{
    std::vector< std::pair< std::string, boost::shared_ptr< BodySettings > > > outputVector;

    std::map< std::string, std::vector< boost::shared_ptr< GravityFieldVariationSettings > > > gravityFieldVariationSettingsList;

    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator bodyIterator = bodySettings.begin( );
         bodyIterator != bodySettings.end( ); bodyIterator ++ )
    {
        if( bodyIterator->second->gravityFieldVariationSettings.size( ) != 0 )
        {
            gravityFieldVariationSettingsList[ bodyIterator->first ] =
                    bodyIterator->second->gravityFieldVariationSettings;
        }
        outputVector.push_back( std::make_pair( bodyIterator->first, bodyIterator->second ) );
    }

    if( gravityFieldVariationSettingsList.size( ) != 0 )
    {
        for( std::map< std::string, std::vector< boost::shared_ptr< GravityFieldVariationSettings > > >::iterator variationIterator =
             gravityFieldVariationSettingsList.begin( ); variationIterator != gravityFieldVariationSettingsList.end( ); variationIterator++ )
        {
            unsigned int currentBodyIndex = -1;

            for( unsigned int k = 0; k < variationIterator->second.size( ); k++ )
            {
                boost::shared_ptr< BasicSolidBodyGravityFieldVariationSettings > variationSettings1 =
                        boost::dynamic_pointer_cast< BasicSolidBodyGravityFieldVariationSettings >( variationIterator->second[ k ] );

                if( variationSettings1 != NULL )
                {

                    for( unsigned int i = 0; i < outputVector.size( ); i++ )
                    {
                        if( outputVector[ i ].first == variationIterator->first )
                        {
                            currentBodyIndex = i;
                        }
                    }

                    std::vector< std::string > deformingBodies = variationSettings1->getDeformingBodies( );
                    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
                    {
                        for( unsigned int j = 0; j < outputVector.size( ); j++ )
                        {
                            if( deformingBodies[ i ] == outputVector[ j ].first )
                            {
                                if( j > currentBodyIndex )
                                {
                                    std::pair< std::string, boost::shared_ptr< BodySettings > > entryToMove = outputVector[ j ];
                                    outputVector.erase( outputVector.begin( ) + j );
                                    outputVector.insert( outputVector.begin( ) + i, entryToMove );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return outputVector;
}


//! Function to create a map of bodies objects.
NamedBodyMap createBodies(
        const std::map< std::string, boost::shared_ptr< BodySettings > >& bodySettings )
{
    std::vector< std::pair< std::string, boost::shared_ptr< BodySettings > > > orderedBodySettings =
            determineBodyCreationOrder( bodySettings );

    // Declare map of bodies that is to be returned.
    NamedBodyMap bodyMap;

    // Create empty body objects.
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        bodyMap[ orderedBodySettings.at( i ).first ] = boost::make_shared< Body >( );
    }

    // Create ephemeris objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->ephemerisSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setEphemeris(
                        createBodyEphemeris( orderedBodySettings.at( i ).second->ephemerisSettings,
                                             orderedBodySettings.at( i ).first ) );
        }
    }

    // Create atmosphere model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->atmosphereSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setAtmosphereModel(
                        createAtmosphereModel( orderedBodySettings.at( i ).second->atmosphereSettings,
                                               orderedBodySettings.at( i ).first ) );
        }
    }

    // Create body shape model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->shapeModelSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setShapeModel(
                        createBodyShapeModel( orderedBodySettings.at( i ).second->shapeModelSettings,
                                              orderedBodySettings.at( i ).first ) );
        }
    }

    // Create rotation model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->rotationModelSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setRotationalEphemeris(
                        createRotationModel( orderedBodySettings.at( i ).second->rotationModelSettings,
                                             orderedBodySettings.at( i ).first ) );
        }
    }

    // Create gravity field model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->gravityFieldSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setGravityFieldModel(
                        createGravityFieldModel( orderedBodySettings.at( i ).second->gravityFieldSettings,
                                                 orderedBodySettings.at( i ).first, bodyMap,
                                                 orderedBodySettings.at( i ).second->gravityFieldVariationSettings ) );
        }
    }

    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->gravityFieldVariationSettings.size( ) > 0 )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setGravityFieldVariationSet(
                        createGravityFieldModelVariationsSet(
                            orderedBodySettings.at( i ).first, bodyMap,
                            orderedBodySettings.at( i ).second->gravityFieldVariationSettings ) );
        }
    }

    // Create aerodynamic coefficient interface objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->aerodynamicCoefficientSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface(
                            orderedBodySettings.at( i ).second->aerodynamicCoefficientSettings,
                            orderedBodySettings.at( i ).first ) );
        }
    }


    // Create radiation pressure coefficient objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        std::map< std::string, boost::shared_ptr< RadiationPressureInterfaceSettings > > radiationPressureSettings =
                orderedBodySettings.at( i ).second->radiationPressureSettings;
        for( std::map< std::string, boost::shared_ptr< RadiationPressureInterfaceSettings > >::iterator
             radiationPressureSettingsIterator = radiationPressureSettings.begin( );
             radiationPressureSettingsIterator != radiationPressureSettings.end( );
             radiationPressureSettingsIterator++ )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setRadiationPressureInterface(
                        radiationPressureSettingsIterator->first,
                        createRadiationPressureInterface(
                            radiationPressureSettingsIterator->second, orderedBodySettings.at( i ).first, bodyMap ) );
        }

    }
    return bodyMap;

}

//! Function to define the global origin and orientation of the reference frame that is to be used in the simulations.
void setGlobalFrameBodyEphemerides( const NamedBodyMap& bodyMap,
                                    const std::string& globalFrameOrigin,
                                    const std::string& globalFrameOrientation )
{
    using namespace tudat::simulation_setup;
    std::string ephemerisFrameOrigin;
    std::string ephemerisFrameOrientation;
    std::string rotationModelFrame;

    // Iterate over all bodies
    for( NamedBodyMap::const_iterator bodyIterator = bodyMap.begin( ); bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        // Check id body contains an ephemeris
        if( bodyIterator->second->getEphemeris( ) != NULL )
        {
            // Retrieve ephemeris origin
            ephemerisFrameOrigin = bodyIterator->second->getEphemeris( )->getReferenceFrameOrigin( );

            // Check if ephemeris origin differs from global origin.
            if( ephemerisFrameOrigin != globalFrameOrigin )
            {
                // Check if correction can be made
                if( bodyMap.count( ephemerisFrameOrigin ) == 0 )
                {
                    throw std::runtime_error(
                            "Error, body " + bodyIterator->first + " has ephemeris in frame " +
                            ephemerisFrameOrigin + ", but no conversion to frame " + globalFrameOrigin +
                            " can be made" );
                }
                else
                {
                    // Retrieve and set frame correction functions.
                    boost::function< basic_mathematics::Vector6d( const double& ) > correctionFunction =
                            boost::bind( &Body::getStateInBaseFrameFromEphemeris, bodyMap.at( ephemerisFrameOrigin ), _1 );
                    bodyIterator->second->setBaseFrameFunction( correctionFunction );

                    boost::function< Eigen::Matrix< long double, 6, 1 >( const double& ) > longCorrectionFunction =
                            boost::bind( &Body::getLongStateInBaseFrameFromEphemeris, bodyMap.at( ephemerisFrameOrigin ), _1 );
                    bodyIterator->second->setBaseFrameLongFunction( longCorrectionFunction );
                }
            }

            // Retrieve ephemeris orientation
            ephemerisFrameOrientation = bodyIterator->second->getEphemeris( )->getReferenceFrameOrientation( );

            // If two are not equal, throw error.
            if( ephemerisFrameOrientation != globalFrameOrientation )
            {
                throw std::runtime_error(
                            "Error, ephemeris orientation of body " + bodyIterator->first +
                            " is not the same as global orientation " + ephemerisFrameOrientation + ", " +
                            globalFrameOrientation );
            }


        }

        // Check if body has rotational ephemeris.
        if( bodyIterator->second->getRotationalEphemeris( ) != NULL )
        {
            // Check if rotational ephemeris base frame orienatation is equal to to global orientation.
            rotationModelFrame = bodyIterator->second->getRotationalEphemeris( )->getBaseFrameOrientation( );

            // Throw error if two frames are not equal.
            if( rotationModelFrame != globalFrameOrientation )
            {
                throw std::runtime_error(
                            "Error, rotation base orientation of body " + bodyIterator->first +
                            " is not the same as global orientation " + rotationModelFrame + ", " +
                            globalFrameOrientation );
            }
        }
    }

}

}

}
