/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrection.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"
#include "Tudat/Astrodynamics/Relativity/metric.h"

namespace tudat
{

namespace observation_models
{

//! Function to create object that computes a single (type of) correction to the light-time
boost::shared_ptr< LightTimeCorrection > createLightTimeCorrections(
        const boost::shared_ptr< LightTimeCorrectionSettings > correctionSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::pair< std::string, std::string >& transmitter,
        const std::pair< std::string, std::string >& receiver )
{

    using namespace tudat::ephemerides;
    using namespace tudat::gravitation;

    boost::shared_ptr< LightTimeCorrection > lightTimeCorrection;

    // Identify type of light time correction to be created.
    switch( correctionSettings->getCorrectionType( ) )
    {
    case first_order_relativistic:
    {
        // Check input consistency
        if( boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings ) != NULL )
        {
            // Retrieve list of bodies causing light time perturbation
            std::vector< std::string > perturbingBodies =
                    boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings )->
                    getPerturbingBodies( );

            std::vector< boost::function< Eigen::Vector6d( const double ) > > perturbingBodyStateFunctions;
            std::vector< boost::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

            // Retrieve mass and state functions for each perturbing body.
            for( unsigned int i = 0; i < perturbingBodies.size( ); i++ )
            {
                if( bodyMap.count( perturbingBodies[ i ] ) == 0 )
                {
                    throw std::runtime_error(
                                "Error when making 1st order relativistic light time correction, could not find body " +
                                perturbingBodies.at( i ) );
                }
                else
                {
                    // Set state function.
                    perturbingBodyStateFunctions.push_back(
                                boost::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, double >,
                                                                         bodyMap.at( perturbingBodies[ i ] ), _1 ) );

                    // Set gravitational parameter function.
                    perturbingBodyGravitationalParameterFunctions.push_back(
                                boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                             bodyMap.at( perturbingBodies[ i ] )->
                                             getGravityFieldModel( ) ) );
                }
            }

            // Create light-time correction function
            lightTimeCorrection = boost::make_shared< FirstOrderLightTimeCorrectionCalculator >(
                        perturbingBodyStateFunctions, perturbingBodyGravitationalParameterFunctions, perturbingBodies,
                        transmitter.first, receiver.first,
                        boost::bind( &relativity::PPNParameterSet::getParameterGamma, relativity::ppnParameterSet ) );

        }
        else
        {
            throw std::runtime_error(
                        "Error, correction settings type (1st order relativistic) does not coincide with data type." );
        }

        break;
    }
    default:
    {
        std::string errorMessage = "Error, light time correction type " +
                std::to_string( correctionSettings->getCorrectionType( ) ) + " not recognized.";
        throw std::runtime_error( errorMessage );

        break;
    }

    }
    return lightTimeCorrection;
}

}

}
