/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/createLightTimeCorrection.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

boost::shared_ptr< LightTimeCorrection > createLightTimeCorrections(
        const boost::shared_ptr< LightTimeCorrectionSettings > correctionSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::pair< std::string, std::string >& transmitter,
        const std::pair< std::string, std::string >& receiver )
{

    using namespace tudat::ephemerides;
    using namespace tudat::gravitation;

    boost::shared_ptr< LightTimeCorrection > lightTimeCorrection;

    switch( correctionSettings->getCorrectionType( ) )
    {

    case first_order_relativistic:
    {
        if( boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings ) != NULL )
        {
            std::vector< std::string > perturbingBodies =
                    boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( correctionSettings )->
                    getPerturbingBodies( );

            std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > perturbingBodyStateFunctions;
            std::vector< boost::function< double( ) > > perturbingBodyGravitationalParameterFunctions;

            for( unsigned int i = 0; i < perturbingBodies.size( ); i++ )
            {
                if( bodyMap.count( perturbingBodies[ i ] ) == 0 )
                {
                    std::cerr<<"Error when making 1st order relativistic light time correction, could not find body "<<perturbingBodies[ i ]<<std::endl;
                }
                else
                {
                    perturbingBodyStateFunctions.push_back( boost::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris,
                                                                         bodyMap.at( perturbingBodies[ i ] ), _1 ) );
                        perturbingBodyGravitationalParameterFunctions.push_back(
                                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                 bodyMap.at( perturbingBodies[ i ] )->
                                                 getGravityFieldModel( ) ) );
                }
            }

            lightTimeCorrection = boost::make_shared< FirstOrderLightTimeCorrectionCalculator >(
                        perturbingBodyStateFunctions, perturbingBodyGravitationalParameterFunctions, perturbingBodies );

        }
        else
        {
            std::cerr<<"Error, correction settings type (1st order relativistic) does not coincide with data type."<<std::endl;
        }

        break;
    }
    default:
    {
        std::cerr<<"Error, light time correction type not recognized."<<std::endl;
        break;
    }

    }
    return lightTimeCorrection;
}

}

}
