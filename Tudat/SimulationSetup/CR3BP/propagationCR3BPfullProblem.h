/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/createStateDerivativeModel.h"
#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"




//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
                                                std::string( "mainTestPropagationCR3BPfullProblem.cpp" ).length( ) );
    std::string outputPath = reducedPath + "SimulationOutput/";
    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

namespace tudat
{

namespace propagators
{


//!Function to transform normalized co-rotating coordinates into cartesian ones
Eigen::Vector6d convertCorotatingNormalizedToCartesianCoordinates(
        const double gravitationalParameterPrimary,
        const double gravitationalParameterSecondary,
        const double distancePrimarySecondary,
        const Eigen::Vector6d& normalizedState,
        const double normalizedTime );

//! Function to transform cartesian coordinates into co-rotating normalized ones
Eigen::Vector6d convertCartesianToCorotatingNormalizedCoordinates(
        const double gravitationalParameterPrimary,
        const double gravitationalParameterSecondary,
        const double distancePrimarySecondary,
        const Eigen::Vector6d& cartesianState,
        const double time );

//! Function to directly setup CR3BP bodyMap
simulation_setup::NamedBodyMap setupBodyMapCR3BPBodyMap(
        const double distancePrimarySecondary,
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& nameBodyToPropagate,
        const double initialTime );

//! Function to directly setup CR3BP acceleration map
basic_astrodynamics::AccelerationMap setupAccelerationMapCR3BP(
        const std::string& namePrimaryBody,
        const std::string& nameSecondaryBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::NamedBodyMap& bodyMap );



//! Function to simultaneously propagate the dynamics in the CR3BP and in the full dynamics problem
//! and compute the difference in state at the end of the propagation
Eigen::Vector6d propagateCR3BPandFullDynamicsProblem(
        const double initialTime,
        const double finalTime,
        const Eigen::Vector6d& initialState,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        simulation_setup::NamedBodyMap& bodyMap,
        const std::vector < std::string >& bodiesCR3BP );

}

}


