/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      Note that the exact implementation of Newton-Raphson as root finder should be updated if
 *      someone would want to use a different root-finding technique.
 *
 *      By default the eccentricity is used as the iteration procedure. This is because in
 *      optimizing a Cassini-like trajectory, the pericenter radius had about 2-4 NaN values in
 *      100000 times the gravity assist calculation. The eccentricity iteration had no NaN values
 *      for a similar run in which 100000 gravity assist calculations were done. Also the
 *      eccentricity seemed to require slightly less iterations (does not necessarily mean it is
 *      faster or more accurate).
 *
 */

#ifndef TUDAT_CREATE_TRANSFER_TRAJECTORY_H
#define TUDAT_CREATE_TRANSFER_TRAJECTORY_H

#include <map>

#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "tudat/astro/mission_segments/transferLeg.h"
#include "tudat/astro/mission_segments/transferNode.h"
#include "tudat/astro/mission_segments/transferTrajectory.h"
#include "tudat/astro/low_thrust/shape_based/sphericalShapingLeg.h"
#include "tudat/astro/low_thrust/shape_based/hodographicShapingLeg.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{
namespace mission_segments
{

//! From  "Problem Description for the 1st ACT Competition on Global Trajectory Optimisation"
const static std::map< std::string, double > DEFAULT_MINIMUM_PERICENTERS =
{   { "Mercury", 2740000.0 },
    { "Venus", 6351800.0 },
    { "Earth" , 6678000.0 },
    { "Mars" , 3689000.0 },
    { "Jupiter" , 600000000.0 },
    { "Saturn" , 70000000.0 }
};


static std::map< TransferLegTypes, bool > legRequiresInputFromPreviousNode =
{
    { unpowered_unperturbed_leg, false },
    { dsm_position_based_leg, false },
    { dsm_velocity_based_leg, true },
    {spherical_shaping_low_thrust_leg, true},
    {hodographic_low_thrust_leg, true}
};

static std::map< TransferLegTypes, bool > legRequiresInputFromFollowingNode =
{
    { unpowered_unperturbed_leg, false },
    { dsm_position_based_leg, false },
    { dsm_velocity_based_leg, false },
    {spherical_shaping_low_thrust_leg, true},
    {hodographic_low_thrust_leg, true}
};


class TransferLegSettings
{
public:
    TransferLegSettings(
            const TransferLegTypes legType ):
        legType_( legType ){ }

    virtual ~TransferLegSettings( ) { }

    TransferLegTypes legType_;
};

class SphericalShapingLegSetting: public TransferLegSettings
{
public:
    SphericalShapingLegSetting (
            const std::shared_ptr<root_finders::RootFinderSettings>& rootFinderSettings,
            const double lowerBoundFreeCoefficient,
            const double upperBoundFreeCoefficient,
            const double initialValueFreeCoefficient,
            const double timeToAzimuthInterpolatorStepSize):
        TransferLegSettings( spherical_shaping_low_thrust_leg ),
        rootFinderSettings_(rootFinderSettings),
        lowerBoundFreeCoefficient_(lowerBoundFreeCoefficient),
        upperBoundFreeCoefficient_(upperBoundFreeCoefficient),
        initialValueFreeCoefficient_(initialValueFreeCoefficient),
        timeToAzimuthInterpolatorStepSize_(timeToAzimuthInterpolatorStepSize) { }

    const std::shared_ptr<root_finders::RootFinderSettings> rootFinderSettings_;
    const double lowerBoundFreeCoefficient_;
    const double upperBoundFreeCoefficient_;
    const double initialValueFreeCoefficient_;
    const double timeToAzimuthInterpolatorStepSize_;
};

class HodographicShapingLegSettings: public TransferLegSettings
{
public:
    HodographicShapingLegSettings(
            const shape_based_methods::HodographicBasisFunctionList& radialVelocityFunctionComponents,
            const shape_based_methods::HodographicBasisFunctionList& normalVelocityFunctionComponents,
            const shape_based_methods::HodographicBasisFunctionList& axialVelocityFunctionComponents)
        : TransferLegSettings( hodographic_low_thrust_leg ),
        radialVelocityFunctionComponents_( radialVelocityFunctionComponents ),
        normalVelocityFunctionComponents_( normalVelocityFunctionComponents ),
        axialVelocityFunctionComponents_( axialVelocityFunctionComponents ),
        numberOfFreeRadialCoefficients_( radialVelocityFunctionComponents.size( ) - 3 ),
        numberOfFreeNormalCoefficients_( normalVelocityFunctionComponents.size( ) - 3 ),
        numberOfFreeAxialCoefficients_( axialVelocityFunctionComponents.size( ) - 3 )
        { }

    const shape_based_methods::HodographicBasisFunctionList radialVelocityFunctionComponents_;
    const shape_based_methods::HodographicBasisFunctionList normalVelocityFunctionComponents_;
    const shape_based_methods::HodographicBasisFunctionList axialVelocityFunctionComponents_;

    const int numberOfFreeRadialCoefficients_;
    const int numberOfFreeNormalCoefficients_;
    const int numberOfFreeAxialCoefficients_;
};

std::shared_ptr< TransferLegSettings > dsmVelocityBasedLeg( );

std::shared_ptr< TransferLegSettings > dsmPositionBasedLeg( );

std::shared_ptr< TransferLegSettings > unpoweredLeg( );

std::shared_ptr< TransferLegSettings > sphericalShapingLeg (
        const std::shared_ptr<root_finders::RootFinderSettings>& rootFinderSettings,
        const double lowerBoundFreeCoefficient = TUDAT_NAN,
        const double upperBoundFreeCoefficient = TUDAT_NAN,
        const double initialValueFreeCoefficient = TUDAT_NAN,
        const double timeToAzimuthInterpolatorStepSize = physical_constants::JULIAN_DAY);

std::shared_ptr< TransferLegSettings > hodographicShapingLeg (
        const shape_based_methods::HodographicBasisFunctionList& radialVelocityFunctionComponents,
        const shape_based_methods::HodographicBasisFunctionList& normalVelocityFunctionComponents,
        const shape_based_methods::HodographicBasisFunctionList& axialVelocityFunctionComponents);

class TransferNodeSettings
{
public:
    TransferNodeSettings(
            const TransferNodeTypes nodeType ):
        nodeType_( nodeType ){ }

    virtual ~TransferNodeSettings( ){ }

    TransferNodeTypes nodeType_;
};

class SwingbyNodeSettings: public TransferNodeSettings
{
public:
    SwingbyNodeSettings(
            const double minimumPeriapsisRadius ):
        TransferNodeSettings( swingby ),
        minimumPeriapsisRadius_( minimumPeriapsisRadius ){ }

    double minimumPeriapsisRadius_;
};

class EscapeAndDepartureNodeSettings: public TransferNodeSettings
{
public:
    EscapeAndDepartureNodeSettings(
            const double departureSemiMajorAxis,
            const double departureEccentricity ):
        TransferNodeSettings( escape_and_departure ),
        departureSemiMajorAxis_( departureSemiMajorAxis ),
        departureEccentricity_( departureEccentricity ){ }

    double departureSemiMajorAxis_;
    double departureEccentricity_;
};

class CaptureAndInsertionNodeSettings: public TransferNodeSettings
{
public:
    CaptureAndInsertionNodeSettings(
            const double captureSemiMajorAxis,
            const double captureEccentricity ):
        TransferNodeSettings( capture_and_insertion ),
        captureSemiMajorAxis_( captureSemiMajorAxis ),
        captureEccentricity_( captureEccentricity ){ }

    TransferNodeTypes nodeType_;
    double captureSemiMajorAxis_;
    double captureEccentricity_;
};

std::shared_ptr< TransferNodeSettings > escapeAndDepartureNode(
        const double departureSemiMajorAxis,
        const double departureEccentricity );

std::shared_ptr< TransferNodeSettings > swingbyNode(
        const double minimumPeriapsisDistance = TUDAT_NAN );

std::shared_ptr< TransferNodeSettings > captureAndInsertionNode(
        const double captureSemiMajorAxis,
        const double captureEccentricity );

std::shared_ptr< TransferLeg > createTransferLeg (const simulation_setup::SystemOfBodies& bodyMap,
                                                  const std::shared_ptr< TransferLegSettings > legSettings,
                                                  const std::string& departureBodyName,
                                                  const std::string& arrivalBodyName,
                                                  const std::string& centralBodyName,
                                                  const std::shared_ptr< TransferNode > departureNode,
                                                  const std::shared_ptr< TransferNode > arrivalNode);

std::shared_ptr<TransferNode> createTransferNode(const simulation_setup::SystemOfBodies &bodyMap,
                                                 const std::shared_ptr<TransferNodeSettings> nodeSettings,
                                                 const std::string &nodeBodyName,
                                                 const std::shared_ptr<TransferLeg> incomingTransferLeg,
                                                 const std::shared_ptr<TransferLeg> outgoingTransferLeg,
                                                 const bool nodeComputesIncomingVelocity,
                                                 const bool nodeComputesOutgoingVelocity);

std::pair< std::vector< std::shared_ptr< TransferLegSettings > >,
std::vector< std::shared_ptr< TransferNodeSettings > > > getMgaTransferTrajectorySettings(
        const std::vector< std::string >& fullBodiesList,
        std::function< std::shared_ptr<TransferLegSettings>( int ) > identicalTransferLegSettingsConstructor,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

void getMgaTransferTrajectorySettings(
        std::vector< std::shared_ptr< TransferLegSettings > >& transferLegSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& transferNodeSettings,
        const std::vector< std::string >& fullBodiesList,
        std::function< std::shared_ptr<TransferLegSettings>( int ) > identicalTransferLegSettingsConstructor,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

void getMgaTransferTrajectorySettingsWithoutDsm(
        std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::string& departureBody,
        const std::string& arrivalBody,
        const std::vector< std::string >& flybyBodies,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

std::pair< std::vector< std::shared_ptr< TransferLegSettings > >,
std::vector< std::shared_ptr< TransferNodeSettings > > > getMgaTransferTrajectorySettingsWithoutDsm(
        const std::vector< std::string >& fullBodiesList,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

void getMgaTransferTrajectorySettingsWithoutDsm(
        std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::vector< std::string >& fullBodiesList,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

void getMgaTransferTrajectorySettingsWithPositionBasedDsm(
        std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::string& departureBody,
        const std::string& arrivalBody,
        const std::vector< std::string >& flybyBodies,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

std::pair< std::vector< std::shared_ptr< TransferLegSettings > >,
std::vector< std::shared_ptr< TransferNodeSettings > > > getMgaTransferTrajectorySettingsWithPositionBasedDsm(
        const std::vector< std::string >& fullBodiesList,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

void getMgaTransferTrajectorySettingsWithPositionBasedDsm(
        std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::vector< std::string >& fullBodiesList,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

void getMgaTransferTrajectorySettingsWithVelocityBasedDsm(
        std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::string& departureBody,
        const std::string& arrivalBody,
        const std::vector< std::string >& flybyBodies,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

std::pair< std::vector< std::shared_ptr< TransferLegSettings > >,
std::vector< std::shared_ptr< TransferNodeSettings > > > getMgaTransferTrajectorySettingsWithVelocityBasedDsm(
        const std::vector< std::string >& fullBodiesList,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

void getMgaTransferTrajectorySettingsWithVelocityBasedDsm(
        std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::vector< std::string >& fullBodiesList,
        const std::pair< double, double > departureOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::pair< double, double > arrivalOrbit = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS );

void getMgaTransferTrajectorySettingsWithSphericalShapingThrust (
        std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::string& departureBody,
        const std::string& arrivalBody,
        const std::vector< std::string >& flybyBodies,
        const std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings,
        const std::pair< double, double > departureOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::pair< double, double > arrivalOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const double lowerBoundFreeCoefficient = TUDAT_NAN,
        const double upperBoundFreeCoefficient = TUDAT_NAN,
        const double initialValueFreeCoefficient = TUDAT_NAN,
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS);

std::pair< std::vector< std::shared_ptr< TransferLegSettings > >,
std::vector< std::shared_ptr< TransferNodeSettings > > > getMgaTransferTrajectorySettingsWithSphericalShapingThrust (
        const std::vector< std::string >& fullBodiesList,
        const std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings,
        const std::pair< double, double > departureOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::pair< double, double > arrivalOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const double lowerBoundFreeCoefficient = TUDAT_NAN,
        const double upperBoundFreeCoefficient = TUDAT_NAN,
        const double initialValueFreeCoefficient = TUDAT_NAN,
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS);

void getMgaTransferTrajectorySettingsWithSphericalShapingThrust (
        std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::vector< std::string >& fullBodiesList,
        const std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings,
        const std::pair< double, double > departureOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::pair< double, double > arrivalOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const double lowerBoundFreeCoefficient = TUDAT_NAN,
        const double upperBoundFreeCoefficient = TUDAT_NAN,
        const double initialValueFreeCoefficient = TUDAT_NAN,
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS);

// Create leg and node settings for transfer consisting of hodographic shaping legs, assuming all legs use the same
// velocity shaping functions (i.e. constituted by the same mathematical functions, but with different coefficients)
void getMgaTransferTrajectorySettingsWithHodographicShapingThrust (
        std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::vector< std::string >& fullBodiesList,
        const std::vector< shape_based_methods::HodographicBasisFunctionList >& radialVelocityFunctionComponentsPerLeg,
        const std::vector< shape_based_methods::HodographicBasisFunctionList >& normalVelocityFunctionComponentsPerLeg,
        const std::vector< shape_based_methods::HodographicBasisFunctionList >& axialVelocityFunctionComponentsPerLeg,
        const std::pair< double, double > departureOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::pair< double, double > arrivalOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS);

void getMgaTransferTrajectorySettingsWithHodographicShapingThrust (
        std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::string& departureBody,
        const std::string& arrivalBody,
        const std::vector< std::string >& flybyBodies,
        const std::vector< shape_based_methods::HodographicBasisFunctionList >& radialVelocityFunctionComponentsPerLeg,
        const std::vector< shape_based_methods::HodographicBasisFunctionList >& normalVelocityFunctionComponentsPerLeg,
        const std::vector< shape_based_methods::HodographicBasisFunctionList >& axialVelocityFunctionComponentsPerLeg,
        const std::pair< double, double > departureOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::pair< double, double > arrivalOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS);

std::pair< std::vector< std::shared_ptr< TransferLegSettings > >,
std::vector< std::shared_ptr< TransferNodeSettings > > > getMgaTransferTrajectorySettingsWithHodographicShapingThrust (
        const std::vector< std::string >& fullBodiesList,
        const std::vector< shape_based_methods::HodographicBasisFunctionList >& radialVelocityFunctionComponentsPerLeg,
        const std::vector< shape_based_methods::HodographicBasisFunctionList >& normalVelocityFunctionComponentsPerLeg,
        const std::vector< shape_based_methods::HodographicBasisFunctionList >& axialVelocityFunctionComponentsPerLeg,
        const std::pair< double, double > departureOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::pair< double, double > arrivalOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS);

std::pair< std::vector< std::shared_ptr< TransferLegSettings > >,
std::vector< std::shared_ptr< TransferNodeSettings > > > getMgaTransferTrajectorySettingsWithHodographicShapingThrust (
        const std::vector< std::string >& fullBodiesList,
        const std::vector< double >& timeOfFlightPerLeg,
        const std::vector< double >& numberOfRevolutionsPerLeg,
        const std::pair< double, double > departureOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::pair< double, double > arrivalOrbit = std::make_pair(TUDAT_NAN, TUDAT_NAN),
        const std::map< std::string, double > minimumPericenterRadii = DEFAULT_MINIMUM_PERICENTERS);

std::shared_ptr< TransferTrajectory > createTransferTrajectory(
        const simulation_setup::SystemOfBodies& bodyMap,
        const std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        const std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::vector< std::string >& nodeIds,
        const std::string& centralBody);

void getParameterVectorDecompositionIndices(
        const std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        const std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        std::vector< std::pair< int, int > >& legParameterIndices,
        std::vector< std::pair< int, int > >& nodeParameterIndices );

void printTransferParameterDefinition(
        const std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        const std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings );

} // namespace mission_segments

} // namespace tudat

#endif // TUDAT_CREATE_TRANSFER_TRAJECTORY_H
