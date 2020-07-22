/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/simulation/simulation.h>
#include <tudat/math/root_finders/secantRootFinder.h>

#include "propagationTargeting.h"

using namespace tudat;

PropagationTargetingProblem::PropagationTargetingProblem(const double altitudeOfPerigee,
        const double altitudeOfApogee, const double altitudeOfTarget, const double longitudeOfTarget,
        const std::shared_ptr< propagators::DependentVariableSaveSettings> dependentVariablesToSave,
        const bool useExtendedDynamics) :
    altitudeOfPerigee_( altitudeOfPerigee ), altitudeOfApogee_( altitudeOfApogee ),
    altitudeOfTarget_( altitudeOfTarget ), longitudeOfTarget_( longitudeOfTarget ), dependentVariablesToSave_( dependentVariablesToSave ),
    useExtendedDynamics_( useExtendedDynamics )
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::input_output;

    // Definition of the orbit
    earthRadius_ = spice_interface::getAverageRadius( "Earth" );
    radiusOfPerigee_ =  earthRadius_ + altitudeOfPerigee_;
    radiusOfApogee_ = earthRadius_ + altitudeOfApogee_;
    earthGravitationalParameter_ = spice_interface::getBodyGravitationalParameter( "Earth" );
    semiMajorAxis_ = (radiusOfApogee_ + radiusOfPerigee_)/2.0;

    //Integration time: half a orbit
    simulationStartEpoch_ = 0.0;
    simulationEndEpoch_ = 1.2 * mathematical_constants::PI *
            std::sqrt(pow(semiMajorAxis_,3)/earthGravitationalParameter_);

    // Create the body Earth from Spice interface
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
    if( useExtendedDynamics_ )
    {
        bodySettings =
                getDefaultBodySettings( {"Earth", "Moon", "Sun"}, simulationStartEpoch_ - 3600.0, simulationEndEpoch_ + 3600.0 );
        bodySettings[ "Moon" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        bodySettings[ "Moon" ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( "J2000" );
    }
    else
    {
        bodySettings =
                getDefaultBodySettings( {"Earth"} );
    }
    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );

    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->ephemerisSettings->resetFrameOrientation( "J2000" );


    //Create bodyMap and add the satellite as an empty body
    bodyMap_ = simulation_setup::createBodies( bodySettings );
    bodyMap_["Satellite"] = std::make_shared<Body>();
    setGlobalFrameBodyEphemerides( bodyMap_, "Earth", "J2000" );
}


std::vector<double> PropagationTargetingProblem::fitness(const std::vector<double> &x) const
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::input_output;

    const double earthRotationRate = 2.0 * mathematical_constants::PI / physical_constants::SIDEREAL_DAY;
    const double fixedStepSize = 30.0;

    //Define position of the target at 35000 km from Earth at 30 deg latitude
    Eigen::Vector3d target;
    target[0] = (earthRadius_ + altitudeOfTarget_) * cos(longitudeOfTarget_*mathematical_constants::PI/180);
    target[1] = (earthRadius_ + altitudeOfTarget_) * sin(longitudeOfTarget_*mathematical_constants::PI/180);
    target[2] = 0.0;

    //Define initial position of satellite at the perigee
    Eigen::Vector6d initialKeplerElements;
    initialKeplerElements[ semiMajorAxisIndex ] = semiMajorAxis_;
    initialKeplerElements[ eccentricityIndex ] = (radiusOfApogee_ - radiusOfPerigee_)/(radiusOfApogee_ + radiusOfPerigee_);
    initialKeplerElements[ inclinationIndex ] = 35.0 * mathematical_constants::PI/180.0;
    initialKeplerElements[ argumentOfPeriapsisIndex ] = x[0] * mathematical_constants::PI/180.0;
    initialKeplerElements[ longitudeOfAscendingNodeIndex ] = x[1] * mathematical_constants::PI/180.0;
    initialKeplerElements[ trueAnomalyIndex ] = 0.0;

    const Eigen::Vector6d systemInitialState = convertKeplerianToCartesianElements(
                initialKeplerElements, earthGravitationalParameter_ );

    //Setup simulation. Simple Keplerian orbit (only central-gravity of Earth)
    std::vector< std::string > bodiesToPropagate = { "Satellite" };
    std::vector< std::string > centralBodies = { "Earth" };
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    if( useExtendedDynamics_ )
    {
        accelerationsOfSatellite[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
                                                           2, 2 ) );
        accelerationsOfSatellite[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                          point_mass_gravity ) );
        accelerationsOfSatellite[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                         point_mass_gravity ) );
        accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    }
    else
    {
        accelerationsOfSatellite[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                           point_mass_gravity ) );
        accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    }
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, bodiesToPropagate, centralBodies );

    //Setup propagator (cowell) and integrator (RK4 fixed stepsize)
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch_, cowell, dependentVariablesToSave_ );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch_, fixedStepSize );

    //Start simulation
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap_, integratorSettings, propagatorSettings, true, false, false );

    //Retrieve results
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    previousStateHistory_ = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    previousFinalState_ = previousStateHistory_.rbegin( )->second;
    previousDependentVariablesHistory_ = dynamicsSimulator.getDependentVariableHistory();
    previousDependentVariablesFinalValues_ = previousDependentVariablesHistory_.rbegin()->second;


    //Find minimum distance from target
    Eigen::Vector3d separationFromTarget =
            Eigen::Quaterniond( Eigen::AngleAxisd(
                                    -earthRotationRate * integrationResult.begin( )->first, Eigen::Vector3d::UnitZ( ) ) ) *
            ( integrationResult.begin( )->second.segment( 0, 3 ) )- target;

    double bestDistanceFromTarget = separationFromTarget.norm( );
    double timeForBestDistanceFromTarget = integrationResult.begin( )->first;

    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
         stateIterator != integrationResult.end( ); stateIterator++ )
    {
        separationFromTarget = Eigen::Quaterniond( Eigen::AngleAxisd(
                                                       -earthRotationRate * stateIterator->first, Eigen::Vector3d::UnitZ( ) ) ) *
                stateIterator->second.segment( 0, 3 ) - target;
        const double distanceFromTarget = separationFromTarget.norm( );

        if( distanceFromTarget < bestDistanceFromTarget )
        {
            bestDistanceFromTarget = distanceFromTarget;
            timeForBestDistanceFromTarget = stateIterator->first;
        }

    }

    std::vector< double > output = {bestDistanceFromTarget} ;

    return output;


}


std::pair<std::vector<double>, std::vector<double>> PropagationTargetingProblem::get_bounds() const
{
    return {{0.0, 0.0}, {360.0, 180.0}};
}
