/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SINGLE_SATELLITE_THRUST_H
#define TUDAT_SINGLE_SATELLITE_THRUST_H

#include "tudat/io/applicationOutput.h"
#include <tudat/simulation/simulation.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::basic_mathematics;
using namespace tudat::basic_astrodynamics;

typedef Eigen::MatrixXd StateType;

class SingleSatelliteThrust {
 public:
  SingleSatelliteThrust(
      const double initial_epoch = 0,
      const double termination_epoch = 14.0 * physical_constants::JULIAN_DAY,
      const double vehicle_mass = 5.0E3,
      const double vehicle_isp = 5000.0,
      const double vehicle_thrust = 25.0,
      const std::string &vehicle_name = "Vehicle",
      const std::string &mission_body = "Earth",
      const std::vector<std::string> &spice_bodies = {"Sun", "Earth", "Moon"}) {

    this->_global_initial_epoch = initial_epoch;
    this->_global_final_epoch = termination_epoch;
    this->_vehicle_isp = vehicle_isp;
    this->_vehicle_mass = vehicle_mass;
    this->_vehicle_name = vehicle_name;
    this->_vehicle_thrust = vehicle_thrust;
    this->_mission_body = mission_body;
    this->_spice_bodies = spice_bodies;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels();

    // Create body objects.
    std::map<std::string, std::shared_ptr<BodySettings>> body_settings =
        getDefaultBodySettings(spice_bodies);

    this->_body_system = createBodies(body_settings);

    // Create vehicle objects.
    this->_body_system["Vehicle"] = std::make_shared<simulation_setup::Body>();
    this->_body_system["Vehicle"]->setConstantBodyMass(vehicle_mass);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector<std::string> bodiesToPropagate;
    std::vector<std::string> centralBodies;
    bodiesToPropagate.push_back(this->_vehicle_name);
    centralBodies.push_back(this->_mission_body);

    this->_central_bodies = centralBodies;
    this->_bodies_to_propagate = bodiesToPropagate;

    // Finalize body creation.
    setGlobalFrameBodyEphemerides(this->_body_system, "SSB", "ECLIPJ2000");

    // Define point mass gravity accelerations of system bodies.
    for (auto body : spice_bodies) {
      this->_base_accelerations_vehicle[body].push_back(
          std::make_shared<AccelerationSettings>(point_mass_gravity));
    }


    // Define propagation termination conditions (stop after 2 weeks).
    std::shared_ptr<PropagationTimeTerminationSettings> terminationSettings =
        std::make_shared<propagators::PropagationTimeTerminationSettings>(14.0 * physical_constants::JULIAN_DAY);

    this->_termination_settings = terminationSettings;
    this->reset();
  }

  void reset() {
    this->_current_initial_epoch = this->_global_initial_epoch;
    this->_current_system_initial_state = this->_global_system_initial_state;
    this->_current_vehicle_mass = this->_vehicle_mass;
  }


  void step(std::vector<double> &action, double time_step) {

    bool done = this->_current_initial_epoch + time_step > this->_global_final_epoch;
    this->_current_final_epoch = done? this->_global_final_epoch: this->_current_initial_epoch + time_step;
    double termination_time = done? this->_global_final_epoch-this->_current_initial_epoch: time_step;


    // Define propagation termination conditions (stop after 2 weeks).
    std::shared_ptr<PropagationTimeTerminationSettings> termination_settings =
        std::make_shared<propagators::PropagationTimeTerminationSettings>(termination_time);


    // actions
    double along_track_thrust_normed = action[0] < -1.0 ? -1.0 : action[0] > 1.0 ? 1.0 : action[0];

    // along-track thrust guidance
    std::shared_ptr<ThrustDirectionGuidanceSettings> along_track_thrust_direction_settings =
        std::make_shared<ThrustDirectionFromStateGuidanceSettings>(this->_mission_body,
                                                                   true,
                                                                   false);

    std::shared_ptr<ThrustMagnitudeSettings> along_track_thrust_magnitude_settings =
        std::make_shared<ConstantThrustMagnitudeSettings>(this->_vehicle_thrust * along_track_thrust_normed,
                                                          this->_vehicle_isp);

    std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings>>> accelerations_of_vehicle(this->_base_accelerations_vehicle);

    accelerations_of_vehicle[this->_vehicle_name].push_back(
        std::make_shared<ThrustAccelerationSettings>(along_track_thrust_direction_settings,
                                                     along_track_thrust_magnitude_settings));

    SelectedAccelerationMap acceleration_map;
    acceleration_map["Vehicle"] = accelerations_of_vehicle;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap acceleration_model_map = createAccelerationModelsMap(
        this->_body_system, acceleration_map, this->_bodies_to_propagate, this->_central_bodies);

    // Define settings for propagation of translational dynamics.
    std::shared_ptr<TranslationalStatePropagatorSettings<double>> translationalPropagatorSettings =
        std::make_shared<TranslationalStatePropagatorSettings<double>>(
            this->_central_bodies,
            acceleration_model_map,
            this->_bodies_to_propagate,
            this->_current_system_initial_state,
            termination_settings);

    // Create settings for propagating the mass of the vehicle.
    std::vector<std::string> bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back(this->_vehicle_name);

    Eigen::VectorXd initialBodyMasses = Eigen::VectorXd(1);
    initialBodyMasses(0) = this->_current_vehicle_mass;

    // Create mass rate models
    std::shared_ptr<MassRateModelSettings> massRateModelSettings =
        std::make_shared<FromThrustMassModelSettings>(true);
    std::map<std::string, std::shared_ptr<basic_astrodynamics::MassRateModel>> massRateModels;
    massRateModels[this->_vehicle_name] = createMassRateModel(
        this->_vehicle_name, massRateModelSettings, this->_body_system, acceleration_model_map);

    std::shared_ptr<SingleArcPropagatorSettings<double>> massPropagatorSettings =
        std::make_shared<MassPropagatorSettings<double>>(
            bodiesWithMassToPropagate, massRateModels, initialBodyMasses, termination_settings);

    // Create list of propagation settings.
    std::vector<std::shared_ptr<SingleArcPropagatorSettings<double>>> propagatorSettingsVector;
    propagatorSettingsVector.push_back(translationalPropagatorSettings);
    propagatorSettingsVector.push_back(massPropagatorSettings);

    // Create propagation settings for mass and translational dynamics concurrently
    std::shared_ptr<PropagatorSettings<double>> propagatorSettings =
        std::make_shared<MultiTypePropagatorSettings<double>>(propagatorSettingsVector, termination_settings);

    // Define integrator settings
    std::shared_ptr<IntegratorSettings<>> integratorSettings =
        std::make_shared<IntegratorSettings<>>(rungeKutta4, 0.0, 30.0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator<> dynamicsSimulator(
        this->_body_system,
        integratorSettings,
        propagatorSettings,
        true, false, false);

    // Retrieve numerical solutions for state and dependent variables
    std::map<double, Eigen::Matrix<double, Eigen::Dynamic, 1>> numericalSolution =
        dynamicsSimulator.getEquationsOfMotionNumericalSolution();

    this->_current_initial_epoch = this->_current_final_epoch;
    this->_current_vehicle_mass = 0;// get final vehicle mass
    this->_current_system_initial_state = 0; // get system final state

    return std::tuple<StateType, double, bool>(state, reward, done)
  }

 private:
  double _global_system_initial_state;
  double _global_initial_epoch;
  double _global_final_epoch;

  double _current_vehicle_mass;
  double _current_system_initial_state;
  double _current_initial_epoch;
  double _current_final_epoch;

  double _vehicle_mass;
  double _vehicle_isp;
  double _vehicle_thrust;

  std::string _vehicle_name;
  std::string _mission_body;

  std::vector<std::string> _spice_bodies;

  std::vector<std::string> _central_bodies;
  std::vector<std::string> _bodies_to_propagate;
  SystemOfBodies _body_system;
  std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings>>> _base_accelerations_vehicle;
};

#endif//TUDAT_SINGLE_SATELLITE_THRUST_H
