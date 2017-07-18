/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "jsonInterface.h"

namespace tudat
{

namespace json_interface
{

//! Keys recognised by json_interface.

/// Simulation
const std::string Keys::simulation = "simulation";

const std::string Keys::Simulation::startEpoch = "startEpoch";
const std::string Keys::Simulation::endEpoch = "endEpoch";
const std::string Keys::Simulation::globalFrameOrigin = "globalFrameOrigin";
const std::string Keys::Simulation::globalFrameOrientation = "globalFrameOrientation";
const std::string Keys::Simulation::spiceKernels = "spiceKernels";
const std::string Keys::Simulation::preloadSpiceData = "preloadSpiceData";


/// Body
const std::string Keys::bodies = "bodies";

////// Body::Aerodynamics
const std::string Keys::Body::aerodynamics = "aerodynamics";

const std::string Keys::Body::Aerodynamics::type = "type";
const std::string Keys::Body::Aerodynamics::referenceArea = "referenceArea";
const std::string Keys::Body::Aerodynamics::dragCoefficient = "dragCoefficient";
const std::string Keys::Body::Aerodynamics::forceCoefficients = "forceCoefficients";
const std::string Keys::Body::Aerodynamics::momentCoefficients = "momentCoefficients";
const std::string Keys::Body::Aerodynamics::areCoefficientsInAerodynamicFrame = "areCoefficientsInAerodynamicFrame";
const std::string Keys::Body::Aerodynamics::areCoefficientsInNegativeAxisDirection = "areCoefficientsInNegativeAxisDirection";

////// Body::Atmosphere
const std::string Keys::Body::atmosphere = "atmosphere";

const std::string Keys::Body::Atmosphere::type = "type";
const std::string Keys::Body::Atmosphere::densityScaleHeight = "densityScaleHeight";
const std::string Keys::Body::Atmosphere::constantTemperature = "constantTemperature";
const std::string Keys::Body::Atmosphere::densityAtZeroAltitude = "densityAtZeroAltitude";
const std::string Keys::Body::Atmosphere::specificGasConstant = "specificGasConstant";
const std::string Keys::Body::Atmosphere::atmosphereFile = "atmosphereFile";
const std::string Keys::Body::Atmosphere::spaceWeatherFile = "spaceWeatherFile";

////// Body::RadiationPressure
const std::string Keys::Body::radiationPressure = "radiationPressure";

const std::string Keys::Body::RadiationPressure::type = "type";
const std::string Keys::Body::RadiationPressure::referenceArea = "referenceArea";
const std::string Keys::Body::RadiationPressure::radiationPressureCoefficient = "radiationPressureCoefficient";
const std::string Keys::Body::RadiationPressure::ocultingBodies = "ocultingBodies";


/// Propagator
const std::string Keys::propagators = "propagators";

const std::string Keys::Propagator::integratedStateType = "integratedStateType";
const std::string Keys::Propagator::type = "type";
const std::string Keys::Propagator::centralBodies = "centralBodies";
const std::string Keys::Propagator::bodiesToPropagate = "bodiesToPropagate";
const std::string Keys::Propagator::initialStates = "initialStates";
const std::string Keys::Propagator::initialStateTypes = "initialStateTypes";

////// Propagator::Termination
const std::string Keys::Propagator::termination = "termination";


////// Propagator::Acceleration
const std::string Keys::Propagator::accelerations = "accelerations";



/// Integrator
const std::string Keys::integrator = "integrator";

const std::string Keys::Integrator::type = "type";
const std::string Keys::Integrator::initialTime = "initialTime";
const std::string Keys::Integrator::initialTimeStep = "initialTimeStep";
const std::string Keys::Integrator::saveFrequency = "saveFrequency";
const std::string Keys::Integrator::rungeKuttaCoefficientSet = "rungeKuttaCoefficientSet";
const std::string Keys::Integrator::minimumStepSize = "minimumStepSize";
const std::string Keys::Integrator::maximumStepSize = "maximumStepSize";
const std::string Keys::Integrator::relativeErrorTolerance = "relativeErrorTolerance";
const std::string Keys::Integrator::absoluteErrorTolerance = "absoluteErrorTolerance";
const std::string Keys::Integrator::safetyFactorForNextStepSize = "safetyFactorForNextStepSize";
const std::string Keys::Integrator::maximumFactorIncreaseForNextStepSize = "maximumFactorIncreaseForNextStepSize";
const std::string Keys::Integrator::minimumFactorDecreaseForNextStepSize = "minimumFactorDecreaseForNextStepSize";


/// Output
const std::string Keys::output = "output";




//! Key trees recognised by `json_interface`.

/// Simulation
const KeyTree KeyTrees::Simulation::startEpoch = { Keys::simulation, Keys::Simulation::startEpoch };
const KeyTree KeyTrees::Simulation::endEpoch = { Keys::simulation, Keys::Simulation::endEpoch };
const KeyTree KeyTrees::Simulation::globalFrameOrigin = { Keys::simulation, Keys::Simulation::globalFrameOrigin };
const KeyTree KeyTrees::Simulation::globalFrameOrientation = { Keys::simulation, Keys::Simulation::globalFrameOrientation };
const KeyTree KeyTrees::Simulation::spiceKernels = { Keys::simulation, Keys::Simulation::spiceKernels };
const KeyTree KeyTrees::Simulation::preloadSpiceData = { Keys::simulation, Keys::Simulation::preloadSpiceData };

/*
/// Integrator
const KeyTree KeyTrees::Integrator::initialTime = { Keys::integrator, Keys::Integrator::initialTime };
*/

} // namespace json_interface

} // namespace tudat
