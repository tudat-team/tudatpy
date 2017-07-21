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

#include "keys.h"

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

const std::string Keys::Body::useDefaultSettings = "useDefaultSettings";
const std::string Keys::Body::mass = "mass";
const std::string Keys::Body::referenceArea = "referenceArea";

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
const std::string Keys::Body::Atmosphere::file = "file";
const std::string Keys::Body::Atmosphere::spaceWeatherFile = "spaceWeatherFile";

////// Body::Ephemeris
const std::string Keys::Body::ephemeris = "ephemeris";
const std::string Keys::Body::Ephemeris::type = "type";
const std::string Keys::Body::Ephemeris::frameOrigin = "frameOrigin";
const std::string Keys::Body::Ephemeris::frameOrientation = "frameOrientation";
const std::string Keys::Body::Ephemeris::makeMultiArc = "makeMultiArc";
const std::string Keys::Body::Ephemeris::correctForStellarAbberation = "correctForStellarAbberation";
const std::string Keys::Body::Ephemeris::correctForLightTimeAbberation = "correctForLightTimeAbberation";
const std::string Keys::Body::Ephemeris::convergeLighTimeAbberation = "convergeLighTimeAbberation";
const std::string Keys::Body::Ephemeris::initialTime = "initialTime";
const std::string Keys::Body::Ephemeris::finalTime = "finalTime";
const std::string Keys::Body::Ephemeris::timeStep = "timeStep";
const std::string Keys::Body::Ephemeris::interpolator = "interpolator";
const std::string Keys::Body::Ephemeris::useLongDoubleStates = "useLongDoubleStates";
const std::string Keys::Body::Ephemeris::bodyIdentifier = "bodyIdentifier";
const std::string Keys::Body::Ephemeris::useCircularCoplanarApproximation = "useCircularCoplanarApproximation";
const std::string Keys::Body::Ephemeris::constantState = "constantState";
// const std::string Keys::Body::Ephemeris::customStateFunction = "customStateFunction";
const std::string Keys::Body::Ephemeris::initialStateInKeplerianElements = "initialStateInKeplerianElements";
const std::string Keys::Body::Ephemeris::epochOfInitialState = "epochOfInitialState";
const std::string Keys::Body::Ephemeris::centralBodyGravitationalParameter = "centralBodyGravitationalParameter";
const std::string Keys::Body::Ephemeris::rootFinderAbsoluteTolerance = "rootFinderAbsoluteTolerance";
const std::string Keys::Body::Ephemeris::rootFinderMaximumNumberOfIterations = "rootFinderMaximumNumberOfIterations";
const std::string Keys::Body::Ephemeris::bodyStateHistory = "bodyStateHistory";

////// Body::GravityField
const std::string Keys::Body::gravityField = "gravityField";
const std::string Keys::Body::GravityField::type = "type";
const std::string Keys::Body::GravityField::gravitationalParameter = "gravitationalParameter";
const std::string Keys::Body::GravityField::referenceRadius = "referenceRadius";
const std::string Keys::Body::GravityField::cosineCoefficients = "cosineCoefficients";
const std::string Keys::Body::GravityField::sineCoefficients = "sineCoefficients";
const std::string Keys::Body::GravityField::associatedReferenceFrame = "associatedReferenceFrame";
const std::string Keys::Body::GravityField::file = "file";
const std::string Keys::Body::GravityField::maximumDegree = "maximumDegree";
const std::string Keys::Body::GravityField::maximumOrder = "maximumOrder";
const std::string Keys::Body::GravityField::gravitationalParameterIndex = "gravitationalParameterIndex";
const std::string Keys::Body::GravityField::referenceRadiusIndex = "referenceRadiusIndex";

////// Body::RadiationPressure
const std::string Keys::Body::radiationPressure = "radiationPressure";
const std::string Keys::Body::RadiationPressure::type = "type";
const std::string Keys::Body::RadiationPressure::referenceArea = "referenceArea";
const std::string Keys::Body::RadiationPressure::radiationPressureCoefficient = "radiationPressureCoefficient";
const std::string Keys::Body::RadiationPressure::ocultingBodies = "ocultingBodies";

////// Body::RotationModel
const std::string Keys::Body::rotationModel = "rotationModel";
const std::string Keys::Body::RotationModel::type = "type";
const std::string Keys::Body::RotationModel::originalFrame = "originalFrame";
const std::string Keys::Body::RotationModel::targetFrame = "targetFrame";
const std::string Keys::Body::RotationModel::initialOrientation = "initialOrientation";
const std::string Keys::Body::RotationModel::initialTime = "initialTime";
const std::string Keys::Body::RotationModel::rotationRate = "rotationRate";

////// Body::ShapeModel
const std::string Keys::Body::shapeModel = "shapeModel";
const std::string Keys::Body::ShapeModel::type = "type";
const std::string Keys::Body::ShapeModel::radius = "radius";
const std::string Keys::Body::ShapeModel::equatorialRadius = "equatorialRadius";
const std::string Keys::Body::ShapeModel::flattening = "flattening";


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
