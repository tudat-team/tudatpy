/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <boost/regex.hpp>

#include "tudat/interface/json/support/utilities.h"
#include "tudat/interface/json/support/keys.h"

namespace tudat
{

namespace json_interface
{

// Special keys (used internally by json_interface, can't be used in JSON files).

const std::string SpecialKeys::root = "~";
const char SpecialKeys::dot = '.';
const std::string SpecialKeys::up = "<-";
const std::string SpecialKeys::rootObject = "#root";
const std::string SpecialKeys::keyPath = "#keypath";

const std::vector< std::string > SpecialKeys::objectContaining =
{
    SpecialKeys::rootObject,
    SpecialKeys::keyPath
};

const std::vector< std::string > SpecialKeys::all =
{
    SpecialKeys::root,
    SpecialKeys::up,
    SpecialKeys::rootObject,
    SpecialKeys::keyPath
};


// Keys recognised by json_interface.

 const std::string Keys::simulationType = "simulationType";
const std::string Keys::initialEpoch = "initialEpoch";
const std::string Keys::finalEpoch = "finalEpoch";
const std::string Keys::globalFrameOrigin = "globalFrameOrigin";
const std::string Keys::globalFrameOrientation = "globalFrameOrientation";


//  Spice
const std::string Keys::spice = "spice";
const std::string Keys::Spice::useStandardKernels = "useStandardKernels";
const std::string Keys::Spice::alternativeKernels = "alternativeKernels";
const std::string Keys::Spice::kernels = "kernels";
const std::string Keys::Spice::preloadEphemeris = "preloadEphemeris";
const std::string Keys::Spice::interpolationOffsets = "interpolationOffsets";
const std::string Keys::Spice::interpolationStep = "interpolationStep";


//  Body
const std::string Keys::bodies = "bodies";

const std::string Keys::Body::useDefaultSettings = "useDefaultSettings";
const std::string Keys::Body::initialState = "initialState";
const std::string Keys::Body::initialStateOrigin = "initialStateOrigin";
const std::string Keys::Body::mass = "mass";
const std::string Keys::Body::rotationalState = "rotationalState";
const std::string Keys::Body::referenceArea = "referenceArea";

// //  Body::State
const std::string Keys::Body::State::type = "type";
// Cartesian
const std::string Keys::Body::State::x = "x";
const std::string Keys::Body::State::y = "y";
const std::string Keys::Body::State::z = "z";
const std::string Keys::Body::State::vx = "vx";
const std::string Keys::Body::State::vy = "vy";
const std::string Keys::Body::State::vz = "vz";
// Keplerian
const std::string Keys::Body::State::centralBodyGravitationalParameter = "centralBodyGravitationalParameter";
const std::string Keys::Body::State::centralBodyAverageRadius = "centralBodyAverageRadius";
const std::string Keys::Body::State::semiMajorAxis = "semiMajorAxis";
const std::string Keys::Body::State::eccentricity = "eccentricity";
const std::string Keys::Body::State::inclination = "inclination";
const std::string Keys::Body::State::argumentOfPeriapsis = "argumentOfPeriapsis";
const std::string Keys::Body::State::longitudeOfAscendingNode = "longitudeOfAscendingNode";
const std::string Keys::Body::State::trueAnomaly = "trueAnomaly";
const std::string Keys::Body::State::meanAnomaly = "meanAnomaly";
const std::string Keys::Body::State::eccentricAnomaly = "eccentricAnomaly";
const std::string Keys::Body::State::semiLatusRectum = "semiLatusRectum";
const std::string Keys::Body::State::meanMotion = "meanMotion";
const std::string Keys::Body::State::period = "period";
const std::string Keys::Body::State::radius = "radius";
const std::string Keys::Body::State::altitude = "altitude";
const std::string Keys::Body::State::periapsisDistance = "periapsisDistance";
const std::string Keys::Body::State::apoapsisDistance = "apoapsisDistance";
const std::string Keys::Body::State::periapsisAltitude = "periapsisAltitude";
const std::string Keys::Body::State::apoapsisAltitude = "apoapsisAltitude";
// Spherical
const std::string Keys::Body::State::epoch = "epoch";
const std::string Keys::Body::State::latitude = "latitude";
const std::string Keys::Body::State::longitude = "longitude";
const std::string Keys::Body::State::speed = "speed";
const std::string Keys::Body::State::flightPathAngle = "flightPathAngle";
const std::string Keys::Body::State::headingAngle = "headingAngle";


// //  Body::aerodynamics
const std::string Keys::Body::aerodynamics = "aerodynamics";
const std::string Keys::Body::Aerodynamics::coefficientsType = "coefficientsType";
const std::string Keys::Body::Aerodynamics::referenceLength = "referenceLength";
const std::string Keys::Body::Aerodynamics::referenceArea = "referenceArea";
const std::string Keys::Body::Aerodynamics::lateralReferenceLength = "lateralReferenceLength";
const std::string Keys::Body::Aerodynamics::momentReferencePoint = "momentReferencePoint";
const std::string Keys::Body::Aerodynamics::independentVariableNames = "independentVariableNames";
const std::string Keys::Body::Aerodynamics::areCoefficientsInAerodynamicFrame = "areCoefficientsInAerodynamicFrame";
const std::string Keys::Body::Aerodynamics::areCoefficientsInNegativeAxisDirection = "areCoefficientsInNegativeAxisDirection";
const std::string Keys::Body::Aerodynamics::controlSurface = "controlSurface";
// Constant
const std::string Keys::Body::Aerodynamics::dragCoefficient = "dragCoefficient";
const std::string Keys::Body::Aerodynamics::forceCoefficients = "forceCoefficients";
const std::string Keys::Body::Aerodynamics::momentCoefficients = "momentCoefficients";
// Tabulated< N >
const std::string Keys::Body::Aerodynamics::independentVariableValues = "independentVariableValues";
// Tabulated< 1 >
const std::string Keys::Body::Aerodynamics::interpolator = "interpolator";


// //  Body::Atmosphere
const std::string Keys::Body::atmosphere = "atmosphere";
const std::string Keys::Body::Atmosphere::type = "type";
const std::string Keys::Body::Atmosphere::densityScaleHeight = "densityScaleHeight";
const std::string Keys::Body::Atmosphere::constantTemperature = "constantTemperature";
const std::string Keys::Body::Atmosphere::densityAtZeroAltitude = "densityAtZeroAltitude";
const std::string Keys::Body::Atmosphere::specificGasConstant = "specificGasConstant";
const std::string Keys::Body::Atmosphere::ratioOfSpecificHeats = "ratioOfSpecificHeats";
const std::string Keys::Body::Atmosphere::file = "file";
const std::string Keys::Body::Atmosphere::independentVariablesNames = "independentVariablesNames";
const std::string Keys::Body::Atmosphere::dependentVariablesNames = "dependentVariablesNames";
const std::string Keys::Body::Atmosphere::boundaryHandling = "boundaryHandling";
const std::string Keys::Body::Atmosphere::spaceWeatherFile = "spaceWeatherFile";

// //  Body::Ephemeris
const std::string Keys::Body::ephemeris = "ephemeris";
const std::string Keys::Body::Ephemeris::type = "type";
const std::string Keys::Body::Ephemeris::frameOrigin = "frameOrigin";
const std::string Keys::Body::Ephemeris::frameOrientation = "frameOrientation";
// const std::string Keys::Body::Ephemeris::makeMultiArc = "makeMultiArc";
const std::string Keys::Body::Ephemeris::correctForStellarAberration = "correctForStellarAberration";
const std::string Keys::Body::Ephemeris::correctForLightTimeAberration = "correctForLightTimeAberration";
const std::string Keys::Body::Ephemeris::convergeLighTimeAberration = "convergeLighTimeAberration";
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

// //  Body::GravityField
const std::string Keys::Body::gravityField = "gravityField";
const std::string Keys::Body::GravityField::type = "type";
const std::string Keys::Body::GravityField::gravitationalParameter = "gravitationalParameter";
const std::string Keys::Body::GravityField::referenceRadius = "referenceRadius";
const std::string Keys::Body::GravityField::cosineCoefficients = "cosineCoefficients";
const std::string Keys::Body::GravityField::sineCoefficients = "sineCoefficients";
const std::string Keys::Body::GravityField::associatedReferenceFrame = "associatedReferenceFrame";
const std::string Keys::Body::GravityField::model = "model";
const std::string Keys::Body::GravityField::file = "file";
const std::string Keys::Body::GravityField::maximumDegree = "maximumDegree";
const std::string Keys::Body::GravityField::maximumOrder = "maximumOrder";
const std::string Keys::Body::GravityField::gravitationalParameterIndex = "gravitationalParameterIndex";
const std::string Keys::Body::GravityField::referenceRadiusIndex = "referenceRadiusIndex";

// //  Body::RadiationPressure
const std::string Keys::Body::radiationPressure = "radiationPressure";
const std::string Keys::Body::RadiationPressure::type = "type";
const std::string Keys::Body::RadiationPressure::sourceBody = "sourceBody";
const std::string Keys::Body::RadiationPressure::occultingBodies = "occultingBodies";
const std::string Keys::Body::RadiationPressure::referenceArea = "referenceArea";
const std::string Keys::Body::RadiationPressure::radiationPressureCoefficient = "radiationPressureCoefficient";

// //  Body::RotationModel
const std::string Keys::Body::rotationModel = "rotationModel";
const std::string Keys::Body::RotationModel::type = "type";
const std::string Keys::Body::RotationModel::originalFrame = "originalFrame";
const std::string Keys::Body::RotationModel::targetFrame = "targetFrame";
const std::string Keys::Body::RotationModel::initialOrientation = "initialOrientation";
const std::string Keys::Body::RotationModel::initialTime = "initialTime";
const std::string Keys::Body::RotationModel::rotationRate = "rotationRate";
const std::string Keys::Body::RotationModel::precessionNutationTheory = "precessionNutationTheory";
const std::string Keys::Body::RotationModel::centralBodyName = "centralBodyName";


// //  Body::ShapeModel
const std::string Keys::Body::shapeModel = "shapeModel";
const std::string Keys::Body::ShapeModel::type = "type";
const std::string Keys::Body::ShapeModel::radius = "radius";
const std::string Keys::Body::ShapeModel::equatorialRadius = "equatorialRadius";
const std::string Keys::Body::ShapeModel::flattening = "flattening";


// //  Body::GravityFieldVariation
const std::string Keys::Body::gravityFieldVariation = "gravityFieldVariation";
const std::string Keys::Body::GravityFieldVariation::bodyDeformationType = "bodyDeformationType";
const std::string Keys::Body::GravityFieldVariation::deformingBodies = "deformingBodies";
const std::string Keys::Body::GravityFieldVariation::loveNumbers = "loveNumbers";
const std::string Keys::Body::GravityFieldVariation::modelInterpolation = "modelInterpolation";
const std::string Keys::Body::GravityFieldVariation::cosineCoefficientCorrections = "cosineCoefficientCorrections";
const std::string Keys::Body::GravityFieldVariation::sineCoefficientCorrections = "sineCoefficientCorrections";
const std::string Keys::Body::GravityFieldVariation::minimumDegree = "minimumDegree";
const std::string Keys::Body::GravityFieldVariation::minimumOrder = "minimumOrder";
const std::string Keys::Body::GravityFieldVariation::interpolator = "interpolator";

// // Body::GroundStation
const std::string Keys::Body::groundStation = "groundStation";
const std::string Keys::Body::GroundStation::stationPosition = "stationPosition";
const std::string Keys::Body::GroundStation::positionElementType = "positionElementType";
const std::string Keys::Body::GroundStation::stationName = "stationName";

//  Variable
const std::string Keys::Variable::type = "type";
const std::string Keys::Variable::dependentVariableType = "dependentVariableType";
const std::string Keys::Variable::body = "body";
const std::string Keys::Variable::componentIndex = "componentIndex";
const std::string Keys::Variable::componentIndices = "componentIndices";
const std::string Keys::Variable::useAccelerationNorm = "useAccelerationNorm";
const std::string Keys::Variable::relativeToBody = "relativeToBody";
const std::string Keys::Variable::accelerationType = "accelerationType";
//const std::string Keys::Variable::bodyUndergoingAcceleration = "bodyUndergoingAcceleration";
const std::string Keys::Variable::bodyExertingAcceleration = "bodyExertingAcceleration";
const std::string Keys::Variable::torqueType = "torqueType";
//const std::string Keys::Variable::bodyUndergoingTorque = "bodyUndergoingTorque";
const std::string Keys::Variable::bodyExertingTorque = "bodyExertingTorque";
const std::string Keys::Variable::baseFrame = "baseFrame";
const std::string Keys::Variable::targetFrame = "targetFrame";
const std::string Keys::Variable::angle = "angle";
const std::string Keys::Variable::deformationType = "deformationType";
const std::string Keys::Variable::identifier = "identifier";
const std::string Keys::Variable::derivativeWrtBody = "derivativeWrtBody";
const std::string Keys::Variable::thirdBody = "thirdBody";


// Parameter
const std::string Keys::parametersToEstimate = "parametersToEstimate";

const std::string Keys::Parameter::parameterType = "parameterType";
const std::string Keys::Parameter::associatedBody = "associatedBody";
const std::string Keys::Parameter::secondaryIdentifier = "secondaryIdentifier";

const std::string Keys::Parameter::initialStateValue = "initialStateValue";
const std::string Keys::Parameter::centralBody = "centralBody";
const std::string Keys::Parameter::frameOrientation = "frameOrientation";
const std::string Keys::Parameter::arcStartTimes = "arcStartTimes";

const std::string Keys::Parameter::coefficientIndices = "coefficientIndices";
const std::string Keys::Parameter::maximumDegree = "maximumDegree";
const std::string Keys::Parameter::minimumDegree = "minimumDegree";
const std::string Keys::Parameter::maximumOrder = "maximumOrder";
const std::string Keys::Parameter::minimumOrder = "minimumOrder";

const std::string Keys::Parameter::deformingBodies = "deformingBodies";

const std::string Keys::Parameter::componentsToEstimate = "componentsToEstimate";

const std::string Keys::Parameter::observableType = "observableType";
const std::string Keys::Parameter::linkEnds = "linkEnds";
const std::string Keys::Parameter::referenceLinkEnd = "referenceLinkEnd";

const std::string Keys::Parameter::degree = "degree";
const std::string Keys::Parameter::orders = "orders";
const std::string Keys::Parameter::useComplexValue = "useComplexValue";

const std::string Keys::estimationSettings = "estimation";

const std::string Keys::Estimation::inverseAprioriCovariance = "inverseAprioriCovariance";
const std::string Keys::Estimation::reintegrateEquationsOnFirstIteration = "reintegrateOnFirstIteration";
const std::string Keys::Estimation::reintegrateVariationalEquations = "reintegrateVariationalEquations";
const std::string Keys::Estimation::saveInformationMatrix = "savePartials";
const std::string Keys::Estimation::printOutput = "printOutput";
const std::string Keys::Estimation::saveResidualsAndParametersFromEachIteration = "saveEstimationHistory";
const std::string Keys::Estimation::saveStateHistoryForEachIteration = "saveStateHistory";

const std::string Keys::Estimation::maximumNumberOfIterations = "maximumIteration";
const std::string Keys::Estimation::minimumResidualChange = "minimumResidualChange";
const std::string Keys::Estimation::minimumResidual = "minimumResidual";
const std::string Keys::Estimation::numberOfIterationsWithoutImprovement = "numberOfUnimprovedIterations";

const std::string Keys::Estimation::dataWeights = "dataWeights";

// Observation
const std::string Keys::observations = "observations";
const std::string Keys::Observation::observableType = "observableType";
const std::string Keys::Observation::lightTimeCorrectionSettingsList = "lightTimeCorrections";
const std::string Keys::Observation::biasSettings = "bias";

const std::string Keys::Observation::transmitterProperTimeRateSettings = "transmitterProperTimeRate";
const std::string Keys::Observation::receiverProperTimeRateSettings = "receiverProperTimeRate";

const std::string Keys::Observation::constantIntegrationTime = "constantIntegrationTime";

const std::string Keys::Observation::oneWayRangeObsevationSettings = "oneWayRangeObsevation";
const std::string Keys::Observation::retransmissionTimes = "retransmissionTimes";

const std::string Keys::Observation::uplinkOneWayDopplerSettings = "uplinkOneWayDoppler";
const std::string Keys::Observation::downlinkOneWayDopplerSettings = "downlinkOneWayDoppler";

const std::string Keys::Observation::properTimeRateType = "properTimeRateType";
const std::string Keys::Observation::centralBody = "centralBody";

const std::string Keys::Observation::lightTimeCorrectionType = "lightTimeCorrectionType";
const std::string Keys::Observation::perturbingBodies = "perturbingBodies";

const std::string Keys::Observation::observationSimulationTimesType = "observationSimulationTimesType";
const std::string Keys::Observation::observationSimulationTimesList = "observationSimulationTimesList";

const std::string Keys::Observation::observableViabilityType = "viabilityType";
const std::string Keys::Observation::associatedLinkEnd = "associatedLinkEnd";
const std::string Keys::Observation::doubleParameter = "doubleParameter";
const std::string Keys::Observation::stringParameter = "stringParameter";


// ObservationBias
const std::string Keys::ObservationBias::biasType = "biasType";
const std::string Keys::ObservationBias::multipleBiasesList = "multipleBiasesList";
const std::string Keys::ObservationBias::constantBias = "constantBias";

const std::string Keys::ObservationBias::arcWiseBiasList = "arcWiseBiasList";
const std::string Keys::ObservationBias::arcStartTimes = "arcStartTimes";
const std::string Keys::ObservationBias::referenceLinkEnd = "referenceLinkEnd";


//  Propagator
const std::string Keys::propagators = "propagators";

const std::string Keys::Propagator::integratedStateType = "integratedStateType";
const std::string Keys::Propagator::initialStates = "initialStates";
const std::string Keys::Propagator::bodiesToPropagate = "bodiesToPropagate";

// Translational
const std::string Keys::Propagator::type = "type";
const std::string Keys::Propagator::centralBodies = "centralBodies";

// //  Acceleration
const std::string Keys::Propagator::accelerations = "accelerations";
const std::string Keys::Propagator::Acceleration::type = "type";
const std::string Keys::Propagator::Acceleration::maximumDegree = "maximumDegree";
const std::string Keys::Propagator::Acceleration::maximumOrder = "maximumOrder";
const std::string Keys::Propagator::Acceleration::maximumDegreeOfBodyExertingAcceleration = "maximumDegreeOfBodyExertingAcceleration";
const std::string Keys::Propagator::Acceleration::maximumOrderOfBodyExertingAcceleration = "maximumOrderOfBodyExertingAcceleration";
const std::string Keys::Propagator::Acceleration::maximumDegreeOfBodyUndergoingAcceleration = "maximumDegreeOfBodyUndergoingAcceleration";
const std::string Keys::Propagator::Acceleration::maximumOrderOfBodyUndergoingAcceleration = "maximumOrderOfBodyUndergoingAcceleration";
const std::string Keys::Propagator::Acceleration::maximumDegreeOfCentralBody = "maximumDegreeOfCentralBody";
const std::string Keys::Propagator::Acceleration::maximumOrderOfCentralBody = "maximumOrderOfCentralBody";
const std::string Keys::Propagator::Acceleration::calculateSchwarzschildCorrection = "calculateSchwarzschildCorrection";
const std::string Keys::Propagator::Acceleration::calculateLenseThirringCorrection = "calculateLenseThirringCorrection";
const std::string Keys::Propagator::Acceleration::calculateDeSitterCorrection = "calculateDeSitterCorrection";
const std::string Keys::Propagator::Acceleration::primaryBody = "primaryBody";
const std::string Keys::Propagator::Acceleration::centralBodyAngularMomentum = "centralBodyAngularMomentum";
const std::string Keys::Propagator::Acceleration::constantAcceleration = "constantAcceleration";
const std::string Keys::Propagator::Acceleration::sineAcceleration = "sineAcceleration";
const std::string Keys::Propagator::Acceleration::cosineAcceleration = "cosineAcceleration";

// // //  Acceleration::Thrust

// // // //  Acceleration::Thrust::Direction
const std::string Keys::Propagator::Acceleration::Thrust::direction = "direction";
const std::string Keys::Propagator::Acceleration::Thrust::Direction::type = "type";
const std::string Keys::Propagator::Acceleration::Thrust::Direction::relativeBody = "relativeBody";
const std::string Keys::Propagator::Acceleration::Thrust::Direction::colinearWithVelocity = "colinearWithVelocity";
const std::string Keys::Propagator::Acceleration::Thrust::Direction::towardsRelativeBody = "towardsRelativeBody";

// // // //  Acceleration::Thrust::Direction
const std::string Keys::Propagator::Acceleration::Thrust::magnitude = "magnitude";
const std::string Keys::Propagator::Acceleration::Thrust::Magnitude::type = "type";
const std::string Keys::Propagator::Acceleration::Thrust::Magnitude::originID = "originID";
const std::string Keys::Propagator::Acceleration::Thrust::Magnitude::constantMagnitude = "constantMagnitude";
const std::string Keys::Propagator::Acceleration::Thrust::Magnitude::specificImpulse = "specificImpulse";
const std::string Keys::Propagator::Acceleration::Thrust::Magnitude::bodyFixedDirection = "bodyFixedDirection";
const std::string Keys::Propagator::Acceleration::Thrust::Magnitude::useAllEngines = "useAllEngines";

const std::string Keys::Propagator::Acceleration::Thrust::dataInterpolation = "dataInterpolation";
const std::string Keys::Propagator::Acceleration::Thrust::specificImpulse = "specificImpulse";
const std::string Keys::Propagator::Acceleration::Thrust::frame = "frame";
const std::string Keys::Propagator::Acceleration::Thrust::centralBody = "centralBody";


// //  Mass rate model
const std::string Keys::Propagator::massRateModels = "massRateModels";
const std::string Keys::Propagator::MassRateModel::type = "type";
const std::string Keys::Propagator::MassRateModel::useAllThrustModels = "useAllThrustModels";
const std::string Keys::Propagator::MassRateModel::associatedThrustSource = "associatedThrustSource";


// //  Torque
const std::string Keys::Propagator::torques = "torques";
const std::string Keys::Propagator::Torque::type = "type";



// Termination
const std::string Keys::termination = "termination";
const std::string Keys::Termination::anyOf = "anyOf";
const std::string Keys::Termination::allOf = "allOf";
const std::string Keys::Termination::variable = "variable";
const std::string Keys::Termination::lowerLimit = "lowerLimit";
const std::string Keys::Termination::upperLimit = "upperLimit";


//  Integrator
const std::string Keys::integrator = "integrator";
const std::string Keys::Integrator::type = "type";
const std::string Keys::Integrator::initialTime = "initialTime";
const std::string Keys::Integrator::stepSize = "stepSize";
const std::string Keys::Integrator::initialStepSize = "initialStepSize";
const std::string Keys::Integrator::saveFrequency = "saveFrequency";
const std::string Keys::Integrator::assessTerminationOnMinorSteps = "assessTerminationOnMinorSteps";
const std::string Keys::Integrator::rungeKuttaCoefficientSet = "rungeKuttaCoefficientSet";
const std::string Keys::Integrator::minimumStepSize = "minimumStepSize";
const std::string Keys::Integrator::maximumStepSize = "maximumStepSize";
const std::string Keys::Integrator::relativeErrorTolerance = "relativeErrorTolerance";
const std::string Keys::Integrator::absoluteErrorTolerance = "absoluteErrorTolerance";
const std::string Keys::Integrator::areTolerancesDefinedAsScalar = "areTolerancesDefinedAsScalar";
const std::string Keys::Integrator::safetyFactorForNextStepSize = "safetyFactorForNextStepSize";
const std::string Keys::Integrator::maximumFactorIncreaseForNextStepSize = "maximumFactorIncreaseForNextStepSize";
const std::string Keys::Integrator::minimumFactorDecreaseForNextStepSize = "minimumFactorDecreaseForNextStepSize";
const std::string Keys::Integrator::bandwidth = "bandwidth";
const std::string Keys::Integrator::extrapolationSequence = "extrapolationSequence";
const std::string Keys::Integrator::maximumNumberOfSteps = "maximumNumberOfSteps";
const std::string Keys::Integrator::maximumOrder = "maximumOrder";
const std::string Keys::Integrator::minimumOrder = "minimumOrder";

//  Interpolation

// //  Interpolation::DataMap
const std::string Keys::Interpolation::DataMap::map = "map";
const std::string Keys::Interpolation::DataMap::file = "file";
const std::string Keys::Interpolation::DataMap::independentVariableValues = "independentVariableValues";
const std::string Keys::Interpolation::DataMap::dependentVariableValues = "dependentVariableValues";
const std::string Keys::Interpolation::DataMap::dependentVariableFirstDerivativeValues = "dependentVariableFirstDerivativeValues";

// //  Interpolation::Interpolator
const std::string Keys::Interpolation::Interpolator::type = "type";
const std::string Keys::Interpolation::Interpolator::lookupScheme = "lookupScheme";
const std::string Keys::Interpolation::Interpolator::useLongDoubleTimeStep = "useLongDoubleTimeStep";
const std::string Keys::Interpolation::Interpolator::order = "order";
const std::string Keys::Interpolation::Interpolator::boundaryHandling = "boundaryHandling";
const std::string Keys::Interpolation::Interpolator::lagrangeBoundaryHandling = "lagrangeBoundaryHandling";

// //  Interpolation::DataInterpolation
const std::string Keys::Interpolation::DataInterpolation::data = "data";
const std::string Keys::Interpolation::DataInterpolation::interpolator = "interpolator";

// //  Interpolation::ModelInterpolation
const std::string Keys::Interpolation::ModelInterpolation::initialTime = "initialTime";
const std::string Keys::Interpolation::ModelInterpolation::finalTime = "finalTime";
const std::string Keys::Interpolation::ModelInterpolation::timeStep = "timeStep";
const std::string Keys::Interpolation::ModelInterpolation::interpolator = "interpolator";


//  Export
const std::string Keys::xport = "export";
const std::string Keys::Export::file = "file";
const std::string Keys::Export::variables = "variables";
const std::string Keys::Export::header = "header";
const std::string Keys::Export::epochsInFirstColumn = "epochsInFirstColumn";
const std::string Keys::Export::onlyInitialStep = "onlyInitialStep";
const std::string Keys::Export::onlyFinalStep = "onlyFinalStep";
const std::string Keys::Export::numericalPrecision = "numericalPrecision";
const std::string Keys::Export::printVariableIndicesToTerminal = "printVariableIndicesToTerminal";

//  Options
const std::string Keys::options = "options";
const std::string Keys::Options::notifyOnPropagationStart = "notifyOnPropagationStart";
const std::string Keys::Options::notifyOnPropagationTermination = "notifyOnPropagationTermination";
const std::string Keys::Options::printInterval = "printInterval";
const std::string Keys::Options::defaultValueUsedForMissingKey = "defaultValueUsedForMissingKey";
const std::string Keys::Options::unusedKey = "unusedKey";
const std::string Keys::Options::fullSettingsFile = "fullSettingsFile";
const std::string Keys::Options::tagOutputFilesIfPropagationFails = "tagOutputFilesIfPropagationFails";


// KEYPATH

//! Get the int-value of an int-convertible key.
int indexFromKey( const std::string& key )
{
    boost::cmatch groups;
    boost::regex_match( key.c_str( ), groups, boost::regex( R"(\@(\d+))" ) );
    if ( groups[ 1 ].matched )
    {
        try
        {
            return std::stoi( groups[ 1 ] );
        }
        catch ( ... ) { }
    }
    return -1;
}

//! Constructor with a single key path string representation.
KeyPath::KeyPath( const std::string& keyPathStringRepresentation ) : std::vector< std::string >( )
{
    const std::vector< std::string > keys = split( keyPathStringRepresentation, SpecialKeys::dot );
    for ( const std::string key : keys )
    {
        boost::cmatch groups;
        boost::regex_match( key.c_str( ), groups, boost::regex( R"((.+?)\[(\d+?)\])" ) );
        if ( groups[ 1 ].matched && groups[ 2 ].matched )
        {
            const std::string arrayKey( groups[ 1 ] );
            const std::string arrayIndex( groups[ 2 ] );
            push_back( arrayKey );
            push_back( "@" + arrayIndex );
        }
        else
        {
            push_back( key );
        }
    }
}

//! String representation for `KeyPath`, as key.subkey.vectorIndex.subsubkey ...
std::ostream& operator << ( std::ostream& stringRepresentation, const KeyPath& keyPath )
{
    bool somethingAdded = false;
    for ( unsigned int i = 0; i < keyPath.size( ); ++i )
    {
        const std::string key = keyPath.at( i );
        if ( key != SpecialKeys::root )
        {
            const int intKey = indexFromKey( key );
            if ( intKey >= 0 )
            {
                stringRepresentation << '[' << intKey << ']';
            }
            else
            {
                if ( somethingAdded )
                {
                    stringRepresentation << SpecialKeys::dot;
                }
                stringRepresentation << key;
            }
            somethingAdded = true;
        }
    }
    return stringRepresentation;
}

//! Get the canonical representation of the key path.
KeyPath KeyPath::canonical( const KeyPath& basePath ) const
{
    // Compound key path
    KeyPath compoundKeyPath;

    // Check absolute paths
    if ( ! this->isAbsolute( ) && basePath.isAbsolute( ) )
    {
        compoundKeyPath = basePath / *this;
    }
    else
    {
        compoundKeyPath = *this;
    }

    // Remove ..
    KeyPath canonicalKeyPath = { };
    for ( std::string key : compoundKeyPath )
    {
        if ( key == SpecialKeys::up )
        {
            if ( canonicalKeyPath.size( ) > 0 )
            {
                canonicalKeyPath.pop_back( );
            }
        }
        else
        {
            canonicalKeyPath.push_back( key );
        }
    }

    return canonicalKeyPath;
}


} // namespace json_interface

} // namespace tudat
