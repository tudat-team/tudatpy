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

//! Special keys (used internally by json_interface, can't be used in JSON files).

const std::string SpecialKeys::root = "~";
const std::string SpecialKeys::up = "..";
const std::string SpecialKeys::rootObject = "#root";
const std::string SpecialKeys::keyPath = "#keypath";

const std::vector< std::string > SpecialKeys::all = { SpecialKeys::root,
                                                      SpecialKeys::up,
                                                      SpecialKeys::rootObject,
                                                      SpecialKeys::keyPath };


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


////// Body::GravityFieldVariation
const std::string Keys::Body::gravityFieldVariations = "gravityFieldVariations";
const std::string Keys::Body::GravityFieldVariation::bodyDeformationType = "bodyDeformationType";
const std::string Keys::Body::GravityFieldVariation::modelInterpolation = "modelInterpolation";
const std::string Keys::Body::GravityFieldVariation::deformingBodies = "deformingBodies";
const std::string Keys::Body::GravityFieldVariation::loveNumbers = "loveNumbers";
const std::string Keys::Body::GravityFieldVariation::referenceRadius = "referenceRadius";
const std::string Keys::Body::GravityFieldVariation::cosineCoefficientCorrections = "cosineCoefficientCorrections";
const std::string Keys::Body::GravityFieldVariation::sineCoefficientCorrections = "sineCoefficientCorrections";
const std::string Keys::Body::GravityFieldVariation::minimumDegree = "minimumDegree";
const std::string Keys::Body::GravityFieldVariation::minimumOrder = "minimumOrder";



/// Propagator
const std::string Keys::propagation = "propagation";
const std::string Keys::Propagator::integratedStateType = "integratedStateType";
const std::string Keys::Propagator::initialStates = "initialStates";

////// Propagator::Termination
const std::string Keys::Propagator::termination = "termination";

////// Propagator::Save
const std::string Keys::Propagator::output = "output";

const std::string Keys::Propagator::printInterval = "printInterval";

// Hybrid
const std::string Keys::Propagator::propagators = "propagators";

// Translational
const std::string Keys::Propagator::type = "type";
const std::string Keys::Propagator::centralBodies = "centralBodies";
const std::string Keys::Propagator::bodiesToPropagate = "bodiesToPropagate";
// const std::string Keys::Propagator::initialStateTypes = "initialStateTypes";

////// Acceleration
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

///////// Acceleration::ThrustDirection
const std::string Keys::Propagator::Acceleration::direction = "direction";
const std::string Keys::Propagator::Acceleration::ThrustDirection::type = "type";
const std::string Keys::Propagator::Acceleration::ThrustDirection::relativeBody = "relativeBody";
const std::string Keys::Propagator::Acceleration::ThrustDirection::colinearWithVelocity = "colinearWithVelocity";
const std::string Keys::Propagator::Acceleration::ThrustDirection::towardsRelativeBody = "towardsRelativeBody";

///////// Acceleration::ThrustDirection
const std::string Keys::Propagator::Acceleration::magnitude = "magnitude";
const std::string Keys::Propagator::Acceleration::ThrustMagnitude::type = "type";
const std::string Keys::Propagator::Acceleration::ThrustMagnitude::originID = "originID";
const std::string Keys::Propagator::Acceleration::ThrustMagnitude::constantMagnitude = "constantMagnitude";
const std::string Keys::Propagator::Acceleration::ThrustMagnitude::specificImpulse = "specificImpulse";
const std::string Keys::Propagator::Acceleration::ThrustMagnitude::bodyFixedDirection = "bodyFixedDirection";
const std::string Keys::Propagator::Acceleration::ThrustMagnitude::useAllEngines = "useAllEngines";


////// Mass rate
const std::string Keys::Propagator::massRates = "massRates";


////// Torque
const std::string Keys::Propagator::torques = "torques";




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


/// Interpolator
const std::string Keys::interpolator = "interpolator";
const std::string Keys::Interpolator::type = "type";
const std::string Keys::Interpolator::lookupScheme = "lookupScheme";
const std::string Keys::Interpolator::useLongDoubleTimeStep = "useLongDoubleTimeStep";
const std::string Keys::Interpolator::order = "order";
const std::string Keys::Interpolator::boundaryHandling = "boundaryHandling";

/// ModelInterpolation
const std::string Keys::ModelInterpolation::initialTime = "initialTime";
const std::string Keys::ModelInterpolation::finalTime = "finalTime";
const std::string Keys::ModelInterpolation::timeStep = "timeStep";
const std::string Keys::ModelInterpolation::interpolator = "interpolator";

/// Output
const std::string Keys::output = "output";




//! Key paths recognised by `json_interface`.

/// Simulation
const KeyPath KeyPaths::Simulation::startEpoch = { Keys::simulation, Keys::Simulation::startEpoch };
const KeyPath KeyPaths::Simulation::endEpoch = { Keys::simulation, Keys::Simulation::endEpoch };
const KeyPath KeyPaths::Simulation::globalFrameOrigin = { Keys::simulation, Keys::Simulation::globalFrameOrigin };
const KeyPath KeyPaths::Simulation::globalFrameOrientation = { Keys::simulation, Keys::Simulation::globalFrameOrientation };
const KeyPath KeyPaths::Simulation::spiceKernels = { Keys::simulation, Keys::Simulation::spiceKernels };
const KeyPath KeyPaths::Simulation::preloadSpiceData = { Keys::simulation, Keys::Simulation::preloadSpiceData };

/*
/// Integrator
const KeyPath KeyPaths::Integrator::initialTime = { Keys::integrator, Keys::Integrator::initialTime };
*/



//! -DOC
KeyPath KeyPath::canonical( const KeyPath& basePath ) const
{
    // Compound key path
    KeyPath compoundKeyPath;

    // Check absolute paths
    if ( this->isAbsolute( ) )
    {
        compoundKeyPath = *this;
    }
    else if ( basePath.isAbsolute( ) )
    {
        compoundKeyPath = basePath / *this;
    }
    else
    {
        throw std::runtime_error( "Could not determine canonical key path because the base path is not absolute." );
    }

    // Remove ..
    KeyPath canonicalKeyPath;
    for ( std::string key : compoundKeyPath )
    {
        if ( key == SpecialKeys::up )
        {
            canonicalKeyPath.pop_back( );
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
