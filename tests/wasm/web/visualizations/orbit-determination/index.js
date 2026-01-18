// Orbit Determination / Differential Correction Visualization
// Demonstrates batch least squares orbit determination with position observations.
// Shows convergence comparison between OMM (secular drift) and Full Force dynamics.
// REQUIRES Tudat WASM bindings - no JavaScript fallbacks per Agents.md

import { configureClockForOrbit, clearOrbitEntities } from '../shared/utils.js';

/**
 * Run orbit determination and visualize results
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @param {Function} log - Logging function
 * @param {string} dynamicsModel - "omm" or "fullforce"
 * @returns {Object} Result data including residuals for charting
 */
export function runOrbitDeterminationVisualization(viewer, orbitEntities, log, dynamicsModel = 'fullforce') {
    clearOrbitEntities(viewer, orbitEntities);

    // Tudat bindings are REQUIRED - no fallbacks
    if (typeof Module === 'undefined' || typeof Module.runOrbitDetermination !== 'function') {
        log('ERROR: Tudat WASM bindings required for Orbit Determination visualization', 'error');
        log('Module.runOrbitDetermination function not available', 'error');
        return null;
    }

    // ISS-like orbit - truth state
    const truthState = {
        semiMajorAxis: 6793,      // km
        eccentricity: 0.0001,
        inclination: 51.6,        // degrees
        raan: 45.0,
        argPeriapsis: 90.0,
        trueAnomaly: 0.0
    };

    // Initial guess - perturbed from truth (represents a priori knowledge)
    const initialGuess = {
        semiMajorAxis: 6793 + 5,  // 5 km error in SMA
        eccentricity: 0.0001,
        inclination: 51.6 + 0.1,  // 0.1 deg error
        raan: 45.0 + 0.2,         // 0.2 deg error
        argPeriapsis: 90.0,
        trueAnomaly: 0.5          // 0.5 deg error
    };

    const period = 5400;          // ~90 min orbit
    const duration = period * 2;  // 2 orbits of observations
    const numObservations = 50;   // Position observations
    const noiseStdDev = 100;      // 100 meter position noise
    const maxIterations = 10;

    log(`Running orbit determination with ${dynamicsModel.toUpperCase()} dynamics...`, 'info');
    log(`Observations: ${numObservations} over ${(duration/3600).toFixed(1)} hours`, 'info');
    log(`Position noise: ${noiseStdDev} m (1-sigma)`, 'info');

    let result;
    try {
        result = Module.runOrbitDetermination(
            JSON.stringify(initialGuess),
            JSON.stringify(truthState),
            duration,
            numObservations,
            noiseStdDev,
            maxIterations,
            dynamicsModel
        );
    } catch (e) {
        log('Orbit determination failed: ' + e.message, 'error');
        return null;
    }

    // Parse result
    const numIterations = Math.round(result[0]);
    const numObs = Math.round(result[1]);
    log(`Converged in ${numIterations} iterations`, 'info');

    const iterationDataSize = 1 + 6 + numObs * 3;  // rms + state + residuals
    const iterations = [];

    let offset = 2;
    for (let iter = 0; iter < numIterations; iter++) {
        const rms = result[offset];
        const state = [];
        for (let j = 0; j < 6; j++) {
            state.push(result[offset + 1 + j]);
        }
        const residuals = [];
        for (let j = 0; j < numObs * 3; j++) {
            residuals.push(result[offset + 7 + j]);
        }
        iterations.push({ rms, state, residuals });
        offset += iterationDataSize;
        log(`Iteration ${iter}: RMS = ${rms.toFixed(2)} m`, 'info');
    }

    // Extract truth trajectory
    const truthTrajectory = [];
    for (let i = 0; i < numObs; i++) {
        truthTrajectory.push({
            x: result[offset + i * 3],
            y: result[offset + i * 3 + 1],
            z: result[offset + i * 3 + 2]
        });
    }
    offset += numObs * 3;

    // Extract observations
    const observations = [];
    for (let i = 0; i < numObs; i++) {
        observations.push({
            x: result[offset + i * 3],
            y: result[offset + i * 3 + 1],
            z: result[offset + i * 3 + 2]
        });
    }
    offset += numObs * 3;

    // Extract estimated trajectory (from final converged state)
    const estimatedTrajectory = [];
    for (let i = 0; i < numObs; i++) {
        estimatedTrajectory.push({
            x: result[offset + i * 3],
            y: result[offset + i * 3 + 1],
            z: result[offset + i * 3 + 2]
        });
    }

    // Visualize in Cesium
    const clock = viewer.clock;
    const startTime = clock.startTime;
    const dt = duration / (numObs - 1);

    // Truth trajectory (white dashed line)
    const truthPositions = truthTrajectory.map(p => new Cesium.Cartesian3(p.x, p.y, p.z));
    const truthOrbit = viewer.entities.add({
        name: 'Truth Trajectory',
        polyline: {
            positions: truthPositions,
            width: 2,
            material: new Cesium.PolylineDashMaterialProperty({
                color: Cesium.Color.WHITE,
                dashLength: 16
            })
        }
    });
    orbitEntities.push(truthOrbit);

    // Observations as points (yellow)
    for (let i = 0; i < numObs; i++) {
        const obs = observations[i];
        const obsEntity = viewer.entities.add({
            name: `Observation ${i + 1}`,
            position: new Cesium.Cartesian3(obs.x, obs.y, obs.z),
            point: {
                pixelSize: 5,
                color: Cesium.Color.YELLOW.withAlpha(0.7),
                outlineColor: Cesium.Color.BLACK,
                outlineWidth: 1
            }
        });
        orbitEntities.push(obsEntity);
    }

    // Estimated orbit as thin white solid line
    const estimatedPositions = estimatedTrajectory.map(p => new Cesium.Cartesian3(p.x, p.y, p.z));
    const estimatedOrbit = viewer.entities.add({
        name: 'Estimated Orbit',
        polyline: {
            positions: estimatedPositions,
            width: 1,
            material: Cesium.Color.WHITE.withAlpha(0.8)
        }
    });
    orbitEntities.push(estimatedOrbit);

    // Create animated satellite following estimated trajectory
    const satelliteColor = dynamicsModel === 'omm' ? Cesium.Color.CYAN : Cesium.Color.LIME;

    // Create sampled position property for animation
    const estimatedSampledPosition = new Cesium.SampledPositionProperty();
    estimatedSampledPosition.setInterpolationOptions({
        interpolationDegree: 5,
        interpolationAlgorithm: Cesium.LagrangePolynomialApproximation
    });

    for (let i = 0; i < numObs; i++) {
        const sampleTime = Cesium.JulianDate.addSeconds(startTime, i * dt, new Cesium.JulianDate());
        const pos = estimatedTrajectory[i];
        estimatedSampledPosition.addSample(sampleTime, new Cesium.Cartesian3(pos.x, pos.y, pos.z));
    }

    // Animated satellite entity
    const satellite = viewer.entities.add({
        name: dynamicsModel === 'omm' ? 'OMM Estimated' : 'Full Force Estimated',
        description: `Estimated orbit from ${dynamicsModel.toUpperCase()} dynamics\nConverged in ${iterations.length} iterations\nFinal RMS: ${iterations[iterations.length - 1].rms.toFixed(1)} m`,
        position: estimatedSampledPosition,
        orientation: new Cesium.VelocityOrientationProperty(estimatedSampledPosition),
        point: {
            pixelSize: 12,
            color: satelliteColor,
            outlineColor: Cesium.Color.WHITE,
            outlineWidth: 2
        },
        path: {
            show: true,
            leadTime: 0,
            trailTime: period * 0.3,
            width: 2,
            material: new Cesium.PolylineGlowMaterialProperty({
                glowPower: 0.3,
                color: satelliteColor
            })
        },
        label: {
            text: dynamicsModel === 'omm' ? 'OMM' : 'Full Force',
            font: '12px monospace',
            fillColor: satelliteColor,
            pixelOffset: new Cesium.Cartesian2(0, -15)
        }
    });
    orbitEntities.push(satellite);

    // Configure clock
    configureClockForOrbit(viewer, duration, null, period / 20);

    // Camera view
    viewer.camera.flyTo({
        destination: Cesium.Cartesian3.fromDegrees(0, 0, 25000000),
        orientation: {
            heading: 0,
            pitch: Cesium.Math.toRadians(-90),
            roll: 0
        },
        duration: 1.0
    });

    // Return data for residuals chart
    return {
        name: 'Orbit Determination',
        description: `${dynamicsModel.toUpperCase()} dynamics - ${numIterations} iterations`,
        iterations: iterations,
        truthTrajectory: truthTrajectory,
        estimatedTrajectory: estimatedTrajectory,
        observations: observations,
        dynamicsModel: dynamicsModel,
        numObservations: numObs,
        noiseStdDev: noiseStdDev,
        period: period
    };
}

/**
 * Show orbit determination visualization with selectable dynamics model
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @param {Function} log - Logging function
 * @param {string} dynamicsModel - "omm" or "fullforce" (default: "fullforce")
 */
export function showOrbitDeterminationVisualization(viewer, orbitEntities, log, dynamicsModel = 'fullforce') {
    const result = runOrbitDeterminationVisualization(viewer, orbitEntities, log, dynamicsModel);

    if (!result) {
        return {
            name: 'Orbit Determination',
            description: 'ERROR: Tudat bindings required',
            chartData: null
        };
    }

    return {
        name: 'Orbit Determination',
        description: result.description,
        chartData: result
    };
}
