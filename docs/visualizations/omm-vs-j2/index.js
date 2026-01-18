// OMM vs Two-Body Model Comparison Visualization
// Shows divergence between OMM mean-element propagation (with J2 secular rates)
// and pure Two-Body (Kepler) propagation. Demonstrates the effect of J2 secular
// perturbations on orbital plane orientation (RAAN drift, argument of perigee drift).
// REQUIRES Tudat WASM bindings - no JavaScript fallbacks per Agents.md

import { configureClockForOrbit, clearOrbitEntities } from '../shared/utils.js';

/**
 * Add OMM vs J2 comparison visualization
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @param {number} period - Orbital period in seconds
 * @param {number} numOrbits - Number of orbits to propagate
 * @param {Function} log - Logging function
 * @param {Object} chartContext - Chart context for separation chart
 */
export function addOMMvsJ2Visualization(viewer, orbitEntities, period, numOrbits, log, chartContext) {
    // Example OMM - ISS-like orbit
    // In a real implementation, this would be parsed from OMM XML/KVN
    const ommElements = {
        semiMajorAxis: 6793,      // km (ISS ~420km altitude)
        eccentricity: 0.0001,     // Nearly circular
        inclination: 51.6,        // degrees
        raan: 45.0,               // degrees
        argPeriapsis: 90.0,       // degrees
        meanAnomaly: 0.0          // degrees (at epoch)
    };

    const totalTime = period * numOrbits;
    const numSamples = 360 * numOrbits;

    // Tudat bindings are REQUIRED - no fallbacks
    if (typeof Module === 'undefined' || typeof Module.propagateOMMvsJ2 !== 'function') {
        log('ERROR: Tudat WASM bindings required for OMM vs J2 visualization', 'error');
        log('Module.propagateOMMvsJ2 function not available', 'error');
        return null;
    }

    log('Using Tudat OMM vs J2 comparison', 'info');

    let ephemerisData;
    try {
        // Convert OMM to JSON string for Tudat
        const ommJson = JSON.stringify(ommElements);
        ephemerisData = Module.propagateOMMvsJ2(ommJson, totalTime, numSamples);
        log(`Got ${ephemerisData.length / 7} samples from Tudat`, 'info');
    } catch (e) {
        log('OMM vs J2 propagation failed: ' + e.message, 'error');
        return null;
    }

    const clock = viewer.clock;
    const startTime = clock.startTime;

    const ommPositions = new Cesium.SampledPositionProperty();
    ommPositions.setInterpolationOptions({
        interpolationDegree: 5,
        interpolationAlgorithm: Cesium.LagrangePolynomialApproximation
    });

    const j2Positions = new Cesium.SampledPositionProperty();
    j2Positions.setInterpolationOptions({
        interpolationDegree: 5,
        interpolationAlgorithm: Cesium.LagrangePolynomialApproximation
    });

    const orbitPositions = [];
    const samplesPerOrbit = Math.floor(numSamples / numOrbits);
    const separationData = [];
    let maxSep = 0;

    for (let i = 0; i < numSamples; i++) {
        const idx = i * 7;
        const t = ephemerisData[idx];
        const sampleTime = Cesium.JulianDate.addSeconds(startTime, t, new Cesium.JulianDate());

        const ommPos = new Cesium.Cartesian3(
            ephemerisData[idx + 1],
            ephemerisData[idx + 2],
            ephemerisData[idx + 3]
        );
        const j2Pos = new Cesium.Cartesian3(
            ephemerisData[idx + 4],
            ephemerisData[idx + 5],
            ephemerisData[idx + 6]
        );

        ommPositions.addSample(sampleTime, ommPos);
        j2Positions.addSample(sampleTime, j2Pos);

        if (i < samplesPerOrbit) {
            orbitPositions.push(ommPos);
        }

        const dx = ephemerisData[idx + 1] - ephemerisData[idx + 4];
        const dy = ephemerisData[idx + 2] - ephemerisData[idx + 5];
        const dz = ephemerisData[idx + 3] - ephemerisData[idx + 6];
        const separation = Math.sqrt(dx*dx + dy*dy + dz*dz);
        separationData.push({ t: t, separation: separation });
        if (separation > maxSep) maxSep = separation;
    }

    log(`Max separation: ${maxSep.toFixed(2)} m over ${numOrbits} orbits`, 'info');

    // Reference orbit (first orbit)
    const refOrbit = viewer.entities.add({
        name: 'Reference Orbit',
        polyline: {
            positions: orbitPositions,
            width: 2,
            material: new Cesium.PolylineDashMaterialProperty({
                color: Cesium.Color.WHITE.withAlpha(0.4),
                dashLength: 16
            })
        }
    });
    orbitEntities.push(refOrbit);

    // OMM satellite (cyan) - mean element propagation
    const ommSat = viewer.entities.add({
        name: 'OMM (Mean Elements)',
        description: `OMM mean element propagation\nIncludes secular drift rates:\n- RAAN drift\n- Argument of perigee drift`,
        position: ommPositions,
        orientation: new Cesium.VelocityOrientationProperty(ommPositions),
        point: {
            pixelSize: 12,
            color: Cesium.Color.CYAN,
            outlineColor: Cesium.Color.WHITE,
            outlineWidth: 2
        },
        path: {
            show: true,
            leadTime: 0,
            trailTime: period * 0.5,
            width: 2,
            material: new Cesium.PolylineGlowMaterialProperty({
                glowPower: 0.3,
                color: Cesium.Color.CYAN
            })
        },
        label: {
            text: 'OMM',
            font: '12px monospace',
            fillColor: Cesium.Color.CYAN,
            pixelOffset: new Cesium.Cartesian2(0, -15)
        },
        viewFrom: new Cesium.Cartesian3(-50000, 0, -20000)
    });
    orbitEntities.push(ommSat);

    // Two-Body satellite (lime) - Keplerian propagation
    const twoBodySat = viewer.entities.add({
        name: 'Two-Body (Kepler)',
        description: `Pure two-body propagation\nFixed orbital plane\nNo perturbations\n${numOrbits} orbits propagated`,
        position: j2Positions,
        orientation: new Cesium.VelocityOrientationProperty(j2Positions),
        point: {
            pixelSize: 12,
            color: Cesium.Color.LIME,
            outlineColor: Cesium.Color.WHITE,
            outlineWidth: 2
        },
        path: {
            show: true,
            leadTime: 0,
            trailTime: period * 0.5,
            width: 2,
            material: new Cesium.PolylineGlowMaterialProperty({
                glowPower: 0.3,
                color: Cesium.Color.LIME
            })
        },
        label: {
            text: 'Kepler',
            font: '12px monospace',
            fillColor: Cesium.Color.LIME,
            pixelOffset: new Cesium.Cartesian2(0, -15)
        }
    });
    orbitEntities.push(twoBodySat);

    // Return separation data for chart rendering
    return { separationData, totalTime, startTime, maxSep };
}

/**
 * Full OMM vs J2 model comparison visualization setup
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @param {Function} log - Logging function
 * @param {Object} chartContext - Context for chart rendering (optional)
 */
export function showOMMvsJ2Visualization(viewer, orbitEntities, log, chartContext = null) {
    clearOrbitEntities(viewer, orbitEntities);

    const period = 5400;  // ISS-like ~90 min orbit
    const numOrbits = 20;

    const chartData = addOMMvsJ2Visualization(viewer, orbitEntities, period, numOrbits, log, chartContext);

    if (!chartData) {
        log('OMM vs J2 visualization requires Tudat WASM bindings', 'error');
        return {
            name: 'OMM vs J2',
            description: 'ERROR: Tudat bindings required',
            chartData: null
        };
    }

    configureClockForOrbit(viewer, period * numOrbits, null, period / 10);

    viewer.camera.flyTo({
        destination: Cesium.Cartesian3.fromDegrees(0, 0, 25000000),
        orientation: {
            heading: 0,
            pitch: Cesium.Math.toRadians(-90),
            roll: 0
        },
        duration: 1.0
    });

    return {
        name: 'OMM vs J2',
        description: 'Mean vs osculating element propagation',
        chartData: chartData
    };
}
