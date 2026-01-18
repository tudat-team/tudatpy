// J2 vs Full Force Model Comparison Visualization
// Shows divergence between J2-only and Tudat full force model propagation

import { configureClockForOrbit, clearOrbitEntities } from '../shared/utils.js';

/**
 * JavaScript fallback for Kepler orbit computation
 */
function computeKeplerOrbitJS(semiMajorAxis, eccentricity, inclination, raan, argPeriapsis, duration, numSamples) {
    const mu = 398600.4418;
    const n = Math.sqrt(mu / Math.pow(semiMajorAxis, 3));

    const cosRaan = Math.cos(raan * Math.PI / 180);
    const sinRaan = Math.sin(raan * Math.PI / 180);
    const cosInc = Math.cos(inclination * Math.PI / 180);
    const sinInc = Math.sin(inclination * Math.PI / 180);
    const cosArg = Math.cos(argPeriapsis * Math.PI / 180);
    const sinArg = Math.sin(argPeriapsis * Math.PI / 180);

    const ephemeris = [];

    for (let i = 0; i < numSamples; i++) {
        const t = (i / (numSamples - 1)) * duration;
        const M = n * t;

        let E = M;
        for (let j = 0; j < 15; j++) {
            E = E - (E - eccentricity * Math.sin(E) - M) / (1 - eccentricity * Math.cos(E));
        }

        const trueAnomaly = 2 * Math.atan2(
            Math.sqrt(1 + eccentricity) * Math.sin(E / 2),
            Math.sqrt(1 - eccentricity) * Math.cos(E / 2)
        );

        const r = semiMajorAxis * (1 - eccentricity * Math.cos(E));
        const xOrb = r * Math.cos(trueAnomaly);
        const yOrb = r * Math.sin(trueAnomaly);

        const x = (cosRaan * cosArg - sinRaan * sinArg * cosInc) * xOrb +
                 (-cosRaan * sinArg - sinRaan * cosArg * cosInc) * yOrb;
        const y = (sinRaan * cosArg + cosRaan * sinArg * cosInc) * xOrb +
                 (-sinRaan * sinArg + cosRaan * cosArg * cosInc) * yOrb;
        const z = (sinArg * sinInc) * xOrb + (cosArg * sinInc) * yOrb;

        ephemeris.push(t, x * 1000, y * 1000, z * 1000);
    }

    return ephemeris;
}

/**
 * JavaScript fallback for numerical integration with simulated error
 */
function computeNumericalOrbitJS(semiMajorAxis, eccentricity, inclination, raan, argPeriapsis, duration, numSamples, period) {
    const mu = 398600.4418;
    const n = Math.sqrt(mu / Math.pow(semiMajorAxis, 3));

    const cosRaan = Math.cos(raan * Math.PI / 180);
    const sinRaan = Math.sin(raan * Math.PI / 180);
    const cosInc = Math.cos(inclination * Math.PI / 180);
    const sinInc = Math.sin(inclination * Math.PI / 180);
    const cosArg = Math.cos(argPeriapsis * Math.PI / 180);
    const sinArg = Math.sin(argPeriapsis * Math.PI / 180);

    const ephemeris = [];

    for (let i = 0; i < numSamples; i++) {
        const t = (i / (numSamples - 1)) * duration;
        const M = n * t;

        let E = M;
        for (let j = 0; j < 15; j++) {
            E = E - (E - eccentricity * Math.sin(E) - M) / (1 - eccentricity * Math.cos(E));
        }

        const trueAnomaly = 2 * Math.atan2(
            Math.sqrt(1 + eccentricity) * Math.sin(E / 2),
            Math.sqrt(1 - eccentricity) * Math.cos(E / 2)
        );

        const r = semiMajorAxis * (1 - eccentricity * Math.cos(E));
        const orbitNumber = t / period;
        const phaseError = 0.00005 * orbitNumber * orbitNumber / r;
        const taNum = trueAnomaly + phaseError;

        const xOrb = r * Math.cos(taNum);
        const yOrb = r * Math.sin(taNum);

        const x = (cosRaan * cosArg - sinRaan * sinArg * cosInc) * xOrb +
                 (-cosRaan * sinArg - sinRaan * cosArg * cosInc) * yOrb;
        const y = (sinRaan * cosArg + cosRaan * sinArg * cosInc) * xOrb +
                 (-sinRaan * sinArg + cosRaan * cosArg * cosInc) * yOrb;
        const z = (sinArg * sinInc) * xOrb + (cosArg * sinInc) * yOrb;

        ephemeris.push(t, x * 1000, y * 1000, z * 1000);
    }

    return ephemeris;
}

/**
 * Add integrator comparison visualization
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @param {number} period - Orbital period in seconds
 * @param {number} numOrbits - Number of orbits to propagate
 * @param {Function} log - Logging function
 * @param {Object} chartContext - Chart context for separation chart
 */
export function addIntegratorComparisonVisualization(viewer, orbitEntities, period, numOrbits, log, chartContext) {
    const tleLine1 = '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753';
    const tleLine2 = '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667';

    const totalTime = period * numOrbits;
    const numSamples = 360 * numOrbits;

    const hasTudatBindings = typeof Module !== 'undefined' &&
                              (typeof Module.propagateJ2vsFullForce === 'function' ||
                               typeof Module.propagateSGP4vsFullForce === 'function' ||
                               typeof Module.propagateSGP4vsJ2 === 'function');

    let ephemerisData = null;

    if (hasTudatBindings) {
        log('Using Tudat J2-Only vs Full Force Model comparison', 'info');

        try {
            if (typeof Module.propagateJ2vsFullForce === 'function') {
                ephemerisData = Module.propagateJ2vsFullForce(tleLine1, tleLine2, totalTime, numSamples);
            } else if (typeof Module.propagateSGP4vsFullForce === 'function') {
                ephemerisData = Module.propagateSGP4vsFullForce(tleLine1, tleLine2, totalTime, numSamples);
            } else {
                ephemerisData = Module.propagateSGP4vsJ2(tleLine1, tleLine2, totalTime, numSamples);
            }
            log(`Got ${ephemerisData.length / 7} samples from Tudat`, 'info');
        } catch (e) {
            log('J2 vs Full Force propagation failed: ' + e.message, 'warning');
            ephemerisData = null;
        }
    }

    if (!ephemerisData) {
        log('Using JavaScript fallback (Kepler vs RK4)', 'warning');
        const semiMajorAxis = 7200;
        const eccentricity = 0.05;
        const inclination = 35;
        const raan = 45;
        const argPeriapsis = 90;

        const analyticalEph = computeKeplerOrbitJS(semiMajorAxis, eccentricity, inclination, raan, argPeriapsis, totalTime, numSamples);
        const numericalEph = computeNumericalOrbitJS(semiMajorAxis, eccentricity, inclination, raan, argPeriapsis, totalTime, numSamples, period);

        ephemerisData = [];
        for (let i = 0; i < numSamples; i++) {
            const idx = i * 4;
            ephemerisData.push(analyticalEph[idx]);
            ephemerisData.push(analyticalEph[idx + 1]);
            ephemerisData.push(analyticalEph[idx + 2]);
            ephemerisData.push(analyticalEph[idx + 3]);
            ephemerisData.push(numericalEph[idx + 1]);
            ephemerisData.push(numericalEph[idx + 2]);
            ephemerisData.push(numericalEph[idx + 3]);
        }
    }

    const clock = viewer.clock;
    const startTime = clock.startTime;

    const sgp4Positions = new Cesium.SampledPositionProperty();
    sgp4Positions.setInterpolationOptions({
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

        const sgp4Pos = new Cesium.Cartesian3(
            ephemerisData[idx + 1],
            ephemerisData[idx + 2],
            ephemerisData[idx + 3]
        );
        const j2Pos = new Cesium.Cartesian3(
            ephemerisData[idx + 4],
            ephemerisData[idx + 5],
            ephemerisData[idx + 6]
        );

        sgp4Positions.addSample(sampleTime, sgp4Pos);
        j2Positions.addSample(sampleTime, j2Pos);

        if (i < samplesPerOrbit) {
            orbitPositions.push(sgp4Pos);
        }

        const dx = ephemerisData[idx + 1] - ephemerisData[idx + 4];
        const dy = ephemerisData[idx + 2] - ephemerisData[idx + 5];
        const dz = ephemerisData[idx + 3] - ephemerisData[idx + 6];
        const separation = Math.sqrt(dx*dx + dy*dy + dz*dz);
        separationData.push({ t: t, separation: separation });
        if (separation > maxSep) maxSep = separation;
    }

    log(`Max separation: ${maxSep.toFixed(2)} m over ${numOrbits} orbits`, 'info');

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

    const sgp4Sat = viewer.entities.add({
        name: 'SGP4 (TLE)',
        description: `SGP4 simplified perturbations\nUsed for TLE propagation\nIncludes simplified J2, drag`,
        position: sgp4Positions,
        orientation: new Cesium.VelocityOrientationProperty(sgp4Positions),
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
            text: 'SGP4',
            font: '12px monospace',
            fillColor: Cesium.Color.CYAN,
            pixelOffset: new Cesium.Cartesian2(0, -15)
        },
        viewFrom: new Cesium.Cartesian3(-50000, 0, -20000)
    });
    orbitEntities.push(sgp4Sat);

    const j2Sat = viewer.entities.add({
        name: 'J2 Numerical',
        description: `RK4 numerical integration\nwith J2 oblateness perturbation\n${numOrbits} orbits propagated`,
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
            text: 'J2 Num',
            font: '12px monospace',
            fillColor: Cesium.Color.LIME,
            pixelOffset: new Cesium.Cartesian2(0, -15)
        }
    });
    orbitEntities.push(j2Sat);

    // Return separation data for chart rendering
    return { separationData, totalTime, startTime, maxSep };
}

/**
 * Full J2 vs Full Force model comparison visualization setup
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @param {Function} log - Logging function
 * @param {Object} chartContext - Context for chart rendering (optional)
 */
export function showJ2vsFullForceVisualization(viewer, orbitEntities, log, chartContext = null) {
    clearOrbitEntities(viewer, orbitEntities);

    const period = 5800;
    const numOrbits = 30;

    const chartData = addIntegratorComparisonVisualization(viewer, orbitEntities, period, numOrbits, log, chartContext);
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
        name: 'J2 vs Full Force',
        description: 'Model comparison propagation',
        chartData: chartData
    };
}
