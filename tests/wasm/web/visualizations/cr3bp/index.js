// CR3BP (Circular Restricted 3-Body Problem) Visualization
// Shows periodic orbits in the Earth-Moon system using Tudat WASM propagation

import { configureClockForOrbit, addMoonOrbit, clearOrbitEntities } from '../shared/utils.js';

/**
 * CR3BP Orbit configurations from JPL Three-Body Periodic Orbits database
 * https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
 */
export function getCR3BPOrbitConfig(orbitType) {
    const mu = 0.01215058560962404;  // Earth-Moon mass ratio
    const L = 389703.264829278;       // km (Earth-Moon distance)
    const TU = 382981.289129055;      // seconds per time unit

    const orbits = {
        'l2-halo': {
            name: 'L2 Northern Halo',
            description: 'Orbit around L2 point beyond Moon\nUsed by: Lunar Gateway, JWST (Sun-Earth)',
            color: '#8b5cf6',
            x0: 1.0827766352133668,
            y0: 0.0,
            z0: 0.20232381728937715,
            vx0: 0.0,
            vy0: -0.20084062472053987,
            vz0: 0.0,
            period: 2.3807980015152594
        },
        'dro': {
            name: 'Distant Retrograde Orbit',
            description: 'Stable orbit around Moon\nUsed by: Artemis missions',
            color: '#10b981',
            x0: 1.09,
            y0: 0.0,
            z0: 0.0,
            vx0: 0.0,
            vy0: 0.18,
            vz0: 0.0,
            period: 3.4
        },
        'l1-halo': {
            name: 'L1 Northern Halo',
            description: 'Orbit around L1 point between Earth and Moon',
            color: '#f59e0b',
            x0: 0.8369,
            y0: 0.0,
            z0: 0.08,
            vx0: 0.0,
            vy0: -0.15,
            vz0: 0.0,
            period: 2.8
        },
        'lyapunov': {
            name: 'L2 Lyapunov',
            description: 'Planar periodic orbit at L2\nIn-plane oscillation around Lagrange point',
            color: '#ec4899',
            x0: 1.180,
            y0: 0.0,
            z0: 0.0,
            vx0: 0.0,
            vy0: -0.08,
            vz0: 0.0,
            period: 3.4
        }
    };

    return { ...orbits[orbitType], mu, L, TU };
}

/**
 * Add CR3BP trajectory visualization using Tudat 3-body propagation
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @param {Object} config - Orbit configuration from getCR3BPOrbitConfig
 * @param {Function} log - Logging function
 */
export function addCR3BPVisualization(viewer, orbitEntities, config, log) {
    if (!config) {
        config = getCR3BPOrbitConfig('l2-halo');
    }

    const { mu, L, TU, x0, y0, z0, vx0, vy0, vz0, period, name, description, color } = config;

    const duration = period * 1.0;
    const numPoints = 500;
    const earthX = -mu;

    // Require Tudat WASM bindings - NO FALLBACK
    if (typeof Module === 'undefined' || typeof Module.propagateCR3BP !== 'function') {
        log('ERROR: Tudat WASM bindings not available for CR3BP', 'error');
        return;
    }

    log(`Computing ${name} trajectory with Tudat...`, 'info');
    const trajectory = Module.propagateCR3BP(mu, x0, y0, z0, vx0, vy0, vz0, duration, numPoints);
    log(`Got ${trajectory.length / 7} trajectory points`, 'info');

    const positions = [];
    const times = [];

    for (let i = 0; i < numPoints && i * 7 < trajectory.length; i++) {
        const idx = i * 7;
        const t = trajectory[idx];
        const x = trajectory[idx + 1];
        const y = trajectory[idx + 2];
        const z = trajectory[idx + 3];

        times.push(t * TU);

        const xEarth = (x - earthX) * L * 1000;
        const yEarth = y * L * 1000;
        const zEarth = z * L * 1000;

        positions.push(new Cesium.Cartesian3(xEarth, yEarth, zEarth));
    }

    if (positions.length > 0) {
        const orbitColor = Cesium.Color.fromCssColorString(color);

        const trajectoryEntity = viewer.entities.add({
            name: name,
            polyline: {
                positions: positions,
                width: 10,
                material: new Cesium.PolylineGlowMaterialProperty({
                    glowPower: 0.2,
                    color: orbitColor
                })
            },
            description: description
        });
        orbitEntities.push(trajectoryEntity);

        const startTime = viewer.clock.startTime;
        const stopTime = viewer.clock.stopTime;
        const viewerDuration = Cesium.JulianDate.secondsDifference(stopTime, startTime);
        const trajectoryDuration = times[times.length - 1] - times[0];
        const timeScale = trajectoryDuration / viewerDuration;

        const property = new Cesium.SampledPositionProperty();

        for (let i = 0; i < positions.length; i++) {
            const trajTime = times[i] - times[0];
            const viewerTime = trajTime / timeScale;
            const time = Cesium.JulianDate.addSeconds(startTime, viewerTime, new Cesium.JulianDate());
            property.addSample(time, positions[i]);
        }

        const spacecraftEntity = viewer.entities.add({
            name: `${name} Spacecraft`,
            position: property,
            point: {
                pixelSize: 12,
                color: orbitColor
            },
            path: {
                show: false
            },
            label: {
                text: name,
                font: '11px monospace',
                fillColor: orbitColor,
                pixelOffset: new Cesium.Cartesian2(12, 0),
                showBackground: true,
                backgroundColor: Cesium.Color.BLACK.withAlpha(0.7)
            }
        });
        orbitEntities.push(spacecraftEntity);

        log(`${name}: ${(trajectoryDuration / 86400).toFixed(1)} days period`, 'info');
    }
}

/**
 * Set camera view for CR3BP visualization
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {boolean} animate - Whether to animate camera movement
 */
export function setCR3BPCameraView(viewer, animate = true) {
    const moonDistance = 384400000;

    viewer.trackedEntity = undefined;

    const cameraDistance = moonDistance * 1.8;

    const cameraOptions = {
        destination: new Cesium.Cartesian3(cameraDistance, 0, 0),
        orientation: {
            direction: new Cesium.Cartesian3(-1, 0, 0),
            up: new Cesium.Cartesian3(0, 0, 1)
        }
    };

    if (animate) {
        cameraOptions.duration = 1.5;
        viewer.camera.flyTo(cameraOptions);
    } else {
        viewer.camera.setView(cameraOptions);
    }
}

/**
 * Full CR3BP visualization setup
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @param {string} orbitType - Type of orbit ('l2-halo', 'dro', 'l1-halo', 'lyapunov')
 * @param {Function} log - Logging function
 * @param {boolean} animate - Whether to animate camera movement
 */
export function showCR3BPVisualization(viewer, orbitEntities, orbitType, log, animate = true) {
    clearOrbitEntities(viewer, orbitEntities);

    addMoonOrbit(viewer, orbitEntities);

    const config = getCR3BPOrbitConfig(orbitType);
    addCR3BPVisualization(viewer, orbitEntities, config, log);

    const periodDays = (config.period * config.TU / 86400).toFixed(1);
    const period = config.period * config.TU;

    configureClockForOrbit(viewer, period, null, period / 360);
    setCR3BPCameraView(viewer, animate);

    return {
        name: config.name,
        periodDays: periodDays,
        description: config.description
    };
}
