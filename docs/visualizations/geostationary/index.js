// Geostationary Orbit Visualization
// Shows GEO period and mean motion

import { configureClockForOrbit, addAnimatedOrbit, clearOrbitEntities } from '../shared/utils.js';

/**
 * Full geostationary orbit visualization setup
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 */
export function showGeostationaryVisualization(viewer, orbitEntities) {
    clearOrbitEntities(viewer, orbitEntities);

    const period = 86164; // Sidereal day
    configureClockForOrbit(viewer, period, null, period / 60);

    addAnimatedOrbit(viewer, orbitEntities, {
        name: 'Geostationary Orbit',
        semiMajorAxis: 42164,
        eccentricity: 0.0,
        inclination: 0.0,
        raan: 0,
        argPeriapsis: 0,
        color: '#ffd93d',
        period: period,
        description: 'Geostationary orbit test\nPeriod = 86,164s (sidereal day)\nKepler\'s 3rd law: T² ∝ a³',
        referenceFrame: 'FIXED'
    });

    return {
        name: 'Geostationary Orbit',
        description: 'GEO period and mean motion'
    };
}
