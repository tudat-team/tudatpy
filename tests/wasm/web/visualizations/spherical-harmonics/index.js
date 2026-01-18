// Spherical Harmonics Gravity Visualization
// Shows J2 gravity field perturbations with multiple orbits

import { configureClockForOrbit, addAnimatedOrbit, clearOrbitEntities } from '../shared/utils.js';

/**
 * Full spherical harmonics visualization setup
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 */
export function showSphericalHarmonicsVisualization(viewer, orbitEntities) {
    clearOrbitEntities(viewer, orbitEntities);

    const period1 = 5400;
    const period2 = 6900;
    configureClockForOrbit(viewer, period1 * 3, null, period1 / 20);

    addAnimatedOrbit(viewer, orbitEntities, {
        name: 'Low LEO (J2 dominated)',
        semiMajorAxis: 6800,
        eccentricity: 0.001,
        inclination: 98.0, // Sun-synchronous
        raan: 0,
        argPeriapsis: 0,
        color: '#ff9f1c',
        period: period1,
        description: 'Sun-synchronous orbit\nJ2 precession ~1°/day\nSpherical harmonics gravity'
    });

    addAnimatedOrbit(viewer, orbitEntities, {
        name: 'Medium LEO',
        semiMajorAxis: 7500,
        eccentricity: 0.01,
        inclination: 45,
        raan: 90,
        argPeriapsis: 0,
        color: '#8b5cf6',
        period: period2,
        description: 'Reference orbit\nfor gravity field comparison'
    });

    return {
        name: 'Spherical Harmonics',
        description: 'J2 gravity field perturbations'
    };
}
