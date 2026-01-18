// Atmospheric Drag Visualization
// Shows NRLMSISE-00 density model effects on LEO orbits

import { configureClockForOrbit, addAnimatedOrbit, addAtmosphereShell, clearOrbitEntities } from '../shared/utils.js';

/**
 * Full atmospheric drag visualization setup
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 */
export function showAtmosphericDragVisualization(viewer, orbitEntities) {
    clearOrbitEntities(viewer, orbitEntities);

    const period = 5400;
    configureClockForOrbit(viewer, period, null, period / 30);

    addAnimatedOrbit(viewer, orbitEntities, {
        name: 'LEO (Drag Region)',
        semiMajorAxis: 6600, // ~230km altitude - significant drag
        eccentricity: 0.001,
        inclination: 51.6, // ISS-like
        raan: 0,
        argPeriapsis: 0,
        color: '#00f0ff',
        period: period,
        description: 'Low Earth orbit in drag region\nAltitude ~230 km\nNRLMSISE-00 / Exponential atmosphere'
    });

    addAtmosphereShell(viewer, orbitEntities);

    return {
        name: 'Atmospheric Drag',
        description: 'NRLMSISE-00 density model'
    };
}
