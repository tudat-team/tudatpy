// Libration Points Visualization
// Shows Earth-Moon L1-L5 Lagrange points

import { configureClockForOrbit, addMoonOrbit, clearOrbitEntities } from '../shared/utils.js';

/**
 * Add Libration points visualization (L1-L5 Lagrange points)
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @param {Function} log - Logging function
 */
export function addLibrationPointsVisualization(viewer, orbitEntities, log) {
    const moonDist = 384400; // km
    const mu = 0.01215; // Earth-Moon mass ratio

    // Check if Tudat bindings are available
    const hasTudatBindings = typeof Module !== 'undefined' &&
                              typeof Module.computeLibrationPoints === 'function';

    let lagrangePoints;

    if (hasTudatBindings) {
        log('Using Tudat library for libration points', 'info');

        try {
            const lpData = Module.computeLibrationPoints(mu);
            // lpData contains [L1x,L1y,L1z, L2x,L2y,L2z, L3x,L3y,L3z, L4x,L4y,L4z, L5x,L5y,L5z]

            lagrangePoints = [
                { name: 'L1', pos: new Cesium.Cartesian3(lpData[0] * moonDist * 1000, lpData[1] * moonDist * 1000, lpData[2] * moonDist * 1000), color: Cesium.Color.RED },
                { name: 'L2', pos: new Cesium.Cartesian3(lpData[3] * moonDist * 1000, lpData[4] * moonDist * 1000, lpData[5] * moonDist * 1000), color: Cesium.Color.RED },
                { name: 'L3', pos: new Cesium.Cartesian3(lpData[6] * moonDist * 1000, lpData[7] * moonDist * 1000, lpData[8] * moonDist * 1000), color: Cesium.Color.RED },
                { name: 'L4', pos: new Cesium.Cartesian3(lpData[9] * moonDist * 1000, lpData[10] * moonDist * 1000, lpData[11] * moonDist * 1000), color: Cesium.Color.LIME },
                { name: 'L5', pos: new Cesium.Cartesian3(lpData[12] * moonDist * 1000, lpData[13] * moonDist * 1000, lpData[14] * moonDist * 1000), color: Cesium.Color.LIME }
            ];
        } catch (e) {
            log('Tudat libration points failed: ' + e.message + ', using JS fallback', 'warning');
            lagrangePoints = null;
        }
    }

    // Fallback to JavaScript approximations
    if (!lagrangePoints) {
        const L1x = moonDist * 0.8369;
        const L2x = moonDist * 1.1562;
        const L3x = -moonDist * 1.0051;
        const L4x = moonDist * 0.5 - moonDist * mu;
        const L4y = moonDist * Math.sqrt(3) / 2;
        const L5y = -L4y;

        lagrangePoints = [
            { name: 'L1', pos: new Cesium.Cartesian3(L1x * 1000, 0, 0), color: Cesium.Color.RED },
            { name: 'L2', pos: new Cesium.Cartesian3(L2x * 1000, 0, 0), color: Cesium.Color.RED },
            { name: 'L3', pos: new Cesium.Cartesian3(L3x * 1000, 0, 0), color: Cesium.Color.RED },
            { name: 'L4', pos: new Cesium.Cartesian3(L4x * 1000, L4y * 1000, 0), color: Cesium.Color.LIME },
            { name: 'L5', pos: new Cesium.Cartesian3(L4x * 1000, L5y * 1000, 0), color: Cesium.Color.LIME }
        ];
    }

    const dataSource = hasTudatBindings ? 'Tudat C++' : 'JS approx';

    lagrangePoints.forEach(lp => {
        const entity = viewer.entities.add({
            name: lp.name,
            description: `${lp.name} Lagrange Point\nData source: ${dataSource}`,
            position: lp.pos,
            point: {
                pixelSize: 10,
                color: lp.color,
                outlineColor: Cesium.Color.WHITE,
                outlineWidth: 1
            },
            label: {
                text: lp.name,
                font: '12px monospace',
                fillColor: lp.color,
                pixelOffset: new Cesium.Cartesian2(12, 0)
            }
        });
        orbitEntities.push(entity);
    });

    // Add Moon orbit and Moon
    addMoonOrbit(viewer, orbitEntities);

    // Add Earth label
    const earthLabel = viewer.entities.add({
        position: Cesium.Cartesian3.ZERO,
        label: {
            text: 'Earth',
            font: '12px monospace',
            fillColor: Cesium.Color.CYAN,
            pixelOffset: new Cesium.Cartesian2(0, 20)
        }
    });
    orbitEntities.push(earthLabel);
}

/**
 * Full libration points visualization setup
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @param {Function} log - Logging function
 */
export function showLibrationPointsVisualization(viewer, orbitEntities, log) {
    clearOrbitEntities(viewer, orbitEntities);

    addLibrationPointsVisualization(viewer, orbitEntities, log);

    const period = 86400 * 27; // Lunar month
    configureClockForOrbit(viewer, period, null, period / 120);

    // Zoom out to see all 5 Lagrange points
    viewer.camera.flyTo({
        destination: Cesium.Cartesian3.fromElements(0, 0, 1200000000),
        orientation: {
            heading: 0,
            pitch: Cesium.Math.toRadians(-90),
            roll: 0
        },
        duration: 1.5
    });

    return {
        name: 'Libration Points',
        description: 'Earth-Moon L1-L5 Lagrange points'
    };
}
