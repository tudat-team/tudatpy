// Lambert Targeting Visualization
// Shows orbital transfer arc computation between two positions

import { configureClockForOrbit, clearOrbitEntities, zoomToFitOrbit } from '../shared/utils.js';

/**
 * Add Lambert transfer arc visualization
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 */
export function addLambertTransferVisualization(viewer, orbitEntities) {
    const earthRadius = 6378.136; // km
    const distanceUnit = earthRadius; // Canonical DU

    // Departure position: (2 R_E, 0, 0)
    const r1 = { x: 2 * distanceUnit, y: 0, z: 0 };
    // Arrival position: (2 R_E, 2√3 R_E, 0)
    const r2 = { x: 2 * distanceUnit, y: 2 * Math.sqrt(3) * distanceUnit, z: 0 };

    // Departure orbit (circular at 2 R_E)
    const depOrbitPositions = [];
    for (let nu = 0; nu <= 360; nu += 5) {
        const rad = nu * Cesium.Math.toRadians(1);
        depOrbitPositions.push(new Cesium.Cartesian3(
            2 * distanceUnit * 1000 * Math.cos(rad),
            2 * distanceUnit * 1000 * Math.sin(rad),
            0
        ));
    }
    const depOrbit = viewer.entities.add({
        name: 'Departure Orbit',
        polyline: {
            positions: depOrbitPositions,
            width: 2,
            material: Cesium.Color.GREEN.withAlpha(0.5)
        }
    });
    orbitEntities.push(depOrbit);

    // Transfer arc (elliptical segment from test case)
    const transferPositions = [];
    const a = 2.5 * distanceUnit;
    const e = 0.2;
    for (let t = 0; t <= 1; t += 0.02) {
        const angle = t * Math.PI / 3;
        const r = a * (1 - e * e) / (1 + e * Math.cos(angle - Math.PI / 6));
        transferPositions.push(new Cesium.Cartesian3(
            r * 1000 * Math.cos(angle),
            r * 1000 * Math.sin(angle),
            0
        ));
    }
    const transferArc = viewer.entities.add({
        name: 'Lambert Transfer Arc',
        polyline: {
            positions: transferPositions,
            width: 4,
            material: new Cesium.PolylineGlowMaterialProperty({
                glowPower: 0.3,
                color: Cesium.Color.ORANGE
            })
        }
    });
    orbitEntities.push(transferArc);

    // Departure point
    const depPoint = viewer.entities.add({
        name: 'Departure',
        position: new Cesium.Cartesian3(r1.x * 1000, r1.y * 1000, r1.z * 1000),
        point: { pixelSize: 12, color: Cesium.Color.LIME },
        label: {
            text: 'Departure\nV = (2736, 6594) m/s',
            font: '10px monospace',
            fillColor: Cesium.Color.LIME,
            pixelOffset: new Cesium.Cartesian2(15, 0)
        }
    });
    orbitEntities.push(depPoint);

    // Arrival point
    const arrPoint = viewer.entities.add({
        name: 'Arrival',
        position: new Cesium.Cartesian3(r2.x * 1000, r2.y * 1000, r2.z * 1000),
        point: { pixelSize: 12, color: Cesium.Color.RED },
        label: {
            text: 'Arrival\nV = (-1368, 4225) m/s',
            font: '10px monospace',
            fillColor: Cesium.Color.RED,
            pixelOffset: new Cesium.Cartesian2(15, 0)
        }
    });
    orbitEntities.push(arrPoint);

    // Velocity vectors (scaled for visibility)
    const vScale = 1000;
    const depVel = { x: 2735.8, y: 6594.3 };
    const arrVel = { x: -1367.9, y: 4225.03 };

    const depVelArrow = viewer.entities.add({
        polyline: {
            positions: [
                new Cesium.Cartesian3(r1.x * 1000, r1.y * 1000, 0),
                new Cesium.Cartesian3(r1.x * 1000 + depVel.x * vScale, r1.y * 1000 + depVel.y * vScale, 0)
            ],
            width: 3,
            material: Cesium.Color.LIME
        }
    });
    orbitEntities.push(depVelArrow);

    const arrVelArrow = viewer.entities.add({
        polyline: {
            positions: [
                new Cesium.Cartesian3(r2.x * 1000, r2.y * 1000, 0),
                new Cesium.Cartesian3(r2.x * 1000 + arrVel.x * vScale, r2.y * 1000 + arrVel.y * vScale, 0)
            ],
            width: 3,
            material: Cesium.Color.RED
        }
    });
    orbitEntities.push(arrVelArrow);
}

/**
 * Full Lambert targeting visualization setup
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 */
export function showLambertVisualization(viewer, orbitEntities) {
    clearOrbitEntities(viewer, orbitEntities);

    addLambertTransferVisualization(viewer, orbitEntities);

    const transferTime = 4034; // From test case
    configureClockForOrbit(viewer, transferTime * 2, null, transferTime / 30);
    zoomToFitOrbit(viewer, 15000, 0.2);

    return {
        name: 'Lambert Targeting',
        description: 'Orbital transfer arc computation'
    };
}
