// Reference Frames Visualization
// Shows J2000 / ECLIPJ2000 rotations and coordinate systems

import { configureClockForOrbit, clearOrbitEntities } from '../shared/utils.js';

/**
 * Add reference frames visualization (J2000 equatorial vs ECLIPJ2000 ecliptic)
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 */
export function addReferenceFrames(viewer, orbitEntities) {
    const axisLength = 15000000; // 15,000 km in meters
    const earthRadius = 6371000;
    const planeRadius = 18000000;
    const obliquity = 23.4393 * Math.PI / 180;

    // ===== J2000 EQUATORIAL FRAME (Cyan) =====
    const xAxisEq = viewer.entities.add({
        name: 'J2000 X (Vernal Equinox)',
        polyline: {
            positions: [new Cesium.Cartesian3(earthRadius, 0, 0), new Cesium.Cartesian3(axisLength, 0, 0)],
            width: 5,
            material: Cesium.Color.CYAN
        }
    });
    orbitEntities.push(xAxisEq);

    const yAxisEq = viewer.entities.add({
        name: 'J2000 Y (Equatorial)',
        polyline: {
            positions: [new Cesium.Cartesian3(0, earthRadius, 0), new Cesium.Cartesian3(0, axisLength, 0)],
            width: 4,
            material: Cesium.Color.CYAN.withAlpha(0.7)
        }
    });
    orbitEntities.push(yAxisEq);

    const zAxisEq = viewer.entities.add({
        name: 'J2000 Z (North Pole)',
        polyline: {
            positions: [new Cesium.Cartesian3(0, 0, earthRadius), new Cesium.Cartesian3(0, 0, axisLength)],
            width: 4,
            material: Cesium.Color.CYAN.withAlpha(0.7)
        }
    });
    orbitEntities.push(zAxisEq);

    // Equatorial plane circle
    const eqPositions = [];
    for (let a = 0; a <= 360; a += 3) {
        const rad = a * Math.PI / 180;
        eqPositions.push(new Cesium.Cartesian3(
            planeRadius * Math.cos(rad),
            planeRadius * Math.sin(rad),
            0
        ));
    }
    const eqPlane = viewer.entities.add({
        name: 'Equatorial Plane',
        polyline: {
            positions: eqPositions,
            width: 2,
            material: Cesium.Color.CYAN.withAlpha(0.5)
        }
    });
    orbitEntities.push(eqPlane);

    // ===== ECLIPJ2000 ECLIPTIC FRAME (Orange) =====
    const yEclEnd = new Cesium.Cartesian3(
        0,
        axisLength * Math.cos(obliquity),
        axisLength * Math.sin(obliquity)
    );
    const yAxisEcl = viewer.entities.add({
        name: 'ECLIPJ2000 Y',
        polyline: {
            positions: [
                new Cesium.Cartesian3(0, earthRadius * Math.cos(obliquity), earthRadius * Math.sin(obliquity)),
                yEclEnd
            ],
            width: 4,
            material: Cesium.Color.ORANGE.withAlpha(0.7)
        }
    });
    orbitEntities.push(yAxisEcl);

    const zEclEnd = new Cesium.Cartesian3(
        0,
        -axisLength * Math.sin(obliquity),
        axisLength * Math.cos(obliquity)
    );
    const zAxisEcl = viewer.entities.add({
        name: 'ECLIPJ2000 Z (Ecliptic Pole)',
        polyline: {
            positions: [
                new Cesium.Cartesian3(0, -earthRadius * Math.sin(obliquity), earthRadius * Math.cos(obliquity)),
                zEclEnd
            ],
            width: 4,
            material: Cesium.Color.ORANGE.withAlpha(0.7)
        }
    });
    orbitEntities.push(zAxisEcl);

    // Ecliptic plane circle
    const eclPositions = [];
    for (let a = 0; a <= 360; a += 3) {
        const rad = a * Math.PI / 180;
        const x = planeRadius * Math.cos(rad);
        const yTemp = planeRadius * Math.sin(rad);
        const y = yTemp * Math.cos(obliquity);
        const z = yTemp * Math.sin(obliquity);
        eclPositions.push(new Cesium.Cartesian3(x, y, z));
    }
    const eclPlane = viewer.entities.add({
        name: 'Ecliptic Plane',
        polyline: {
            positions: eclPositions,
            width: 2,
            material: Cesium.Color.ORANGE.withAlpha(0.5)
        }
    });
    orbitEntities.push(eclPlane);

    // ===== LABELS =====
    const j2000Label = viewer.entities.add({
        position: new Cesium.Cartesian3(0, 0, axisLength * 1.1),
        label: {
            text: 'J2000 Z\n(Celestial Pole)',
            font: '12px monospace',
            fillColor: Cesium.Color.CYAN,
            heightReference: Cesium.HeightReference.NONE
        }
    });
    orbitEntities.push(j2000Label);

    const eclPoleLabel = viewer.entities.add({
        position: zEclEnd,
        label: {
            text: 'ECLIP Z\n(Ecliptic Pole)',
            font: '12px monospace',
            fillColor: Cesium.Color.ORANGE,
            heightReference: Cesium.HeightReference.NONE,
            pixelOffset: new Cesium.Cartesian2(10, 0)
        }
    });
    orbitEntities.push(eclPoleLabel);

    const vernalLabel = viewer.entities.add({
        position: new Cesium.Cartesian3(axisLength * 1.05, 0, 0),
        label: {
            text: 'X (Vernal Equinox)\n♈ First Point of Aries',
            font: '12px monospace',
            fillColor: Cesium.Color.WHITE,
            heightReference: Cesium.HeightReference.NONE
        }
    });
    orbitEntities.push(vernalLabel);

    const obliqLabel = viewer.entities.add({
        position: new Cesium.Cartesian3(0, axisLength * 0.6, axisLength * 0.3),
        label: {
            text: '← 23.4° obliquity →',
            font: '11px monospace',
            fillColor: Cesium.Color.YELLOW,
            heightReference: Cesium.HeightReference.NONE
        }
    });
    orbitEntities.push(obliqLabel);

    const eqPlaneLabel = viewer.entities.add({
        position: new Cesium.Cartesian3(planeRadius * 0.7, planeRadius * 0.7, 0),
        label: {
            text: 'Equatorial Plane',
            font: '11px monospace',
            fillColor: Cesium.Color.CYAN,
            heightReference: Cesium.HeightReference.NONE
        }
    });
    orbitEntities.push(eqPlaneLabel);

    const eclPlaneLabel = viewer.entities.add({
        position: new Cesium.Cartesian3(planeRadius * 0.7, planeRadius * 0.5, planeRadius * 0.35),
        label: {
            text: 'Ecliptic Plane',
            font: '11px monospace',
            fillColor: Cesium.Color.ORANGE,
            heightReference: Cesium.HeightReference.NONE
        }
    });
    orbitEntities.push(eclPlaneLabel);
}

/**
 * Full reference frames visualization setup
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 */
export function showReferenceFramesVisualization(viewer, orbitEntities) {
    clearOrbitEntities(viewer, orbitEntities);

    addReferenceFrames(viewer, orbitEntities);

    configureClockForOrbit(viewer, 86400, null, 100);

    return {
        name: 'Reference Frames',
        description: 'J2000 / ECLIPJ2000 rotations'
    };
}
