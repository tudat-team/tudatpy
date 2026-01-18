// Shared visualization utilities for Tudat WASM test runner
// These utilities are used by multiple visualization modules

/**
 * Configure the Cesium clock for a specific orbital period
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {number} periodSeconds - Orbital period in seconds
 * @param {Cesium.JulianDate} epochDate - Optional epoch date (defaults to J2000)
 * @param {number} multiplier - Optional time multiplier
 */
export function configureClockForOrbit(viewer, periodSeconds, epochDate = null, multiplier = null) {
    const clock = viewer.clock;
    const start = epochDate || Cesium.JulianDate.fromIso8601('2000-01-01T12:00:00Z');
    const stop = Cesium.JulianDate.addSeconds(start, periodSeconds, new Cesium.JulianDate());

    clock.startTime = start;
    clock.currentTime = Cesium.JulianDate.clone(start);
    clock.stopTime = stop;
    clock.clockRange = Cesium.ClockRange.LOOP_STOP;
    clock.multiplier = multiplier || Math.max(1, periodSeconds / 60);
    clock.shouldAnimate = true;

    viewer.timeline.zoomTo(start, stop);

    if (viewer.clockViewModel) {
        viewer.clockViewModel.shouldAnimate = true;
    }
}

/**
 * Zoom camera to fit an orbit based on semi-major axis
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {number} semiMajorAxis - Semi-major axis in km
 * @param {number} eccentricity - Orbital eccentricity (default 0)
 */
export function zoomToFitOrbit(viewer, semiMajorAxis, eccentricity = 0) {
    const apoapsis = semiMajorAxis * (1 + eccentricity);
    const cameraDistance = apoapsis * 3.5 * 1000;

    viewer.trackedEntity = undefined;

    const direction = new Cesium.Cartesian3(-0.5, -0.3, -0.7);
    Cesium.Cartesian3.normalize(direction, direction);

    viewer.camera.flyTo({
        destination: new Cesium.Cartesian3(
            cameraDistance * 0.5,
            cameraDistance * 0.3,
            cameraDistance * 0.7
        ),
        orientation: {
            direction: direction,
            up: new Cesium.Cartesian3(0, 0, 1)
        },
        duration: 1.0
    });
}

/**
 * Add an animated orbit with time-varying satellite position
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities for cleanup
 * @param {Object} params - Orbit parameters
 */
export function addAnimatedOrbit(viewer, orbitEntities, params) {
    const {
        name,
        semiMajorAxis,
        eccentricity,
        inclination,
        raan = 0,
        argPeriapsis = 0,
        trueAnomaly = 0,
        color,
        period,
        description,
        referenceFrame = 'INERTIAL'
    } = params;

    const useFixedFrame = referenceFrame === 'FIXED' ||
        (Math.abs(inclination) < 0.1 && eccentricity < 0.01 && Math.abs(period - 86164) < 100);

    const a = semiMajorAxis;
    const e = eccentricity;
    const i = inclination * Cesium.Math.toRadians(1);
    const omega = argPeriapsis * Cesium.Math.toRadians(1);
    const Omega = raan * Cesium.Math.toRadians(1);

    const p = a * (1 - e * e);

    if (useFixedFrame) {
        const geoPositions = [];
        for (let lon = 0; lon <= 360; lon += 2) {
            geoPositions.push(Cesium.Cartesian3.fromDegrees(lon, 0, (semiMajorAxis - 6371) * 1000));
        }
        const orbitEntity = viewer.entities.add({
            name: name + ' Orbit',
            polyline: {
                positions: geoPositions,
                width: 10,
                material: new Cesium.PolylineGlowMaterialProperty({
                    glowPower: 0.3,
                    color: Cesium.Color.fromCssColorString(color).withAlpha(0.9)
                })
            }
        });
        orbitEntities.push(orbitEntity);
    }

    zoomToFitOrbit(viewer, semiMajorAxis, eccentricity);

    const clock = viewer.clock;
    const startTime = clock.startTime;
    let positionProperty;

    const numSamplesPerPeriod = 360;
    const numPeriods = 10;

    if (useFixedFrame) {
        const geoLongitude = raan;
        const geoAltitude = (semiMajorAxis - 6371) * 1000;
        const fixedPosition = Cesium.Cartesian3.fromDegrees(geoLongitude, 0, geoAltitude);
        positionProperty = new Cesium.ConstantPositionProperty(fixedPosition, Cesium.ReferenceFrame.FIXED);
    } else {
        positionProperty = new Cesium.SampledPositionProperty(Cesium.ReferenceFrame.INERTIAL);

        for (let orbitNum = 0; orbitNum < numPeriods; orbitNum++) {
            for (let sample = 0; sample <= numSamplesPerPeriod; sample++) {
                const fraction = sample / numSamplesPerPeriod;
                const totalTime = (orbitNum + fraction) * period;
                const time = Cesium.JulianDate.addSeconds(startTime, totalTime, new Cesium.JulianDate());

                const M = (2 * Math.PI * fraction) + (trueAnomaly * Cesium.Math.toRadians(1));

                let E = M;
                for (let iter = 0; iter < 10; iter++) {
                    E = M + e * Math.sin(E);
                }

                const nu = 2 * Math.atan2(
                    Math.sqrt(1 + e) * Math.sin(E / 2),
                    Math.sqrt(1 - e) * Math.cos(E / 2)
                );

                const r = p / (1 + e * Math.cos(nu));

                const xPQW = r * Math.cos(nu);
                const yPQW = r * Math.sin(nu);

                const cosO = Math.cos(Omega), sinO = Math.sin(Omega);
                const cosi = Math.cos(i), sini = Math.sin(i);
                const cosw = Math.cos(omega), sinw = Math.sin(omega);

                const x = (cosO * cosw - sinO * sinw * cosi) * xPQW + (-cosO * sinw - sinO * cosw * cosi) * yPQW;
                const y = (sinO * cosw + cosO * sinw * cosi) * xPQW + (-sinO * sinw + cosO * cosw * cosi) * yPQW;
                const zPos = (sinw * sini) * xPQW + (cosw * sini) * yPQW;

                positionProperty.addSample(time, new Cesium.Cartesian3(x * 1000, y * 1000, zPos * 1000));
            }
        }
    }

    let orientationProperty;
    if (!useFixedFrame) {
        orientationProperty = new Cesium.VelocityOrientationProperty(positionProperty);
    }

    const satEntity = viewer.entities.add({
        name: name,
        description: description,
        position: positionProperty,
        orientation: orientationProperty,
        point: {
            pixelSize: 12,
            color: Cesium.Color.fromCssColorString(color),
            outlineColor: Cesium.Color.WHITE,
            outlineWidth: 2
        },
        path: useFixedFrame ? undefined : {
            show: true,
            leadTime: 0,
            trailTime: period,
            width: 3,
            material: Cesium.Color.fromCssColorString(color)
        },
        label: {
            text: name,
            font: '12px monospace',
            style: Cesium.LabelStyle.FILL_AND_OUTLINE,
            outlineWidth: 2,
            verticalOrigin: Cesium.VerticalOrigin.BOTTOM,
            pixelOffset: new Cesium.Cartesian2(0, -15),
            fillColor: Cesium.Color.fromCssColorString(color),
            outlineColor: Cesium.Color.BLACK
        },
        viewFrom: new Cesium.Cartesian3(-50000, 0, 20000)
    });
    orbitEntities.push(satEntity);
}

/**
 * Add Moon's orbit for cislunar visualizations
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 * @returns {Cesium.Entity} Moon entity reference
 */
export function addMoonOrbit(viewer, orbitEntities) {
    const moonDistance = 384400;
    const moonInclination = 5.145;

    const positions = [];
    for (let nu = 0; nu <= 360; nu += 5) {
        const nuRad = nu * Cesium.Math.toRadians(1);
        const iRad = moonInclination * Cesium.Math.toRadians(1);

        const x = moonDistance * Math.cos(nuRad);
        const y = moonDistance * Math.sin(nuRad) * Math.cos(iRad);
        const z = moonDistance * Math.sin(nuRad) * Math.sin(iRad);

        positions.push(new Cesium.Cartesian3(x * 1000, y * 1000, z * 1000));
    }

    const moonOrbitEntity = viewer.entities.add({
        name: 'Moon Orbit',
        polyline: {
            positions: positions,
            width: 2,
            material: Cesium.Color.GRAY.withAlpha(0.4)
        }
    });
    orbitEntities.push(moonOrbitEntity);

    const moonEntity = viewer.entities.add({
        name: 'Moon',
        position: new Cesium.Cartesian3(moonDistance * 1000, 0, 0),
        point: {
            pixelSize: 15,
            color: Cesium.Color.LIGHTGRAY,
            outlineColor: Cesium.Color.WHITE,
            outlineWidth: 2
        },
        label: {
            text: 'Moon',
            font: '12px monospace',
            fillColor: Cesium.Color.WHITE,
            pixelOffset: new Cesium.Cartesian2(0, -18)
        }
    });
    orbitEntities.push(moonEntity);

    return moonEntity;
}

/**
 * Clear all orbit entities from the viewer
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array of entities to remove
 */
export function clearOrbitEntities(viewer, orbitEntities) {
    orbitEntities.forEach(entity => {
        viewer.entities.remove(entity);
    });
    orbitEntities.length = 0;
}

/**
 * Add atmosphere shell visualization for drag demonstrations
 * @param {Cesium.Viewer} viewer - Cesium viewer instance
 * @param {Array} orbitEntities - Array to track created entities
 */
export function addAtmosphereShell(viewer, orbitEntities) {
    const earthRadius = 6371;
    const layers = [
        { altitude: 100, color: '#ff6b6b', alpha: 0.15, name: 'Thermosphere (100km)' },
        { altitude: 200, color: '#ffd93d', alpha: 0.12, name: 'Upper Atmosphere (200km)' },
        { altitude: 400, color: '#00f0ff', alpha: 0.08, name: 'LEO Region (400km)' }
    ];

    layers.forEach(layer => {
        const positions = [];
        const r = (earthRadius + layer.altitude) * 1000;

        for (let lat = -80; lat <= 80; lat += 20) {
            const latRad = lat * Cesium.Math.toRadians(1);
            const ringR = r * Math.cos(latRad);
            const z = r * Math.sin(latRad);

            for (let lon = 0; lon <= 360; lon += 10) {
                const lonRad = lon * Cesium.Math.toRadians(1);
                positions.push(new Cesium.Cartesian3(
                    ringR * Math.cos(lonRad),
                    ringR * Math.sin(lonRad),
                    z
                ));
            }
        }

        const shellEntity = viewer.entities.add({
            name: layer.name,
            polyline: {
                positions: positions,
                width: 1,
                material: Cesium.Color.fromCssColorString(layer.color).withAlpha(layer.alpha)
            }
        });
        orbitEntities.push(shellEntity);
    });
}
