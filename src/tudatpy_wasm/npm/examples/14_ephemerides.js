/**
 * Example 14: Ephemerides (Kepler and Tabulated)
 *
 * This example demonstrates Keplerian ephemeris propagation and
 * tabulated ephemeris interpolation for spacecraft and planetary positions.
 *
 * Run with: node 14_ephemerides.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Ephemerides Example ===\n');

    const tudat = await createTudatModule();

    // ========================================
    // Kepler Ephemeris
    // ========================================
    console.log('1. Keplerian Ephemeris');
    console.log('======================\n');

    console.log('A Keplerian ephemeris propagates orbits using the two-body problem.');
    console.log('Given initial Keplerian elements, position/velocity at any time');
    console.log('can be computed analytically.\n');

    const earthGM = 3.986004418e14;  // m^3/s^2

    // ISS-like orbit
    const semiMajorAxis = 6778137.0;  // ~400 km altitude
    const eccentricity = 0.0001;       // Nearly circular
    const inclination = 51.6 * Math.PI / 180;  // ISS inclination
    const raan = 0.0;
    const argPeriapsis = 0.0;
    const trueAnomaly0 = 0.0;

    const orbitalPeriod = 2 * Math.PI * Math.sqrt(Math.pow(semiMajorAxis, 3) / earthGM);
    const meanMotion = 2 * Math.PI / orbitalPeriod;

    console.log('ISS-like orbit parameters:');
    console.log(`  Semi-major axis: ${(semiMajorAxis / 1000).toFixed(2)} km`);
    console.log(`  Eccentricity: ${eccentricity}`);
    console.log(`  Inclination: ${(inclination * 180 / Math.PI).toFixed(1)} degrees`);
    console.log(`  Orbital period: ${(orbitalPeriod / 60).toFixed(2)} minutes`);
    console.log(`  Mean motion: ${(meanMotion * 86400 * 180 / Math.PI).toFixed(2)} deg/day\n`);

    // Propagate Kepler orbit at different times
    console.log('Keplerian propagation (quarter-orbit intervals):');
    console.log('Time (min)    True Anomaly (deg)    Altitude (km)');
    console.log('----------    ------------------    -------------');

    for (let i = 0; i <= 4; i++) {
        const time = i * orbitalPeriod / 4;

        // Convert time to mean anomaly
        const meanAnomaly = meanMotion * time;

        // Solve Kepler's equation (simple iteration for nearly circular orbit)
        let E = meanAnomaly;
        for (let iter = 0; iter < 10; iter++) {
            E = meanAnomaly + eccentricity * Math.sin(E);
        }

        // True anomaly from eccentric anomaly
        const trueAnomaly = 2 * Math.atan2(
            Math.sqrt(1 + eccentricity) * Math.sin(E / 2),
            Math.sqrt(1 - eccentricity) * Math.cos(E / 2)
        );

        // Radius at this true anomaly
        const radius = semiMajorAxis * (1 - eccentricity * eccentricity) /
                       (1 + eccentricity * Math.cos(trueAnomaly));
        const altitude = radius - 6378137.0;

        console.log(`${(time / 60).toFixed(1).padStart(10)}    ${(trueAnomaly * 180 / Math.PI).toFixed(1).padStart(18)}    ${(altitude / 1000).toFixed(2).padStart(13)}`);
    }

    // ========================================
    // Cartesian State from Keplerian Elements
    // ========================================
    console.log('\n\n2. State Vector Computation');
    console.log('===========================\n');

    const kepler = new tudat.Vector6d();
    kepler.set(0, semiMajorAxis);
    kepler.set(1, eccentricity);
    kepler.set(2, inclination);
    kepler.set(3, argPeriapsis);
    kepler.set(4, raan);
    kepler.set(5, trueAnomaly0);

    const cartesian = tudat.astro.element_conversion.keplerian_to_cartesian(kepler, earthGM);

    console.log('Initial state (t=0):');
    console.log(`  Position: [${(cartesian.get(0)/1000).toFixed(3)}, ${(cartesian.get(1)/1000).toFixed(3)}, ${(cartesian.get(2)/1000).toFixed(3)}] km`);
    console.log(`  Velocity: [${(cartesian.get(3)/1000).toFixed(4)}, ${(cartesian.get(4)/1000).toFixed(4)}, ${(cartesian.get(5)/1000).toFixed(4)}] km/s`);

    // State at quarter orbit (true anomaly = 90 deg)
    kepler.set(5, Math.PI / 2);
    const cartesian90 = tudat.astro.element_conversion.keplerian_to_cartesian(kepler, earthGM);

    console.log('\nState at true anomaly = 90°:');
    console.log(`  Position: [${(cartesian90.get(0)/1000).toFixed(3)}, ${(cartesian90.get(1)/1000).toFixed(3)}, ${(cartesian90.get(2)/1000).toFixed(3)}] km`);
    console.log(`  Velocity: [${(cartesian90.get(3)/1000).toFixed(4)}, ${(cartesian90.get(4)/1000).toFixed(4)}, ${(cartesian90.get(5)/1000).toFixed(4)}] km/s`);

    // ========================================
    // Tabulated Ephemeris
    // ========================================
    console.log('\n\n3. Tabulated Ephemeris');
    console.log('======================\n');

    console.log('Tabulated ephemerides store discrete state vectors and interpolate');
    console.log('between them. This is used for complex trajectories or loaded from');
    console.log('external sources (SPICE, numerical propagation results, etc.).\n');

    // Generate a tabulated ephemeris from numerical propagation
    const numPoints = 20;
    const totalTime = orbitalPeriod;
    const dt = totalTime / (numPoints - 1);

    const ephemerisTable = [];

    console.log('Generating tabulated ephemeris (20 points over one orbit)...');

    for (let i = 0; i < numPoints; i++) {
        const time = i * dt;
        const meanAnomaly = meanMotion * time;

        // Solve Kepler's equation
        let E = meanAnomaly;
        for (let iter = 0; iter < 10; iter++) {
            E = meanAnomaly + eccentricity * Math.sin(E);
        }

        const trueAnomaly = 2 * Math.atan2(
            Math.sqrt(1 + eccentricity) * Math.sin(E / 2),
            Math.sqrt(1 - eccentricity) * Math.cos(E / 2)
        );

        // Compute state at this true anomaly
        kepler.set(5, trueAnomaly);
        const state = tudat.astro.element_conversion.keplerian_to_cartesian(kepler, earthGM);

        ephemerisTable.push({
            time: time,
            position: [state.get(0), state.get(1), state.get(2)],
            velocity: [state.get(3), state.get(4), state.get(5)]
        });

        state.delete();
    }

    console.log(`  Created ${ephemerisTable.length} ephemeris points`);
    console.log(`  Time span: 0 to ${(totalTime / 60).toFixed(2)} minutes`);
    console.log(`  Time step: ${(dt / 60).toFixed(2)} minutes\n`);

    // Interpolate at arbitrary times
    console.log('Interpolation (linear) at intermediate times:');
    console.log('Time (min)    X (km)         Y (km)         Z (km)');
    console.log('----------    ------         ------         ------');

    function interpolateEphemeris(table, time) {
        // Find bracketing indices
        let i = 0;
        while (i < table.length - 1 && table[i + 1].time < time) {
            i++;
        }

        if (i >= table.length - 1) {
            return table[table.length - 1];
        }

        // Linear interpolation factor
        const t0 = table[i].time;
        const t1 = table[i + 1].time;
        const f = (time - t0) / (t1 - t0);

        return {
            position: [
                table[i].position[0] + f * (table[i + 1].position[0] - table[i].position[0]),
                table[i].position[1] + f * (table[i + 1].position[1] - table[i].position[1]),
                table[i].position[2] + f * (table[i + 1].position[2] - table[i].position[2])
            ],
            velocity: [
                table[i].velocity[0] + f * (table[i + 1].velocity[0] - table[i].velocity[0]),
                table[i].velocity[1] + f * (table[i + 1].velocity[1] - table[i].velocity[1]),
                table[i].velocity[2] + f * (table[i + 1].velocity[2] - table[i].velocity[2])
            ]
        };
    }

    const queryTimes = [0, 11.5, 23.1, 45.7, 69.3, 92.0];
    for (const t of queryTimes) {
        const tSeconds = t * 60;
        const interp = interpolateEphemeris(ephemerisTable, tSeconds);
        console.log(`${t.toFixed(1).padStart(10)}    ${(interp.position[0]/1000).toFixed(2).padStart(10)}    ${(interp.position[1]/1000).toFixed(2).padStart(10)}    ${(interp.position[2]/1000).toFixed(2).padStart(10)}`);
    }

    // ========================================
    // Higher-Order Interpolation
    // ========================================
    console.log('\n\n4. Interpolation Methods');
    console.log('========================\n');

    console.log('Different interpolation methods for ephemeris data:\n');
    console.log('  Linear:       Simple, fast, discontinuous velocity');
    console.log('  Lagrange:     Polynomial interpolation, smooth');
    console.log('  Hermite:      Uses both position and velocity data');
    console.log('  Cubic spline: Smooth second derivatives\n');

    // Compare interpolation at midpoint
    const midTime = totalTime / 2;

    // True value (Keplerian)
    const meanAnomalyMid = meanMotion * midTime;
    let Emid = meanAnomalyMid;
    for (let iter = 0; iter < 10; iter++) {
        Emid = meanAnomalyMid + eccentricity * Math.sin(Emid);
    }
    const trueAnomalyMid = 2 * Math.atan2(
        Math.sqrt(1 + eccentricity) * Math.sin(Emid / 2),
        Math.sqrt(1 - eccentricity) * Math.cos(Emid / 2)
    );
    kepler.set(5, trueAnomalyMid);
    const trueState = tudat.astro.element_conversion.keplerian_to_cartesian(kepler, earthGM);

    // Linear interpolation
    const linearInterp = interpolateEphemeris(ephemerisTable, midTime);

    // Compute error
    const posError = Math.sqrt(
        Math.pow(linearInterp.position[0] - trueState.get(0), 2) +
        Math.pow(linearInterp.position[1] - trueState.get(1), 2) +
        Math.pow(linearInterp.position[2] - trueState.get(2), 2)
    );

    console.log(`Interpolation at t = ${(midTime / 60).toFixed(2)} minutes:`);
    console.log(`  True position:        [${(trueState.get(0)/1000).toFixed(3)}, ${(trueState.get(1)/1000).toFixed(3)}, ${(trueState.get(2)/1000).toFixed(3)}] km`);
    console.log(`  Linear interpolated:  [${(linearInterp.position[0]/1000).toFixed(3)}, ${(linearInterp.position[1]/1000).toFixed(3)}, ${(linearInterp.position[2]/1000).toFixed(3)}] km`);
    console.log(`  Position error:       ${(posError).toFixed(2)} m`);

    trueState.delete();

    // ========================================
    // Planetary Ephemerides
    // ========================================
    console.log('\n\n5. Planetary Ephemerides');
    console.log('========================\n');

    console.log('For solar system bodies, ephemerides typically come from:');
    console.log('  - SPICE kernels (JPL DE430, DE440)');
    console.log('  - Analytical theories (VSOP87, ELP2000)');
    console.log('  - Numerical integration\n');

    // Simplified planetary ephemeris (circular approximation)
    const planets = {
        Mercury: { sma: 57.91e9, period: 87.97 },
        Venus:   { sma: 108.2e9, period: 224.7 },
        Earth:   { sma: 149.6e9, period: 365.25 },
        Mars:    { sma: 227.9e9, period: 687.0 },
        Jupiter: { sma: 778.5e9, period: 4333 },
        Saturn:  { sma: 1432e9, period: 10759 }
    };

    console.log('Approximate heliocentric state at J2000.0:');
    console.log('Planet       X (AU)      Y (AU)      Z (AU)');
    console.log('------       ------      ------      ------');

    const AU = 149597870700;  // m
    const sunGM = 1.32712440018e20;

    // Approximate initial positions (based on J2000 epochs)
    const initialLongitudes = {
        Mercury: 252.3 * Math.PI / 180,
        Venus:   181.9 * Math.PI / 180,
        Earth:   100.5 * Math.PI / 180,
        Mars:    355.4 * Math.PI / 180,
        Jupiter: 34.4 * Math.PI / 180,
        Saturn:  49.9 * Math.PI / 180
    };

    for (const [name, params] of Object.entries(planets)) {
        const lon = initialLongitudes[name];
        const x = params.sma * Math.cos(lon) / AU;
        const y = params.sma * Math.sin(lon) / AU;
        console.log(`${name.padEnd(12)} ${x.toFixed(4).padStart(8)}    ${y.toFixed(4).padStart(8)}    ${(0).toFixed(4).padStart(8)}`);
    }

    // ========================================
    // Rotational Ephemeris
    // ========================================
    console.log('\n\n6. Rotational Ephemeris');
    console.log('=======================\n');

    console.log('Rotational ephemerides describe body orientation over time.\n');

    // Earth rotation
    const earthRotationRate = 7.292115e-5;  // rad/s (sidereal)
    const siderealDay = 2 * Math.PI / earthRotationRate;

    console.log('Earth rotation:');
    console.log(`  Sidereal rotation rate: ${(earthRotationRate * 86400 * 180 / Math.PI).toFixed(4)} deg/day`);
    console.log(`  Sidereal day: ${(siderealDay / 3600).toFixed(4)} hours`);
    console.log(`  Precession rate: ~50.3 arcsec/year`);
    console.log(`  Nutation: ±9.2 arcsec (18.6 year period)\n`);

    // Moon rotation (synchronous)
    const moonOrbitalPeriod = 27.32 * 86400;  // seconds
    const moonRotationRate = 2 * Math.PI / moonOrbitalPeriod;

    console.log('Moon rotation (synchronous):');
    console.log(`  Rotation period: ${(moonOrbitalPeriod / 86400).toFixed(2)} days`);
    console.log(`  Rotation rate: ${(moonRotationRate * 86400 * 180 / Math.PI).toFixed(4)} deg/day`);
    console.log(`  Libration amplitude: ~7° (longitude), ~6° (latitude)`);

    // Clean up
    kepler.delete();
    cartesian.delete();
    cartesian90.delete();

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
