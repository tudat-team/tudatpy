/**
 * Example 16: Two-Line Element (TLE) Propagation
 *
 * This example demonstrates TLE parsing, SGP4/SDP4 propagation,
 * and satellite tracking using TLE data.
 *
 * Run with: node 16_tle_propagation.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== TLE Propagation Example ===\n');

    const tudat = await createTudatModule();

    // ========================================
    // TLE Format Overview
    // ========================================
    console.log('1. TLE Format Overview');
    console.log('======================\n');

    console.log('A Two-Line Element (TLE) set contains orbital parameters');
    console.log('in a standardized format for satellite tracking.\n');

    // Example TLE (ISS - typical format)
    const tleLine0 = 'ISS (ZARYA)';
    const tleLine1 = '1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9993';
    const tleLine2 = '2 25544  51.6400 247.4627 0007417  80.6015 279.5714 15.50055555123456';

    console.log('Example TLE (ISS):');
    console.log(`  Line 0: ${tleLine0}`);
    console.log(`  Line 1: ${tleLine1}`);
    console.log(`  Line 2: ${tleLine2}\n`);

    // ========================================
    // TLE Parsing
    // ========================================
    console.log('\n2. TLE Element Breakdown');
    console.log('========================\n');

    console.log('Line 1 fields:');
    console.log('  01-01: Line number (1)');
    console.log('  03-07: Satellite catalog number (25544)');
    console.log('  08-08: Classification (U = unclassified)');
    console.log('  10-11: International designator year (98)');
    console.log('  12-14: International designator launch (067)');
    console.log('  15-17: International designator piece (A)');
    console.log('  19-20: Epoch year (24)');
    console.log('  21-32: Epoch day fraction (001.50000000)');
    console.log('  34-43: Mean motion derivative (d²n/dt² / 2)');
    console.log('  45-52: Mean motion 2nd derivative (d³n/dt³ / 6)');
    console.log('  54-61: BSTAR drag term');
    console.log('  63-63: Ephemeris type (0)');
    console.log('  65-68: Element set number');
    console.log('  69-69: Checksum\n');

    console.log('Line 2 fields:');
    console.log('  01-01: Line number (2)');
    console.log('  03-07: Satellite catalog number');
    console.log('  09-16: Inclination (degrees)');
    console.log('  18-25: RAAN (degrees)');
    console.log('  27-33: Eccentricity (decimal assumed)');
    console.log('  35-42: Argument of perigee (degrees)');
    console.log('  44-51: Mean anomaly (degrees)');
    console.log('  53-63: Mean motion (revs/day)');
    console.log('  64-68: Revolution number at epoch');
    console.log('  69-69: Checksum\n');

    // Parse TLE elements
    function parseTLE(line1, line2) {
        return {
            catalogNumber: parseInt(line1.substring(2, 7)),
            classification: line1.charAt(7),
            intlDesig: line1.substring(9, 17).trim(),
            epochYear: parseInt(line1.substring(18, 20)),
            epochDay: parseFloat(line1.substring(20, 32)),
            meanMotionDot: parseFloat(line1.substring(33, 43)),
            meanMotionDDot: parseFloat((line1.substring(44, 50) + 'e' + line1.substring(50, 52)).replace(' ', '')),
            bstar: parseFloat((line1.substring(53, 59) + 'e' + line1.substring(59, 61)).replace(' ', '')),
            ephemerisType: parseInt(line1.charAt(62)),
            elementSetNum: parseInt(line1.substring(64, 68)),

            inclination: parseFloat(line2.substring(8, 16)),
            raan: parseFloat(line2.substring(17, 25)),
            eccentricity: parseFloat('0.' + line2.substring(26, 33)),
            argPerigee: parseFloat(line2.substring(34, 42)),
            meanAnomaly: parseFloat(line2.substring(43, 51)),
            meanMotion: parseFloat(line2.substring(52, 63)),
            revNumber: parseInt(line2.substring(63, 68))
        };
    }

    const tle = parseTLE(tleLine1, tleLine2);

    console.log('Parsed TLE elements:');
    console.log(`  Satellite: ${tleLine0} (NORAD ID: ${tle.catalogNumber})`);
    console.log(`  Epoch: 20${tle.epochYear} day ${tle.epochDay.toFixed(8)}`);
    console.log(`  Inclination: ${tle.inclination.toFixed(4)}°`);
    console.log(`  RAAN: ${tle.raan.toFixed(4)}°`);
    console.log(`  Eccentricity: ${tle.eccentricity.toFixed(7)}`);
    console.log(`  Arg of perigee: ${tle.argPerigee.toFixed(4)}°`);
    console.log(`  Mean anomaly: ${tle.meanAnomaly.toFixed(4)}°`);
    console.log(`  Mean motion: ${tle.meanMotion.toFixed(8)} rev/day`);
    console.log(`  BSTAR: ${tle.bstar.toExponential(4)}`);

    // ========================================
    // Orbital Parameters from TLE
    // ========================================
    console.log('\n\n3. Derived Orbital Parameters');
    console.log('==============================\n');

    const earthGM = 3.986004418e14;  // m³/s²
    const earthRadius = 6378.137;    // km

    // Mean motion to semi-major axis
    const meanMotionRadSec = tle.meanMotion * 2 * Math.PI / 86400;
    const semiMajorAxis = Math.pow(earthGM / (meanMotionRadSec * meanMotionRadSec), 1/3);

    // Orbital period
    const orbitalPeriod = 86400 / tle.meanMotion;

    // Apogee and perigee
    const apogeeRadius = semiMajorAxis * (1 + tle.eccentricity);
    const perigeeRadius = semiMajorAxis * (1 - tle.eccentricity);
    const apogeeAltitude = apogeeRadius / 1000 - earthRadius;
    const perigeeAltitude = perigeeRadius / 1000 - earthRadius;

    console.log('Derived parameters:');
    console.log(`  Semi-major axis: ${(semiMajorAxis / 1000).toFixed(3)} km`);
    console.log(`  Orbital period: ${(orbitalPeriod / 60).toFixed(2)} minutes`);
    console.log(`  Apogee altitude: ${apogeeAltitude.toFixed(2)} km`);
    console.log(`  Perigee altitude: ${perigeeAltitude.toFixed(2)} km`);
    console.log(`  Revolutions per day: ${tle.meanMotion.toFixed(4)}`);

    // Velocity at perigee and apogee
    const velocityPerigee = Math.sqrt(earthGM * (2/perigeeRadius - 1/semiMajorAxis));
    const velocityApogee = Math.sqrt(earthGM * (2/apogeeRadius - 1/semiMajorAxis));

    console.log(`  Velocity at perigee: ${(velocityPerigee / 1000).toFixed(3)} km/s`);
    console.log(`  Velocity at apogee: ${(velocityApogee / 1000).toFixed(3)} km/s`);

    // ========================================
    // SGP4 Propagation Model
    // ========================================
    console.log('\n\n4. SGP4/SDP4 Propagation');
    console.log('========================\n');

    console.log('The SGP4 (Simplified General Perturbations 4) model includes:');
    console.log('  - Earth oblateness (J2, J3, J4)');
    console.log('  - Atmospheric drag (using BSTAR)');
    console.log('  - Solar/lunar gravitational effects (in SDP4 for deep space)');
    console.log('  - Secular and periodic variations\n');

    console.log('SGP4 is used for near-Earth satellites (period < 225 minutes)');
    console.log('SDP4 is used for deep-space satellites (period > 225 minutes)\n');

    // Simplified SGP4-like propagation (not full implementation)
    // This demonstrates the concept without full SGP4 complexity

    console.log('Propagating ISS orbit for 24 hours (simplified model):\n');

    // Convert TLE elements to initial Keplerian elements
    const keplerElements = new tudat.Vector6d();
    keplerElements.set(0, semiMajorAxis);
    keplerElements.set(1, tle.eccentricity);
    keplerElements.set(2, tle.inclination * Math.PI / 180);
    keplerElements.set(3, tle.argPerigee * Math.PI / 180);
    keplerElements.set(4, tle.raan * Math.PI / 180);
    keplerElements.set(5, tle.meanAnomaly * Math.PI / 180);

    // Convert mean anomaly to true anomaly (simplified)
    let E = keplerElements.get(5);
    for (let i = 0; i < 10; i++) {
        E = keplerElements.get(5) + tle.eccentricity * Math.sin(E);
    }
    const trueAnomaly = 2 * Math.atan2(
        Math.sqrt(1 + tle.eccentricity) * Math.sin(E / 2),
        Math.sqrt(1 - tle.eccentricity) * Math.cos(E / 2)
    );
    keplerElements.set(5, trueAnomaly);

    // Get initial Cartesian state
    const initialState = tudat.astro.element_conversion.keplerian_to_cartesian(keplerElements, earthGM);

    console.log('Time (hr)    Altitude (km)    Lat (deg)    Lon (deg)');
    console.log('---------    -------------    ---------    ---------');

    // Propagate using two-body with simple J2 approximation
    const J2 = 1.08263e-3;
    const propagationDuration = 24 * 3600;  // 24 hours
    const numSteps = 96;  // Every 15 minutes
    const dt = propagationDuration / numSteps;

    let position = [initialState.get(0), initialState.get(1), initialState.get(2)];
    let velocity = [initialState.get(3), initialState.get(4), initialState.get(5)];

    for (let step = 0; step <= numSteps; step += 4) {  // Print every hour
        const time = step * dt;

        // Current position characteristics
        const r = Math.sqrt(position[0]**2 + position[1]**2 + position[2]**2);
        const altitude = r / 1000 - earthRadius;

        // Latitude (geocentric)
        const lat = Math.asin(position[2] / r) * 180 / Math.PI;

        // Longitude (simplified, ignoring Earth rotation for display)
        let lon = Math.atan2(position[1], position[0]) * 180 / Math.PI;

        // Adjust for Earth rotation
        const earthRotationRate = 360 / 86164.1;  // deg/sec (sidereal)
        lon -= earthRotationRate * time;
        lon = ((lon + 180) % 360 + 360) % 360 - 180;

        if (step <= numSteps) {
            console.log(`${(time / 3600).toFixed(1).padStart(9)}    ${altitude.toFixed(1).padStart(13)}    ${lat.toFixed(1).padStart(9)}    ${lon.toFixed(1).padStart(9)}`);
        }

        // Simple propagation step (Euler for demonstration)
        if (step < numSteps) {
            for (let i = 0; i < 4; i++) {  // 4 substeps per output
                const substep_r = Math.sqrt(position[0]**2 + position[1]**2 + position[2]**2);

                // Two-body acceleration
                const accel = [
                    -earthGM * position[0] / Math.pow(substep_r, 3),
                    -earthGM * position[1] / Math.pow(substep_r, 3),
                    -earthGM * position[2] / Math.pow(substep_r, 3)
                ];

                // Simplified J2 perturbation
                const z_r = position[2] / substep_r;
                const J2_factor = 1.5 * J2 * Math.pow(earthRadius * 1000 / substep_r, 2);
                accel[0] += earthGM * position[0] / Math.pow(substep_r, 3) * J2_factor * (1 - 5 * z_r * z_r);
                accel[1] += earthGM * position[1] / Math.pow(substep_r, 3) * J2_factor * (1 - 5 * z_r * z_r);
                accel[2] += earthGM * position[2] / Math.pow(substep_r, 3) * J2_factor * (3 - 5 * z_r * z_r);

                // Update state
                position[0] += velocity[0] * dt / 4;
                position[1] += velocity[1] * dt / 4;
                position[2] += velocity[2] * dt / 4;
                velocity[0] += accel[0] * dt / 4;
                velocity[1] += accel[1] * dt / 4;
                velocity[2] += accel[2] * dt / 4;
            }
        }
    }

    // ========================================
    // Ground Track
    // ========================================
    console.log('\n\n5. Ground Track Characteristics');
    console.log('================================\n');

    // Ground track repeat
    const revsPerDay = tle.meanMotion;
    const groundTrackRepeat = 86400 / orbitalPeriod - Math.floor(86400 / orbitalPeriod);

    console.log(`ISS ground track characteristics:`);
    console.log(`  Revolutions per day: ${revsPerDay.toFixed(2)}`);
    console.log(`  Westward shift per orbit: ${(360 / revsPerDay).toFixed(2)}°`);
    console.log(`  Latitude range: ±${tle.inclination.toFixed(1)}°`);

    // Time between ascending node crossings
    const nodalPeriod = 86400 / revsPerDay;
    console.log(`  Nodal period: ${(nodalPeriod / 60).toFixed(2)} minutes`);

    // ========================================
    // Multiple Satellites
    // ========================================
    console.log('\n\n6. Satellite Constellation Example');
    console.log('===================================\n');

    // Example TLEs for different satellite types
    const satellites = [
        {
            name: 'ISS',
            meanMotion: 15.50,
            inclination: 51.6,
            altitude: 420
        },
        {
            name: 'Hubble Space Telescope',
            meanMotion: 15.09,
            inclination: 28.5,
            altitude: 540
        },
        {
            name: 'Starlink (typical)',
            meanMotion: 15.06,
            inclination: 53.0,
            altitude: 550
        },
        {
            name: 'GPS (typical)',
            meanMotion: 2.00,
            inclination: 55.0,
            altitude: 20200
        },
        {
            name: 'GEO (geostationary)',
            meanMotion: 1.00,
            inclination: 0.0,
            altitude: 35786
        }
    ];

    console.log('Satellite              Rev/day    Inc (°)    Alt (km)    Period (hr)');
    console.log('=========              =======    =======    ========    ===========');

    for (const sat of satellites) {
        const period = 24 / sat.meanMotion;
        console.log(`${sat.name.padEnd(22)} ${sat.meanMotion.toFixed(2).padStart(7)}    ${sat.inclination.toFixed(1).padStart(7)}    ${sat.altitude.toString().padStart(8)}    ${period.toFixed(2).padStart(11)}`);
    }

    // ========================================
    // TLE Accuracy and Updates
    // ========================================
    console.log('\n\n7. TLE Accuracy Considerations');
    console.log('==============================\n');

    console.log('TLE accuracy degrades over time:');
    console.log('  - Fresh TLE (< 1 day): ~1-5 km accuracy');
    console.log('  - 1 week old: ~10-50 km accuracy');
    console.log('  - 1 month old: Could be 100+ km off\n');

    console.log('Factors affecting TLE accuracy:');
    console.log('  - Atmospheric drag variability');
    console.log('  - Solar activity (affects drag)');
    console.log('  - Maneuvers (TLE becomes invalid)');
    console.log('  - Orbital altitude (lower = more drag)\n');

    console.log('Recommendations:');
    console.log('  - Use TLEs within ±7 days of epoch');
    console.log('  - Update TLEs regularly from Space-Track.org');
    console.log('  - For precision work, use SP3 ephemerides');

    // Clean up
    keplerElements.delete();
    initialState.delete();

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
