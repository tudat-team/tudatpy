/**
 * Example 46: TLE-Based Ephemeris Setup
 *
 * Ported from: examples/tudatpy/estimation/set_ephemeris_from_tle.py
 *
 * This example demonstrates how to use Two-Line Element (TLE) data
 * to set up ephemeris for satellite tracking and orbit prediction.
 *
 * Key concepts:
 * - TLE format and parsing
 * - Mean elements vs osculating elements
 * - SGP4/SDP4 propagation
 * - Ephemeris generation from TLE
 * - TLE accuracy and limitations
 *
 * Run with: node 46_tle_ephemeris.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== TLE-Based Ephemeris Setup ===\n');

    const tudat = await createTudatModule();

    // Constants
    const GM_EARTH = 3.986004418e14;
    const R_EARTH = 6.378137e6;
    const J2 = 1.08263e-3;
    const DEG_TO_RAD = Math.PI / 180;
    const MINUTES_PER_DAY = 1440;

    /**
     * TLE Format Overview:
     *
     * Line 0: Satellite name (optional)
     * Line 1: 1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
     * Line 2: 2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN
     *
     * Line 1 fields:
     * - Column 1: Line number
     * - Columns 3-7: Satellite catalog number
     * - Column 8: Classification (U=unclassified)
     * - Columns 10-17: International designator (launch year, launch number, piece)
     * - Columns 19-32: Epoch (year and day fraction)
     * - Columns 34-43: First derivative of mean motion (ballistic coefficient)
     * - Columns 45-52: Second derivative of mean motion
     * - Columns 54-61: B* drag term
     * - Column 63: Ephemeris type
     * - Columns 65-68: Element set number
     * - Column 69: Checksum
     *
     * Line 2 fields:
     * - Column 1: Line number
     * - Columns 3-7: Satellite catalog number
     * - Columns 9-16: Inclination (degrees)
     * - Columns 18-25: RAAN (degrees)
     * - Columns 27-33: Eccentricity (decimal point assumed)
     * - Columns 35-42: Argument of perigee (degrees)
     * - Columns 44-51: Mean anomaly (degrees)
     * - Columns 53-63: Mean motion (revolutions per day)
     * - Columns 64-68: Revolution number at epoch
     */

    // Example TLEs for various satellites
    const tleData = [
        {
            name: 'ISS (ZARYA)',
            line1: '1 25544U 98067A   24015.50000000  .00016717  00000-0  30000-3 0  9993',
            line2: '2 25544  51.6400 247.4627 0006303 130.5360 325.0288 15.49815856434567'
        },
        {
            name: 'STARLINK-1007',
            line1: '1 44713U 19074A   24015.50000000  .00001264  00000-0  90000-4 0  9999',
            line2: '2 44713  53.0539 123.9846 0001234  89.3452 270.7789 15.06391987234567'
        },
        {
            name: 'HUBBLE SPACE TELESCOPE',
            line1: '1 20580U 90037B   24015.50000000  .00001580  00000-0  75000-4 0  9995',
            line2: '2 20580  28.4700  30.6879 0002695 123.4567 236.5432 15.09436987345678'
        }
    ];

    /**
     * Parse TLE lines and extract orbital elements
     */
    function parseTLE(line1, line2) {
        // Line 1 parsing
        const catalogNumber = parseInt(line1.substring(2, 7));
        const classification = line1.charAt(7);
        const intlDesignator = line1.substring(9, 17).trim();

        // Epoch
        const epochYear = parseInt(line1.substring(18, 20));
        const epochDay = parseFloat(line1.substring(20, 32));
        const fullYear = epochYear < 57 ? 2000 + epochYear : 1900 + epochYear;

        // Mean motion derivatives
        const ndot = parseFloat(line1.substring(33, 43));
        const nddot = parseFloat(line1.substring(44, 50).replace(' ', '').replace('-', 'e-').replace('+', 'e+') || '0');
        const bstar = parseFloat(line1.substring(53, 59).replace(' ', '') + 'e' + line1.substring(59, 61));

        // Line 2 parsing
        const inclination = parseFloat(line2.substring(8, 16)) * DEG_TO_RAD;
        const raan = parseFloat(line2.substring(17, 25)) * DEG_TO_RAD;
        const eccentricity = parseFloat('0.' + line2.substring(26, 33));
        const argPerigee = parseFloat(line2.substring(34, 42)) * DEG_TO_RAD;
        const meanAnomaly = parseFloat(line2.substring(43, 51)) * DEG_TO_RAD;
        const meanMotion = parseFloat(line2.substring(52, 63));  // rev/day
        const revNumber = parseInt(line2.substring(63, 68));

        // Compute semi-major axis from mean motion
        // n = sqrt(GM/a^3), n in rad/s
        const n_radPerSec = meanMotion * 2 * Math.PI / 86400;
        const a = Math.pow(GM_EARTH / (n_radPerSec * n_radPerSec), 1/3);

        // Compute period
        const period = 86400 / meanMotion;  // seconds

        return {
            catalogNumber,
            classification,
            intlDesignator,
            epoch: { year: fullYear, dayOfYear: epochDay },
            ndot,
            nddot,
            bstar,
            inclination,
            raan,
            eccentricity,
            argPerigee,
            meanAnomaly,
            meanMotion,
            revNumber,
            semiMajorAxis: a,
            period
        };
    }

    /**
     * Compute checksum for TLE line
     */
    function computeChecksum(line) {
        let sum = 0;
        for (let i = 0; i < 68; i++) {
            const char = line.charAt(i);
            if (char >= '0' && char <= '9') {
                sum += parseInt(char);
            } else if (char === '-') {
                sum += 1;
            }
        }
        return sum % 10;
    }

    // Parse and display TLE data
    console.log('Parsed TLE Data:');
    console.log('═'.repeat(70));

    const satellites = tleData.map(tle => {
        const elements = parseTLE(tle.line1, tle.line2);
        return { name: tle.name, ...elements };
    });

    satellites.forEach(sat => {
        console.log(`\n${sat.name}`);
        console.log('─'.repeat(50));
        console.log(`  Catalog #: ${sat.catalogNumber} (${sat.classification})`);
        console.log(`  Int'l ID:  ${sat.intlDesignator}`);
        console.log(`  Epoch:     ${sat.epoch.year} day ${sat.epoch.dayOfYear.toFixed(8)}`);
        console.log(`  Elements:`);
        console.log(`    a = ${(sat.semiMajorAxis / 1000).toFixed(3)} km`);
        console.log(`    e = ${sat.eccentricity.toFixed(7)}`);
        console.log(`    i = ${(sat.inclination / DEG_TO_RAD).toFixed(4)}°`);
        console.log(`    Ω = ${(sat.raan / DEG_TO_RAD).toFixed(4)}°`);
        console.log(`    ω = ${(sat.argPerigee / DEG_TO_RAD).toFixed(4)}°`);
        console.log(`    M = ${(sat.meanAnomaly / DEG_TO_RAD).toFixed(4)}°`);
        console.log(`  Mean motion: ${sat.meanMotion.toFixed(8)} rev/day`);
        console.log(`  Period: ${(sat.period / 60).toFixed(2)} min`);
        console.log(`  B*: ${sat.bstar.toExponential(4)}`);

        // Derived quantities
        const altitude = sat.semiMajorAxis - R_EARTH;
        const perigee = sat.semiMajorAxis * (1 - sat.eccentricity) - R_EARTH;
        const apogee = sat.semiMajorAxis * (1 + sat.eccentricity) - R_EARTH;

        console.log(`  Altitude: ${(altitude / 1000).toFixed(1)} km (mean)`);
        console.log(`  Perigee:  ${(perigee / 1000).toFixed(1)} km`);
        console.log(`  Apogee:   ${(apogee / 1000).toFixed(1)} km`);
    });

    // SGP4 Propagation (simplified)
    console.log('\n\n--- Simplified SGP4 Propagation ---\n');

    /**
     * Simplified SGP4-like propagation
     * (Full SGP4 is more complex - this demonstrates the concept)
     */
    function propagateSGP4(elements, dtMinutes) {
        // Time in minutes since epoch
        const dt = dtMinutes * 60;  // seconds

        // Mean motion in rad/s
        const n0 = elements.meanMotion * 2 * Math.PI / 86400;

        // J2 secular perturbations
        const a0 = elements.semiMajorAxis;
        const e0 = elements.eccentricity;
        const i0 = elements.inclination;
        const p = a0 * (1 - e0 * e0);
        const cosI = Math.cos(i0);
        const sinI = Math.sin(i0);

        // Secular rates
        const RAANdot = -1.5 * n0 * J2 * (R_EARTH / p) ** 2 * cosI;
        const argPdot = 0.75 * n0 * J2 * (R_EARTH / p) ** 2 * (4 - 5 * sinI * sinI);
        const Mdot = n0 + 0.75 * n0 * J2 * (R_EARTH / p) ** 2 * Math.sqrt(1 - e0*e0) * (2 - 3 * sinI * sinI);

        // Updated elements
        const RAAN = elements.raan + RAANdot * dt;
        const argP = elements.argPerigee + argPdot * dt;
        const M = elements.meanAnomaly + Mdot * dt;

        // Solve Kepler's equation
        let E = M;
        for (let iter = 0; iter < 10; iter++) {
            E = E - (E - e0 * Math.sin(E) - M) / (1 - e0 * Math.cos(E));
        }

        // True anomaly
        const nu = 2 * Math.atan2(
            Math.sqrt(1 + e0) * Math.sin(E / 2),
            Math.sqrt(1 - e0) * Math.cos(E / 2)
        );

        // Position in orbital plane
        const r = a0 * (1 - e0 * Math.cos(E));
        const xOrb = r * Math.cos(nu);
        const yOrb = r * Math.sin(nu);

        // Rotation to ECI
        const cosO = Math.cos(RAAN), sinO = Math.sin(RAAN);
        const cosw = Math.cos(argP), sinw = Math.sin(argP);
        const cosi = Math.cos(i0), sini = Math.sin(i0);

        const x = (cosO*cosw - sinO*sinw*cosi) * xOrb + (-cosO*sinw - sinO*cosw*cosi) * yOrb;
        const y = (sinO*cosw + cosO*sinw*cosi) * xOrb + (-sinO*sinw + cosO*cosw*cosi) * yOrb;
        const z = (sinw*sini) * xOrb + (cosw*sini) * yOrb;

        return { x, y, z, r };
    }

    // Generate ephemeris for ISS
    const iss = satellites[0];
    console.log(`Generating ephemeris for ${iss.name}:`);
    console.log('Time [min] | X [km]      | Y [km]      | Z [km]      | R [km]');
    console.log('───────────┼─────────────┼─────────────┼─────────────┼───────────');

    for (let t = 0; t <= 90; t += 15) {
        const pos = propagateSGP4(iss, t);
        console.log(`${t.toString().padStart(10)} | ${(pos.x/1000).toFixed(1).padStart(11)} | ${(pos.y/1000).toFixed(1).padStart(11)} | ${(pos.z/1000).toFixed(1).padStart(11)} | ${(pos.r/1000).toFixed(1).padStart(9)}`);
    }

    // TLE accuracy analysis
    console.log('\n--- TLE Accuracy Considerations ---\n');

    console.log('Typical TLE accuracy:');
    console.log('  - Position: 1-5 km at epoch');
    console.log('  - Along-track error: ~1 km/day growth');
    console.log('  - Cross-track: relatively stable');
    console.log('  - Radial: ~0.1-0.5 km');

    console.log('\nFactors affecting TLE accuracy:');
    console.log('  1. Age of TLE (error grows with time)');
    console.log('  2. Atmospheric drag (LEO satellites)');
    console.log('  3. Solar activity (affects drag)');
    console.log('  4. Orbital altitude and eccentricity');
    console.log('  5. Maneuver history');

    // Best practices
    console.log('\n--- Best Practices for TLE Usage ---\n');

    console.log('1. Use recent TLEs (< 1-2 days old for LEO)');
    console.log('2. Propagate forward only, not backward');
    console.log('3. Use SGP4/SDP4 propagator (not Keplerian)');
    console.log('4. Consider TLE update frequency:');
    console.log('   - ISS: Multiple times per day');
    console.log('   - LEO: Daily');
    console.log('   - MEO/GEO: Weekly');
    console.log('5. Account for conjunction analysis margins');

    // Comparison: Mean vs Osculating elements
    console.log('\n--- Mean vs Osculating Elements ---\n');

    console.log('TLE elements are MEAN elements (averaged over perturbations)');
    console.log('For precise work, convert to osculating elements:\n');

    const oscCorrections = {
        'Semi-major axis': '~1-10 km difference',
        'Eccentricity': '~0.001-0.01 difference',
        'Inclination': '~0.01-0.1° difference',
        'RAAN': '~0.1-1° difference',
        'Arg. perigee': '~0.1-1° difference',
        'Mean anomaly': '~0.1-1° difference'
    };

    for (const [element, correction] of Object.entries(oscCorrections)) {
        console.log(`  ${element.padEnd(16)}: ${correction}`);
    }

    console.log('\n=== TLE Ephemeris Setup Complete ===');
}

main().catch(console.error);
