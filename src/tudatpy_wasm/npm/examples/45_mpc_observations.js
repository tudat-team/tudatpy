/**
 * Example 45: Minor Planet Center (MPC) Observation Handling
 *
 * Ported from: examples/tudatpy/estimation/estimation_with_mpc.py
 *             examples/tudatpy/estimation/retrieving_mpc_observation_data.py
 *
 * This example demonstrates how to work with Minor Planet Center
 * observation data format, which is the standard for asteroid and
 * comet astrometry.
 *
 * Key concepts:
 * - MPC observation format (80-column)
 * - Right ascension and declination
 * - Observatory codes
 * - Observation weighting
 * - Astrometric reduction
 *
 * Run with: node 45_mpc_observations.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Minor Planet Center (MPC) Observation Handling ===\n');

    const tudat = await createTudatModule();

    // Constants
    const DEG_TO_RAD = Math.PI / 180;
    const RAD_TO_DEG = 180 / Math.PI;
    const ARCSEC_TO_RAD = Math.PI / 648000;
    const AU = 1.496e11;

    /**
     * MPC 80-column format specification:
     * Columns 1-5: Packed provisional/permanent designation
     * Column 6: Discovery asterisk (if applicable)
     * Column 7: Note 1
     * Column 8: Note 2
     * Columns 15-25: Date of observation (YYYY MM DD.ddddd)
     * Columns 33-44: Right ascension (HH MM SS.ss)
     * Columns 45-56: Declination (sDD MM SS.s)
     * Columns 66-71: Magnitude
     * Column 72: Band
     * Columns 78-80: Observatory code
     */

    // Example MPC observations (simulated for asteroid 433 Eros)
    const mpcObservations = [
        '00433         C2024 01 15.54321 05 23 45.67 +15 34 56.7          16.2 V      703',
        '00433         C2024 01 16.55432 05 24 12.34 +15 36 12.3          16.1 V      703',
        '00433         C2024 01 17.56543 05 24 39.01 +15 37 28.9          16.0 V      703',
        '00433         C2024 01 18.57654 05 25 05.68 +15 38 45.5          15.9 V      F51',
        '00433         C2024 01 19.58765 05 25 32.35 +15 40 02.1          15.8 V      F51',
        '00433         C2024 01 20.59876 05 25 59.02 +15 41 18.7          15.7 V      G96'
    ];

    console.log('MPC 80-Column Format Overview:');
    console.log('─'.repeat(80));
    console.log('Cols  1-5:  Object designation');
    console.log('Cols 15-25: Date (YYYY MM DD.ddddd)');
    console.log('Cols 33-44: Right Ascension (HH MM SS.ss)');
    console.log('Cols 45-56: Declination (±DD MM SS.s)');
    console.log('Cols 66-71: Magnitude');
    console.log('Col  72:    Photometric band');
    console.log('Cols 78-80: Observatory code');
    console.log('─'.repeat(80));

    /**
     * Parse MPC observation line
     */
    function parseMPCLine(line) {
        // Ensure line is at least 80 characters
        const paddedLine = line.padEnd(80);

        // Extract fields
        const designation = paddedLine.substring(0, 5).trim();
        const discovery = paddedLine.charAt(5) === '*';
        const note1 = paddedLine.charAt(6);
        const note2 = paddedLine.charAt(7);

        // Date parsing
        const year = parseInt(paddedLine.substring(14, 18));
        const month = parseInt(paddedLine.substring(19, 21));
        const dayFrac = parseFloat(paddedLine.substring(22, 31));
        const day = Math.floor(dayFrac);
        const hourFrac = (dayFrac - day) * 24;
        const hour = Math.floor(hourFrac);
        const minFrac = (hourFrac - hour) * 60;
        const minute = Math.floor(minFrac);
        const second = (minFrac - minute) * 60;

        // Right ascension parsing (HH MM SS.ss)
        const raHours = parseInt(paddedLine.substring(32, 34));
        const raMinutes = parseInt(paddedLine.substring(35, 37));
        const raSeconds = parseFloat(paddedLine.substring(38, 43));
        const ra = (raHours + raMinutes / 60 + raSeconds / 3600) * 15 * DEG_TO_RAD;  // Convert to radians

        // Declination parsing (±DD MM SS.s)
        const decSign = paddedLine.charAt(44) === '-' ? -1 : 1;
        const decDegrees = parseInt(paddedLine.substring(45, 47));
        const decMinutes = parseInt(paddedLine.substring(48, 50));
        const decSeconds = parseFloat(paddedLine.substring(51, 55));
        const dec = decSign * (decDegrees + decMinutes / 60 + decSeconds / 3600) * DEG_TO_RAD;

        // Magnitude
        const magStr = paddedLine.substring(65, 70).trim();
        const magnitude = magStr ? parseFloat(magStr) : null;
        const band = paddedLine.charAt(71).trim() || 'V';

        // Observatory code
        const obsCode = paddedLine.substring(77, 80).trim();

        return {
            designation,
            discovery,
            note1,
            note2,
            date: { year, month, day, hour, minute, second },
            ra,  // radians
            dec, // radians
            magnitude,
            band,
            obsCode
        };
    }

    /**
     * Format RA for display
     */
    function formatRA(ra) {
        const raDeg = ra * RAD_TO_DEG;
        const raHours = raDeg / 15;
        const h = Math.floor(raHours);
        const m = Math.floor((raHours - h) * 60);
        const s = ((raHours - h) * 60 - m) * 60;
        return `${h.toString().padStart(2, '0')}h ${m.toString().padStart(2, '0')}m ${s.toFixed(2).padStart(5, '0')}s`;
    }

    /**
     * Format Dec for display
     */
    function formatDec(dec) {
        const decDeg = dec * RAD_TO_DEG;
        const sign = decDeg >= 0 ? '+' : '-';
        const absDec = Math.abs(decDeg);
        const d = Math.floor(absDec);
        const m = Math.floor((absDec - d) * 60);
        const s = ((absDec - d) * 60 - m) * 60;
        return `${sign}${d.toString().padStart(2, '0')}° ${m.toString().padStart(2, '0')}' ${s.toFixed(1).padStart(4, '0')}"`;
    }

    // Parse all observations
    console.log('\nParsed Observations:');
    console.log('═'.repeat(80));

    const observations = mpcObservations.map(line => parseMPCLine(line));

    observations.forEach((obs, i) => {
        const dateStr = `${obs.date.year}-${obs.date.month.toString().padStart(2, '0')}-${obs.date.day.toString().padStart(2, '0')}`;
        console.log(`Obs ${i + 1}: ${obs.designation}`);
        console.log(`  Date: ${dateStr} ${obs.date.hour.toString().padStart(2, '0')}:${obs.date.minute.toString().padStart(2, '0')}:${obs.date.second.toFixed(1).padStart(4, '0')}`);
        console.log(`  RA:   ${formatRA(obs.ra)} (${(obs.ra * RAD_TO_DEG).toFixed(4)}°)`);
        console.log(`  Dec:  ${formatDec(obs.dec)} (${(obs.dec * RAD_TO_DEG).toFixed(4)}°)`);
        console.log(`  Mag:  ${obs.magnitude} ${obs.band}`);
        console.log(`  Obs:  ${obs.obsCode}`);
    });

    // Observatory codes
    console.log('\n--- Observatory Codes ---\n');

    const observatories = {
        '703': { name: 'Catalina Sky Survey', lat: 32.4167, lon: -110.7317, alt: 2510 },
        'F51': { name: 'Pan-STARRS 1', lat: 20.7072, lon: -156.2575, alt: 3052 },
        'G96': { name: 'Mt. Lemmon Survey', lat: 32.4433, lon: -110.7883, alt: 2776 },
        '950': { name: 'La Palma', lat: 28.7606, lon: -17.8819, alt: 2326 },
        '568': { name: 'Mauna Kea', lat: 19.8260, lon: -155.4720, alt: 4205 }
    };

    console.log('Observatories used in this example:');
    console.log('Code | Name                  | Latitude  | Longitude  | Altitude');
    console.log('─────┼───────────────────────┼───────────┼────────────┼──────────');

    const usedObs = new Set(observations.map(o => o.obsCode));
    for (const code of usedObs) {
        if (observatories[code]) {
            const obs = observatories[code];
            console.log(`${code}  | ${obs.name.padEnd(21)} | ${obs.lat.toFixed(4).padStart(9)}° | ${obs.lon.toFixed(4).padStart(10)}° | ${obs.alt.toString().padStart(6)} m`);
        }
    }

    // Observation quality and weighting
    console.log('\n--- Observation Weighting ---\n');

    // Typical astrometric uncertainties by observatory class
    const uncertainties = {
        'professional': 0.2,  // arcsec
        'amateur': 0.5,       // arcsec
        'satellite': 1.0,     // arcsec
        'historical': 2.0     // arcsec
    };

    console.log('Typical astrometric uncertainties:');
    for (const [type, sigma] of Object.entries(uncertainties)) {
        console.log(`  ${type.padEnd(12)}: ${sigma.toFixed(1)} arcsec (${(sigma * ARCSEC_TO_RAD * 1e6).toFixed(1)} µrad)`);
    }

    // Weight calculation
    const refSigma = 0.5;  // Reference uncertainty [arcsec]
    console.log(`\nWeight calculation (σ_ref = ${refSigma} arcsec):`);
    console.log('Weight = (σ_ref / σ_obs)²\n');

    observations.forEach((obs, i) => {
        // Assume professional observatory
        const sigma = 0.3;  // arcsec
        const weight = Math.pow(refSigma / sigma, 2);
        console.log(`  Obs ${i + 1} (${obs.obsCode}): σ = ${sigma.toFixed(1)}" -> w = ${weight.toFixed(2)}`);
    });

    // Motion analysis
    console.log('\n--- Apparent Motion Analysis ---\n');

    // Compute angular motion between observations
    for (let i = 1; i < observations.length; i++) {
        const obs1 = observations[i - 1];
        const obs2 = observations[i];

        // Time difference in hours
        const dt = 24;  // Approximately 1 day between observations

        // Angular separation
        const dRA = obs2.ra - obs1.ra;
        const dDec = obs2.dec - obs1.dec;
        const cosDec = Math.cos((obs1.dec + obs2.dec) / 2);
        const sep = Math.sqrt((dRA * cosDec) ** 2 + dDec ** 2);

        // Motion rate
        const rateArcsecPerHour = sep * RAD_TO_DEG * 3600 / dt;

        console.log(`Obs ${i} -> ${i + 1}:`);
        console.log(`  ΔRA: ${(dRA * RAD_TO_DEG * 3600).toFixed(1)}" ΔDec: ${(dDec * RAD_TO_DEG * 3600).toFixed(1)}"`);
        console.log(`  Separation: ${(sep * RAD_TO_DEG * 3600).toFixed(1)}" Rate: ${rateArcsecPerHour.toFixed(2)}"/hr`);
    }

    // Residual computation (simulated)
    console.log('\n--- Simulated Residuals ---\n');

    // Assume we have computed positions from orbit determination
    // These would come from numerical integration of the orbit
    const computedPositions = observations.map((obs, i) => ({
        ra: obs.ra + (Math.random() - 0.5) * 0.5 * ARCSEC_TO_RAD,  // Add small offset
        dec: obs.dec + (Math.random() - 0.5) * 0.5 * ARCSEC_TO_RAD
    }));

    console.log('Observation Residuals (O-C):');
    console.log('Obs | ΔRA ["]  | ΔDec ["] | RSS ["]');
    console.log('────┼──────────┼──────────┼─────────');

    let sumSqRA = 0, sumSqDec = 0;
    observations.forEach((obs, i) => {
        const comp = computedPositions[i];
        const cosDec = Math.cos(obs.dec);

        const dRA = (obs.ra - comp.ra) * cosDec * RAD_TO_DEG * 3600;
        const dDec = (obs.dec - comp.dec) * RAD_TO_DEG * 3600;
        const rss = Math.sqrt(dRA * dRA + dDec * dDec);

        sumSqRA += dRA * dRA;
        sumSqDec += dDec * dDec;

        console.log(`${(i + 1).toString().padStart(3)} | ${dRA.toFixed(3).padStart(8)} | ${dDec.toFixed(3).padStart(8)} | ${rss.toFixed(3).padStart(7)}`);
    });

    const rmsRA = Math.sqrt(sumSqRA / observations.length);
    const rmsDec = Math.sqrt(sumSqDec / observations.length);

    console.log('────┴──────────┴──────────┴─────────');
    console.log(`RMS: ${rmsRA.toFixed(3).padStart(12)} | ${rmsDec.toFixed(3).padStart(8)}`);

    // MPC submission format
    console.log('\n--- MPC Submission Notes ---\n');
    console.log('For submitting observations to MPC:');
    console.log('  1. Use standard 80-column format');
    console.log('  2. Include all required fields');
    console.log('  3. Use correct observatory code');
    console.log('  4. Report magnitude in appropriate band');
    console.log('  5. Include discovery asterisk if applicable');
    console.log('  6. Submit via WAMO (Web-based Astrometric Measurement Output)');

    console.log('\n=== MPC Observation Handling Complete ===');
}

main().catch(console.error);
