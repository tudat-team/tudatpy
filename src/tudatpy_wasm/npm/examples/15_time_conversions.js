/**
 * Example 15: Time Conversions and Earth Orientation
 *
 * This example demonstrates time scale conversions, Julian dates,
 * and Earth orientation parameters.
 *
 * Run with: node 15_time_conversions.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Time Conversions and Earth Orientation ===\n');

    const tudat = await createTudatModule();

    // ========================================
    // Time Scales
    // ========================================
    console.log('1. Astronomical Time Scales');
    console.log('===========================\n');

    console.log('Common time scales used in astrodynamics:\n');
    console.log('  UTC  - Coordinated Universal Time (civil time)');
    console.log('  TAI  - International Atomic Time (continuous)');
    console.log('  TT   - Terrestrial Time (for Earth surface)');
    console.log('  TDB  - Barycentric Dynamical Time (solar system dynamics)');
    console.log('  UT1  - Universal Time (Earth rotation angle)');
    console.log('  GPS  - GPS Time (GPS system epoch)\n');

    // Time scale offsets
    console.log('Time scale relationships (as of 2024):');
    console.log('  TAI - UTC = 37 seconds (leap seconds)');
    console.log('  TT - TAI  = 32.184 seconds (fixed)');
    console.log('  TT - UTC  = 69.184 seconds');
    console.log('  GPS - UTC = 18 seconds');
    console.log('  UT1 - UTC ≈ 0 ± 0.9 seconds (varies)\n');

    const leapSeconds = 37;  // Current (2024)
    const ttMinusTai = 32.184;

    // ========================================
    // Julian Date Conversions
    // ========================================
    console.log('\n2. Julian Date Conversions');
    console.log('==========================\n');

    console.log('Julian Date (JD): Continuous count of days since');
    console.log('  noon on January 1, 4713 BC (proleptic Julian calendar)\n');

    // Reference epochs
    const J2000_JD = 2451545.0;  // January 1, 2000, 12:00 TT
    const MJD_OFFSET = 2400000.5;  // Modified Julian Date offset
    const UNIX_EPOCH_JD = 2440587.5;  // January 1, 1970, 00:00 UTC

    console.log('Reference epochs:');
    console.log(`  J2000.0: JD ${J2000_JD} (Jan 1, 2000, 12:00 TT)`);
    console.log(`  MJD epoch: JD ${MJD_OFFSET} (Nov 17, 1858)`);
    console.log(`  Unix epoch: JD ${UNIX_EPOCH_JD} (Jan 1, 1970, 00:00 UTC)\n`);

    // Convert calendar date to Julian Date
    function calendarToJD(year, month, day, hour = 12) {
        // Algorithm from Meeus, "Astronomical Algorithms"
        if (month <= 2) {
            year -= 1;
            month += 12;
        }
        const A = Math.floor(year / 100);
        const B = 2 - A + Math.floor(A / 4);
        const JD = Math.floor(365.25 * (year + 4716)) +
                   Math.floor(30.6001 * (month + 1)) +
                   day + hour / 24 + B - 1524.5;
        return JD;
    }

    // Convert Julian Date to calendar date
    function jdToCalendar(jd) {
        const Z = Math.floor(jd + 0.5);
        const F = jd + 0.5 - Z;
        let A;
        if (Z < 2299161) {
            A = Z;
        } else {
            const alpha = Math.floor((Z - 1867216.25) / 36524.25);
            A = Z + 1 + alpha - Math.floor(alpha / 4);
        }
        const B = A + 1524;
        const C = Math.floor((B - 122.1) / 365.25);
        const D = Math.floor(365.25 * C);
        const E = Math.floor((B - D) / 30.6001);

        const day = B - D - Math.floor(30.6001 * E) + F;
        const month = E < 14 ? E - 1 : E - 13;
        const year = month > 2 ? C - 4716 : C - 4715;

        return { year, month, day: Math.floor(day), hour: (day % 1) * 24 };
    }

    // Example conversions
    console.log('Calendar to Julian Date conversions:');
    const dates = [
        { year: 2000, month: 1, day: 1, hour: 12 },   // J2000.0
        { year: 2024, month: 7, day: 4, hour: 0 },    // Recent date
        { year: 1969, month: 7, day: 20, hour: 20 },  // Moon landing
    ];

    for (const d of dates) {
        const jd = calendarToJD(d.year, d.month, d.day, d.hour);
        const mjd = jd - MJD_OFFSET;
        console.log(`  ${d.year}-${String(d.month).padStart(2, '0')}-${String(d.day).padStart(2, '0')} ${String(d.hour).padStart(2, '0')}:00 → JD ${jd.toFixed(5)}, MJD ${mjd.toFixed(5)}`);
    }

    // ========================================
    // Epoch Conversions
    // ========================================
    console.log('\n\n3. Epoch Conversions');
    console.log('====================\n');

    // Convert between different epoch representations
    function jdToSecondsFromJ2000(jd) {
        return (jd - J2000_JD) * 86400;
    }

    function secondsFromJ2000ToJD(seconds) {
        return J2000_JD + seconds / 86400;
    }

    function jdToUnixTime(jd) {
        return (jd - UNIX_EPOCH_JD) * 86400;
    }

    function unixTimeToJD(unix) {
        return UNIX_EPOCH_JD + unix / 86400;
    }

    console.log('Epoch conversion examples:');
    console.log('  J2000.0 in seconds from J2000: 0');
    console.log('  J2000.0 as Unix timestamp: ' + jdToUnixTime(J2000_JD).toFixed(0));

    const now = new Date();
    const unixNow = now.getTime() / 1000;
    const jdNow = unixTimeToJD(unixNow);
    const secFromJ2000 = jdToSecondsFromJ2000(jdNow);

    console.log(`\n  Current time (Unix): ${unixNow.toFixed(0)}`);
    console.log(`  Current time (JD): ${jdNow.toFixed(5)}`);
    console.log(`  Current time (sec from J2000): ${secFromJ2000.toFixed(0)}`);
    console.log(`  Current time (days from J2000): ${(secFromJ2000 / 86400).toFixed(3)}`);

    // ========================================
    // Sidereal Time
    // ========================================
    console.log('\n\n4. Sidereal Time');
    console.log('================\n');

    console.log('Greenwich Mean Sidereal Time (GMST) measures Earth\'s rotation');
    console.log('relative to the vernal equinox.\n');

    // GMST calculation (simplified)
    function computeGMST(jd) {
        // Julian centuries from J2000.0
        const T = (jd - J2000_JD) / 36525;

        // GMST at 0h UT (in degrees)
        let gmst = 280.46061837 +
                   360.98564736629 * (jd - J2000_JD) +
                   0.000387933 * T * T -
                   T * T * T / 38710000;

        // Normalize to 0-360
        gmst = ((gmst % 360) + 360) % 360;

        return gmst;
    }

    // GMST at J2000.0
    const gmstJ2000 = computeGMST(J2000_JD);
    console.log(`GMST at J2000.0: ${gmstJ2000.toFixed(4)}°`);

    // GMST progression over a day
    console.log('\nGMST progression over one solar day:');
    console.log('Hour (UT)    GMST (degrees)    GMST (hours)');
    console.log('---------    --------------    ------------');

    const testJD = calendarToJD(2024, 3, 20, 0);  // Vernal equinox
    for (let h = 0; h <= 24; h += 6) {
        const jd = testJD + h / 24;
        const gmst = computeGMST(jd);
        const gmstHours = gmst / 15;
        console.log(`${String(h).padStart(9)}    ${gmst.toFixed(2).padStart(14)}    ${gmstHours.toFixed(2).padStart(12)}`);
    }

    console.log('\nNote: GMST advances ~361° per solar day (sidereal day is shorter)');

    // ========================================
    // Earth Orientation Parameters
    // ========================================
    console.log('\n\n5. Earth Orientation Parameters (EOP)');
    console.log('=====================================\n');

    console.log('EOPs describe Earth\'s actual orientation vs. the ideal model:\n');
    console.log('  xp, yp    - Polar motion (position of rotation pole)');
    console.log('  UT1-UTC   - Earth rotation angle correction');
    console.log('  dX, dY    - Celestial pole offsets (nutation corrections)');
    console.log('  LOD       - Length of day excess\n');

    // Typical EOP values
    console.log('Typical EOP magnitudes:');
    console.log('  Polar motion xp, yp: ±0.5 arcsec (varies over ~430 days)');
    console.log('  UT1-UTC: ±0.9 seconds (kept within bounds by leap seconds)');
    console.log('  dX, dY: ±0.001 arcsec');
    console.log('  LOD excess: ~1-3 milliseconds\n');

    // Impact on positioning
    const polarMotionArcsec = 0.3;  // typical value
    const earthRadius = 6378137;  // meters
    const polarMotionRad = polarMotionArcsec * Math.PI / (180 * 3600);
    const surfaceDisplacement = earthRadius * polarMotionRad;

    console.log('Impact on surface positioning:');
    console.log(`  Polar motion of ${polarMotionArcsec} arcsec → ${surfaceDisplacement.toFixed(2)} m at equator`);
    console.log(`  UT1-UTC error of 1 second → ~465 m at equator (Earth rotation)`);

    // ========================================
    // Frame Transformations
    // ========================================
    console.log('\n\n6. Reference Frame Transformations');
    console.log('===================================\n');

    console.log('Key reference frames in astrodynamics:\n');
    console.log('  GCRF  - Geocentric Celestial Reference Frame (inertial)');
    console.log('  ITRF  - International Terrestrial Reference Frame (Earth-fixed)');
    console.log('  ICRF  - International Celestial Reference Frame (J2000-based)');
    console.log('  EME2000 - Earth Mean Equator and Equinox of J2000\n');

    // GCRF to ITRF transformation involves:
    console.log('GCRF to ITRF transformation chain:');
    console.log('  1. Precession (P) - Slow drift of equinox');
    console.log('  2. Nutation (N) - Periodic wobble of Earth axis');
    console.log('  3. Earth Rotation (R) - Sidereal rotation');
    console.log('  4. Polar Motion (W) - Wobble of rotation pole\n');
    console.log('  ITRF = W · R · N · P · GCRF\n');

    // Earth rotation rate
    const omegaEarth = 7.292115e-5;  // rad/s
    console.log('Earth rotation parameters:');
    console.log(`  Angular velocity: ${omegaEarth.toExponential(6)} rad/s`);
    console.log(`  Sidereal day: ${(2 * Math.PI / omegaEarth / 3600).toFixed(6)} hours`);
    console.log(`  Surface velocity at equator: ${(omegaEarth * earthRadius).toFixed(2)} m/s`);

    // ========================================
    // Light Time and Aberration
    // ========================================
    console.log('\n\n7. Light Time Corrections');
    console.log('=========================\n');

    console.log('When observing distant objects, light travel time matters:\n');

    const c = 299792458;  // m/s
    const AU = 149597870700;  // m

    const distances = [
        { name: 'Moon', dist: 384400e3 },
        { name: 'Sun', dist: AU },
        { name: 'Mars (closest)', dist: 54.6e9 },
        { name: 'Jupiter', dist: 588e9 },
        { name: 'Voyager 1', dist: 24e12 }
    ];

    console.log('Object           Distance        Light time');
    console.log('------           --------        ----------');

    for (const obj of distances) {
        const lightTime = obj.dist / c;
        let timeStr;
        if (lightTime < 60) {
            timeStr = `${lightTime.toFixed(2)} sec`;
        } else if (lightTime < 3600) {
            timeStr = `${(lightTime / 60).toFixed(2)} min`;
        } else {
            timeStr = `${(lightTime / 3600).toFixed(2)} hr`;
        }
        console.log(`${obj.name.padEnd(16)} ${(obj.dist / 1e9).toFixed(2).padStart(10)} km    ${timeStr.padStart(10)}`);
    }

    console.log('\nAberration: Apparent shift in object position due to observer motion');
    console.log(`  Earth orbital velocity: ~30 km/s`);
    console.log(`  Maximum stellar aberration: ~20.5 arcsec`);

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
