/**
 * Example 17: Escape and Capture Maneuvers
 *
 * This example demonstrates escape and capture trajectories,
 * hyperbolic excess velocity, and sphere of influence calculations.
 *
 * Run with: node 17_escape_capture.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Escape and Capture Maneuvers ===\n');

    const tudat = await createTudatModule();

    // Constants
    const sunGM = 1.32712440018e20;    // m³/s²
    const earthGM = 3.986004418e14;    // m³/s²
    const marsGM = 4.282837e13;        // m³/s²
    const moonGM = 4.9028695e12;       // m³/s²

    const earthRadius = 6378137;       // m
    const marsRadius = 3396200;        // m
    const AU = 149597870700;           // m

    // ========================================
    // Sphere of Influence
    // ========================================
    console.log('1. Sphere of Influence (SOI)');
    console.log('============================\n');

    console.log('The SOI defines the region where a body\'s gravity dominates.');
    console.log('Beyond the SOI, the Sun\'s gravity becomes dominant.\n');

    // SOI radius: r_SOI = a * (m/M)^(2/5)
    function sphereOfInfluence(sma, planetMass, sunMass) {
        return sma * Math.pow(planetMass / sunMass, 2/5);
    }

    const bodies = [
        { name: 'Mercury', sma: 0.387 * AU, mass: 3.302e23 },
        { name: 'Venus', sma: 0.723 * AU, mass: 4.869e24 },
        { name: 'Earth', sma: 1.000 * AU, mass: 5.972e24 },
        { name: 'Moon (from Earth)', sma: 384400e3, mass: 7.342e22, parent: earthGM },
        { name: 'Mars', sma: 1.524 * AU, mass: 6.417e23 },
        { name: 'Jupiter', sma: 5.203 * AU, mass: 1.898e27 },
        { name: 'Saturn', sma: 9.537 * AU, mass: 5.683e26 }
    ];

    const sunMass = 1.989e30;

    console.log('Body          SMA (AU)    SOI Radius');
    console.log('----          --------    ----------');

    for (const body of bodies) {
        let soi;
        if (body.parent) {
            // Moon's SOI relative to Earth
            soi = sphereOfInfluence(body.sma, body.mass, 5.972e24);
            console.log(`${body.name.padEnd(13)} ${(body.sma / 1e6).toFixed(0).padStart(8)} km    ${(soi / 1e3).toFixed(0)} km`);
        } else {
            soi = sphereOfInfluence(body.sma, body.mass, sunMass);
            const soiAU = soi / AU;
            if (soiAU > 0.1) {
                console.log(`${body.name.padEnd(13)} ${(body.sma / AU).toFixed(3).padStart(8)}       ${soiAU.toFixed(3)} AU (${(soi / 1e9).toFixed(1)} million km)`);
            } else {
                console.log(`${body.name.padEnd(13)} ${(body.sma / AU).toFixed(3).padStart(8)}       ${(soi / 1e6).toFixed(0)} km`);
            }
        }
    }

    const earthSOI = sphereOfInfluence(AU, 5.972e24, sunMass);
    console.log(`\nEarth's SOI: ${(earthSOI / 1e6).toFixed(0)} km = ${(earthSOI / AU).toFixed(4)} AU`);
    console.log(`This is where patched conic approximation switches reference frame.`);

    // ========================================
    // Escape Velocity
    // ========================================
    console.log('\n\n2. Escape Velocity');
    console.log('==================\n');

    console.log('Escape velocity is the minimum velocity to leave a body\'s');
    console.log('gravitational influence (reach zero velocity at infinity).\n');

    // v_esc = sqrt(2 * GM / r)
    function escapeVelocity(gm, radius) {
        return Math.sqrt(2 * gm / radius);
    }

    const escapeData = [
        { name: 'Earth (surface)', gm: earthGM, r: earthRadius },
        { name: 'Earth (LEO 400km)', gm: earthGM, r: earthRadius + 400e3 },
        { name: 'Earth (GEO)', gm: earthGM, r: earthRadius + 35786e3 },
        { name: 'Moon (surface)', gm: moonGM, r: 1737400 },
        { name: 'Mars (surface)', gm: marsGM, r: marsRadius },
        { name: 'Sun (Earth orbit)', gm: sunGM, r: AU }
    ];

    console.log('Location             Escape Velocity');
    console.log('--------             ---------------');

    for (const data of escapeData) {
        const vEsc = escapeVelocity(data.gm, data.r);
        console.log(`${data.name.padEnd(20)} ${(vEsc / 1000).toFixed(3)} km/s`);
    }

    // ========================================
    // Hyperbolic Excess Velocity (C3)
    // ========================================
    console.log('\n\n3. Hyperbolic Excess Velocity (V∞) and C3');
    console.log('==========================================\n');

    console.log('C3 = V∞² is the characteristic energy for interplanetary travel.');
    console.log('It represents the kinetic energy per unit mass at infinity.\n');

    // V_infinity for various missions
    const missions = [
        { name: 'Earth escape (minimum)', c3: 0 },
        { name: 'Lunar mission', c3: 0.5 },
        { name: 'Mars mission (low)', c3: 8 },
        { name: 'Mars mission (typical)', c3: 12 },
        { name: 'Venus mission', c3: 7 },
        { name: 'Jupiter mission', c3: 80 },
        { name: 'Solar escape (from Earth)', c3: 152 },
        { name: 'Interstellar (Voyager)', c3: 100 }
    ];

    console.log('Mission                     C3 (km²/s²)    V∞ (km/s)');
    console.log('-------                     -----------    ---------');

    for (const mission of missions) {
        const vInf = Math.sqrt(mission.c3);
        console.log(`${mission.name.padEnd(27)} ${mission.c3.toFixed(1).padStart(11)}    ${vInf.toFixed(2).padStart(9)}`);
    }

    // ========================================
    // Escape Maneuver from Parking Orbit
    // ========================================
    console.log('\n\n4. Escape from Parking Orbit');
    console.log('============================\n');

    const parkingAltitude = 400e3;  // 400 km
    const parkingRadius = earthRadius + parkingAltitude;

    // Circular velocity in parking orbit
    const vCircular = Math.sqrt(earthGM / parkingRadius);

    // Escape velocity at parking orbit
    const vEscape = escapeVelocity(earthGM, parkingRadius);

    console.log(`Parking orbit: ${parkingAltitude / 1000} km altitude (circular)`);
    console.log(`  Orbital velocity: ${(vCircular / 1000).toFixed(3)} km/s`);
    console.log(`  Escape velocity: ${(vEscape / 1000).toFixed(3)} km/s`);
    console.log(`  Minimum ΔV for escape: ${((vEscape - vCircular) / 1000).toFixed(3)} km/s\n`);

    // Required velocity for different C3 values
    // v_periapsis = sqrt(v_esc² + C3)
    console.log('Required periapsis velocity for different C3 values:');
    console.log('C3 (km²/s²)    V_periapsis (km/s)    ΔV from LEO (km/s)');
    console.log('-----------    ------------------    ------------------');

    const c3Values = [0, 5, 10, 15, 20, 30, 50, 100];

    for (const c3 of c3Values) {
        const c3_si = c3 * 1e6;  // Convert to m²/s²
        const vPeriapsis = Math.sqrt(vEscape * vEscape + c3_si);
        const deltaV = vPeriapsis - vCircular;
        console.log(`${c3.toString().padStart(11)}    ${(vPeriapsis / 1000).toFixed(3).padStart(18)}    ${(deltaV / 1000).toFixed(3).padStart(18)}`);
    }

    // ========================================
    // Hyperbolic Trajectory
    // ========================================
    console.log('\n\n5. Hyperbolic Departure Trajectory');
    console.log('===================================\n');

    // Example: Mars mission with C3 = 10 km²/s²
    const targetC3 = 10e6;  // m²/s²
    const vInfinity = Math.sqrt(targetC3);

    console.log(`Target C3: ${targetC3 / 1e6} km²/s² (V∞ = ${(vInfinity / 1000).toFixed(3)} km/s)`);

    // Hyperbolic orbit parameters
    // Energy: E = C3 / 2 = -GM/2a  =>  a = -GM / C3
    const a_hyp = -earthGM / targetC3;  // Negative for hyperbola
    const e_hyp = 1 - parkingRadius / a_hyp;

    // Asymptotic angle
    const asymptoticAngle = Math.acos(-1 / e_hyp);

    // Turn angle
    const turnAngle = Math.PI - 2 * asymptoticAngle;

    console.log(`\nHyperbolic orbit parameters:`);
    console.log(`  Semi-major axis: ${(a_hyp / 1000).toFixed(0)} km (negative for hyperbola)`);
    console.log(`  Eccentricity: ${e_hyp.toFixed(4)}`);
    console.log(`  Asymptotic angle: ${(asymptoticAngle * 180 / Math.PI).toFixed(2)}°`);
    console.log(`  Turn angle: ${(turnAngle * 180 / Math.PI).toFixed(2)}°`);

    // Velocity at periapsis
    const vPeriapsis = Math.sqrt(earthGM * (2 / parkingRadius - 1 / a_hyp));
    console.log(`  Velocity at periapsis: ${(vPeriapsis / 1000).toFixed(3)} km/s`);
    console.log(`  Required ΔV from LEO: ${((vPeriapsis - vCircular) / 1000).toFixed(3)} km/s`);

    // ========================================
    // Capture at Target Planet
    // ========================================
    console.log('\n\n6. Capture at Mars');
    console.log('==================\n');

    // Arriving at Mars with V∞ = 2.5 km/s (typical)
    const marsVInf = 2.5e3;  // m/s
    const captureAltitude = 400e3;  // 400 km
    const captureRadius = marsRadius + captureAltitude;

    console.log(`Approach conditions:`);
    console.log(`  V∞ at Mars: ${(marsVInf / 1000).toFixed(2)} km/s`);
    console.log(`  Target periapsis: ${captureAltitude / 1000} km altitude\n`);

    // Velocity at periapsis on hyperbolic approach
    const vMarsEscape = escapeVelocity(marsGM, captureRadius);
    const vApproachPeriapsis = Math.sqrt(marsVInf * marsVInf + vMarsEscape * vMarsEscape);

    // Circular orbit velocity at capture altitude
    const vMarsCircular = Math.sqrt(marsGM / captureRadius);

    // Direct capture to circular orbit
    const deltaVCapture = vApproachPeriapsis - vMarsCircular;

    console.log(`Capture maneuver (direct to circular):`);
    console.log(`  Mars escape velocity at periapsis: ${(vMarsEscape / 1000).toFixed(3)} km/s`);
    console.log(`  Approach velocity at periapsis: ${(vApproachPeriapsis / 1000).toFixed(3)} km/s`);
    console.log(`  Circular orbit velocity: ${(vMarsCircular / 1000).toFixed(3)} km/s`);
    console.log(`  ΔV for capture: ${(deltaVCapture / 1000).toFixed(3)} km/s`);

    // Elliptical capture (more efficient)
    console.log(`\nElliptical capture (apoapsis at 33000 km altitude):`);
    const ellipseApoapsis = marsRadius + 33000e3;
    const vEllipsePeriapsis = Math.sqrt(marsGM * (2 / captureRadius - 2 / (captureRadius + ellipseApoapsis)));
    const deltaVElliptic = vApproachPeriapsis - vEllipsePeriapsis;

    console.log(`  ΔV for elliptical capture: ${(deltaVElliptic / 1000).toFixed(3)} km/s`);
    console.log(`  ΔV savings: ${((deltaVCapture - deltaVElliptic) / 1000).toFixed(3)} km/s`);

    // ========================================
    // Oberth Effect
    // ========================================
    console.log('\n\n7. Oberth Effect');
    console.log('================\n');

    console.log('The Oberth effect: maneuvers are more efficient at higher speeds.');
    console.log('ΔV applied at periapsis produces more C3 than at apoapsis.\n');

    // Compare: same ΔV at periapsis vs at higher altitude
    const deltaV_maneuver = 1.0e3;  // 1 km/s

    // At periapsis (LEO)
    const v1_low = vCircular + deltaV_maneuver;
    const energy1 = 0.5 * v1_low * v1_low - earthGM / parkingRadius;
    const c3_low = 2 * energy1;

    // At higher altitude (say 10000 km)
    const highAltitude = 10000e3;
    const highRadius = earthRadius + highAltitude;
    const vCircularHigh = Math.sqrt(earthGM / highRadius);
    const v1_high = vCircularHigh + deltaV_maneuver;
    const energy2 = 0.5 * v1_high * v1_high - earthGM / highRadius;
    const c3_high = 2 * energy2;

    console.log(`Effect of 1 km/s ΔV at different altitudes:`);
    console.log(`  At LEO (${parkingAltitude / 1000} km):`);
    console.log(`    Initial velocity: ${(vCircular / 1000).toFixed(3)} km/s`);
    console.log(`    Final velocity: ${(v1_low / 1000).toFixed(3)} km/s`);
    console.log(`    Resulting C3: ${(c3_low / 1e6).toFixed(2)} km²/s²`);
    console.log(`  At ${highAltitude / 1000} km altitude:`);
    console.log(`    Initial velocity: ${(vCircularHigh / 1000).toFixed(3)} km/s`);
    console.log(`    Final velocity: ${(v1_high / 1000).toFixed(3)} km/s`);
    console.log(`    Resulting C3: ${(c3_high / 1e6).toFixed(2)} km²/s²`);
    console.log(`\n  C3 gain at low altitude: ${((c3_low - c3_high) / 1e6).toFixed(2)} km²/s² more`);
    console.log(`  This is why interplanetary burns happen at periapsis!`);

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
