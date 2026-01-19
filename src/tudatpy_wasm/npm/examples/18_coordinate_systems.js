/**
 * Example 18: Coordinate Systems and Transformations
 *
 * This example demonstrates various coordinate systems used in
 * astrodynamics and transformations between them.
 *
 * Run with: node 18_coordinate_systems.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Coordinate Systems and Transformations ===\n');

    const tudat = await createTudatModule();

    // ========================================
    // Cartesian Coordinates
    // ========================================
    console.log('1. Cartesian Coordinates');
    console.log('========================\n');

    console.log('Cartesian (X, Y, Z) coordinates are the fundamental representation');
    console.log('for position and velocity in inertial reference frames.\n');

    const earthGM = 3.986004418e14;
    const earthRadius = 6378137;

    // Example position (ISS-like orbit)
    const position = [6778000, 0, 0];  // m
    const velocity = [0, 7672, 0];      // m/s

    console.log('Example state (Cartesian):');
    console.log(`  Position: [${position[0]/1000}, ${position[1]/1000}, ${position[2]/1000}] km`);
    console.log(`  Velocity: [${velocity[0]/1000}, ${velocity[1]/1000}, ${velocity[2]/1000}] km/s`);

    const r = Math.sqrt(position[0]**2 + position[1]**2 + position[2]**2);
    const v = Math.sqrt(velocity[0]**2 + velocity[1]**2 + velocity[2]**2);
    console.log(`  |r| = ${(r/1000).toFixed(3)} km`);
    console.log(`  |v| = ${(v/1000).toFixed(4)} km/s`);

    // ========================================
    // Spherical Coordinates
    // ========================================
    console.log('\n\n2. Spherical Coordinates');
    console.log('========================\n');

    console.log('Spherical coordinates (r, θ, φ) are natural for radial problems.');
    console.log('  r: radial distance');
    console.log('  θ: azimuth angle (in x-y plane from x-axis)');
    console.log('  φ: elevation angle (from x-y plane)\n');

    function cartesianToSpherical(x, y, z) {
        const r = Math.sqrt(x*x + y*y + z*z);
        const theta = Math.atan2(y, x);
        const phi = Math.asin(z / r);
        return { r, theta, phi };
    }

    function sphericalToCartesian(r, theta, phi) {
        const x = r * Math.cos(phi) * Math.cos(theta);
        const y = r * Math.cos(phi) * Math.sin(theta);
        const z = r * Math.sin(phi);
        return { x, y, z };
    }

    // Convert example position
    const spherical = cartesianToSpherical(position[0], position[1], position[2]);
    console.log('Example position in spherical:');
    console.log(`  r = ${(spherical.r / 1000).toFixed(3)} km`);
    console.log(`  θ = ${(spherical.theta * 180 / Math.PI).toFixed(4)}°`);
    console.log(`  φ = ${(spherical.phi * 180 / Math.PI).toFixed(4)}°`);

    // Verify round-trip
    const backToCart = sphericalToCartesian(spherical.r, spherical.theta, spherical.phi);
    console.log(`  Round-trip: [${(backToCart.x/1000).toFixed(3)}, ${(backToCart.y/1000).toFixed(3)}, ${(backToCart.z/1000).toFixed(3)}] km ✓`);

    // ========================================
    // Geodetic Coordinates
    // ========================================
    console.log('\n\n3. Geodetic Coordinates');
    console.log('=======================\n');

    console.log('Geodetic coordinates (latitude, longitude, altitude) are used');
    console.log('for positions relative to Earth\'s surface.\n');

    const earthFlattening = 1 / 298.257223563;
    const earthSemiMajor = 6378137.0;
    const earthSemiMinor = earthSemiMajor * (1 - earthFlattening);
    const earthEccentricitySq = 1 - (earthSemiMinor / earthSemiMajor) ** 2;

    function cartesianToGeodetic(x, y, z) {
        // Iterative algorithm
        const p = Math.sqrt(x*x + y*y);
        const lon = Math.atan2(y, x);

        // Initial estimate
        let lat = Math.atan2(z, p * (1 - earthEccentricitySq));
        let N, alt;

        // Iterate
        for (let i = 0; i < 10; i++) {
            const sinLat = Math.sin(lat);
            N = earthSemiMajor / Math.sqrt(1 - earthEccentricitySq * sinLat * sinLat);
            alt = p / Math.cos(lat) - N;
            lat = Math.atan2(z, p * (1 - earthEccentricitySq * N / (N + alt)));
        }

        return {
            lat: lat * 180 / Math.PI,
            lon: lon * 180 / Math.PI,
            alt: alt
        };
    }

    // Example: Convert spacecraft position to geodetic
    const geodetic = cartesianToGeodetic(position[0], position[1], position[2]);
    console.log('Example position in geodetic:');
    console.log(`  Latitude: ${geodetic.lat.toFixed(4)}°`);
    console.log(`  Longitude: ${geodetic.lon.toFixed(4)}°`);
    console.log(`  Altitude: ${(geodetic.alt / 1000).toFixed(3)} km`);

    // Ground track point example
    const groundPositions = [
        { name: 'Kennedy Space Center', lat: 28.5729, lon: -80.6490 },
        { name: 'Baikonur Cosmodrome', lat: 45.9650, lon: 63.3050 },
        { name: 'Guiana Space Centre', lat: 5.2320, lon: -52.7756 },
        { name: 'Tanegashima', lat: 30.4000, lon: 130.9667 }
    ];

    console.log('\nLaunch site coordinates:');
    for (const site of groundPositions) {
        console.log(`  ${site.name}: ${site.lat.toFixed(2)}°N, ${Math.abs(site.lon).toFixed(2)}°${site.lon >= 0 ? 'E' : 'W'}`);
    }

    // ========================================
    // RSW (Radial-Along-track-Cross-track)
    // ========================================
    console.log('\n\n4. RSW Coordinate System');
    console.log('========================\n');

    console.log('RSW (also called RIC or RTN) is an orbit-centered frame:');
    console.log('  R: Radial (outward from central body)');
    console.log('  S: Along-track (in velocity direction)');
    console.log('  W: Cross-track (orbit normal)\n');

    // Compute RSW basis vectors
    const rHat = [position[0]/r, position[1]/r, position[2]/r];

    // Angular momentum vector (r × v)
    const h = [
        position[1]*velocity[2] - position[2]*velocity[1],
        position[2]*velocity[0] - position[0]*velocity[2],
        position[0]*velocity[1] - position[1]*velocity[0]
    ];
    const hMag = Math.sqrt(h[0]**2 + h[1]**2 + h[2]**2);
    const wHat = [h[0]/hMag, h[1]/hMag, h[2]/hMag];

    // S = W × R
    const sHat = [
        wHat[1]*rHat[2] - wHat[2]*rHat[1],
        wHat[2]*rHat[0] - wHat[0]*rHat[2],
        wHat[0]*rHat[1] - wHat[1]*rHat[0]
    ];

    console.log('RSW basis vectors (inertial components):');
    console.log(`  R: [${rHat[0].toFixed(4)}, ${rHat[1].toFixed(4)}, ${rHat[2].toFixed(4)}]`);
    console.log(`  S: [${sHat[0].toFixed(4)}, ${sHat[1].toFixed(4)}, ${sHat[2].toFixed(4)}]`);
    console.log(`  W: [${wHat[0].toFixed(4)}, ${wHat[1].toFixed(4)}, ${wHat[2].toFixed(4)}]`);

    // Express velocity in RSW
    const vR = velocity[0]*rHat[0] + velocity[1]*rHat[1] + velocity[2]*rHat[2];
    const vS = velocity[0]*sHat[0] + velocity[1]*sHat[1] + velocity[2]*sHat[2];
    const vW = velocity[0]*wHat[0] + velocity[1]*wHat[1] + velocity[2]*wHat[2];

    console.log(`\nVelocity in RSW frame:`);
    console.log(`  v_R (radial): ${(vR / 1000).toFixed(4)} km/s`);
    console.log(`  v_S (along-track): ${(vS / 1000).toFixed(4)} km/s`);
    console.log(`  v_W (cross-track): ${(vW / 1000).toFixed(4)} km/s`);

    // ========================================
    // Perifocal (PQW) Coordinates
    // ========================================
    console.log('\n\n5. Perifocal (PQW) Coordinates');
    console.log('==============================\n');

    console.log('The perifocal frame is aligned with the orbit:');
    console.log('  P: Points to periapsis');
    console.log('  Q: In orbital plane, 90° ahead of P');
    console.log('  W: Orbit normal (same as RSW)\n');

    // For a circular orbit at true anomaly = 0, P ≈ R and Q ≈ S
    console.log('For the example circular orbit:');
    console.log('  P ≈ R (at periapsis/true anomaly = 0)');
    console.log('  Q ≈ S');
    console.log('  W = orbit normal\n');

    // Position in perifocal frame: r = p / (1 + e*cos(ν)) in P-Q plane
    // For circular orbit: r = a (constant)
    console.log('Position in perifocal frame (true anomaly = 0°):');
    console.log(`  p = ${(r / 1000).toFixed(3)} km, q = 0 km, w = 0 km`);

    // ========================================
    // Topocentric (ENU) Coordinates
    // ========================================
    console.log('\n\n6. Topocentric (ENU) Coordinates');
    console.log('================================\n');

    console.log('ENU (East-North-Up) is used for ground-based observations:');
    console.log('  E: East (local horizon)');
    console.log('  N: North (local horizon)');
    console.log('  U: Up (zenith)\n');

    // Observer at Cape Canaveral
    const observerLat = 28.5 * Math.PI / 180;
    const observerLon = -80.5 * Math.PI / 180;
    const observerAlt = 0;

    // Rotation from ECEF to ENU
    // This is a simplified calculation
    console.log(`Observer at Cape Canaveral (${28.5}°N, ${80.5}°W):`);
    console.log('  ENU frame rotated relative to ECEF by:');
    console.log(`    Longitude rotation: ${(-80.5).toFixed(2)}° about Z`);
    console.log(`    Latitude rotation: ${(90 - 28.5).toFixed(2)}° about Y`);

    // Azimuth-Elevation from position
    // (Simplified - assuming observer at origin and satellite at position)
    const elevation = Math.asin(position[2] / r);
    const azimuth = Math.atan2(position[1], position[0]);

    console.log(`\nSatellite as seen from observer (simplified):`);
    console.log(`  Azimuth: ${(azimuth * 180 / Math.PI).toFixed(2)}°`);
    console.log(`  Elevation: ${(elevation * 180 / Math.PI).toFixed(2)}°`);
    console.log(`  Range: ${(r / 1000).toFixed(2)} km`);

    // ========================================
    // Body-Fixed vs Inertial
    // ========================================
    console.log('\n\n7. Body-Fixed vs Inertial Frames');
    console.log('================================\n');

    console.log('Key distinction:');
    console.log('  Inertial: Fixed relative to distant stars (GCRF, J2000)');
    console.log('  Body-fixed: Rotates with the body (ITRF, planet-fixed)\n');

    const earthRotRate = 7.292115e-5;  // rad/s

    console.log('Earth rotation effects:');
    console.log(`  Angular velocity: ${(earthRotRate * 1e6).toFixed(3)} μrad/s`);
    console.log(`  Period: ${(2 * Math.PI / earthRotRate / 3600).toFixed(4)} hours (sidereal day)`);
    console.log(`  Surface velocity at equator: ${(earthRotRate * earthRadius).toFixed(2)} m/s\n`);

    // Velocity transformation
    console.log('Velocity transformation from inertial to body-fixed:');
    console.log('  v_fixed = v_inertial - ω × r');
    console.log('  where ω = [0, 0, 7.292115e-5] rad/s\n');

    // Example: velocity correction for a satellite
    const omegaCrossR = [
        -earthRotRate * position[1],
        earthRotRate * position[0],
        0
    ];
    const vFixed = [
        velocity[0] - omegaCrossR[0],
        velocity[1] - omegaCrossR[1],
        velocity[2] - omegaCrossR[2]
    ];

    console.log('Example satellite velocity:');
    console.log(`  Inertial: [${(velocity[0]/1000).toFixed(4)}, ${(velocity[1]/1000).toFixed(4)}, ${(velocity[2]/1000).toFixed(4)}] km/s`);
    console.log(`  Body-fixed: [${(vFixed[0]/1000).toFixed(4)}, ${(vFixed[1]/1000).toFixed(4)}, ${(vFixed[2]/1000).toFixed(4)}] km/s`);

    // ========================================
    // Keplerian Elements (Using Tudat)
    // ========================================
    console.log('\n\n8. Keplerian Elements Conversion');
    console.log('=================================\n');

    // Create Cartesian state vector
    const cartesian = new tudat.Vector6d();
    cartesian.set(0, position[0]);
    cartesian.set(1, position[1]);
    cartesian.set(2, position[2]);
    cartesian.set(3, velocity[0]);
    cartesian.set(4, velocity[1]);
    cartesian.set(5, velocity[2]);

    // Convert to Keplerian
    const keplerian = tudat.astro.element_conversion.cartesian_to_keplerian(cartesian, earthGM);

    console.log('Cartesian to Keplerian conversion (using Tudat):');
    console.log(`  Semi-major axis: ${(keplerian.get(0) / 1000).toFixed(3)} km`);
    console.log(`  Eccentricity: ${keplerian.get(1).toFixed(6)}`);
    console.log(`  Inclination: ${(keplerian.get(2) * 180 / Math.PI).toFixed(4)}°`);
    console.log(`  Arg of periapsis: ${(keplerian.get(3) * 180 / Math.PI).toFixed(4)}°`);
    console.log(`  RAAN: ${(keplerian.get(4) * 180 / Math.PI).toFixed(4)}°`);
    console.log(`  True anomaly: ${(keplerian.get(5) * 180 / Math.PI).toFixed(4)}°`);

    // Clean up
    cartesian.delete();
    keplerian.delete();

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
