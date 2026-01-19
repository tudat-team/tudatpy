/**
 * Example 40: Asteroid Orbit Optimization
 *
 * Ported from: examples/tudatpy/pygmo/asteroid_orbit_optimization
 *
 * This example demonstrates optimization of an asteroid observation
 * mission orbit using evolutionary algorithms.
 *
 * Key concepts:
 * - Mission orbit optimization
 * - Multi-objective optimization
 * - Observation geometry
 * - Coverage analysis
 * - Constraint handling
 *
 * Run with: node 40_asteroid_optimization.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Asteroid Orbit Optimization ===\n');

    const tudat = await createTudatModule();

    // Asteroid parameters (Eros-like)
    const asteroid = {
        name: 'Eros',
        GM: 4.463e5,           // m^3/s^2
        radius: 8400,          // m (mean radius)
        rotationPeriod: 5.27 * 3600,  // s
        a: 1.458 * 1.496e11,   // AU -> m (heliocentric)
        e: 0.223,
        dimensions: [34400, 11200, 11200]  // m (elongated shape)
    };

    console.log('Target Asteroid:');
    console.log(`  Name: ${asteroid.name}`);
    console.log(`  GM: ${asteroid.GM.toExponential(3)} m³/s²`);
    console.log(`  Mean radius: ${asteroid.radius} m`);
    console.log(`  Rotation period: ${(asteroid.rotationPeriod / 3600).toFixed(2)} hours`);
    console.log(`  Dimensions: ${asteroid.dimensions.join(' × ')} m`);

    // Mission requirements
    const mission = {
        minAltitude: 5000,      // m (above surface)
        maxAltitude: 50000,     // m
        minSunAngle: 30,        // deg (illumination constraint)
        coverageTarget: 0.90,   // 90% surface coverage goal
        observationDuration: 30 * 86400,  // 30 days
        groundResolution: 1     // m/pixel target
    };

    console.log('\nMission Requirements:');
    console.log(`  Altitude range: ${mission.minAltitude/1000} - ${mission.maxAltitude/1000} km`);
    console.log(`  Min sun angle: ${mission.minSunAngle}°`);
    console.log(`  Coverage target: ${mission.coverageTarget * 100}%`);
    console.log(`  Duration: ${mission.observationDuration / 86400} days`);

    // Orbit design variables
    // x = [sma, ecc, inc, raan, aop]
    const bounds = {
        sma: [asteroid.radius + mission.minAltitude, asteroid.radius + mission.maxAltitude],
        ecc: [0, 0.3],
        inc: [0, Math.PI],
        raan: [0, 2 * Math.PI],
        aop: [0, 2 * Math.PI]
    };

    // Orbital mechanics around asteroid
    function computeOrbitPeriod(sma) {
        return 2 * Math.PI * Math.sqrt(sma**3 / asteroid.GM);
    }

    function computeGroundTrack(sma, inc, duration) {
        const orbitalPeriod = computeOrbitPeriod(sma);
        const numOrbits = duration / orbitalPeriod;
        const asteroidRotations = duration / asteroid.rotationPeriod;

        // Ground track repeat ratio
        const repeatRatio = numOrbits / asteroidRotations;

        // Latitude coverage from inclination
        const maxLat = Math.min(inc, Math.PI - inc) * 180 / Math.PI;

        return { numOrbits, asteroidRotations, repeatRatio, maxLat };
    }

    // Coverage estimation (simplified)
    function estimateCoverage(sma, ecc, inc) {
        // Swath width approximation
        const altitude = sma * (1 - ecc) - asteroid.radius;  // Periapsis altitude
        const fov = 30 * Math.PI / 180;  // 30 deg FOV
        const swathWidth = 2 * altitude * Math.tan(fov / 2);

        // Ground track analysis
        const gt = computeGroundTrack(sma, inc, mission.observationDuration);

        // Simplified coverage model
        // Coverage increases with number of ground tracks and inclination
        const lonCoverage = Math.min(1, gt.repeatRatio * swathWidth / (2 * Math.PI * asteroid.radius));
        const latCoverage = Math.min(1, Math.sin(inc) * 2);  // Higher inc covers more latitudes

        return lonCoverage * latCoverage;
    }

    // Delta-V for orbit maintenance (simplified)
    function estimateDeltaV(sma, ecc) {
        // Solar radiation pressure perturbation (very simplified)
        const altitude = sma - asteroid.radius;
        const srpAccel = 1e-8 * (1.496e11 / asteroid.a)**2;  // ~1e-8 m/s² at 1 AU

        // Accumulation over mission
        const deltaV = srpAccel * mission.observationDuration;

        // Add margin for eccentricity maintenance
        const eccCorrection = ecc * Math.sqrt(asteroid.GM / sma) * 0.1;

        return deltaV + eccCorrection;
    }

    // Resolution at periapsis
    function computeResolution(sma, ecc) {
        const periapsis = sma * (1 - ecc);
        const altitude = periapsis - asteroid.radius;
        // Assuming 1 mrad/pixel camera
        return altitude * 1e-3;  // m/pixel
    }

    // Fitness function (multi-objective)
    function fitness(x) {
        const sma = x[0];
        const ecc = x[1];
        const inc = x[2];

        // Objective 1: Maximize coverage (minimize negative coverage)
        const coverage = estimateCoverage(sma, ecc, inc);
        const f1 = -coverage;

        // Objective 2: Minimize delta-V
        const deltaV = estimateDeltaV(sma, ecc);
        const f2 = deltaV;

        // Objective 3: Maximize resolution (minimize periapsis altitude)
        const resolution = computeResolution(sma, ecc);
        const f3 = resolution;

        // Constraints
        const periapsis = sma * (1 - ecc);
        const apoapsis = sma * (1 + ecc);

        // Altitude constraints
        const c1 = Math.max(0, (asteroid.radius + mission.minAltitude) - periapsis);
        const c2 = Math.max(0, apoapsis - (asteroid.radius + mission.maxAltitude));

        // Stability constraint (Hill sphere approximation)
        const hillRadius = asteroid.a * Math.pow(asteroid.GM / (3 * 1.989e30), 1/3);
        const c3 = Math.max(0, apoapsis - hillRadius * 0.5);

        return {
            objectives: [f1, f2, f3],
            constraints: [c1, c2, c3],
            coverage, deltaV, resolution
        };
    }

    // Simple genetic algorithm optimization
    console.log('\n--- Running Orbit Optimization ---\n');

    const populationSize = 50;
    const generations = 100;

    // Initialize population
    function randomIndividual() {
        return [
            bounds.sma[0] + Math.random() * (bounds.sma[1] - bounds.sma[0]),
            bounds.ecc[0] + Math.random() * (bounds.ecc[1] - bounds.ecc[0]),
            bounds.inc[0] + Math.random() * (bounds.inc[1] - bounds.inc[0])
        ];
    }

    let population = [];
    for (let i = 0; i < populationSize; i++) {
        const x = randomIndividual();
        const f = fitness(x);
        population.push({ x, f });
    }

    // Evolution
    for (let gen = 0; gen < generations; gen++) {
        // Sort by fitness (using weighted sum for simplicity)
        population.sort((a, b) => {
            const penaltyA = a.f.constraints.reduce((s, c) => s + c * 1000, 0);
            const penaltyB = b.f.constraints.reduce((s, c) => s + c * 1000, 0);
            const scoreA = a.f.objectives[0] + a.f.objectives[1] / 100 + a.f.objectives[2] / 10 + penaltyA;
            const scoreB = b.f.objectives[0] + b.f.objectives[1] / 100 + b.f.objectives[2] / 10 + penaltyB;
            return scoreA - scoreB;
        });

        // Keep best half
        const survivors = population.slice(0, populationSize / 2);

        // Generate offspring
        const offspring = [];
        while (offspring.length < populationSize / 2) {
            // Select parents
            const p1 = survivors[Math.floor(Math.random() * survivors.length)];
            const p2 = survivors[Math.floor(Math.random() * survivors.length)];

            // Crossover
            const child = [
                (p1.x[0] + p2.x[0]) / 2 + (Math.random() - 0.5) * 1000,
                (p1.x[1] + p2.x[1]) / 2 + (Math.random() - 0.5) * 0.05,
                (p1.x[2] + p2.x[2]) / 2 + (Math.random() - 0.5) * 0.1
            ];

            // Clamp to bounds
            child[0] = Math.max(bounds.sma[0], Math.min(bounds.sma[1], child[0]));
            child[1] = Math.max(bounds.ecc[0], Math.min(bounds.ecc[1], child[1]));
            child[2] = Math.max(bounds.inc[0], Math.min(bounds.inc[1], child[2]));

            const f = fitness(child);
            offspring.push({ x: child, f });
        }

        population = [...survivors, ...offspring];

        // Progress report
        if (gen % 20 === 0 || gen === generations - 1) {
            const best = population[0];
            console.log(`Gen ${gen}: Coverage=${(best.f.coverage * 100).toFixed(1)}%, ΔV=${best.f.deltaV.toFixed(2)} m/s, Res=${best.f.resolution.toFixed(2)} m/px`);
        }
    }

    // Best solution
    const best = population[0];
    const sma = best.x[0];
    const ecc = best.x[1];
    const inc = best.x[2];

    console.log('\n=== Optimal Orbit ===');
    console.log(`Semi-major axis: ${(sma / 1000).toFixed(2)} km`);
    console.log(`Eccentricity: ${ecc.toFixed(4)}`);
    console.log(`Inclination: ${(inc * 180/Math.PI).toFixed(1)}°`);

    const periapsis = sma * (1 - ecc);
    const apoapsis = sma * (1 + ecc);
    const period = computeOrbitPeriod(sma);

    console.log(`\nDerived Parameters:`);
    console.log(`  Periapsis altitude: ${((periapsis - asteroid.radius) / 1000).toFixed(2)} km`);
    console.log(`  Apoapsis altitude: ${((apoapsis - asteroid.radius) / 1000).toFixed(2)} km`);
    console.log(`  Orbital period: ${(period / 3600).toFixed(2)} hours`);

    const gt = computeGroundTrack(sma, inc, mission.observationDuration);
    console.log(`  Orbits during mission: ${gt.numOrbits.toFixed(1)}`);
    console.log(`  Asteroid rotations: ${gt.asteroidRotations.toFixed(1)}`);
    console.log(`  Ground track repeat ratio: ${gt.repeatRatio.toFixed(3)}`);

    console.log(`\nMission Performance:`);
    console.log(`  Estimated coverage: ${(best.f.coverage * 100).toFixed(1)}%`);
    console.log(`  Delta-V budget: ${best.f.deltaV.toFixed(2)} m/s`);
    console.log(`  Best resolution: ${best.f.resolution.toFixed(2)} m/pixel`);

    // Compare with other orbit types
    console.log('\n--- Comparison with Standard Orbits ---');

    const orbitTypes = [
        { name: 'Circular polar', sma: 20000, ecc: 0, inc: Math.PI/2 },
        { name: 'Elliptical equatorial', sma: 25000, ecc: 0.2, inc: 0.1 },
        { name: 'Sun-synchronous-like', sma: 15000, ecc: 0.1, inc: 1.4 }
    ];

    console.log('\nOrbit Type          | Coverage | ΔV [m/s] | Res [m/px]');
    console.log('--------------------|----------|----------|------------');

    for (const orbit of orbitTypes) {
        const f = fitness([orbit.sma, orbit.ecc, orbit.inc]);
        console.log(`${orbit.name.padEnd(19)} | ${(f.coverage * 100).toFixed(1).padStart(7)}% | ${f.deltaV.toFixed(2).padStart(8)} | ${f.resolution.toFixed(2).padStart(10)}`);
    }
    console.log(`${'Optimized'.padEnd(19)} | ${(best.f.coverage * 100).toFixed(1).padStart(7)}% | ${best.f.deltaV.toFixed(2).padStart(8)} | ${best.f.resolution.toFixed(2).padStart(10)}`);

    // Mission design notes
    console.log('\n--- Mission Design Notes ---');
    console.log('Key considerations for asteroid proximity operations:');
    console.log('  1. Irregular gravity field (non-spherical shape)');
    console.log('  2. Solar radiation pressure dominates at low altitudes');
    console.log('  3. Third-body perturbations from Sun');
    console.log('  4. Outgassing effects (for comets)');
    console.log('  5. Navigation challenges (weak gravity, irregular shape)');

    console.log('\n=== Asteroid Orbit Optimization Complete ===');
}

main().catch(console.error);
