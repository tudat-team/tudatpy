/**
 * Example 28: Two-Stage Rocket Ascent on Mars
 *
 * Ported from: examples/tudatpy/propagation/two_stage_rocket_ascent.py
 *
 * This example demonstrates simulation of a two-stage rocket ascent
 * trajectory on Mars, with dual-thrust solid motor modeling and
 * stage separation at apogee.
 *
 * Key concepts:
 * - Two-stage rocket dynamics
 * - Dual-thrust solid propellant modeling (boost phase)
 * - Stage separation at apogee
 * - Mars atmospheric drag
 * - Custom termination conditions
 * - Multi-body propagation (translational + mass)
 *
 * Run with: node 28_two_stage_rocket.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Two-Stage Rocket Ascent on Mars ===\n');

    const tudat = await createTudatModule();

    // Mars physical constants
    const GM_MARS = 4.282837e13;        // Mars gravitational parameter [m^3/s^2]
    const R_MARS = 3.3895e6;            // Mars equatorial radius [m]
    const g0 = 9.80665;                 // Standard gravity for Isp [m/s^2]

    // Mars exponential atmosphere parameters
    const rho0_mars = 0.020;            // Surface density [kg/m^3]
    const H_mars = 11100;               // Scale height [m]

    // Stage 1 parameters (first stage + second stage attached)
    const stage1 = {
        wetMass: 370,           // Total wet mass [kg]
        dryMass: 185,           // Dry mass after propellant burned [kg]
        thrust: 4250,           // Nominal thrust [N]
        boostFactor: 1.75,      // Boost phase multiplier
        boostDuration: 15,      // Boost phase duration [s]
        Isp: 275,               // Specific impulse [s]
        thrustAngle: 40,        // Thrust angle from vertical [deg]
        CD: 0.85,               // Drag coefficient
        CL: 0.40,               // Lift coefficient
        refArea: 0.25           // Reference area [m^2]
    };

    // Stage 2 parameters (second stage only, after separation)
    const stage2 = {
        wetMass: 85,            // Wet mass [kg]
        dryMass: 38.25,         // Dry mass [kg]
        thrust: 2250,           // Thrust [N]
        boostFactor: 1.75,      // Boost phase multiplier
        boostDuration: 15,      // Boost phase duration [s]
        Isp: 280,               // Specific impulse [s]
        thrustAngle: 90,        // Thrust angle from vertical (horizontal) [deg]
        CD: 0.55,               // Drag coefficient
        CL: 0.25,               // Lift coefficient
        refArea: 0.25           // Reference area [m^2]
    };

    console.log('Stage 1 Configuration:');
    console.log(`  Wet mass: ${stage1.wetMass} kg, Dry mass: ${stage1.dryMass} kg`);
    console.log(`  Thrust: ${stage1.thrust} N (boost: ${stage1.thrust * stage1.boostFactor} N)`);
    console.log(`  Isp: ${stage1.Isp} s, Thrust angle: ${stage1.thrustAngle} deg`);

    console.log('\nStage 2 Configuration:');
    console.log(`  Wet mass: ${stage2.wetMass} kg, Dry mass: ${stage2.dryMass} kg`);
    console.log(`  Thrust: ${stage2.thrust} N (boost: ${stage2.thrust * stage2.boostFactor} N)`);
    console.log(`  Isp: ${stage2.Isp} s, Thrust angle: ${stage2.thrustAngle} deg`);

    // Initial conditions (launch from Mars surface)
    const launchAltitude = 500;         // Launch altitude [m]
    const launchLatitude = 18.85;       // Latitude [deg]
    const launchSpeed = 1.0;            // Initial speed [m/s]
    const flightPathAngle = 40;         // Initial flight path angle [deg]

    console.log('\nLaunch Conditions:');
    console.log(`  Altitude: ${launchAltitude} m`);
    console.log(`  Latitude: ${launchLatitude} deg`);
    console.log(`  Initial speed: ${launchSpeed} m/s`);
    console.log(`  Flight path angle: ${flightPathAngle} deg`);

    // Convert launch conditions to Cartesian state
    const r0 = R_MARS + launchAltitude;
    const lat = launchLatitude * Math.PI / 180;
    const fpa = flightPathAngle * Math.PI / 180;

    // Position in body-fixed frame (simplified, ignoring longitude/rotation)
    const x0 = r0 * Math.cos(lat);
    const y0 = 0;
    const z0 = r0 * Math.sin(lat);

    // Velocity direction (radial and tangential components)
    const vr = launchSpeed * Math.sin(fpa);  // Radial velocity
    const vt = launchSpeed * Math.cos(fpa);  // Tangential velocity

    // Velocity in inertial frame (simplified)
    const vx0 = vr * Math.cos(lat) - vt * Math.sin(lat);
    const vy0 = 0;
    const vz0 = vr * Math.sin(lat) + vt * Math.cos(lat);

    // Mars atmospheric density model
    function atmosphericDensity(altitude) {
        if (altitude < 0) altitude = 0;
        return rho0_mars * Math.exp(-altitude / H_mars);
    }

    // Compute thrust magnitude with boost phase
    function getThrustMagnitude(stage, timeSinceBurn) {
        if (timeSinceBurn < stage.boostDuration) {
            return stage.thrust * stage.boostFactor;
        }
        return stage.thrust;
    }

    // Compute mass flow rate
    function getMassFlowRate(stage, timeSinceBurn) {
        const thrust = getThrustMagnitude(stage, timeSinceBurn);
        return thrust / (stage.Isp * g0);
    }

    // State vector: [x, y, z, vx, vy, vz, mass]
    // Simulation parameters
    const dt = 0.1;  // Time step [s] - small for accuracy during high-thrust phase
    const maxTime = 225 * 60;  // Max simulation time [s]

    // Derivatives function for a given stage
    function derivatives(t, state, stage, burnStartTime) {
        const pos = [state[0], state[1], state[2]];
        const vel = [state[3], state[4], state[5]];
        const mass = state[6];

        // Distance from Mars center
        const r = Math.sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        const altitude = r - R_MARS;

        // Gravitational acceleration (point mass)
        const gMag = GM_MARS / (r * r);
        const aGrav = [
            -gMag * pos[0] / r,
            -gMag * pos[1] / r,
            -gMag * pos[2] / r
        ];

        // Velocity magnitude
        const v = Math.sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);

        // Aerodynamic acceleration (drag only, simplified)
        let aDrag = [0, 0, 0];
        if (altitude < 200000 && v > 0.1) {  // Only below 200 km
            const rho = atmosphericDensity(altitude);
            const dragForce = 0.5 * rho * v * v * stage.CD * stage.refArea;
            const dragAccel = dragForce / mass;
            // Drag opposes velocity
            aDrag = [
                -dragAccel * vel[0] / v,
                -dragAccel * vel[1] / v,
                -dragAccel * vel[2] / v
            ];
        }

        // Thrust acceleration
        let aThrust = [0, 0, 0];
        let mdot = 0;

        // Only thrust if propellant remains
        if (mass > stage.dryMass) {
            const timeSinceBurn = t - burnStartTime;
            const thrustMag = getThrustMagnitude(stage, timeSinceBurn);
            mdot = -getMassFlowRate(stage, timeSinceBurn);  // Negative for mass loss

            // Thrust direction: angle from vertical (radial direction)
            const thrustAngleRad = stage.thrustAngle * Math.PI / 180;

            // Vertical (radial outward) direction
            const radial = [pos[0] / r, pos[1] / r, pos[2] / r];

            // Tangential direction (in direction of motion, simplified)
            let tangent;
            if (v > 0.1) {
                tangent = [vel[0] / v, vel[1] / v, vel[2] / v];
            } else {
                // If nearly stationary, use perpendicular to radial
                tangent = [-radial[1], radial[0], 0];
                const tangMag = Math.sqrt(tangent[0]*tangent[0] + tangent[1]*tangent[1]);
                if (tangMag > 0.01) {
                    tangent = [tangent[0] / tangMag, tangent[1] / tangMag, 0];
                }
            }

            // Thrust direction: combination of radial and tangential
            const thrustDir = [
                Math.cos(thrustAngleRad) * radial[0] + Math.sin(thrustAngleRad) * tangent[0],
                Math.cos(thrustAngleRad) * radial[1] + Math.sin(thrustAngleRad) * tangent[1],
                Math.cos(thrustAngleRad) * radial[2] + Math.sin(thrustAngleRad) * tangent[2]
            ];

            const thrustAccel = thrustMag / mass;
            aThrust = [
                thrustAccel * thrustDir[0],
                thrustAccel * thrustDir[1],
                thrustAccel * thrustDir[2]
            ];
        }

        // Total acceleration
        const ax = aGrav[0] + aDrag[0] + aThrust[0];
        const ay = aGrav[1] + aDrag[1] + aThrust[1];
        const az = aGrav[2] + aDrag[2] + aThrust[2];

        return [vel[0], vel[1], vel[2], ax, ay, az, mdot];
    }

    // RK4 integrator
    function rk4Step(t, state, stage, burnStartTime, h) {
        const k1 = derivatives(t, state, stage, burnStartTime);
        const s1 = state.map((v, i) => v + 0.5 * h * k1[i]);
        const k2 = derivatives(t + 0.5 * h, s1, stage, burnStartTime);
        const s2 = state.map((v, i) => v + 0.5 * h * k2[i]);
        const k3 = derivatives(t + 0.5 * h, s2, stage, burnStartTime);
        const s3 = state.map((v, i) => v + h * k3[i]);
        const k4 = derivatives(t + h, s3, stage, burnStartTime);
        return state.map((v, i) => v + (h / 6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
    }

    console.log('\n--- Stage 1 Ascent (to apogee) ---\n');

    // Initialize state for Stage 1
    let state = [x0, y0, z0, vx0, vy0, vz0, stage1.wetMass];
    let t = 0;
    const stage1BurnStart = 0;

    // Track maximum altitude for apogee detection
    let lastAltitude = 0;
    let stage1Data = [];

    // Stage 1 propagation - until apogee or propellant exhausted
    let reachedApogee = false;
    while (t < maxTime && !reachedApogee) {
        const r = Math.sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
        const altitude = r - R_MARS;
        const v = Math.sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);

        // Store data
        if (Math.floor(t) % 1 === 0 || stage1Data.length === 0) {
            stage1Data.push({
                t: t,
                altitude: altitude / 1000,
                velocity: v,
                mass: state[6]
            });
        }

        // Check for apogee (altitude starts decreasing) after initial ascent
        if (t > 1 && altitude < lastAltitude && state[6] >= stage1.dryMass) {
            // Propellant exhausted and descending - apogee passed
            reachedApogee = true;
            console.log(`Stage 1 apogee reached at t = ${t.toFixed(1)} s`);
            console.log(`  Altitude: ${(altitude / 1000).toFixed(2)} km`);
            console.log(`  Velocity: ${v.toFixed(1)} m/s`);
            console.log(`  Mass: ${state[6].toFixed(2)} kg (dry: ${stage1.dryMass} kg)`);
            break;
        }

        // Check for surface impact
        if (altitude < 0) {
            console.log('Stage 1 impacted surface!');
            break;
        }

        lastAltitude = altitude;

        // Integrate
        state = rk4Step(t, state, stage1, stage1BurnStart, dt);
        t += dt;
    }

    // Stage separation
    console.log('\n--- Stage Separation ---');
    const separationTime = t;
    const separationAltitude = (Math.sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]) - R_MARS) / 1000;
    console.log(`Separation at t = ${separationTime.toFixed(1)} s, altitude = ${separationAltitude.toFixed(2)} km`);

    // Stage 2 initial state (same position/velocity, new mass)
    state = [state[0], state[1], state[2], state[3], state[4], state[5], stage2.wetMass];
    const stage2BurnStart = t;

    console.log('\n--- Stage 2 Burn (horizontal) ---\n');

    // Stage 2 propagation
    let stage2Data = [];
    const stage2MaxTime = t + 300 * 60;  // 300 minutes max

    while (t < stage2MaxTime) {
        const r = Math.sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
        const altitude = r - R_MARS;
        const v = Math.sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);

        // Store data periodically
        if (Math.floor(t - separationTime) % 10 === 0 || stage2Data.length === 0) {
            stage2Data.push({
                t: t,
                altitude: altitude / 1000,
                velocity: v,
                mass: state[6]
            });
        }

        // Check for surface impact
        if (altitude < 500) {
            console.log(`Stage 2 descended below 500m at t = ${t.toFixed(1)} s`);
            break;
        }

        // Check if orbit achieved (velocity > orbital velocity at current altitude)
        const orbitalVelocity = Math.sqrt(GM_MARS / r);
        if (v > orbitalVelocity * 0.95 && altitude > 50000) {
            console.log(`Stage 2 achieved near-orbital velocity at t = ${t.toFixed(1)} s`);
            console.log(`  Altitude: ${(altitude / 1000).toFixed(2)} km`);
            console.log(`  Velocity: ${v.toFixed(1)} m/s (orbital: ${orbitalVelocity.toFixed(1)} m/s)`);
            break;
        }

        // Integrate
        state = rk4Step(t, state, stage2, stage2BurnStart, dt);
        t += dt;
    }

    // Final state analysis
    const finalR = Math.sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
    const finalAltitude = finalR - R_MARS;
    const finalV = Math.sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);

    console.log('\n=== Final State ===');
    console.log(`Total flight time: ${(t / 60).toFixed(2)} minutes`);
    console.log(`Final altitude: ${(finalAltitude / 1000).toFixed(2)} km`);
    console.log(`Final velocity: ${(finalV / 1000).toFixed(3)} km/s`);
    console.log(`Final mass: ${state[6].toFixed(2)} kg`);

    // Orbital analysis
    const specificEnergy = 0.5 * finalV * finalV - GM_MARS / finalR;
    const circularVelocity = Math.sqrt(GM_MARS / finalR);

    console.log('\n--- Orbital Analysis ---');
    console.log(`Circular velocity at altitude: ${(circularVelocity / 1000).toFixed(3)} km/s`);
    console.log(`Current velocity: ${(finalV / 1000).toFixed(3)} km/s`);
    console.log(`Velocity ratio: ${(finalV / circularVelocity * 100).toFixed(1)}%`);

    if (specificEnergy < 0) {
        const sma = -GM_MARS / (2 * specificEnergy);

        // Angular momentum for eccentricity
        const h = Math.sqrt(
            Math.pow(state[1]*state[5] - state[2]*state[4], 2) +
            Math.pow(state[2]*state[3] - state[0]*state[5], 2) +
            Math.pow(state[0]*state[4] - state[1]*state[3], 2)
        );
        const e = Math.sqrt(1 + 2 * specificEnergy * h * h / (GM_MARS * GM_MARS));

        const periapsis = sma * (1 - e);
        const apoapsis = sma * (1 + e);

        console.log(`\nOrbit type: Elliptical (captured)`);
        console.log(`Semi-major axis: ${(sma / 1000).toFixed(0)} km`);
        console.log(`Eccentricity: ${e.toFixed(4)}`);
        console.log(`Periapsis altitude: ${((periapsis - R_MARS) / 1000).toFixed(0)} km`);
        console.log(`Apoapsis altitude: ${((apoapsis - R_MARS) / 1000).toFixed(0)} km`);

        if (periapsis < R_MARS) {
            console.log('WARNING: Periapsis below surface - orbit not stable!');
        }
    } else {
        console.log(`\nOrbit type: Hyperbolic (escape trajectory)`);
        const vInf = Math.sqrt(2 * specificEnergy);
        console.log(`Excess velocity: ${(vInf / 1000).toFixed(3)} km/s`);
    }

    // Delta-V analysis
    const stage1DeltaV = stage1.Isp * g0 * Math.log(stage1.wetMass / stage1.dryMass);
    const stage2DeltaV = stage2.Isp * g0 * Math.log(stage2.wetMass / stage2.dryMass);
    const totalDeltaV = stage1DeltaV + stage2DeltaV;

    console.log('\n--- Delta-V Budget ---');
    console.log(`Stage 1 delta-V: ${(stage1DeltaV / 1000).toFixed(3)} km/s`);
    console.log(`Stage 2 delta-V: ${(stage2DeltaV / 1000).toFixed(3)} km/s`);
    console.log(`Total delta-V: ${(totalDeltaV / 1000).toFixed(3)} km/s`);

    // Mars orbital velocity at low orbit
    const marsOrbitV = Math.sqrt(GM_MARS / (R_MARS + 200000));  // 200 km altitude
    console.log(`\nRequired velocity for 200 km circular orbit: ${(marsOrbitV / 1000).toFixed(3)} km/s`);
    console.log(`Theoretical delta-V margin: ${((totalDeltaV - marsOrbitV) / 1000).toFixed(3)} km/s`);

    // Trajectory summary
    console.log('\n--- Trajectory Summary ---');
    console.log(`Stage 1 data points: ${stage1Data.length}`);
    console.log(`Stage 2 data points: ${stage2Data.length}`);

    const maxAlt1 = Math.max(...stage1Data.map(d => d.altitude));
    const maxAlt2 = Math.max(...stage2Data.map(d => d.altitude));
    console.log(`Stage 1 max altitude: ${maxAlt1.toFixed(1)} km`);
    console.log(`Stage 2 max altitude: ${maxAlt2.toFixed(1)} km`);

    console.log('\n=== Two-Stage Rocket Ascent Complete ===');
}

main().catch(console.error);
