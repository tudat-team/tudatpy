/**
 * Example 31: Coupled Translational-Rotational Dynamics
 *
 * Ported from: examples/tudatpy/propagation/coupled_translational_rotational_dynamics.py
 *
 * This example demonstrates coupled propagation of both translational
 * (orbital) and rotational (attitude) dynamics of a spacecraft.
 *
 * Key concepts:
 * - Coupled state propagation
 * - Attitude dynamics with gravity gradient torque
 * - Quaternion attitude representation
 * - Inertia tensor effects
 * - Spin stabilization
 *
 * Run with: node 31_coupled_dynamics.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Coupled Translational-Rotational Dynamics ===\n');

    const tudat = await createTudatModule();

    // Physical constants
    const GM_EARTH = 3.986004418e14;  // Earth gravitational parameter [m^3/s^2]
    const R_EARTH = 6.378137e6;       // Earth radius [m]

    // Spacecraft parameters
    const spacecraft = {
        mass: 500,  // kg
        // Inertia tensor [kg*m^2] - slightly asymmetric for realistic tumbling
        Ixx: 100,
        Iyy: 150,
        Izz: 200,
        // Products of inertia (assumed zero for principal axes)
        Ixy: 0,
        Ixz: 0,
        Iyz: 0
    };

    console.log('Spacecraft Configuration:');
    console.log(`  Mass: ${spacecraft.mass} kg`);
    console.log(`  Inertia: Ixx=${spacecraft.Ixx}, Iyy=${spacecraft.Iyy}, Izz=${spacecraft.Izz} kg*m^2`);

    // Orbit parameters (LEO)
    const altitude = 400e3;  // 400 km
    const r0 = R_EARTH + altitude;
    const v0 = Math.sqrt(GM_EARTH / r0);  // Circular velocity
    const orbitalPeriod = 2 * Math.PI * Math.sqrt(r0 * r0 * r0 / GM_EARTH);

    console.log(`\nOrbit Configuration:`);
    console.log(`  Altitude: ${altitude / 1000} km`);
    console.log(`  Orbital velocity: ${(v0 / 1000).toFixed(3)} km/s`);
    console.log(`  Orbital period: ${(orbitalPeriod / 60).toFixed(1)} minutes`);

    // Initial attitude (quaternion: [q0, q1, q2, q3] where q0 is scalar)
    // Start with small rotation from body-fixed aligned with orbital frame
    const theta0 = 5 * Math.PI / 180;  // 5 degree initial tilt
    const q0 = [Math.cos(theta0/2), Math.sin(theta0/2), 0, 0];  // Rotation about x-axis

    // Initial angular velocity [rad/s]
    // Include orbital rate for nadir-pointing + small perturbation
    const n = 2 * Math.PI / orbitalPeriod;  // Mean motion
    const omega0 = [0.001, -n + 0.002, 0.001];  // Small wobble around nadir-pointing

    console.log(`\nInitial Attitude:`);
    console.log(`  Quaternion: [${q0.map(x => x.toFixed(4)).join(', ')}]`);
    console.log(`  Angular velocity: [${omega0.map(x => (x * 180 / Math.PI).toFixed(4)).join(', ')}] deg/s`);

    // State vector: [x, y, z, vx, vy, vz, q0, q1, q2, q3, wx, wy, wz]
    // Position in orbital plane, velocity perpendicular
    let state = [
        r0, 0, 0,           // Position
        0, v0, 0,           // Velocity
        q0[0], q0[1], q0[2], q0[3],  // Quaternion
        omega0[0], omega0[1], omega0[2]  // Angular velocity
    ];

    // Simulation parameters
    const dt = 1.0;  // Time step [s]
    const numOrbits = 2;
    const totalTime = numOrbits * orbitalPeriod;

    // Compute gravity gradient torque
    function gravityGradientTorque(pos, quat, I) {
        const r = Math.sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        const r3 = r * r * r;

        // Unit vector from spacecraft to Earth center (nadir) in inertial frame
        const nadir_i = [-pos[0]/r, -pos[1]/r, -pos[2]/r];

        // Transform nadir vector to body frame using quaternion
        // q * v * q^(-1) rotation
        const q = quat;
        const nadir_b = quaternionRotate(q, nadir_i);

        // Gravity gradient coefficient
        const k = 3 * GM_EARTH / r3;

        // Gravity gradient torque: T = k * (n x I*n) where n is nadir in body frame
        const In = [
            I.Ixx * nadir_b[0] + I.Ixy * nadir_b[1] + I.Ixz * nadir_b[2],
            I.Ixy * nadir_b[0] + I.Iyy * nadir_b[1] + I.Iyz * nadir_b[2],
            I.Ixz * nadir_b[0] + I.Iyz * nadir_b[1] + I.Izz * nadir_b[2]
        ];

        // Cross product: nadir_b x In
        const torque = [
            k * (nadir_b[1] * In[2] - nadir_b[2] * In[1]),
            k * (nadir_b[2] * In[0] - nadir_b[0] * In[2]),
            k * (nadir_b[0] * In[1] - nadir_b[1] * In[0])
        ];

        return torque;
    }

    // Quaternion rotation of a vector
    function quaternionRotate(q, v) {
        // q = [w, x, y, z], v = [vx, vy, vz]
        // Result = q * [0, v] * q^(-1)
        const w = q[0], x = q[1], y = q[2], z = q[3];

        // Rotation matrix from quaternion (for efficiency)
        const R = [
            [1 - 2*(y*y + z*z), 2*(x*y - w*z), 2*(x*z + w*y)],
            [2*(x*y + w*z), 1 - 2*(x*x + z*z), 2*(y*z - w*x)],
            [2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x*x + y*y)]
        ];

        return [
            R[0][0]*v[0] + R[0][1]*v[1] + R[0][2]*v[2],
            R[1][0]*v[0] + R[1][1]*v[1] + R[1][2]*v[2],
            R[2][0]*v[0] + R[2][1]*v[1] + R[2][2]*v[2]
        ];
    }

    // Normalize quaternion
    function normalizeQuaternion(q) {
        const norm = Math.sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
        return [q[0]/norm, q[1]/norm, q[2]/norm, q[3]/norm];
    }

    // Euler's equations for rotational dynamics
    function eulerEquations(omega, torque, I) {
        // T = I * omega_dot + omega x (I * omega)
        // omega_dot = I^(-1) * (T - omega x (I * omega))

        const Iw = [
            I.Ixx * omega[0],
            I.Iyy * omega[1],
            I.Izz * omega[2]
        ];

        // omega x (I * omega)
        const wxIw = [
            omega[1] * Iw[2] - omega[2] * Iw[1],
            omega[2] * Iw[0] - omega[0] * Iw[2],
            omega[0] * Iw[1] - omega[1] * Iw[0]
        ];

        // omega_dot = I^(-1) * (torque - wxIw)
        return [
            (torque[0] - wxIw[0]) / I.Ixx,
            (torque[1] - wxIw[1]) / I.Iyy,
            (torque[2] - wxIw[2]) / I.Izz
        ];
    }

    // Quaternion derivative
    function quaternionDerivative(q, omega) {
        // q_dot = 0.5 * q * omega_quat
        // where omega_quat = [0, omega]
        const w = q[0], x = q[1], y = q[2], z = q[3];
        const wx = omega[0], wy = omega[1], wz = omega[2];

        return [
            0.5 * (-x*wx - y*wy - z*wz),
            0.5 * (w*wx + y*wz - z*wy),
            0.5 * (w*wy + z*wx - x*wz),
            0.5 * (w*wz + x*wy - y*wx)
        ];
    }

    // Coupled derivatives
    function derivatives(t, s) {
        const pos = [s[0], s[1], s[2]];
        const vel = [s[3], s[4], s[5]];
        const quat = [s[6], s[7], s[8], s[9]];
        const omega = [s[10], s[11], s[12]];

        // Translational dynamics (two-body)
        const r = Math.sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        const r3 = r * r * r;
        const acc = [
            -GM_EARTH * pos[0] / r3,
            -GM_EARTH * pos[1] / r3,
            -GM_EARTH * pos[2] / r3
        ];

        // Rotational dynamics
        const torque = gravityGradientTorque(pos, quat, spacecraft);
        const omegaDot = eulerEquations(omega, torque, spacecraft);
        const quatDot = quaternionDerivative(quat, omega);

        return [
            vel[0], vel[1], vel[2],
            acc[0], acc[1], acc[2],
            quatDot[0], quatDot[1], quatDot[2], quatDot[3],
            omegaDot[0], omegaDot[1], omegaDot[2]
        ];
    }

    // RK4 integrator
    function rk4Step(t, s, h) {
        const k1 = derivatives(t, s);
        const s1 = s.map((v, i) => v + 0.5 * h * k1[i]);
        const k2 = derivatives(t + 0.5 * h, s1);
        const s2 = s.map((v, i) => v + 0.5 * h * k2[i]);
        const k3 = derivatives(t + 0.5 * h, s2);
        const s3 = s.map((v, i) => v + h * k3[i]);
        const k4 = derivatives(t + h, s3);
        const result = s.map((v, i) => v + (h / 6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));

        // Normalize quaternion after integration
        const qNorm = normalizeQuaternion([result[6], result[7], result[8], result[9]]);
        result[6] = qNorm[0];
        result[7] = qNorm[1];
        result[8] = qNorm[2];
        result[9] = qNorm[3];

        return result;
    }

    console.log('\n--- Running Coupled Propagation ---\n');

    // Storage
    const trajectory = [];
    let t = 0;

    while (t <= totalTime) {
        // Store data periodically
        if (Math.floor(t) % 60 === 0 || trajectory.length === 0) {
            const r = Math.sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
            const omegaMag = Math.sqrt(state[10]*state[10] + state[11]*state[11] + state[12]*state[12]);

            // Convert quaternion to Euler angles for display
            const q = [state[6], state[7], state[8], state[9]];
            const euler = quaternionToEuler(q);

            trajectory.push({
                t: t / 60,  // minutes
                altitude: (r - R_EARTH) / 1000,  // km
                roll: euler.roll * 180 / Math.PI,
                pitch: euler.pitch * 180 / Math.PI,
                yaw: euler.yaw * 180 / Math.PI,
                omegaMag: omegaMag * 180 / Math.PI  // deg/s
            });
        }

        // Integrate
        state = rk4Step(t, state, dt);
        t += dt;
    }

    // Convert quaternion to Euler angles (roll, pitch, yaw)
    function quaternionToEuler(q) {
        const w = q[0], x = q[1], y = q[2], z = q[3];

        // Roll (x-axis rotation)
        const sinr_cosp = 2 * (w * x + y * z);
        const cosr_cosp = 1 - 2 * (x * x + y * y);
        const roll = Math.atan2(sinr_cosp, cosr_cosp);

        // Pitch (y-axis rotation)
        const sinp = 2 * (w * y - z * x);
        const pitch = Math.abs(sinp) >= 1 ? Math.sign(sinp) * Math.PI / 2 : Math.asin(sinp);

        // Yaw (z-axis rotation)
        const siny_cosp = 2 * (w * z + x * y);
        const cosy_cosp = 1 - 2 * (y * y + z * z);
        const yaw = Math.atan2(siny_cosp, cosy_cosp);

        return { roll, pitch, yaw };
    }

    console.log(`Propagated ${numOrbits} orbits (${(totalTime / 60).toFixed(1)} minutes)`);
    console.log(`Data points: ${trajectory.length}`);

    // Analyze attitude behavior
    console.log('\n=== Attitude Analysis ===');

    const rolls = trajectory.map(d => d.roll);
    const pitches = trajectory.map(d => d.pitch);
    const yaws = trajectory.map(d => d.yaw);

    console.log('\nRoll angle:');
    console.log(`  Min: ${Math.min(...rolls).toFixed(2)} deg`);
    console.log(`  Max: ${Math.max(...rolls).toFixed(2)} deg`);
    console.log(`  Range: ${(Math.max(...rolls) - Math.min(...rolls)).toFixed(2)} deg`);

    console.log('\nPitch angle:');
    console.log(`  Min: ${Math.min(...pitches).toFixed(2)} deg`);
    console.log(`  Max: ${Math.max(...pitches).toFixed(2)} deg`);
    console.log(`  Range: ${(Math.max(...pitches) - Math.min(...pitches)).toFixed(2)} deg`);

    console.log('\nYaw angle:');
    console.log(`  Min: ${Math.min(...yaws).toFixed(2)} deg`);
    console.log(`  Max: ${Math.max(...yaws).toFixed(2)} deg`);
    console.log(`  Range: ${(Math.max(...yaws) - Math.min(...yaws)).toFixed(2)} deg`);

    // Angular velocity analysis
    const omegaMags = trajectory.map(d => d.omegaMag);
    console.log('\nAngular velocity magnitude:');
    console.log(`  Min: ${Math.min(...omegaMags).toFixed(4)} deg/s`);
    console.log(`  Max: ${Math.max(...omegaMags).toFixed(4)} deg/s`);

    // Gravity gradient stability analysis
    console.log('\n--- Gravity Gradient Stability ---');
    const Ix = spacecraft.Ixx, Iy = spacecraft.Iyy, Iz = spacecraft.Izz;

    // Stability conditions for gravity gradient stabilization
    // Pitch stable if: Iz > Iy and Iz > Ix
    // Roll/yaw stable if: (Iz - Iy)(Iz - Ix) > 0 and Iy > Ix
    const pitchStable = Iz > Iy && Iz > Ix;
    const rollYawStable = (Iz - Iy) * (Iz - Ix) > 0 && Iy > Ix;

    console.log(`Pitch axis (maximum inertia): ${pitchStable ? 'STABLE' : 'UNSTABLE'}`);
    console.log(`Roll/Yaw coupling: ${rollYawStable ? 'STABLE' : 'UNSTABLE'}`);

    // Natural frequencies
    const omegaN = Math.sqrt(3 * GM_EARTH / Math.pow(r0, 3));  // Orbital rate
    const omegaPitch = omegaN * Math.sqrt(3 * (Iz - Iy) / Ix);
    const omegaRoll = omegaN * Math.sqrt((Iz - Iy) * (Iz - Ix) / (Ix * Iy));

    console.log(`\nNatural frequencies:`);
    console.log(`  Orbital rate: ${(omegaN * 180 / Math.PI * 60).toFixed(4)} deg/min`);
    if (!isNaN(omegaPitch) && omegaPitch > 0) {
        console.log(`  Pitch libration: ${(omegaPitch * 180 / Math.PI * 60).toFixed(4)} deg/min`);
        console.log(`  Pitch period: ${(2 * Math.PI / omegaPitch / 60).toFixed(1)} minutes`);
    }

    // Sample trajectory output
    console.log('\n--- Sample Trajectory Points ---');
    console.log('Time[min] | Alt[km] | Roll[deg] | Pitch[deg] | Yaw[deg]');
    console.log('----------|---------|-----------|------------|----------');
    for (let i = 0; i < trajectory.length; i += Math.floor(trajectory.length / 10)) {
        const p = trajectory[i];
        console.log(`${p.t.toFixed(1).padStart(9)} | ${p.altitude.toFixed(1).padStart(7)} | ${p.roll.toFixed(2).padStart(9)} | ${p.pitch.toFixed(2).padStart(10)} | ${p.yaw.toFixed(2).padStart(9)}`);
    }

    console.log('\n=== Coupled Dynamics Simulation Complete ===');
}

main().catch(console.error);
