/**
 * Example 24: Rotational Dynamics
 *
 * Ported from: examples/tudatpy/propagation/coupled_translational_rotational_dynamics.py
 *
 * This example demonstrates satellite attitude/rotational dynamics
 * including torques and angular momentum.
 *
 * Run with: node 24_rotational_dynamics.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Rotational Dynamics Example ===\n');

    const tudat = await createTudatModule();

    // Spacecraft inertia tensor (kg*m^2) - typical small satellite
    const Ixx = 100.0;  // Roll moment of inertia
    const Iyy = 120.0;  // Pitch moment of inertia
    const Izz = 80.0;   // Yaw moment of inertia (spin axis)

    console.log('Spacecraft Inertia Tensor (principal axes):');
    console.log(`  Ixx (roll):  ${Ixx} kg*m^2`);
    console.log(`  Iyy (pitch): ${Iyy} kg*m^2`);
    console.log(`  Izz (yaw):   ${Izz} kg*m^2`);

    // Initial angular velocity (rad/s)
    const omega0_x = 0.01;   // Small roll rate
    const omega0_y = 0.005;  // Small pitch rate
    const omega0_z = 0.1;    // Main spin around z-axis (6 deg/s = 1 rpm)

    console.log('\nInitial Angular Velocity:');
    console.log(`  omega_x: ${(omega0_x * 180 / Math.PI).toFixed(3)} deg/s`);
    console.log(`  omega_y: ${(omega0_y * 180 / Math.PI).toFixed(3)} deg/s`);
    console.log(`  omega_z: ${(omega0_z * 180 / Math.PI).toFixed(3)} deg/s`);

    // Initial quaternion (identity = aligned with inertial frame)
    let q0 = 1.0, q1 = 0.0, q2 = 0.0, q3 = 0.0;

    // Angular momentum (constant in torque-free motion)
    const Hx = Ixx * omega0_x;
    const Hy = Iyy * omega0_y;
    const Hz = Izz * omega0_z;
    const H = Math.sqrt(Hx*Hx + Hy*Hy + Hz*Hz);

    console.log('\nAngular Momentum:');
    console.log(`  Hx: ${Hx.toFixed(4)} kg*m^2/s`);
    console.log(`  Hy: ${Hy.toFixed(4)} kg*m^2/s`);
    console.log(`  Hz: ${Hz.toFixed(4)} kg*m^2/s`);
    console.log(`  |H|: ${H.toFixed(4)} kg*m^2/s`);

    // Rotational kinetic energy (constant in torque-free motion)
    const T = 0.5 * (Ixx * omega0_x*omega0_x + Iyy * omega0_y*omega0_y + Izz * omega0_z*omega0_z);
    console.log(`\nRotational Kinetic Energy: ${T.toFixed(4)} J`);

    // Nutation frequency for a symmetric top (Izz < Ixx = Iyy)
    // omega_nutation = omega_spin * (Izz - Ixx) / Ixx
    // For asymmetric body, behavior is more complex
    const lambda = (Izz - Ixx) / Ixx * omega0_z;
    console.log(`\nNutation rate estimate: ${(lambda * 180 / Math.PI).toFixed(4)} deg/s`);

    // Propagate rotational dynamics using Euler's equations
    // I * d(omega)/dt = -omega x (I * omega) + torque
    // For torque-free: d(omega)/dt = I^-1 * (-omega x (I * omega))

    const dt = 0.1;  // 100 ms time step
    const totalTime = 60.0;  // 1 minute propagation
    const numSteps = Math.floor(totalTime / dt);

    let omega_x = omega0_x;
    let omega_y = omega0_y;
    let omega_z = omega0_z;

    console.log('\nPropagating torque-free rotational motion...\n');
    console.log('Time (s) | omega_x (deg/s) | omega_y (deg/s) | omega_z (deg/s) | |omega| (deg/s)');
    console.log('---------|-----------------|-----------------|-----------------|----------------');

    for (let step = 0; step <= numSteps; step++) {
        const t = step * dt;

        // Print every 10 seconds
        if (step % 100 === 0) {
            const omega_mag = Math.sqrt(omega_x*omega_x + omega_y*omega_y + omega_z*omega_z);
            console.log(`${t.toFixed(1).padStart(8)} | ${(omega_x * 180 / Math.PI).toFixed(4).padStart(15)} | ${(omega_y * 180 / Math.PI).toFixed(4).padStart(15)} | ${(omega_z * 180 / Math.PI).toFixed(4).padStart(15)} | ${(omega_mag * 180 / Math.PI).toFixed(4).padStart(14)}`);
        }

        // Euler's equations (torque-free)
        // d(omega_x)/dt = (Iyy - Izz) / Ixx * omega_y * omega_z
        // d(omega_y)/dt = (Izz - Ixx) / Iyy * omega_z * omega_x
        // d(omega_z)/dt = (Ixx - Iyy) / Izz * omega_x * omega_y

        const domega_x = (Iyy - Izz) / Ixx * omega_y * omega_z;
        const domega_y = (Izz - Ixx) / Iyy * omega_z * omega_x;
        const domega_z = (Ixx - Iyy) / Izz * omega_x * omega_y;

        // RK4 integration
        const k1_x = domega_x;
        const k1_y = domega_y;
        const k1_z = domega_z;

        const ox2 = omega_x + 0.5 * dt * k1_x;
        const oy2 = omega_y + 0.5 * dt * k1_y;
        const oz2 = omega_z + 0.5 * dt * k1_z;
        const k2_x = (Iyy - Izz) / Ixx * oy2 * oz2;
        const k2_y = (Izz - Ixx) / Iyy * oz2 * ox2;
        const k2_z = (Ixx - Iyy) / Izz * ox2 * oy2;

        const ox3 = omega_x + 0.5 * dt * k2_x;
        const oy3 = omega_y + 0.5 * dt * k2_y;
        const oz3 = omega_z + 0.5 * dt * k2_z;
        const k3_x = (Iyy - Izz) / Ixx * oy3 * oz3;
        const k3_y = (Izz - Ixx) / Iyy * oz3 * ox3;
        const k3_z = (Ixx - Iyy) / Izz * ox3 * oy3;

        const ox4 = omega_x + dt * k3_x;
        const oy4 = omega_y + dt * k3_y;
        const oz4 = omega_z + dt * k3_z;
        const k4_x = (Iyy - Izz) / Ixx * oy4 * oz4;
        const k4_y = (Izz - Ixx) / Iyy * oz4 * ox4;
        const k4_z = (Ixx - Iyy) / Izz * ox4 * oy4;

        omega_x += dt / 6 * (k1_x + 2*k2_x + 2*k3_x + k4_x);
        omega_y += dt / 6 * (k1_y + 2*k2_y + 2*k3_y + k4_y);
        omega_z += dt / 6 * (k1_z + 2*k2_z + 2*k3_z + k4_z);

        // Also integrate quaternion for attitude
        // dq/dt = 0.5 * omega_quat * q
        const dq0 = 0.5 * (-q1*omega_x - q2*omega_y - q3*omega_z);
        const dq1 = 0.5 * (q0*omega_x + q2*omega_z - q3*omega_y);
        const dq2 = 0.5 * (q0*omega_y + q3*omega_x - q1*omega_z);
        const dq3 = 0.5 * (q0*omega_z + q1*omega_y - q2*omega_x);

        q0 += dq0 * dt;
        q1 += dq1 * dt;
        q2 += dq2 * dt;
        q3 += dq3 * dt;

        // Normalize quaternion
        const qnorm = Math.sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
        q0 /= qnorm;
        q1 /= qnorm;
        q2 /= qnorm;
        q3 /= qnorm;
    }

    // Verify conservation laws
    const Hx_final = Ixx * omega_x;
    const Hy_final = Iyy * omega_y;
    const Hz_final = Izz * omega_z;
    const H_final = Math.sqrt(Hx_final*Hx_final + Hy_final*Hy_final + Hz_final*Hz_final);
    const T_final = 0.5 * (Ixx * omega_x*omega_x + Iyy * omega_y*omega_y + Izz * omega_z*omega_z);

    console.log('\nConservation Check:');
    console.log(`  Angular momentum: ${H.toFixed(6)} -> ${H_final.toFixed(6)} (error: ${((H_final - H) / H * 100).toFixed(4)}%)`);
    console.log(`  Kinetic energy:   ${T.toFixed(6)} -> ${T_final.toFixed(6)} (error: ${((T_final - T) / T * 100).toFixed(4)}%)`);

    // Final attitude (quaternion to Euler angles)
    const roll = Math.atan2(2*(q0*q1 + q2*q3), 1 - 2*(q1*q1 + q2*q2));
    const pitch = Math.asin(2*(q0*q2 - q3*q1));
    const yaw = Math.atan2(2*(q0*q3 + q1*q2), 1 - 2*(q2*q2 + q3*q3));

    console.log('\nFinal Attitude (Euler angles):');
    console.log(`  Roll:  ${(roll * 180 / Math.PI).toFixed(2)} deg`);
    console.log(`  Pitch: ${(pitch * 180 / Math.PI).toFixed(2)} deg`);
    console.log(`  Yaw:   ${(yaw * 180 / Math.PI).toFixed(2)} deg`);

    console.log('\n=== Rotational dynamics complete ===');
}

main().catch(console.error);
