/**
 * Estimation Dynamical Models Example
 * Port of: tudatpy/examples/estimation/estimation_dynamical_models.py
 *
 * Demonstrates the use of different dynamical models for observation
 * simulation vs. estimation, using the Mars Express spacecraft.
 * Shows how model discrepancies affect estimation accuracy.
 */

export function showEstimationDynamicalModelsExample(chartContainer, log, params = {}) {
    const config = {
        duration: params.duration ?? 10,           // days
        numOrbits: params.numOrbits ?? 15,
        observationInterval: params.observationInterval ?? 60,  // seconds
        noiseLevel: params.noiseLevel ?? 1.0       // m/s for Doppler
    };

    log('Running Estimation with Different Dynamical Models...', 'info');
    log(`Simulation duration: ${config.duration} days`, 'info');
    log(`Observation interval: ${config.observationInterval} seconds`, 'info');

    const startTime = performance.now();
    const result = simulateMarsExpressEstimation(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Truth model: Mars SH(4,4) + Phobos + Deimos + Sun + Jupiter + Earth`, 'info');
    log(`Estimation model: Mars SH(4,4) + Sun + Jupiter + Earth (NO moons)`, 'info');
    log(`Observations simulated: ${result.observations.length}`, 'info');
    log(`Visible observations: ${result.visibleObs.length}`, 'info');
    log(`Final position error: ${(result.finalError.position / 1000).toFixed(2)} km`, 'info');

    renderEstimationDynamicalFigures(chartContainer, result, config);

    return {
        name: 'Estimation Dynamical Models',
        description: 'Mars Express orbit estimation with model mismatch',
        ...result,
        config
    };
}

function simulateMarsExpressEstimation(config) {
    // Mars system parameters
    const GM_MARS = 4.282837e13;  // Mars GM
    const R_MARS = 3.3895e6;      // Mars radius
    const J2_MARS = 1.96045e-3;   // Mars J2

    // Moon parameters (perturbation source in "truth" model)
    const GM_PHOBOS = 7.087e5;    // Phobos GM
    const GM_DEIMOS = 9.8e4;      // Deimos GM
    const a_PHOBOS = 9.376e6;     // Phobos semi-major axis
    const a_DEIMOS = 2.346e7;     // Deimos semi-major axis
    const T_PHOBOS = 27553;       // Phobos period (s)
    const T_DEIMOS = 109075;      // Deimos period (s)

    // Third body parameters
    const GM_SUN = 1.327e20;
    const r_SUN = 2.28e11;        // Mars-Sun distance (average)
    const GM_JUPITER = 1.267e17;
    const r_JUPITER = 6.3e11;     // Mars-Jupiter average distance

    // MEX initial orbit (highly elliptical)
    const periapsis = 350e3 + R_MARS;  // 350 km altitude
    const apoapsis = 10000e3 + R_MARS; // 10,000 km altitude
    const sma = (periapsis + apoapsis) / 2;
    const ecc = (apoapsis - periapsis) / (apoapsis + periapsis);
    const period = 2 * Math.PI * Math.sqrt(sma ** 3 / GM_MARS);

    // Initial state (periapsis, in orbital plane)
    const v0 = Math.sqrt(GM_MARS * (2 / periapsis - 1 / sma));
    let state_truth = [periapsis, 0, 0, 0, v0, 0.05 * v0];  // Slight out-of-plane

    // Ground station (simulating New Norcia tracking)
    const stationLat = -31.0 * Math.PI / 180;  // New Norcia latitude
    const earthMarsDist = 1.5e11;  // Average Earth-Mars distance (simplified)

    // Simulation parameters
    const dt = config.observationInterval;
    const totalTime = config.duration * 86400;
    const numSteps = Math.floor(totalTime / dt);

    // Acceleration function for truth model (includes moons)
    function accelTruth(pos, t) {
        const r = Math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2);

        // Central body (Mars)
        let ax = -GM_MARS * pos[0] / r**3;
        let ay = -GM_MARS * pos[1] / r**3;
        let az = -GM_MARS * pos[2] / r**3;

        // J2 perturbation
        const j2_factor = 1.5 * J2_MARS * (R_MARS / r) ** 2;
        const z_factor = 5 * (pos[2] / r) ** 2 - 1;
        ax += GM_MARS / r**3 * j2_factor * z_factor * pos[0];
        ay += GM_MARS / r**3 * j2_factor * z_factor * pos[1];
        az += GM_MARS / r**3 * j2_factor * (5 * (pos[2] / r) ** 2 - 3) * pos[2];

        // Phobos (circular orbit)
        const phobos_angle = 2 * Math.PI * t / T_PHOBOS;
        const phobos_pos = [a_PHOBOS * Math.cos(phobos_angle), a_PHOBOS * Math.sin(phobos_angle), 0];
        const dp = [pos[0] - phobos_pos[0], pos[1] - phobos_pos[1], pos[2] - phobos_pos[2]];
        const rp = Math.sqrt(dp[0]**2 + dp[1]**2 + dp[2]**2);
        const rp0 = a_PHOBOS;
        ax -= GM_PHOBOS * (dp[0] / rp**3 + phobos_pos[0] / rp0**3);
        ay -= GM_PHOBOS * (dp[1] / rp**3 + phobos_pos[1] / rp0**3);
        az -= GM_PHOBOS * (dp[2] / rp**3 + phobos_pos[2] / rp0**3);

        // Deimos (circular orbit)
        const deimos_angle = 2 * Math.PI * t / T_DEIMOS;
        const deimos_pos = [a_DEIMOS * Math.cos(deimos_angle), a_DEIMOS * Math.sin(deimos_angle), 0];
        const dd = [pos[0] - deimos_pos[0], pos[1] - deimos_pos[1], pos[2] - deimos_pos[2]];
        const rd = Math.sqrt(dd[0]**2 + dd[1]**2 + dd[2]**2);
        const rd0 = a_DEIMOS;
        ax -= GM_DEIMOS * (dd[0] / rd**3 + deimos_pos[0] / rd0**3);
        ay -= GM_DEIMOS * (dd[1] / rd**3 + deimos_pos[1] / rd0**3);
        az -= GM_DEIMOS * (dd[2] / rd**3 + deimos_pos[2] / rd0**3);

        // Sun (simplified third body)
        const sun_pos = [r_SUN * Math.cos(2e-7 * t), r_SUN * Math.sin(2e-7 * t), 0];
        const ds = [pos[0] - sun_pos[0], pos[1] - sun_pos[1], pos[2]];
        const rs = Math.sqrt(ds[0]**2 + ds[1]**2 + ds[2]**2);
        ax -= GM_SUN * ds[0] / rs**3;
        ay -= GM_SUN * ds[1] / rs**3;
        az -= GM_SUN * ds[2] / rs**3;

        return [ax, ay, az];
    }

    // Acceleration function for estimation model (NO moons)
    function accelEstimation(pos, t) {
        const r = Math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2);

        // Central body (Mars)
        let ax = -GM_MARS * pos[0] / r**3;
        let ay = -GM_MARS * pos[1] / r**3;
        let az = -GM_MARS * pos[2] / r**3;

        // J2 perturbation
        const j2_factor = 1.5 * J2_MARS * (R_MARS / r) ** 2;
        const z_factor = 5 * (pos[2] / r) ** 2 - 1;
        ax += GM_MARS / r**3 * j2_factor * z_factor * pos[0];
        ay += GM_MARS / r**3 * j2_factor * z_factor * pos[1];
        az += GM_MARS / r**3 * j2_factor * (5 * (pos[2] / r) ** 2 - 3) * pos[2];

        // Sun (simplified third body) - same as truth
        const sun_pos = [r_SUN * Math.cos(2e-7 * t), r_SUN * Math.sin(2e-7 * t), 0];
        const ds = [pos[0] - sun_pos[0], pos[1] - sun_pos[1], pos[2]];
        const rs = Math.sqrt(ds[0]**2 + ds[1]**2 + ds[2]**2);
        ax -= GM_SUN * ds[0] / rs**3;
        ay -= GM_SUN * ds[1] / rs**3;
        az -= GM_SUN * ds[2] / rs**3;

        return [ax, ay, az];
    }

    // RK4 integrator
    function rk4Step(state, t, dt, accelFn) {
        const pos = [state[0], state[1], state[2]];
        const vel = [state[3], state[4], state[5]];

        const a1 = accelFn(pos, t);
        const k1 = [vel[0], vel[1], vel[2], a1[0], a1[1], a1[2]];

        const p2 = pos.map((p, i) => p + 0.5 * dt * k1[i]);
        const v2 = vel.map((v, i) => v + 0.5 * dt * k1[i + 3]);
        const a2 = accelFn(p2, t + 0.5 * dt);
        const k2 = [v2[0], v2[1], v2[2], a2[0], a2[1], a2[2]];

        const p3 = pos.map((p, i) => p + 0.5 * dt * k2[i]);
        const v3 = vel.map((v, i) => v + 0.5 * dt * k2[i + 3]);
        const a3 = accelFn(p3, t + 0.5 * dt);
        const k3 = [v3[0], v3[1], v3[2], a3[0], a3[1], a3[2]];

        const p4 = pos.map((p, i) => p + dt * k3[i]);
        const v4 = vel.map((v, i) => v + dt * k3[i + 3]);
        const a4 = accelFn(p4, t + dt);
        const k4 = [v4[0], v4[1], v4[2], a4[0], a4[1], a4[2]];

        return state.map((s, i) => s + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
    }

    // Propagate truth trajectory
    const truthTrajectory = [];
    const observations = [];
    let state = [...state_truth];

    for (let i = 0; i <= numSteps; i++) {
        const t = i * dt;
        const day = t / 86400;
        const r = Math.sqrt(state[0]**2 + state[1]**2 + state[2]**2);
        const altitude = (r - R_MARS) / 1000;  // km

        truthTrajectory.push({
            t: day,
            x: state[0] / 1000,  // km
            y: state[1] / 1000,
            z: state[2] / 1000,
            r: r / 1000,
            altitude: altitude
        });

        // Simulate range-rate observation (with noise)
        const rangeRate = (state[0] * state[3] + state[1] * state[4] + state[2] * state[5]) / r;
        const noise = (Math.random() - 0.5) * 2 * config.noiseLevel;
        observations.push({
            t: day,
            value: rangeRate + noise,
            trueValue: rangeRate
        });

        if (i < numSteps) {
            state = rk4Step(state, t, dt, accelTruth);
        }
    }

    // Propagate estimation trajectory (without moon perturbations)
    const estimatedTrajectory = [];
    state = [...state_truth];  // Same initial condition

    for (let i = 0; i <= numSteps; i++) {
        const t = i * dt;
        const day = t / 86400;
        const r = Math.sqrt(state[0]**2 + state[1]**2 + state[2]**2);

        estimatedTrajectory.push({
            t: day,
            x: state[0] / 1000,
            y: state[1] / 1000,
            z: state[2] / 1000,
            r: r / 1000
        });

        if (i < numSteps) {
            state = rk4Step(state, t, dt, accelEstimation);
        }
    }

    // Compute position difference between truth and estimation
    const positionDifference = [];
    for (let i = 0; i < truthTrajectory.length; i++) {
        const dx = truthTrajectory[i].x - estimatedTrajectory[i].x;
        const dy = truthTrajectory[i].y - estimatedTrajectory[i].y;
        const dz = truthTrajectory[i].z - estimatedTrajectory[i].z;
        positionDifference.push({
            t: truthTrajectory[i].t,
            dx: dx,
            dy: dy,
            dz: dz,
            total: Math.sqrt(dx**2 + dy**2 + dz**2)
        });
    }

    // Compute residuals (difference between observed and predicted)
    const residuals = [];
    for (let i = 0; i < observations.length; i++) {
        const truthR = Math.sqrt(
            truthTrajectory[i].x**2 + truthTrajectory[i].y**2 + truthTrajectory[i].z**2
        ) * 1000;  // back to m
        const estR = Math.sqrt(
            estimatedTrajectory[i].x**2 + estimatedTrajectory[i].y**2 + estimatedTrajectory[i].z**2
        ) * 1000;

        residuals.push({
            t: observations[i].t,
            value: (truthR - estR) / 1000  // km
        });
    }

    // Simulate observation visibility (elevation angle constraint)
    const visibleObs = observations.filter((obs, i) => {
        const alt = truthTrajectory[i].altitude;
        return alt > 500 && alt < 8000;  // Simplified visibility
    });

    // Final error
    const lastIdx = positionDifference.length - 1;
    const finalError = {
        position: positionDifference[lastIdx].total * 1000,  // m
        dx: positionDifference[lastIdx].dx * 1000,
        dy: positionDifference[lastIdx].dy * 1000,
        dz: positionDifference[lastIdx].dz * 1000
    };

    return {
        truthTrajectory,
        estimatedTrajectory,
        positionDifference,
        observations,
        visibleObs,
        residuals,
        finalError,
        marsRadius: R_MARS / 1000
    };
}

function renderEstimationDynamicalFigures(container, result, config) {
    container.innerHTML = '';

    const containerRect = container.getBoundingClientRect();
    const containerWidth = containerRect.width || 600;
    const containerHeight = containerRect.height || 500;

    const wrapper = document.createElement('div');
    wrapper.style.cssText = `
        display: flex;
        flex-direction: column;
        gap: 20px;
        padding: 15px;
        width: 100%;
        height: 100%;
        overflow-y: auto;
        box-sizing: border-box;
    `;
    container.appendChild(wrapper);

    const chartWidth = containerWidth - 50;
    const chartHeight = Math.max(200, (containerHeight - 100) / 3);

    // Figure 1: MEX Trajectory (Y-Z plane)
    createFigure(wrapper, 'Mars Express Trajectory (Y-Z Plane)', chartWidth, chartHeight, (svg, w, h) => {
        renderTrajectoryYZ(svg, result, w, h);
    });

    // Figure 2: Position Error Over Time
    createFigure(wrapper, 'Position Error: Truth vs Estimation Model', chartWidth, chartHeight, (svg, w, h) => {
        renderPositionError(svg, result, w, h);
    });

    // Figure 3: Range Residuals
    createFigure(wrapper, 'Final Residuals (Model Mismatch)', chartWidth, chartHeight, (svg, w, h) => {
        renderResiduals(svg, result, w, h);
    });
}

function createFigure(parent, title, width, height, renderFn) {
    const container = document.createElement('div');
    container.style.cssText = `
        background: var(--bg-card);
        border: 1px solid var(--border-glow);
        padding: 10px;
        flex-shrink: 0;
    `;

    const titleEl = document.createElement('div');
    titleEl.style.cssText = 'font-family: "Orbitron", sans-serif; font-size: 12px; color: var(--cyan); margin-bottom: 8px; text-align: center;';
    titleEl.textContent = title;
    container.appendChild(titleEl);

    const chartDiv = document.createElement('div');
    container.appendChild(chartDiv);

    const svg = d3.select(chartDiv)
        .append('svg')
        .attr('width', width)
        .attr('height', height);

    renderFn(svg, width, height);
    parent.appendChild(container);
}

function renderTrajectoryYZ(svg, result, width, height) {
    const margin = { top: 15, right: 20, bottom: 35, left: 55 };

    const data = result.truthTrajectory;
    const maxExtent = Math.max(
        d3.max(data, d => Math.abs(d.y)),
        d3.max(data, d => Math.abs(d.z)),
        result.marsRadius * 1.2
    );

    const x = d3.scaleLinear()
        .domain([-maxExtent * 1.1, maxExtent * 1.1])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-maxExtent * 1.1, maxExtent * 1.1])
        .range([height - margin.bottom, margin.top]);

    // Draw Mars
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', Math.max(3, (x(result.marsRadius) - x(0))))
        .attr('fill', '#c1440e');

    // Draw truth trajectory
    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 1.5)
        .attr('d', d3.line().x(d => x(d.y)).y(d => y(d.z)));

    // Draw estimated trajectory
    svg.append('path')
        .datum(result.estimatedTrajectory)
        .attr('fill', 'none')
        .attr('stroke', '#ff6666')
        .attr('stroke-width', 1)
        .attr('stroke-dasharray', '4,2')
        .attr('d', d3.line().x(d => x(d.y)).y(d => y(d.z)));

    // Legend
    svg.append('line').attr('x1', width - 130).attr('x2', width - 110).attr('y1', 25).attr('y2', 25)
        .attr('stroke', 'var(--cyan)').attr('stroke-width', 2);
    svg.append('text').attr('x', width - 105).attr('y', 28)
        .attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('Truth');

    svg.append('line').attr('x1', width - 130).attr('x2', width - 110).attr('y1', 40).attr('y2', 40)
        .attr('stroke', '#ff6666').attr('stroke-width', 1).attr('stroke-dasharray', '4,2');
    svg.append('text').attr('x', width - 105).attr('y', 43)
        .attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('Estimated');

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(5))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(5))
        .attr('color', 'var(--text-dim)');

    svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 3)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Y [km]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Z [km]');
}

function renderPositionError(svg, result, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const data = result.positionDifference;

    const x = d3.scaleLinear()
        .domain(d3.extent(data, d => d.t))
        .range([margin.left, width - margin.right]);

    const maxError = d3.max(data, d => Math.max(Math.abs(d.dx), Math.abs(d.dy), Math.abs(d.dz), d.total));

    const y = d3.scaleLinear()
        .domain([-maxError * 0.1, maxError * 1.1])
        .range([height - margin.bottom, margin.top]);

    // Draw error components
    const colors = { dx: '#ff6666', dy: '#66ff66', dz: '#6666ff', total: 'var(--cyan)' };

    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', colors.dx)
        .attr('stroke-width', 1)
        .attr('d', d3.line().x(d => x(d.t)).y(d => y(d.dx)));

    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', colors.dy)
        .attr('stroke-width', 1)
        .attr('d', d3.line().x(d => x(d.t)).y(d => y(d.dy)));

    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', colors.dz)
        .attr('stroke-width', 1)
        .attr('d', d3.line().x(d => x(d.t)).y(d => y(d.dz)));

    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', colors.total)
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.t)).y(d => y(d.total)));

    // Legend
    const legendItems = [
        { label: 'Δx', color: colors.dx },
        { label: 'Δy', color: colors.dy },
        { label: 'Δz', color: colors.dz },
        { label: '||ΔX||', color: colors.total }
    ];
    legendItems.forEach((item, i) => {
        svg.append('line')
            .attr('x1', width - 80)
            .attr('x2', width - 60)
            .attr('y1', 15 + i * 13)
            .attr('y2', 15 + i * 13)
            .attr('stroke', item.color)
            .attr('stroke-width', item.label === '||ΔX||' ? 2 : 1);
        svg.append('text')
            .attr('x', width - 55)
            .attr('y', 18 + i * 13)
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '9px')
            .text(item.label);
    });

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(6))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(4))
        .attr('color', 'var(--text-dim)');

    svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 3)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Time [days]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('ΔX [km]');
}

function renderResiduals(svg, result, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    // Subsample for cleaner display
    const step = Math.max(1, Math.floor(result.residuals.length / 200));
    const data = result.residuals.filter((_, i) => i % step === 0);

    const x = d3.scaleLinear()
        .domain(d3.extent(data, d => d.t))
        .range([margin.left, width - margin.right]);

    const maxRes = d3.max(data, d => Math.abs(d.value));

    const y = d3.scaleLinear()
        .domain([-maxRes * 1.1, maxRes * 1.1])
        .range([height - margin.bottom, margin.top]);

    // Zero line
    svg.append('line')
        .attr('x1', margin.left)
        .attr('x2', width - margin.right)
        .attr('y1', y(0))
        .attr('y2', y(0))
        .attr('stroke', 'var(--text-dim)')
        .attr('stroke-width', 0.5)
        .attr('stroke-dasharray', '3,3');

    // Draw residuals as scatter
    svg.selectAll('.residual')
        .data(data)
        .enter()
        .append('circle')
        .attr('cx', d => x(d.t))
        .attr('cy', d => y(d.value))
        .attr('r', 2)
        .attr('fill', 'var(--cyan)')
        .attr('opacity', 0.6);

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(6))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(4))
        .attr('color', 'var(--text-dim)');

    svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 3)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Time [days]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Residual [km]');
}
