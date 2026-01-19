/**
 * Coupled Translational-Rotational Dynamics Example
 * Port of: tudatpy/examples/propagation/coupled_translational_rotational_dynamics.py
 *
 * Demonstrates coupled orbit and attitude propagation with gravity gradient torque.
 */

export function showCoupledDynamicsExample(chartContainer, log, params = {}) {
    const config = {
        altitude: params.altitude ?? 400,           // km
        duration: params.duration ?? 5400,          // 1.5 orbits
        numPoints: params.numPoints ?? 300,
        initialSpin: params.initialSpin ?? 0.01     // rad/s
    };

    log('Running Coupled Dynamics Example...', 'info');
    log(`Altitude: ${config.altitude} km`, 'info');
    log(`Duration: ${(config.duration / 60).toFixed(0)} minutes`, 'info');
    log(`Initial spin rate: ${(config.initialSpin * 180 / Math.PI).toFixed(2)} deg/s`, 'info');

    const startTime = performance.now();
    const result = propagateCoupledDynamics(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Max gravity gradient torque: ${(result.maxTorque * 1e6).toFixed(3)} µN·m`, 'info');
    log(`Attitude oscillation amplitude: ${result.oscillationAmplitude.toFixed(2)} deg`, 'info');

    renderCoupledDynamicsFigures(chartContainer, result, config);

    return {
        name: 'Coupled Dynamics',
        description: 'Coupled translational-rotational propagation',
        ...result,
        config
    };
}

function propagateCoupledDynamics(config) {
    const RE = 6378.137e3;
    const mu = 3.986004418e14;

    const r0 = RE + config.altitude * 1000;
    const v0 = Math.sqrt(mu / r0);
    const period = 2 * Math.PI * Math.sqrt(r0 ** 3 / mu);
    const n = 2 * Math.PI / period;

    // Spacecraft inertia (elongated body)
    const Ix = 100, Iy = 500, Iz = 500;  // kg·m²

    // Initial conditions
    let theta = 0;        // Pitch angle (rotation about normal)
    let omega = config.initialSpin;  // Angular velocity

    const dt = config.duration / config.numPoints;
    const history = [];
    let maxTorque = 0;

    for (let i = 0; i < config.numPoints; i++) {
        const t = i * dt;

        // Orbital position (circular orbit in xy plane)
        const trueAnomaly = n * t;
        const r = r0;

        // Gravity gradient torque (simplified pitch dynamics)
        // T_gg = (3/2) * (mu/r³) * (Iz - Iy) * sin(2θ)
        const ggTorque = 1.5 * (mu / (r ** 3)) * (Iz - Iy) * Math.sin(2 * theta);
        maxTorque = Math.max(maxTorque, Math.abs(ggTorque));

        // Euler's equation for pitch: Ix * dω/dt = T_gg
        const alpha = ggTorque / Ix;

        history.push({
            t: t / 60,  // minutes
            theta: theta * 180 / Math.PI,  // degrees
            omega: omega * 180 / Math.PI,  // deg/s
            torque: ggTorque * 1e6,        // µN·m
            trueAnomaly: trueAnomaly * 180 / Math.PI,
            orbits: t / period
        });

        // Integrate attitude dynamics (RK4)
        const k1_theta = omega;
        const k1_omega = alpha;

        const theta_mid = theta + 0.5 * dt * k1_theta;
        const omega_mid = omega + 0.5 * dt * k1_omega;
        const torque_mid = 1.5 * (mu / (r ** 3)) * (Iz - Iy) * Math.sin(2 * theta_mid);
        const k2_omega = torque_mid / Ix;

        const k3_omega = k2_omega;  // Simplified

        const theta_end = theta + dt * omega_mid;
        const torque_end = 1.5 * (mu / (r ** 3)) * (Iz - Iy) * Math.sin(2 * theta_end);
        const k4_omega = torque_end / Ix;

        omega += dt / 6 * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega);
        theta += dt * omega;
    }

    // Compute oscillation amplitude
    const thetas = history.map(h => h.theta);
    const oscillationAmplitude = (Math.max(...thetas) - Math.min(...thetas)) / 2;

    return {
        history,
        period: period / 60,
        maxTorque,
        oscillationAmplitude
    };
}

function renderCoupledDynamicsFigures(container, result, config) {
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
    const chartHeight = Math.max(150, (containerHeight - 100) / 3);

    // Figure 1: Pitch angle vs time
    createFigure(wrapper, 'Pitch Angle', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.history, 't', 'theta',
            'Time [min]', 'Pitch [deg]', 'var(--cyan)', w, h);
    });

    // Figure 2: Angular velocity
    createFigure(wrapper, 'Angular Velocity', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.history, 't', 'omega',
            'Time [min]', 'ω [deg/s]', 'var(--green)', w, h);
    });

    // Figure 3: Gravity gradient torque
    createFigure(wrapper, 'Gravity Gradient Torque', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.history, 't', 'torque',
            'Time [min]', 'Torque [µN·m]', 'var(--orange)', w, h);
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

function renderLineChart(svg, data, xKey, yKey, xLabel, yLabel, color, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain(d3.extent(data, d => d[xKey]))
        .range([margin.left, width - margin.right]);

    const yExtent = d3.extent(data, d => d[yKey]);
    const yPadding = (yExtent[1] - yExtent[0]) * 0.1 || 1;

    const y = d3.scaleLinear()
        .domain([yExtent[0] - yPadding, yExtent[1] + yPadding])
        .range([height - margin.bottom, margin.top]);

    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', color)
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d[xKey])).y(d => y(d[yKey])));

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
        .attr('y', height - 5)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text(xLabel);

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text(yLabel);
}
