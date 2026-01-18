/**
 * Re-entry Trajectory Example
 * Port of: tudatpy/examples/propagation/reentry_trajectory.py
 *
 * Demonstrates Space Shuttle re-entry trajectory with aerodynamic guidance.
 * Shows altitude, velocity, flight path angle, and heating rate evolution.
 */

export function showReentryTrajectoryExample(chartContainer, log, params = {}) {
    const config = {
        initialAltitude: params.initialAltitude ?? 120,  // km
        initialVelocity: params.initialVelocity ?? 7500, // m/s
        initialFlightPathAngle: params.initialFlightPathAngle ?? -1.5,  // degrees
        targetAltitude: params.targetAltitude ?? 25,  // km (end condition)
        duration: params.duration ?? 2400,  // 40 minutes max
        numPoints: params.numPoints ?? 500
    };

    log('Running Re-entry Trajectory Example...', 'info');
    log(`Initial altitude: ${config.initialAltitude} km`, 'info');
    log(`Initial velocity: ${config.initialVelocity} m/s`, 'info');
    log(`Flight path angle: ${config.initialFlightPathAngle} deg`, 'info');

    // Simulate re-entry trajectory (simplified model until WASM binding is added)
    const startTime = performance.now();
    const result = simulateReentry(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Final altitude: ${result.trajectory[result.trajectory.length - 1].altitude.toFixed(1)} km`, 'info');
    log(`Max heating rate: ${result.maxHeatingRate.toFixed(1)} W/cm²`, 'info');
    log(`Max deceleration: ${result.maxDeceleration.toFixed(1)} g`, 'info');

    renderReentryFigures(chartContainer, result, config);

    return {
        name: 'Re-entry Trajectory',
        description: 'Space Shuttle atmospheric re-entry simulation',
        ...result,
        config
    };
}

function simulateReentry(config) {
    const trajectory = [];
    const RE = 6378.137;  // Earth radius km
    const g0 = 9.80665;   // Standard gravity m/s²
    const mu = 3.986004418e14;  // m³/s²

    // Simple exponential atmosphere model
    const rho0 = 1.225;  // kg/m³ at sea level
    const H = 8500;  // Scale height meters

    // Spacecraft parameters (Space Shuttle-like)
    const mass = 80000;  // kg
    const Sref = 250;    // Reference area m²
    const CL = 0.8;      // Lift coefficient (typical re-entry)
    const CD = 1.2;      // Drag coefficient

    // Initial state
    let altitude = config.initialAltitude;  // km
    let velocity = config.initialVelocity;  // m/s
    let gamma = config.initialFlightPathAngle * Math.PI / 180;  // radians
    let range = 0;  // km

    const dt = config.duration / config.numPoints;
    let maxHeatingRate = 0;
    let maxDeceleration = 0;

    for (let i = 0; i < config.numPoints; i++) {
        const t = i * dt;
        const h = altitude * 1000;  // Convert to meters

        // Atmospheric density
        const rho = rho0 * Math.exp(-h / H);

        // Dynamic pressure
        const q = 0.5 * rho * velocity * velocity;

        // Aerodynamic forces (per unit mass)
        const drag = q * Sref * CD / mass;
        const lift = q * Sref * CL / mass;

        // Gravity at altitude
        const r = (RE + altitude) * 1000;  // meters
        const g = mu / (r * r);

        // Stagnation point heating (Sutton-Graves approximation)
        const rn = 1.0;  // Nose radius meters
        const heatingRate = 1.83e-4 * Math.sqrt(rho / rn) * Math.pow(velocity, 3) / 10000;  // W/cm²

        // Deceleration in g's
        const deceleration = drag / g0;

        maxHeatingRate = Math.max(maxHeatingRate, heatingRate);
        maxDeceleration = Math.max(maxDeceleration, deceleration);

        // Bank angle profile (simplified)
        const bankAngle = Math.max(0, 70 - (120 - altitude) * 0.7);  // degrees

        // Store trajectory point
        trajectory.push({
            t: t / 60,  // Convert to minutes
            altitude,
            velocity: velocity / 1000,  // km/s
            gamma: gamma * 180 / Math.PI,  // degrees
            mach: velocity / 340,  // Approximate Mach number
            heatingRate,
            deceleration,
            bankAngle,
            range,
            dynamicPressure: q / 1000  // kPa
        });

        // Check termination
        if (altitude <= config.targetAltitude) break;

        // Equations of motion (simplified planar entry)
        const dv = -drag - g * Math.sin(gamma);
        const dgamma = (lift * Math.cos(bankAngle * Math.PI / 180) - g * Math.cos(gamma) + velocity * velocity * Math.cos(gamma) / r) / velocity;
        const dh = velocity * Math.sin(gamma);
        const dr = velocity * Math.cos(gamma) / 1000;  // km/s

        // Euler integration
        velocity += dv * dt;
        gamma += dgamma * dt;
        altitude += dh * dt / 1000;  // km
        range += dr * dt;

        // Prevent numerical issues
        velocity = Math.max(velocity, 100);
        altitude = Math.max(altitude, 0);
    }

    return { trajectory, maxHeatingRate, maxDeceleration };
}

function renderReentryFigures(container, result, config) {
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
    const chartHeight = Math.max(160, (containerHeight - 120) / 4);

    // Figure 1: Altitude vs Time
    createFigure(wrapper, 'Altitude Profile', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.trajectory, 't', 'altitude', 'Time [min]', 'Altitude [km]', 'var(--cyan)', w, h);
    });

    // Figure 2: Velocity vs Time
    createFigure(wrapper, 'Velocity Profile', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.trajectory, 't', 'velocity', 'Time [min]', 'Velocity [km/s]', 'var(--orange)', w, h);
    });

    // Figure 3: Heating Rate vs Time
    createFigure(wrapper, 'Stagnation Point Heating Rate', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.trajectory, 't', 'heatingRate', 'Time [min]', 'Heating [W/cm²]', 'var(--red)', w, h);
    });

    // Figure 4: Deceleration vs Time
    createFigure(wrapper, 'Deceleration Load', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.trajectory, 't', 'deceleration', 'Time [min]', 'Deceleration [g]', 'var(--purple)', w, h);
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
        .domain([Math.max(0, yExtent[0] - yPadding), yExtent[1] + yPadding])
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
