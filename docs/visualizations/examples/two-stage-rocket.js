/**
 * Two-Stage Rocket Ascent Example
 * Port of: tudatpy/examples/propagation/two_stage_rocket_ascent.py
 *
 * Demonstrates rocket ascent trajectory with staging.
 * Shows altitude, velocity, acceleration, and dynamic pressure.
 */

export function showTwoStageRocketExample(chartContainer, log, params = {}) {
    const config = {
        // Stage 1 parameters
        stage1Mass: params.stage1Mass ?? 500000,      // kg (dry + propellant)
        stage1Propellant: params.stage1Propellant ?? 400000,  // kg
        stage1Thrust: params.stage1Thrust ?? 7000000, // N
        stage1Isp: params.stage1Isp ?? 280,           // s

        // Stage 2 parameters
        stage2Mass: params.stage2Mass ?? 50000,       // kg (dry + propellant)
        stage2Propellant: params.stage2Propellant ?? 40000,   // kg
        stage2Thrust: params.stage2Thrust ?? 1000000, // N
        stage2Isp: params.stage2Isp ?? 350,           // s

        payloadMass: params.payloadMass ?? 5000,      // kg
        duration: params.duration ?? 600,              // 10 minutes
        numPoints: params.numPoints ?? 1000
    };

    log('Running Two-Stage Rocket Ascent...', 'info');
    log(`Stage 1: ${config.stage1Thrust/1e6} MN thrust, ${config.stage1Isp}s Isp`, 'info');
    log(`Stage 2: ${config.stage2Thrust/1e6} MN thrust, ${config.stage2Isp}s Isp`, 'info');

    const startTime = performance.now();
    const result = simulateRocketAscent(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Stage separation: ${result.stageSeparationTime.toFixed(1)} s at ${result.stageSeparationAltitude.toFixed(1)} km`, 'info');
    log(`Final altitude: ${result.finalAltitude.toFixed(1)} km`, 'info');
    log(`Final velocity: ${result.finalVelocity.toFixed(0)} m/s`, 'info');
    log(`Max Q: ${result.maxQ.toFixed(0)} kPa at ${result.maxQAltitude.toFixed(1)} km`, 'info');

    renderRocketFigures(chartContainer, result, config);

    return {
        name: 'Two-Stage Rocket',
        description: 'Rocket ascent with staging simulation',
        ...result,
        config
    };
}

function simulateRocketAscent(config) {
    const trajectory = [];
    const RE = 6378.137e3;  // Earth radius m
    const g0 = 9.80665;     // m/s²
    const mu = 3.986004418e14;  // m³/s²

    // Atmosphere model
    const rho0 = 1.225;  // kg/m³
    const H = 8500;      // Scale height m

    // Rocket parameters
    const Cd = 0.3;      // Drag coefficient
    const Aref = 10;     // Reference area m²

    // Initial state
    let altitude = 0;     // m
    let velocity = 0;     // m/s
    let gamma = 89 * Math.PI / 180;  // Flight path angle (near vertical)
    let mass = config.stage1Mass + config.stage2Mass + config.payloadMass;
    let downrange = 0;

    // Stage tracking
    let currentStage = 1;
    let propellantRemaining = config.stage1Propellant;
    let thrust = config.stage1Thrust;
    let isp = config.stage1Isp;
    let mdot = thrust / (isp * g0);

    let stageSeparationTime = 0;
    let stageSeparationAltitude = 0;
    let maxQ = 0;
    let maxQAltitude = 0;

    const dt = config.duration / config.numPoints;

    for (let i = 0; i < config.numPoints; i++) {
        const t = i * dt;

        // Atmospheric density
        const rho = altitude < 100000 ? rho0 * Math.exp(-altitude / H) : 0;

        // Dynamic pressure
        const q = 0.5 * rho * velocity * velocity;

        // Track max Q
        if (q > maxQ) {
            maxQ = q;
            maxQAltitude = altitude / 1000;
        }

        // Gravity at altitude
        const r = RE + altitude;
        const g = mu / (r * r);

        // Drag force (per unit mass)
        const drag = q * Aref * Cd / mass;

        // Thrust acceleration
        let thrustAccel = 0;
        if (propellantRemaining > 0) {
            thrustAccel = thrust / mass;
            propellantRemaining -= mdot * dt;
            mass -= mdot * dt;

            // Check for staging
            if (currentStage === 1 && propellantRemaining <= 0) {
                // Stage separation
                stageSeparationTime = t;
                stageSeparationAltitude = altitude / 1000;

                // Drop stage 1
                mass = config.stage2Mass + config.payloadMass;
                propellantRemaining = config.stage2Propellant;
                thrust = config.stage2Thrust;
                isp = config.stage2Isp;
                mdot = thrust / (isp * g0);
                currentStage = 2;
            }
        }

        // Store trajectory point
        trajectory.push({
            t,
            altitude: altitude / 1000,  // km
            velocity,
            acceleration: (thrustAccel - drag - g * Math.sin(gamma)) / g0,  // g's
            gamma: gamma * 180 / Math.PI,
            mass,
            dynamicPressure: q / 1000,  // kPa
            downrange: downrange / 1000,  // km
            stage: currentStage,
            mach: velocity / 340  // Approximate
        });

        // Equations of motion
        const dv = thrustAccel - drag - g * Math.sin(gamma);

        // Gravity turn rate (avoid division by zero when velocity is small)
        let dgamma = 0;
        if (velocity > 10) {
            dgamma = (g / velocity) * (velocity * velocity / (g * r) - 1) * Math.cos(gamma);
        }

        // Gravity turn (pitch over after ~10s)
        if (t > 10 && t < 60 && gamma > 45 * Math.PI / 180) {
            // Gentle pitch over
            gamma -= 0.5 * Math.PI / 180 * dt / 10;
        }

        // Update state
        velocity += dv * dt;
        gamma += dgamma * dt;
        altitude += velocity * Math.sin(gamma) * dt;
        downrange += velocity * Math.cos(gamma) * dt;

        // Clamp values
        velocity = Math.max(velocity, 1);  // Keep minimum velocity to avoid div by zero
        altitude = Math.max(altitude, 0);
        gamma = Math.max(0.01, Math.min(gamma, Math.PI / 2));  // Keep gamma slightly above 0

        // Stop if crashed
        if (altitude < 0 && i > 10) break;
    }

    const finalAltitude = trajectory[trajectory.length - 1].altitude;
    const finalVelocity = trajectory[trajectory.length - 1].velocity;

    return {
        trajectory,
        stageSeparationTime,
        stageSeparationAltitude,
        finalAltitude,
        finalVelocity,
        maxQ: maxQ / 1000,
        maxQAltitude
    };
}

function renderRocketFigures(container, result, config) {
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
    const chartHeight = Math.max(150, (containerHeight - 120) / 4);

    // Figure 1: Altitude vs Time
    createFigure(wrapper, 'Altitude Profile', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChartWithStaging(svg, result.trajectory, 't', 'altitude', 'Time [s]', 'Altitude [km]', 'var(--cyan)', w, h, result.stageSeparationTime);
    });

    // Figure 2: Velocity vs Time
    createFigure(wrapper, 'Velocity Profile', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChartWithStaging(svg, result.trajectory, 't', 'velocity', 'Time [s]', 'Velocity [m/s]', 'var(--green)', w, h, result.stageSeparationTime);
    });

    // Figure 3: Acceleration vs Time
    createFigure(wrapper, 'Acceleration Profile', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChartWithStaging(svg, result.trajectory, 't', 'acceleration', 'Time [s]', 'Acceleration [g]', 'var(--orange)', w, h, result.stageSeparationTime);
    });

    // Figure 4: Dynamic Pressure vs Time
    createFigure(wrapper, 'Dynamic Pressure (Max-Q)', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChartWithStaging(svg, result.trajectory, 't', 'dynamicPressure', 'Time [s]', 'q [kPa]', 'var(--red)', w, h, result.stageSeparationTime);
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

function renderLineChartWithStaging(svg, data, xKey, yKey, xLabel, yLabel, color, width, height, stageSepTime) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain(d3.extent(data, d => d[xKey]))
        .range([margin.left, width - margin.right]);

    const yExtent = d3.extent(data, d => d[yKey]);
    const yPadding = (yExtent[1] - yExtent[0]) * 0.1 || 1;
    const y = d3.scaleLinear()
        .domain([Math.min(0, yExtent[0] - yPadding), yExtent[1] + yPadding])
        .range([height - margin.bottom, margin.top]);

    // Stage separation line
    if (stageSepTime > 0) {
        svg.append('line')
            .attr('x1', x(stageSepTime))
            .attr('x2', x(stageSepTime))
            .attr('y1', margin.top)
            .attr('y2', height - margin.bottom)
            .attr('stroke', 'var(--yellow)')
            .attr('stroke-width', 1)
            .attr('stroke-dasharray', '4,4');

        svg.append('text')
            .attr('x', x(stageSepTime) + 5)
            .attr('y', margin.top + 15)
            .attr('fill', 'var(--yellow)')
            .attr('font-size', '9px')
            .text('Stage Sep');
    }

    // Main line
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
