/**
 * Thrust Satellite Engine Example
 * Port of: tudatpy/examples/propagation/thrust_satellite_engine.py
 *
 * Demonstrates low-thrust orbit transfer with mass depletion.
 * Shows orbital element evolution during thrust arcs.
 */

export function showThrustSatelliteExample(chartContainer, log, params = {}) {
    const config = {
        initialAltitude: params.initialAltitude ?? 300,  // km
        targetAltitude: params.targetAltitude ?? 800,    // km
        thrust: params.thrust ?? 0.5,                     // N
        isp: params.isp ?? 3000,                          // s
        initialMass: params.initialMass ?? 500,           // kg
        duration: params.duration ?? 30 * 86400,          // 30 days
        numPoints: params.numPoints ?? 1000
    };

    log('Running Thrust Satellite Example...', 'info');
    log(`Initial altitude: ${config.initialAltitude} km`, 'info');
    log(`Target altitude: ${config.targetAltitude} km`, 'info');
    log(`Thrust: ${config.thrust} N, Isp: ${config.isp} s`, 'info');

    const startTime = performance.now();
    const result = simulateThrustTransfer(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Final altitude: ${result.finalAltitude.toFixed(1)} km`, 'info');
    log(`Propellant used: ${result.propellantUsed.toFixed(2)} kg`, 'info');
    log(`Delta-V achieved: ${result.deltaV.toFixed(1)} m/s`, 'info');

    renderThrustFigures(chartContainer, result, config);

    return {
        name: 'Thrust Satellite',
        description: 'Low-thrust orbit transfer simulation',
        ...result,
        config
    };
}

function simulateThrustTransfer(config) {
    const trajectory = [];
    const RE = 6378.137;  // Earth radius km
    const mu = 3.986004418e5;  // km³/s²
    const g0 = 9.80665e-3;  // km/s²

    // Initial orbital elements
    let a = RE + config.initialAltitude;  // km
    let e = 0.001;
    let mass = config.initialMass;

    // Thrust parameters
    const thrust = config.thrust / 1000;  // Convert N to kN
    const mdot = thrust / (config.isp * g0);  // kg/s

    const dt = config.duration / config.numPoints;
    let totalDeltaV = 0;
    let thrusting = true;

    for (let i = 0; i < config.numPoints; i++) {
        const t = i * dt;

        // Current state
        const altitude = a - RE;
        const period = 2 * Math.PI * Math.sqrt(a * a * a / mu);
        const velocity = Math.sqrt(mu / a);

        // Check if we've reached target (stop thrust)
        if (altitude >= config.targetAltitude) {
            thrusting = false;
        }

        // Store trajectory point
        trajectory.push({
            t: t / 86400,  // days
            altitude,
            semiMajorAxis: a,
            eccentricity: e,
            period: period / 60,  // minutes
            velocity,
            mass,
            thrusting: thrusting ? 1 : 0
        });

        if (thrusting && mass > config.initialMass * 0.1) {
            // Apply thrust (tangential, prograde)
            const accel = thrust / mass;  // km/s²

            // Simplified Gauss variational equations for circular orbit
            // da/dt = 2 * a² * v * aT / mu
            const dadr = 2 * a * a * velocity * accel / mu;

            // Update orbital elements
            a += dadr * dt;

            // Update mass
            mass -= mdot * dt;

            // Accumulate delta-V
            totalDeltaV += accel * dt * 1000;  // m/s
        }
    }

    const finalAltitude = trajectory[trajectory.length - 1].altitude;
    const propellantUsed = config.initialMass - trajectory[trajectory.length - 1].mass;

    return { trajectory, finalAltitude, propellantUsed, deltaV: totalDeltaV };
}

function renderThrustFigures(container, result, config) {
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

    // Figure 1: Semi-major axis evolution
    createFigure(wrapper, 'Semi-major Axis Evolution', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.trajectory, 't', 'semiMajorAxis', 'Time [days]', 'Semi-major axis [km]', 'var(--cyan)', w, h);
    });

    // Figure 2: Altitude evolution
    createFigure(wrapper, 'Altitude Evolution', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.trajectory, 't', 'altitude', 'Time [days]', 'Altitude [km]', 'var(--green)', w, h);
    });

    // Figure 3: Mass evolution
    createFigure(wrapper, 'Spacecraft Mass', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.trajectory, 't', 'mass', 'Time [days]', 'Mass [kg]', 'var(--orange)', w, h);
    });

    // Figure 4: Orbital period
    createFigure(wrapper, 'Orbital Period', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.trajectory, 't', 'period', 'Time [days]', 'Period [min]', 'var(--purple)', w, h);
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
