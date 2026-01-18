/**
 * Hohmann Transfer Example
 * Basic orbital mechanics demonstration.
 *
 * Demonstrates classical Hohmann transfer between two circular orbits.
 * Shows delta-V requirements and transfer trajectory.
 */

export function showHohmannTransferExample(chartContainer, log, params = {}) {
    const config = {
        initialAltitude: params.initialAltitude ?? 300,   // km (LEO)
        targetAltitude: params.targetAltitude ?? 35786,   // km (GEO)
        numPoints: params.numPoints ?? 200
    };

    log('Running Hohmann Transfer Example...', 'info');
    log(`Initial altitude: ${config.initialAltitude} km (LEO)`, 'info');
    log(`Target altitude: ${config.targetAltitude} km (GEO)`, 'info');

    const startTime = performance.now();
    const result = computeHohmannTransfer(config);
    const elapsed = performance.now() - startTime;
    log(`Computation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`ΔV₁ (departure): ${result.deltaV1.toFixed(0)} m/s`, 'info');
    log(`ΔV₂ (arrival): ${result.deltaV2.toFixed(0)} m/s`, 'info');
    log(`Total ΔV: ${result.totalDeltaV.toFixed(0)} m/s`, 'info');
    log(`Transfer time: ${(result.transferTime / 3600).toFixed(2)} hours`, 'info');

    renderHohmannFigures(chartContainer, result, config);

    return {
        name: 'Hohmann Transfer',
        description: 'LEO to GEO transfer trajectory',
        ...result,
        config
    };
}

function computeHohmannTransfer(config) {
    const RE = 6378.137;  // Earth radius km
    const mu = 3.986004418e5;  // km³/s²

    const r1 = RE + config.initialAltitude;  // Initial orbit radius
    const r2 = RE + config.targetAltitude;   // Target orbit radius

    // Circular orbit velocities
    const v1_circ = Math.sqrt(mu / r1);
    const v2_circ = Math.sqrt(mu / r2);

    // Transfer orbit parameters
    const a_transfer = (r1 + r2) / 2;
    const e_transfer = (r2 - r1) / (r2 + r1);

    // Velocities at periapsis and apoapsis of transfer orbit
    const v1_transfer = Math.sqrt(mu * (2/r1 - 1/a_transfer));
    const v2_transfer = Math.sqrt(mu * (2/r2 - 1/a_transfer));

    // Delta-V requirements
    const deltaV1 = v1_transfer - v1_circ;  // First burn (prograde)
    const deltaV2 = v2_circ - v2_transfer;  // Second burn (prograde)
    const totalDeltaV = deltaV1 + deltaV2;

    // Transfer time (half the transfer orbit period)
    const transferTime = Math.PI * Math.sqrt(a_transfer * a_transfer * a_transfer / mu);

    // Generate trajectory points
    const initialOrbit = [];
    const transferOrbit = [];
    const targetOrbit = [];

    // Initial circular orbit
    for (let i = 0; i <= config.numPoints; i++) {
        const theta = (i / config.numPoints) * 2 * Math.PI;
        initialOrbit.push({
            x: r1 * Math.cos(theta),
            y: r1 * Math.sin(theta)
        });
    }

    // Transfer orbit (from periapsis to apoapsis)
    for (let i = 0; i <= config.numPoints / 2; i++) {
        const theta = (i / (config.numPoints / 2)) * Math.PI;  // 0 to π
        const r = a_transfer * (1 - e_transfer * e_transfer) / (1 + e_transfer * Math.cos(theta));
        transferOrbit.push({
            x: r * Math.cos(theta),
            y: r * Math.sin(theta),
            theta,
            r
        });
    }

    // Target circular orbit
    for (let i = 0; i <= config.numPoints; i++) {
        const theta = (i / config.numPoints) * 2 * Math.PI;
        targetOrbit.push({
            x: r2 * Math.cos(theta),
            y: r2 * Math.sin(theta)
        });
    }

    // Velocity profile during transfer
    const velocityProfile = [];
    for (let i = 0; i <= config.numPoints / 2; i++) {
        const t = (i / (config.numPoints / 2)) * transferTime;
        const theta = (i / (config.numPoints / 2)) * Math.PI;
        const r = a_transfer * (1 - e_transfer * e_transfer) / (1 + e_transfer * Math.cos(theta));
        const v = Math.sqrt(mu * (2/r - 1/a_transfer));

        velocityProfile.push({
            t: t / 3600,  // hours
            r: r - RE,    // altitude
            v,
            theta: theta * 180 / Math.PI
        });
    }

    return {
        r1, r2,
        v1_circ, v2_circ,
        a_transfer, e_transfer,
        v1_transfer, v2_transfer,
        deltaV1: deltaV1 * 1000,  // Convert to m/s
        deltaV2: deltaV2 * 1000,
        totalDeltaV: totalDeltaV * 1000,
        transferTime,
        initialOrbit,
        transferOrbit,
        targetOrbit,
        velocityProfile
    };
}

function renderHohmannFigures(container, result, config) {
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
    const chartHeight = Math.max(250, (containerHeight - 100) / 2);

    // Figure 1: Orbit view
    createFigure(wrapper, 'Hohmann Transfer Trajectory', chartWidth, chartHeight, (svg, w, h) => {
        renderOrbitView(svg, result, w, h);
    });

    // Figure 2: Velocity profile
    createFigure(wrapper, 'Velocity and Altitude During Transfer', chartWidth, chartHeight * 0.7, (svg, w, h) => {
        renderDualAxisChart(svg, result.velocityProfile, w, h);
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

function renderOrbitView(svg, result, width, height) {
    const margin = { top: 20, right: 20, bottom: 40, left: 50 };
    const RE = 6378.137;

    const maxR = result.r2 * 1.1;

    const x = d3.scaleLinear()
        .domain([-maxR, maxR])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-maxR, maxR])
        .range([height - margin.bottom, margin.top]);

    // Draw Earth
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', Math.abs(x(RE) - x(0)))
        .attr('fill', 'url(#earthGrad)')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 1);

    // Earth gradient
    const earthGrad = svg.append('defs')
        .append('radialGradient')
        .attr('id', 'earthGrad')
        .attr('cx', '35%')
        .attr('cy', '35%');
    earthGrad.append('stop').attr('offset', '0%').attr('stop-color', '#4a9fff');
    earthGrad.append('stop').attr('offset', '100%').attr('stop-color', '#0066cc');

    // Draw initial orbit (dashed)
    svg.append('path')
        .datum(result.initialOrbit)
        .attr('fill', 'none')
        .attr('stroke', 'var(--green)')
        .attr('stroke-width', 1.5)
        .attr('stroke-dasharray', '4,4')
        .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

    // Draw target orbit (dashed)
    svg.append('path')
        .datum(result.targetOrbit)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 1.5)
        .attr('stroke-dasharray', '4,4')
        .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

    // Draw transfer orbit (solid, highlighted)
    svg.append('path')
        .datum(result.transferOrbit)
        .attr('fill', 'none')
        .attr('stroke', 'var(--orange)')
        .attr('stroke-width', 3)
        .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

    // Burn markers
    // First burn (periapsis)
    svg.append('circle')
        .attr('cx', x(result.r1))
        .attr('cy', y(0))
        .attr('r', 8)
        .attr('fill', 'var(--red)');

    svg.append('text')
        .attr('x', x(result.r1) + 12)
        .attr('y', y(0) + 4)
        .attr('fill', 'var(--red)')
        .attr('font-size', '10px')
        .text(`ΔV₁: ${(result.deltaV1/1000).toFixed(2)} km/s`);

    // Second burn (apoapsis)
    svg.append('circle')
        .attr('cx', x(-result.r2))
        .attr('cy', y(0))
        .attr('r', 8)
        .attr('fill', 'var(--yellow)');

    svg.append('text')
        .attr('x', x(-result.r2) + 12)
        .attr('y', y(0) + 4)
        .attr('fill', 'var(--yellow)')
        .attr('font-size', '10px')
        .text(`ΔV₂: ${(result.deltaV2/1000).toFixed(2)} km/s`);

    // Legend
    const legend = svg.append('g')
        .attr('transform', `translate(${margin.left + 10}, ${margin.top + 10})`);

    const items = [
        { color: 'var(--green)', label: 'Initial (LEO)', dash: '4,4' },
        { color: 'var(--orange)', label: 'Transfer', dash: null },
        { color: 'var(--cyan)', label: 'Target (GEO)', dash: '4,4' }
    ];

    items.forEach((item, i) => {
        legend.append('line')
            .attr('x1', 0).attr('y1', i * 16)
            .attr('x2', 20).attr('y2', i * 16)
            .attr('stroke', item.color)
            .attr('stroke-width', 2)
            .attr('stroke-dasharray', item.dash);

        legend.append('text')
            .attr('x', 25).attr('y', i * 16 + 4)
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '9px')
            .text(item.label);
    });

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(5).tickFormat(d => `${(d/1000).toFixed(0)}k`))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(5).tickFormat(d => `${(d/1000).toFixed(0)}k`))
        .attr('color', 'var(--text-dim)');

    svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 5)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('X [km]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Y [km]');
}

function renderDualAxisChart(svg, data, width, height) {
    const margin = { top: 10, right: 60, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain(d3.extent(data, d => d.t))
        .range([margin.left, width - margin.right]);

    const yLeft = d3.scaleLinear()
        .domain(d3.extent(data, d => d.r))
        .range([height - margin.bottom, margin.top]);

    const yRight = d3.scaleLinear()
        .domain(d3.extent(data, d => d.v))
        .range([height - margin.bottom, margin.top]);

    // Altitude line
    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.t)).y(d => yLeft(d.r)));

    // Velocity line
    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', 'var(--orange)')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.t)).y(d => yRight(d.v)));

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(5))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(yLeft).ticks(4).tickFormat(d => `${(d/1000).toFixed(0)}k`))
        .attr('color', 'var(--cyan)');

    svg.append('g')
        .attr('transform', `translate(${width - margin.right},0)`)
        .call(d3.axisRight(yRight).ticks(4))
        .attr('color', 'var(--orange)');

    svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 5)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Time [hours]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--cyan)')
        .attr('font-size', '10px')
        .text('Altitude [km]');

    svg.append('text')
        .attr('transform', 'rotate(90)')
        .attr('x', height / 2)
        .attr('y', -width + 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--orange)')
        .attr('font-size', '10px')
        .text('Velocity [km/s]');
}
