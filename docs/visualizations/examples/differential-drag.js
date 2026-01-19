/**
 * Differential Drag Separation Example
 * Port of: tudatpy/examples/propagation/separation_satellites_diff_drag.py
 *
 * Demonstrates satellite separation using differential drag in LEO.
 */

export function showDifferentialDragExample(chartContainer, log, params = {}) {
    const config = {
        altitude: params.altitude ?? 350,           // km
        duration: params.duration ?? 7,             // days
        numPoints: params.numPoints ?? 1000,
        bc1: params.bc1 ?? 50,                       // Ballistic coefficient 1 [kg/m²]
        bc2: params.bc2 ?? 100                       // Ballistic coefficient 2 [kg/m²]
    };

    log('Running Differential Drag Example...', 'info');
    log(`Altitude: ${config.altitude} km`, 'info');
    log(`Ballistic coefficients: ${config.bc1} and ${config.bc2} kg/m²`, 'info');
    log(`Duration: ${config.duration} days`, 'info');

    const startTime = performance.now();
    const result = simulateDifferentialDrag(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Final along-track separation: ${result.finalSeparation.toFixed(1)} km`, 'info');
    log(`Altitude difference: ${(result.altitudeDiff * 1000).toFixed(0)} m`, 'info');

    renderDiffDragFigures(chartContainer, result, config);

    return {
        name: 'Differential Drag',
        description: 'Satellite separation using differential drag',
        ...result,
        config
    };
}

function simulateDifferentialDrag(config) {
    const RE = 6378.137e3;
    const mu = 3.986004418e14;
    const J2 = 1.08263e-3;
    const rho0 = 1.225;           // Sea level density [kg/m³]
    const H = 8500;               // Scale height [m]

    const r0 = RE + config.altitude * 1000;
    const v0 = Math.sqrt(mu / r0);
    const period = 2 * Math.PI * Math.sqrt(r0 ** 3 / mu);
    const n = 2 * Math.PI / period;

    // Initial states (both satellites start together)
    const inclination = 51.6 * Math.PI / 180;  // ISS-like

    // Satellite states: [a, e, i, omega, RAAN, M]
    let sat1 = { a: r0, e: 0.0001, i: inclination, omega: 0, RAAN: 0, M: 0 };
    let sat2 = { a: r0, e: 0.0001, i: inclination, omega: 0, RAAN: 0, M: 0 };

    const dt = config.duration * 86400 / config.numPoints;
    const history = [];

    // Drag coefficients
    const Cd = 2.2;
    const area1 = 1 / config.bc1;  // Effective area/mass ratio
    const area2 = 1 / config.bc2;

    for (let i = 0; i < config.numPoints; i++) {
        const t = i * dt;
        const day = t / 86400;

        // Compute positions
        const r1 = sat1.a * (1 - sat1.e * Math.cos(sat1.M));
        const r2 = sat2.a * (1 - sat2.e * Math.cos(sat2.M));

        const alt1 = (r1 - RE) / 1000;  // km
        const alt2 = (r2 - RE) / 1000;

        // Along-track separation (from mean anomaly difference)
        const dM = sat1.M - sat2.M;
        const separation = dM * sat1.a / 1000;  // km

        history.push({
            t: day,
            alt1,
            alt2,
            separation: Math.abs(separation),
            dAlt: (r1 - r2) / 1000,
            period1: 2 * Math.PI * Math.sqrt(sat1.a ** 3 / mu) / 60,
            period2: 2 * Math.PI * Math.sqrt(sat2.a ** 3 / mu) / 60
        });

        // Atmospheric density at current altitudes
        const rho1 = rho0 * Math.exp(-(r1 - RE) / H);
        const rho2 = rho0 * Math.exp(-(r2 - RE) / H);

        // Drag acceleration
        const v1 = Math.sqrt(mu / r1);
        const v2 = Math.sqrt(mu / r2);

        const aDrag1 = 0.5 * rho1 * v1 ** 2 * Cd * area1;
        const aDrag2 = 0.5 * rho2 * v2 ** 2 * Cd * area2;

        // Semi-major axis decay (secular effect)
        // da/dt = -2 * a² * ρ * v * Cd * A/m / sqrt(μ/a)
        const dadtFactor1 = -2 * sat1.a ** 2 * aDrag1 / Math.sqrt(mu * sat1.a);
        const dadtFactor2 = -2 * sat2.a ** 2 * aDrag2 / Math.sqrt(mu * sat2.a);

        // Update semi-major axes
        sat1.a += dadtFactor1 * dt;
        sat2.a += dadtFactor2 * dt;

        // Update mean anomalies
        const n1 = Math.sqrt(mu / sat1.a ** 3);
        const n2 = Math.sqrt(mu / sat2.a ** 3);
        sat1.M += n1 * dt;
        sat2.M += n2 * dt;

        // Keep M in [0, 2π]
        sat1.M = sat1.M % (2 * Math.PI);
        sat2.M = sat2.M % (2 * Math.PI);
    }

    const finalPoint = history[history.length - 1];

    return {
        history,
        finalSeparation: finalPoint.separation,
        altitudeDiff: finalPoint.dAlt,
        initialAltitude: config.altitude,
        finalAlt1: finalPoint.alt1,
        finalAlt2: finalPoint.alt2
    };
}

function renderDiffDragFigures(container, result, config) {
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

    // Figure 1: Along-track separation
    createFigure(wrapper, 'Along-track Separation', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.history, 't', 'separation',
            'Time [days]', 'Separation [km]', 'var(--cyan)', w, h);
    });

    // Figure 2: Altitudes
    createFigure(wrapper, 'Satellite Altitudes', chartWidth, chartHeight, (svg, w, h) => {
        renderMultiLine(svg, result.history, 't',
            ['alt1', 'alt2'],
            ['var(--green)', 'var(--orange)'],
            [`BC=${config.bc1}`, `BC=${config.bc2}`],
            'Time [days]', 'Altitude [km]', w, h);
    });

    // Figure 3: Altitude difference
    createFigure(wrapper, 'Altitude Difference', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.history, 't', 'dAlt',
            'Time [days]', 'ΔAlt [km]', 'var(--purple)', w, h);
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
        .domain([0, yExtent[1] + yPadding])
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

function renderMultiLine(svg, data, xKey, yKeys, colors, labels, xLabel, yLabel, width, height) {
    const margin = { top: 10, right: 100, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain(d3.extent(data, d => d[xKey]))
        .range([margin.left, width - margin.right]);

    let yMin = Infinity, yMax = -Infinity;
    yKeys.forEach(key => {
        data.forEach(d => {
            yMin = Math.min(yMin, d[key]);
            yMax = Math.max(yMax, d[key]);
        });
    });
    const yPadding = (yMax - yMin) * 0.1 || 1;

    const y = d3.scaleLinear()
        .domain([yMin - yPadding, yMax + yPadding])
        .range([height - margin.bottom, margin.top]);

    yKeys.forEach((key, i) => {
        svg.append('path')
            .datum(data)
            .attr('fill', 'none')
            .attr('stroke', colors[i])
            .attr('stroke-width', 1.5)
            .attr('d', d3.line().x(d => x(d[xKey])).y(d => y(d[key])));
    });

    const legend = svg.append('g')
        .attr('transform', `translate(${width - margin.right + 10}, ${margin.top})`);

    labels.forEach((label, i) => {
        legend.append('line')
            .attr('x1', 0).attr('y1', i * 16)
            .attr('x2', 18).attr('y2', i * 16)
            .attr('stroke', colors[i])
            .attr('stroke-width', 2);

        legend.append('text')
            .attr('x', 22).attr('y', i * 16 + 4)
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '9px')
            .text(label);
    });

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
