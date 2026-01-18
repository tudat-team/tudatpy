/**
 * Covariance Propagation Example
 * Port of: tudatpy/examples/estimation/covariance_propagation_example.py
 *
 * Demonstrates how uncertainty in initial state propagates over time.
 * Shows error ellipsoid growth due to orbital dynamics.
 */

export function showCovariancePropagationExample(chartContainer, log, params = {}) {
    const config = {
        initialAltitude: params.initialAltitude ?? 400,  // km
        duration: params.duration ?? 86400,               // 1 day
        numPoints: params.numPoints ?? 200,
        initialPosStd: params.initialPosStd ?? 50,        // meters
        initialVelStd: params.initialVelStd ?? 0.05       // m/s
    };

    log('Running Covariance Propagation Example...', 'info');
    log(`Initial altitude: ${config.initialAltitude} km`, 'info');
    log(`Initial position σ: ${config.initialPosStd} m`, 'info');
    log(`Initial velocity σ: ${config.initialVelStd} m/s`, 'info');

    const startTime = performance.now();
    const result = propagateCovariance(config);
    const elapsed = performance.now() - startTime;
    log(`Analysis completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Final position σ: ${result.finalPosStd.toFixed(0)} m`, 'info');
    log(`Final velocity σ: ${result.finalVelStd.toFixed(3)} m/s`, 'info');
    log(`Position growth factor: ${(result.finalPosStd / config.initialPosStd).toFixed(1)}x`, 'info');

    renderCovarianceFigures(chartContainer, result, config);

    return {
        name: 'Covariance Propagation',
        description: 'Uncertainty propagation analysis',
        ...result,
        config
    };
}

function propagateCovariance(config) {
    const RE = 6378.137e3;  // Earth radius m
    const mu = 3.986004418e14;  // m³/s²

    const r0 = RE + config.initialAltitude * 1000;  // meters
    const v0 = Math.sqrt(mu / r0);
    const period = 2 * Math.PI * Math.sqrt(r0 * r0 * r0 / mu);
    const n = Math.sqrt(mu / (r0 * r0 * r0));  // Mean motion

    // Initial covariance matrix (diagonal, in local orbital frame)
    const sigmaR = config.initialPosStd;
    const sigmaT = config.initialPosStd;
    const sigmaN = config.initialPosStd;
    const sigmaVr = config.initialVelStd;
    const sigmaVt = config.initialVelStd;
    const sigmaVn = config.initialVelStd;

    // Covariance history
    const covHistory = [];
    const dt = config.duration / config.numPoints;

    for (let i = 0; i < config.numPoints; i++) {
        const t = i * dt;
        const nt = n * t;

        // Simplified Clohessy-Wiltshire (Hill) state transition matrix
        // For circular orbit relative motion
        const phi11 = 4 - 3 * Math.cos(nt);
        const phi12 = 0;
        const phi13 = 0;
        const phi14 = Math.sin(nt) / n;
        const phi15 = 2 * (1 - Math.cos(nt)) / n;
        const phi16 = 0;

        const phi21 = 6 * (Math.sin(nt) - nt);
        const phi22 = 1;
        const phi23 = 0;
        const phi24 = 2 * (Math.cos(nt) - 1) / n;
        const phi25 = (4 * Math.sin(nt) - 3 * nt) / n;
        const phi26 = 0;

        const phi31 = 0;
        const phi32 = 0;
        const phi33 = Math.cos(nt);
        const phi34 = 0;
        const phi35 = 0;
        const phi36 = Math.sin(nt) / n;

        // Propagated position variances (simplified diagonal approximation)
        const sigmaR_t = Math.sqrt(
            (phi11 * sigmaR) ** 2 +
            (phi14 * sigmaVr) ** 2 +
            (phi15 * sigmaVt) ** 2
        );

        const sigmaT_t = Math.sqrt(
            (phi21 * sigmaR) ** 2 +
            (phi22 * sigmaT) ** 2 +
            (phi24 * sigmaVr) ** 2 +
            (phi25 * sigmaVt) ** 2
        );

        const sigmaN_t = Math.sqrt(
            (phi33 * sigmaN) ** 2 +
            (phi36 * sigmaVn) ** 2
        );

        // Propagated velocity variances
        const sigmaVr_t = Math.sqrt(
            (3 * n * Math.sin(nt) * sigmaR) ** 2 +
            (Math.cos(nt) * sigmaVr) ** 2 +
            (2 * Math.sin(nt) * sigmaVt) ** 2
        );

        const sigmaVt_t = Math.sqrt(
            (6 * n * (Math.cos(nt) - 1) * sigmaR) ** 2 +
            (-2 * Math.sin(nt) * sigmaVr) ** 2 +
            ((4 * Math.cos(nt) - 3) * sigmaVt) ** 2
        );

        const sigmaVn_t = Math.sqrt(
            (-n * Math.sin(nt) * sigmaN) ** 2 +
            (Math.cos(nt) * sigmaVn) ** 2
        );

        // RSS values
        const sigmaPos = Math.sqrt(sigmaR_t ** 2 + sigmaT_t ** 2 + sigmaN_t ** 2);
        const sigmaVel = Math.sqrt(sigmaVr_t ** 2 + sigmaVt_t ** 2 + sigmaVn_t ** 2);

        covHistory.push({
            t: t / 3600,  // hours
            sigmaR: sigmaR_t,
            sigmaT: sigmaT_t,
            sigmaN: sigmaN_t,
            sigmaVr: sigmaVr_t * 1000,  // mm/s
            sigmaVt: sigmaVt_t * 1000,
            sigmaVn: sigmaVn_t * 1000,
            sigmaPos,
            sigmaVel,
            orbits: t / period
        });
    }

    const finalPoint = covHistory[covHistory.length - 1];

    return {
        covHistory,
        period: period / 3600,
        finalPosStd: finalPoint.sigmaPos,
        finalVelStd: finalPoint.sigmaVel
    };
}

function renderCovarianceFigures(container, result, config) {
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
    const chartHeight = Math.max(180, (containerHeight - 100) / 3);

    // Figure 1: Position uncertainty components
    createFigure(wrapper, 'Position Uncertainty (RTN Components)', chartWidth, chartHeight, (svg, w, h) => {
        renderMultiLine(svg, result.covHistory, 't',
            ['sigmaR', 'sigmaT', 'sigmaN'],
            ['var(--cyan)', 'var(--green)', 'var(--orange)'],
            ['Radial', 'Transverse', 'Normal'],
            'Time [hours]', 'σ [m]', w, h);
    });

    // Figure 2: Total position uncertainty
    createFigure(wrapper, 'Total Position Uncertainty (3D RSS)', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.covHistory, 't', 'sigmaPos',
            'Time [hours]', 'σ_pos [m]', 'var(--purple)', w, h);
    });

    // Figure 3: Velocity uncertainty
    createFigure(wrapper, 'Velocity Uncertainty (RTN Components)', chartWidth, chartHeight, (svg, w, h) => {
        renderMultiLine(svg, result.covHistory, 't',
            ['sigmaVr', 'sigmaVt', 'sigmaVn'],
            ['var(--red)', 'var(--yellow)', 'var(--cyan)'],
            ['Radial', 'Transverse', 'Normal'],
            'Time [hours]', 'σ [mm/s]', w, h);
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
        .domain([0, yMax + yPadding])
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
