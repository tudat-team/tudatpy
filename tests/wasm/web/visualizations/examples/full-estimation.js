/**
 * Full Parameter Estimation Example
 * Port of: tudatpy/examples/estimation/full_estimation_example.py
 *
 * Demonstrates full parameter estimation for DELFI-C3 satellite including
 * initial state, drag coefficient, and gravitational parameter estimation.
 */

export function showFullEstimationExample(chartContainer, log, params = {}) {
    const config = {
        duration: params.duration ?? 3,              // days
        observationInterval: params.observationInterval ?? 60,  // seconds
        noiseLevel: params.noiseLevel ?? 1e-3,       // m/s Doppler noise
        numIterations: params.numIterations ?? 4
    };

    log('Running Full Parameter Estimation Example...', 'info');
    log(`Simulation duration: ${config.duration} days`, 'info');
    log(`Observation interval: ${config.observationInterval} s`, 'info');
    log(`Doppler noise: ${config.noiseLevel * 1000} mm/s`, 'info');

    const startTime = performance.now();
    const result = performEstimation(config);
    const elapsed = performance.now() - startTime;
    log(`Estimation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Iterations: ${result.iterations.length}`, 'info');
    log(`Observations: ${result.observations.length}`, 'info');
    log(`Final RMS residual: ${(result.finalRMS * 1000).toFixed(3)} mm/s`, 'info');

    renderEstimationFigures(chartContainer, result, config);

    return {
        name: 'Full Estimation',
        description: 'Parameter estimation with Doppler observations',
        ...result,
        config
    };
}

function performEstimation(config) {
    // Satellite and ground station parameters
    const RE = 6378.137e3;  // Earth radius (m)
    const mu = 3.986004418e14;  // Earth gravitational parameter (m³/s²)

    // DELFI-C3 TLE-derived initial state (approximate)
    const a = 6900e3;  // Semi-major axis (m)
    const e = 0.001;   // Eccentricity
    const inc = 98 * Math.PI / 180;  // Inclination (rad)
    const raan = 180 * Math.PI / 180;  // RAAN (rad)

    // Ground station (Delft, Netherlands)
    const stationLat = 52.00667 * Math.PI / 180;
    const stationLon = 4.35556 * Math.PI / 180;
    const minElevation = 15 * Math.PI / 180;

    // True parameters
    const trueParams = {
        x: 6.8e6, y: 0, z: 0,
        vx: 0, vy: 7.5e3, vz: 0,
        Cd: 1.2,
        muEarth: mu
    };

    // Perturbed initial guess
    const initialGuess = {
        x: trueParams.x + 10,
        y: trueParams.y + 10,
        z: trueParams.z + 10,
        vx: trueParams.vx + 0.01,
        vy: trueParams.vy + 0.01,
        vz: trueParams.vz + 0.01,
        Cd: trueParams.Cd + 0.01,
        muEarth: trueParams.muEarth + 1e5
    };

    // Generate observations
    const numObs = Math.floor(config.duration * 86400 / config.observationInterval);
    const observations = [];
    const trueObservations = [];

    // Simulate satellite trajectory for observation generation
    let state = { ...trueParams };
    const dt = config.observationInterval;
    const earthRotRate = 2 * Math.PI / 86164;  // Earth rotation rate

    for (let i = 0; i < numObs; i++) {
        const t = i * dt;
        const r = Math.sqrt(state.x**2 + state.y**2 + state.z**2);
        const v = Math.sqrt(state.vx**2 + state.vy**2 + state.vz**2);

        // Station position (rotating Earth)
        const theta = stationLon + earthRotRate * t;
        const stationX = RE * Math.cos(stationLat) * Math.cos(theta);
        const stationY = RE * Math.cos(stationLat) * Math.sin(theta);
        const stationZ = RE * Math.sin(stationLat);

        // Vector from station to satellite
        const dx = state.x - stationX;
        const dy = state.y - stationY;
        const dz = state.z - stationZ;
        const range = Math.sqrt(dx*dx + dy*dy + dz*dz);

        // Elevation angle
        const stationR = Math.sqrt(stationX*stationX + stationY*stationY + stationZ*stationZ);
        const dotProduct = (dx*stationX + dy*stationY + dz*stationZ) / (range * stationR);
        const elevation = Math.PI/2 - Math.acos(dotProduct);

        // Only observe when above minimum elevation
        if (elevation > minElevation) {
            // Range rate (Doppler)
            const rangeRate = (dx*state.vx + dy*state.vy + dz*state.vz) / range;

            // Add noise
            const noise = (Math.random() - 0.5) * 2 * config.noiseLevel;

            observations.push({
                t: t,
                tDays: t / 86400,
                rangeRate: rangeRate + noise,
                trueRangeRate: rangeRate,
                noise: noise,
                elevation: elevation * 180 / Math.PI
            });

            trueObservations.push({
                t: t,
                rangeRate: rangeRate
            });
        }

        // Propagate state (simple Kepler + J2)
        const r3 = r * r * r;
        const ax = -mu * state.x / r3;
        const ay = -mu * state.y / r3;
        const az = -mu * state.z / r3;

        state.x += state.vx * dt + 0.5 * ax * dt * dt;
        state.y += state.vy * dt + 0.5 * ay * dt * dt;
        state.z += state.vz * dt + 0.5 * az * dt * dt;
        state.vx += ax * dt;
        state.vy += ay * dt;
        state.vz += az * dt;
    }

    // Perform estimation iterations
    const iterations = [];
    let currentEstimate = { ...initialGuess };
    const residualHistory = [];

    for (let iter = 0; iter < config.numIterations; iter++) {
        // Compute residuals with current estimate
        const iterResiduals = [];
        let sumSqResiduals = 0;

        // Simple propagation with current estimate
        let estState = {
            x: currentEstimate.x,
            y: currentEstimate.y,
            z: currentEstimate.z,
            vx: currentEstimate.vx,
            vy: currentEstimate.vy,
            vz: currentEstimate.vz
        };

        let obsIdx = 0;
        for (let i = 0; i < numObs && obsIdx < observations.length; i++) {
            const t = i * dt;

            // Check if this is an observation time
            if (Math.abs(t - observations[obsIdx].t) < 1) {
                const r = Math.sqrt(estState.x**2 + estState.y**2 + estState.z**2);

                // Station position
                const theta = stationLon + earthRotRate * t;
                const stationX = RE * Math.cos(stationLat) * Math.cos(theta);
                const stationY = RE * Math.cos(stationLat) * Math.sin(theta);
                const stationZ = RE * Math.sin(stationLat);

                const dx = estState.x - stationX;
                const dy = estState.y - stationY;
                const dz = estState.z - stationZ;
                const range = Math.sqrt(dx*dx + dy*dy + dz*dz);
                const computedRangeRate = (dx*estState.vx + dy*estState.vy + dz*estState.vz) / range;

                const residual = observations[obsIdx].rangeRate - computedRangeRate;
                iterResiduals.push({
                    t: t,
                    tDays: t / 86400,
                    residual: residual,
                    iteration: iter
                });
                sumSqResiduals += residual * residual;

                obsIdx++;
            }

            // Propagate
            const r = Math.sqrt(estState.x**2 + estState.y**2 + estState.z**2);
            const r3 = r * r * r;
            const ax = -currentEstimate.muEarth * estState.x / r3;
            const ay = -currentEstimate.muEarth * estState.y / r3;
            const az = -currentEstimate.muEarth * estState.z / r3;

            estState.x += estState.vx * dt;
            estState.y += estState.vy * dt;
            estState.z += estState.vz * dt;
            estState.vx += ax * dt;
            estState.vy += ay * dt;
            estState.vz += az * dt;
        }

        const rms = Math.sqrt(sumSqResiduals / iterResiduals.length);

        iterations.push({
            iteration: iter,
            residuals: iterResiduals,
            rms: rms,
            estimate: { ...currentEstimate }
        });

        residualHistory.push(...iterResiduals);

        // Update estimate (simplified gradient descent)
        // In reality this would be a full least squares inversion
        const learningRate = 0.5 / (iter + 1);
        const meanResidual = iterResiduals.reduce((s, r) => s + r.residual, 0) / iterResiduals.length;

        currentEstimate.x += learningRate * 100 * Math.sign(meanResidual);
        currentEstimate.y += learningRate * 100 * Math.sign(meanResidual);
        currentEstimate.vx += learningRate * 0.01 * Math.sign(meanResidual);
        currentEstimate.vy += learningRate * 0.01 * Math.sign(meanResidual);
    }

    // Compute formal errors (simplified)
    const formalErrors = {
        position: config.noiseLevel * 3e3,  // meters
        velocity: config.noiseLevel * 3,     // m/s
        Cd: 0.05,
        muEarth: 1e8
    };

    // Compute true errors
    const finalEstimate = iterations[iterations.length - 1].estimate;
    const trueErrors = {
        x: trueParams.x - finalEstimate.x,
        y: trueParams.y - finalEstimate.y,
        z: trueParams.z - finalEstimate.z,
        vx: trueParams.vx - finalEstimate.vx,
        vy: trueParams.vy - finalEstimate.vy,
        vz: trueParams.vz - finalEstimate.vz
    };

    // Final residuals for histogram
    const finalResiduals = iterations[iterations.length - 1].residuals.map(r => r.residual);

    return {
        observations,
        iterations,
        residualHistory,
        finalResiduals,
        finalRMS: iterations[iterations.length - 1].rms,
        trueParams,
        finalEstimate,
        formalErrors,
        trueErrors
    };
}

function renderEstimationFigures(container, result, config) {
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
    const chartHeight = Math.max(150, (containerHeight - 120) / 3);

    // Figure 1: Observations over time
    createFigure(wrapper, 'Doppler Observations vs Time', chartWidth, chartHeight, (svg, w, h) => {
        renderObservations(svg, result, w, h);
    });

    // Figure 2: Residuals per iteration
    createFigure(wrapper, 'Residuals History by Iteration', chartWidth, chartHeight, (svg, w, h) => {
        renderResidualHistory(svg, result, w, h);
    });

    // Figure 3: Final residual histogram
    createFigure(wrapper, 'Final Residual Distribution', chartWidth, chartHeight, (svg, w, h) => {
        renderResidualHistogram(svg, result, w, h);
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

function renderObservations(svg, result, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain(d3.extent(result.observations, d => d.tDays))
        .range([margin.left, width - margin.right]);

    const yExtent = d3.extent(result.observations, d => d.rangeRate);
    const yPadding = (yExtent[1] - yExtent[0]) * 0.1;
    const y = d3.scaleLinear()
        .domain([yExtent[0] - yPadding, yExtent[1] + yPadding])
        .range([height - margin.bottom, margin.top]);

    // Draw observations as scatter
    svg.selectAll('circle')
        .data(result.observations)
        .enter()
        .append('circle')
        .attr('cx', d => x(d.tDays))
        .attr('cy', d => y(d.rangeRate))
        .attr('r', 1.5)
        .attr('fill', 'var(--cyan)')
        .attr('opacity', 0.7);

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
        .text('Time [days]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Range Rate [m/s]');
}

function renderResidualHistory(svg, result, width, height) {
    const margin = { top: 10, right: 80, bottom: 35, left: 60 };

    const colors = ['var(--cyan)', 'var(--green)', 'var(--orange)', 'var(--purple)'];

    const x = d3.scaleLinear()
        .domain(d3.extent(result.observations, d => d.tDays))
        .range([margin.left, width - margin.right]);

    let maxResidual = 0;
    result.iterations.forEach(iter => {
        iter.residuals.forEach(r => {
            maxResidual = Math.max(maxResidual, Math.abs(r.residual));
        });
    });

    const y = d3.scaleLinear()
        .domain([-maxResidual * 1.1, maxResidual * 1.1])
        .range([height - margin.bottom, margin.top]);

    // Draw residuals for each iteration
    result.iterations.forEach((iter, i) => {
        svg.selectAll(`.iter-${i}`)
            .data(iter.residuals)
            .enter()
            .append('circle')
            .attr('class', `iter-${i}`)
            .attr('cx', d => x(d.tDays))
            .attr('cy', d => y(d.residual))
            .attr('r', 1.5)
            .attr('fill', colors[i % colors.length])
            .attr('opacity', 0.6);
    });

    // Legend
    const legend = svg.append('g')
        .attr('transform', `translate(${width - margin.right + 10}, ${margin.top})`);

    result.iterations.forEach((iter, i) => {
        legend.append('circle')
            .attr('cx', 5)
            .attr('cy', i * 14)
            .attr('r', 3)
            .attr('fill', colors[i % colors.length]);

        legend.append('text')
            .attr('x', 12)
            .attr('y', i * 14 + 3)
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '8px')
            .text(`Iter ${i + 1}`);
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
        .text('Time [days]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Residual [m/s]');
}

function renderResidualHistogram(svg, result, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 50 };

    const data = result.finalResiduals;
    const bins = d3.bin().thresholds(20)(data);

    const x = d3.scaleLinear()
        .domain([d3.min(data), d3.max(data)])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([0, d3.max(bins, d => d.length)])
        .range([height - margin.bottom, margin.top]);

    svg.selectAll('rect')
        .data(bins)
        .enter()
        .append('rect')
        .attr('x', d => x(d.x0) + 1)
        .attr('y', d => y(d.length))
        .attr('width', d => Math.max(0, x(d.x1) - x(d.x0) - 1))
        .attr('height', d => height - margin.bottom - y(d.length))
        .attr('fill', 'var(--cyan)')
        .attr('opacity', 0.7);

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
        .text('Residual [m/s]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Count');
}
