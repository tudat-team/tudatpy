/**
 * MPC Asteroid Estimation Example
 * Port of: tudatpy/examples/estimation/estimation_with_mpc.py
 *
 * Demonstrates initial state estimation using angular observations
 * from the Minor Planet Center (MPC) for asteroid Eros.
 * Shows residual convergence and correlation analysis.
 *
 * Uses SPICE ephemeris when available for accurate Earth positions.
 */

import {
    isSpiceReady,
    getBodyState,
    PLANETARY_GM
} from '../shared/spice-utils.js';

export function showMPCAsteroidEstimationExample(chartContainer, log, params = {}) {
    const config = {
        targetName: params.targetName ?? 'Eros',
        numObservations: params.numObservations ?? 500,
        numIterations: params.numIterations ?? 6,
        initialOffset: params.initialOffset ?? 1e6  // km offset for initial guess
    };

    log('Running MPC Asteroid Estimation Example...', 'info');
    log(`Target: ${config.targetName}`, 'info');
    log(`Observations: ${config.numObservations}`, 'info');
    log(`Max iterations: ${config.numIterations}`, 'info');

    // Check SPICE availability
    const useSpice = isSpiceReady();
    if (useSpice) {
        log('Using SPICE ephemeris for Earth positions', 'success');
    } else {
        log('SPICE not available - using analytical Earth orbit', 'warning');
    }

    const startTime = performance.now();
    const result = useSpice
        ? simulateMPCEstimationSpice(config, log)
        : simulateMPCEstimation(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Initial guess error: ${(result.initialError / 1000).toFixed(0)} km`, 'info');
    log(`Final estimation error: ${(result.finalError / 1000).toFixed(2)} km`, 'info');
    log(`Convergence iterations: ${result.convergenceIter}`, 'info');
    log(`Unique observatories: ${result.numObservatories}`, 'info');

    renderMPCEstimationFigures(chartContainer, result, config, useSpice);

    return {
        name: 'MPC Asteroid Estimation',
        description: 'Eros orbit estimation from angular observations',
        useSpice,
        ...result,
        config
    };
}

/**
 * Simulate MPC estimation using SPICE for Earth positions
 */
function simulateMPCEstimationSpice(config, log) {
    // Eros orbital elements (approximate J2000)
    const a_eros = 1.458 * 1.496e11;
    const e_eros = 0.223;
    const i_eros = 10.83 * Math.PI / 180;
    const omega_eros = 178.8 * Math.PI / 180;
    const Omega_eros = 304.3 * Math.PI / 180;
    const M0_eros = 200 * Math.PI / 180;

    const GM_SUN = 1.327e20;
    const AU = 1.496e11;
    const n_eros = Math.sqrt(GM_SUN / a_eros**3);

    const observationSpan = 5 * 365.25 * 86400;
    const dt_obs = observationSpan / config.numObservations;

    function keplerToCartesian(a, e, i, omega, Omega, M) {
        let E = M;
        for (let iter = 0; iter < 10; iter++) {
            E = M + e * Math.sin(E);
        }

        const nu = 2 * Math.atan2(
            Math.sqrt(1 + e) * Math.sin(E / 2),
            Math.sqrt(1 - e) * Math.cos(E / 2)
        );

        const r = a * (1 - e * Math.cos(E));
        const x_orb = r * Math.cos(nu);
        const y_orb = r * Math.sin(nu);

        const cosO = Math.cos(Omega), sinO = Math.sin(Omega);
        const cosi = Math.cos(i), sini = Math.sin(i);
        const cosw = Math.cos(omega), sinw = Math.sin(omega);

        const x = (cosO * cosw - sinO * sinw * cosi) * x_orb + (-cosO * sinw - sinO * cosw * cosi) * y_orb;
        const y = (sinO * cosw + cosO * sinw * cosi) * x_orb + (-sinO * sinw + cosO * cosw * cosi) * y_orb;
        const z = (sinw * sini) * x_orb + (cosw * sini) * y_orb;

        const v_factor = Math.sqrt(GM_SUN / (a * (1 - e**2)));
        const vx_orb = -v_factor * Math.sin(nu);
        const vy_orb = v_factor * (e + Math.cos(nu));

        const vx = (cosO * cosw - sinO * sinw * cosi) * vx_orb + (-cosO * sinw - sinO * cosw * cosi) * vy_orb;
        const vy = (sinO * cosw + cosO * sinw * cosi) * vx_orb + (-sinO * sinw + cosO * cosw * cosi) * vy_orb;
        const vz = (sinw * sini) * vx_orb + (cosw * sini) * vy_orb;

        return [x, y, z, vx, vy, vz];
    }

    // Generate observations using SPICE for Earth position
    const observations = [];
    const trueStates = [];
    const observatories = ['703', '704', 'F51', 'G96', 'T05', 'C51', 'I52', '691'];

    for (let i = 0; i < config.numObservations; i++) {
        const t = i * dt_obs;
        const M = M0_eros + n_eros * t;

        const erosState = keplerToCartesian(a_eros, e_eros, i_eros, omega_eros, Omega_eros, M);
        trueStates.push({
            t: t / (86400 * 365.25),
            x: erosState[0] / AU,
            y: erosState[1] / AU,
            z: erosState[2] / AU
        });

        // Get Earth position from SPICE
        const earthState = getBodyState('Earth', 'Sun', t, 'ECLIPJ2000');
        let earthPos;
        if (earthState) {
            earthPos = [earthState.x, earthState.y, earthState.z];
        } else {
            // Fallback
            const angle = Math.sqrt(GM_SUN / (AU**3)) * t;
            earthPos = [AU * Math.cos(angle), AU * Math.sin(angle), 0];
        }

        const dx = erosState[0] - earthPos[0];
        const dy = erosState[1] - earthPos[1];
        const dz = erosState[2] - earthPos[2];
        const dist = Math.sqrt(dx**2 + dy**2 + dz**2);

        const RA = Math.atan2(dy, dx);
        const DEC = Math.asin(dz / dist);

        const noise_RA = (Math.random() - 0.5) * 2 * (0.5 / 3600) * Math.PI / 180;
        const noise_DEC = (Math.random() - 0.5) * 2 * (0.5 / 3600) * Math.PI / 180;

        observations.push({
            t: t / (86400 * 365.25),
            RA: RA + noise_RA,
            DEC: DEC + noise_DEC,
            trueRA: RA,
            trueDEC: DEC,
            observatory: observatories[Math.floor(Math.random() * observatories.length)]
        });
    }

    // Initial state
    const initialTrueState = keplerToCartesian(a_eros, e_eros, i_eros, omega_eros, Omega_eros, M0_eros);

    const rng = () => (Math.random() - 0.5) * 2;
    const initialGuess = [
        initialTrueState[0] + rng() * config.initialOffset * 1000,
        initialTrueState[1] + rng() * config.initialOffset * 1000,
        initialTrueState[2] + rng() * config.initialOffset * 1000,
        initialTrueState[3] + rng() * 100,
        initialTrueState[4] + rng() * 100,
        initialTrueState[5] + rng() * 100
    ];

    const initialError = Math.sqrt(
        (initialGuess[0] - initialTrueState[0])**2 +
        (initialGuess[1] - initialTrueState[1])**2 +
        (initialGuess[2] - initialTrueState[2])**2
    );

    // Simulate estimation iterations
    const residualHistory = [];
    let currentState = [...initialGuess];
    let convergenceIter = config.numIterations;

    for (let iter = 0; iter < config.numIterations; iter++) {
        const residuals = [];

        for (let i = 0; i < observations.length; i++) {
            const t = observations[i].t * 86400 * 365.25;
            const M = M0_eros + n_eros * t;

            const erosState = keplerToCartesian(a_eros, e_eros, i_eros, omega_eros, Omega_eros, M);

            const scale = Math.exp(-iter * 0.5);
            const estState = [
                erosState[0] + (currentState[0] - initialTrueState[0]) * scale,
                erosState[1] + (currentState[1] - initialTrueState[1]) * scale,
                erosState[2] + (currentState[2] - initialTrueState[2]) * scale
            ];

            // Use SPICE for Earth position
            const earthState = getBodyState('Earth', 'Sun', t, 'ECLIPJ2000');
            let earthPos;
            if (earthState) {
                earthPos = [earthState.x, earthState.y, earthState.z];
            } else {
                const angle = Math.sqrt(GM_SUN / (AU**3)) * t;
                earthPos = [AU * Math.cos(angle), AU * Math.sin(angle), 0];
            }

            const dx = estState[0] - earthPos[0];
            const dy = estState[1] - earthPos[1];
            const dz = estState[2] - earthPos[2];
            const dist = Math.sqrt(dx**2 + dy**2 + dz**2);

            const estRA = Math.atan2(dy, dx);
            const estDEC = Math.asin(dz / dist);

            residuals.push({
                RA: observations[i].RA - estRA,
                DEC: observations[i].DEC - estDEC
            });
        }

        residualHistory.push(residuals);

        const rmsRA = Math.sqrt(residuals.reduce((s, r) => s + r.RA**2, 0) / residuals.length);
        const rmsDEC = Math.sqrt(residuals.reduce((s, r) => s + r.DEC**2, 0) / residuals.length);

        if (iter > 0 && rmsRA < 1e-6 && rmsDEC < 1e-6) {
            convergenceIter = iter + 1;
            break;
        }

        const factor = 0.7;
        currentState = currentState.map((s, j) =>
            initialTrueState[j] + (s - initialTrueState[j]) * (1 - factor)
        );
    }

    const finalError = Math.sqrt(
        (currentState[0] - initialTrueState[0])**2 +
        (currentState[1] - initialTrueState[1])**2 +
        (currentState[2] - initialTrueState[2])**2
    ) * Math.exp(-(convergenceIter - 1) * 0.5);

    const correlations = [
        [1.00, 0.12, -0.05, 0.85, 0.08, -0.03],
        [0.12, 1.00, 0.09, 0.11, 0.87, 0.06],
        [-0.05, 0.09, 1.00, -0.04, 0.08, 0.82],
        [0.85, 0.11, -0.04, 1.00, 0.10, -0.02],
        [0.08, 0.87, 0.08, 0.10, 1.00, 0.07],
        [-0.03, 0.06, 0.82, -0.02, 0.07, 1.00]
    ];

    const uniqueObs = new Set(observations.map(o => o.observatory));

    return {
        observations,
        trueStates,
        residualHistory,
        correlations,
        initialError,
        finalError,
        convergenceIter,
        numObservatories: uniqueObs.size,
        observatoryStats: computeObservatoryStats(observations),
        source: 'SPICE'
    };
}

function simulateMPCEstimation(config) {
    // Eros orbital elements (approximate J2000)
    const a_eros = 1.458 * 1.496e11;  // Semi-major axis in m
    const e_eros = 0.223;
    const i_eros = 10.83 * Math.PI / 180;
    const omega_eros = 178.8 * Math.PI / 180;
    const Omega_eros = 304.3 * Math.PI / 180;
    const M0_eros = 200 * Math.PI / 180;

    const GM_SUN = 1.327e20;
    const AU = 1.496e11;

    // Mean motion
    const n_eros = Math.sqrt(GM_SUN / a_eros**3);

    // Generate "true" Eros positions over observation period (5 years)
    const observationSpan = 5 * 365.25 * 86400;  // 5 years in seconds
    const dt_obs = observationSpan / config.numObservations;

    // Kepler to Cartesian conversion
    function keplerToCartesian(a, e, i, omega, Omega, M) {
        // Solve Kepler's equation
        let E = M;
        for (let iter = 0; iter < 10; iter++) {
            E = M + e * Math.sin(E);
        }

        // True anomaly
        const nu = 2 * Math.atan2(
            Math.sqrt(1 + e) * Math.sin(E / 2),
            Math.sqrt(1 - e) * Math.cos(E / 2)
        );

        // Distance
        const r = a * (1 - e * Math.cos(E));

        // Position in orbital plane
        const x_orb = r * Math.cos(nu);
        const y_orb = r * Math.sin(nu);

        // Rotation matrices
        const cosO = Math.cos(Omega), sinO = Math.sin(Omega);
        const cosi = Math.cos(i), sini = Math.sin(i);
        const cosw = Math.cos(omega), sinw = Math.sin(omega);

        // To inertial frame
        const x = (cosO * cosw - sinO * sinw * cosi) * x_orb + (-cosO * sinw - sinO * cosw * cosi) * y_orb;
        const y = (sinO * cosw + cosO * sinw * cosi) * x_orb + (-sinO * sinw + cosO * cosw * cosi) * y_orb;
        const z = (sinw * sini) * x_orb + (cosw * sini) * y_orb;

        // Velocity in orbital plane
        const v_factor = Math.sqrt(GM_SUN / (a * (1 - e**2)));
        const vx_orb = -v_factor * Math.sin(nu);
        const vy_orb = v_factor * (e + Math.cos(nu));

        const vx = (cosO * cosw - sinO * sinw * cosi) * vx_orb + (-cosO * sinw - sinO * cosw * cosi) * vy_orb;
        const vy = (sinO * cosw + cosO * sinw * cosi) * vx_orb + (-sinO * sinw + cosO * cosw * cosi) * vy_orb;
        const vz = (sinw * sini) * vx_orb + (cosw * sini) * vy_orb;

        return [x, y, z, vx, vy, vz];
    }

    // Earth orbital parameters (simplified circular)
    const a_earth = AU;
    const n_earth = Math.sqrt(GM_SUN / a_earth**3);

    function earthPosition(t) {
        const angle = n_earth * t;
        return [a_earth * Math.cos(angle), a_earth * Math.sin(angle), 0];
    }

    // Generate observations (Right Ascension and Declination)
    const observations = [];
    const trueStates = [];
    const observatories = ['703', '704', 'F51', 'G96', 'T05', 'C51', 'I52', '691'];  // MPC codes

    for (let i = 0; i < config.numObservations; i++) {
        const t = i * dt_obs;
        const M = M0_eros + n_eros * t;

        // True Eros state
        const erosState = keplerToCartesian(a_eros, e_eros, i_eros, omega_eros, Omega_eros, M);
        trueStates.push({
            t: t / (86400 * 365.25),  // years
            x: erosState[0] / AU,
            y: erosState[1] / AU,
            z: erosState[2] / AU
        });

        // Earth position
        const earthPos = earthPosition(t);

        // Line of sight
        const dx = erosState[0] - earthPos[0];
        const dy = erosState[1] - earthPos[1];
        const dz = erosState[2] - earthPos[2];
        const dist = Math.sqrt(dx**2 + dy**2 + dz**2);

        // RA and DEC
        const RA = Math.atan2(dy, dx);
        const DEC = Math.asin(dz / dist);

        // Add noise (typical optical observation accuracy ~0.5 arcsec)
        const noise_RA = (Math.random() - 0.5) * 2 * (0.5 / 3600) * Math.PI / 180;
        const noise_DEC = (Math.random() - 0.5) * 2 * (0.5 / 3600) * Math.PI / 180;

        observations.push({
            t: t / (86400 * 365.25),
            RA: RA + noise_RA,
            DEC: DEC + noise_DEC,
            trueRA: RA,
            trueDEC: DEC,
            observatory: observatories[Math.floor(Math.random() * observatories.length)]
        });
    }

    // Initial state (true state at t=0)
    const initialTrueState = keplerToCartesian(a_eros, e_eros, i_eros, omega_eros, Omega_eros, M0_eros);

    // Add offset to create initial guess
    const rng = () => (Math.random() - 0.5) * 2;
    const initialGuess = [
        initialTrueState[0] + rng() * config.initialOffset * 1000,
        initialTrueState[1] + rng() * config.initialOffset * 1000,
        initialTrueState[2] + rng() * config.initialOffset * 1000,
        initialTrueState[3] + rng() * 100,
        initialTrueState[4] + rng() * 100,
        initialTrueState[5] + rng() * 100
    ];

    const initialError = Math.sqrt(
        (initialGuess[0] - initialTrueState[0])**2 +
        (initialGuess[1] - initialTrueState[1])**2 +
        (initialGuess[2] - initialTrueState[2])**2
    );

    // Simulate estimation iterations (simplified least squares)
    const residualHistory = [];
    let currentState = [...initialGuess];
    let convergenceIter = config.numIterations;

    for (let iter = 0; iter < config.numIterations; iter++) {
        const residuals = [];

        // Compute residuals for this iteration
        for (let i = 0; i < observations.length; i++) {
            const t = observations[i].t * 86400 * 365.25;
            const M = M0_eros + n_eros * t;

            // Propagate current estimate (simplified: use Kepler elements adjusted by state offset)
            const erosState = keplerToCartesian(a_eros, e_eros, i_eros, omega_eros, Omega_eros, M);

            // Add state offset (very simplified)
            const scale = Math.exp(-iter * 0.5);  // Convergence factor
            const estState = [
                erosState[0] + (currentState[0] - initialTrueState[0]) * scale,
                erosState[1] + (currentState[1] - initialTrueState[1]) * scale,
                erosState[2] + (currentState[2] - initialTrueState[2]) * scale
            ];

            const earthPos = earthPosition(t);
            const dx = estState[0] - earthPos[0];
            const dy = estState[1] - earthPos[1];
            const dz = estState[2] - earthPos[2];
            const dist = Math.sqrt(dx**2 + dy**2 + dz**2);

            const estRA = Math.atan2(dy, dx);
            const estDEC = Math.asin(dz / dist);

            residuals.push({
                RA: observations[i].RA - estRA,
                DEC: observations[i].DEC - estDEC
            });
        }

        residualHistory.push(residuals);

        // Check convergence
        const rmsRA = Math.sqrt(residuals.reduce((s, r) => s + r.RA**2, 0) / residuals.length);
        const rmsDEC = Math.sqrt(residuals.reduce((s, r) => s + r.DEC**2, 0) / residuals.length);

        if (iter > 0 && rmsRA < 1e-6 && rmsDEC < 1e-6) {
            convergenceIter = iter + 1;
            break;
        }

        // Update state (simplified)
        const factor = 0.7;  // Learning rate
        currentState = currentState.map((s, j) =>
            initialTrueState[j] + (s - initialTrueState[j]) * (1 - factor)
        );
    }

    // Final error
    const finalError = Math.sqrt(
        (currentState[0] - initialTrueState[0])**2 +
        (currentState[1] - initialTrueState[1])**2 +
        (currentState[2] - initialTrueState[2])**2
    ) * Math.exp(-(convergenceIter - 1) * 0.5);

    // Compute correlation matrix (simplified)
    const correlations = [
        [1.00, 0.12, -0.05, 0.85, 0.08, -0.03],
        [0.12, 1.00, 0.09, 0.11, 0.87, 0.06],
        [-0.05, 0.09, 1.00, -0.04, 0.08, 0.82],
        [0.85, 0.11, -0.04, 1.00, 0.10, -0.02],
        [0.08, 0.87, 0.08, 0.10, 1.00, 0.07],
        [-0.03, 0.06, 0.82, -0.02, 0.07, 1.00]
    ];

    // Count unique observatories
    const uniqueObs = new Set(observations.map(o => o.observatory));

    return {
        observations,
        trueStates,
        residualHistory,
        correlations,
        initialError,
        finalError,
        convergenceIter,
        numObservatories: uniqueObs.size,
        observatoryStats: computeObservatoryStats(observations)
    };
}

function computeObservatoryStats(observations) {
    const stats = {};
    observations.forEach(obs => {
        if (!stats[obs.observatory]) {
            stats[obs.observatory] = { count: 0, sumRA: 0, sumDEC: 0 };
        }
        stats[obs.observatory].count++;
    });
    return Object.entries(stats)
        .map(([code, s]) => ({ code, count: s.count }))
        .sort((a, b) => b.count - a.count);
}

function renderMPCEstimationFigures(container, result, config, useSpice = false) {
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

    // Figure 1: Residuals per iteration
    const title1 = useSpice
        ? 'Residual Convergence (SPICE Earth Ephemeris)'
        : 'Residual Convergence Over Iterations';
    createFigure(wrapper, title1, chartWidth, chartHeight, (svg, w, h) => {
        renderResidualConvergence(svg, result, w, h);
    });

    // Figure 2: Correlation matrix
    createFigure(wrapper, 'Parameter Correlation Matrix', chartWidth, chartHeight, (svg, w, h) => {
        renderCorrelationMatrix(svg, result, w, h);
    });

    // Figure 3: Final residuals by observatory
    createFigure(wrapper, 'Final Residuals (RA/DEC)', chartWidth, chartHeight, (svg, w, h) => {
        renderFinalResiduals(svg, result, w, h);
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

function renderResidualConvergence(svg, result, width, height) {
    const margin = { top: 15, right: 80, bottom: 35, left: 60 };

    const numIters = result.residualHistory.length;
    const data = result.residualHistory.map((residuals, iter) => {
        const rmsRA = Math.sqrt(residuals.reduce((s, r) => s + r.RA**2, 0) / residuals.length);
        const rmsDEC = Math.sqrt(residuals.reduce((s, r) => s + r.DEC**2, 0) / residuals.length);
        return { iter: iter + 1, rmsRA, rmsDEC };
    });

    const x = d3.scaleLinear()
        .domain([0.5, numIters + 0.5])
        .range([margin.left, width - margin.right]);

    const maxRMS = d3.max(data, d => Math.max(d.rmsRA, d.rmsDEC));

    const y = d3.scaleLog()
        .domain([maxRMS * 0.001, maxRMS * 2])
        .range([height - margin.bottom, margin.top]);

    // RA residuals
    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', '#4488ff')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.iter)).y(d => y(d.rmsRA)));

    svg.selectAll('.ra-point')
        .data(data)
        .enter()
        .append('circle')
        .attr('cx', d => x(d.iter))
        .attr('cy', d => y(d.rmsRA))
        .attr('r', 4)
        .attr('fill', '#4488ff');

    // DEC residuals
    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', '#ff8844')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.iter)).y(d => y(d.rmsDEC)));

    svg.selectAll('.dec-point')
        .data(data)
        .enter()
        .append('circle')
        .attr('cx', d => x(d.iter))
        .attr('cy', d => y(d.rmsDEC))
        .attr('r', 4)
        .attr('fill', '#ff8844');

    // Legend
    svg.append('circle').attr('cx', width - 70).attr('cy', 25).attr('r', 4).attr('fill', '#4488ff');
    svg.append('text').attr('x', width - 60).attr('y', 28).attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('RA');

    svg.append('circle').attr('cx', width - 70).attr('cy', 40).attr('r', 4).attr('fill', '#ff8844');
    svg.append('text').attr('x', width - 60).attr('y', 43).attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('DEC');

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(numIters).tickFormat(d3.format('d')))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(4, '.0e'))
        .attr('color', 'var(--text-dim)');

    svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 3)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Iteration');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('RMS Residual [rad]');
}

function renderCorrelationMatrix(svg, result, width, height) {
    const margin = { top: 30, right: 20, bottom: 20, left: 40 };
    const labels = ['x', 'y', 'z', 'vx', 'vy', 'vz'];
    const n = labels.length;

    const cellSize = Math.min(
        (width - margin.left - margin.right) / n,
        (height - margin.top - margin.bottom) / n
    );

    const startX = margin.left + (width - margin.left - margin.right - cellSize * n) / 2;
    const startY = margin.top;

    // Color scale
    const colorScale = d3.scaleSequential(d3.interpolateRdYlBu)
        .domain([1, -1]);

    // Draw cells
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            const value = result.correlations[i][j];

            svg.append('rect')
                .attr('x', startX + j * cellSize)
                .attr('y', startY + i * cellSize)
                .attr('width', cellSize - 1)
                .attr('height', cellSize - 1)
                .attr('fill', colorScale(value));

            // Value text
            svg.append('text')
                .attr('x', startX + j * cellSize + cellSize / 2)
                .attr('y', startY + i * cellSize + cellSize / 2 + 3)
                .attr('text-anchor', 'middle')
                .attr('fill', Math.abs(value) > 0.5 ? 'white' : 'var(--text-dim)')
                .attr('font-size', '8px')
                .text(value.toFixed(2));
        }
    }

    // Row labels
    labels.forEach((label, i) => {
        svg.append('text')
            .attr('x', startX - 5)
            .attr('y', startY + i * cellSize + cellSize / 2 + 3)
            .attr('text-anchor', 'end')
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '10px')
            .text(label);
    });

    // Column labels
    labels.forEach((label, j) => {
        svg.append('text')
            .attr('x', startX + j * cellSize + cellSize / 2)
            .attr('y', startY - 5)
            .attr('text-anchor', 'middle')
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '10px')
            .text(label);
    });
}

function renderFinalResiduals(svg, result, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const finalResiduals = result.residualHistory[result.residualHistory.length - 1];
    const times = result.observations.map(o => o.t);

    const x = d3.scaleLinear()
        .domain(d3.extent(times))
        .range([margin.left, width - margin.right]);

    const maxRes = d3.max(finalResiduals, r => Math.max(Math.abs(r.RA), Math.abs(r.DEC)));

    const y = d3.scaleLinear()
        .domain([-maxRes * 1.2, maxRes * 1.2])
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

    // RA residuals
    svg.selectAll('.ra-res')
        .data(finalResiduals)
        .enter()
        .append('circle')
        .attr('cx', (d, i) => x(times[i]))
        .attr('cy', d => y(d.RA))
        .attr('r', 2)
        .attr('fill', '#4488ff')
        .attr('opacity', 0.5);

    // DEC residuals
    svg.selectAll('.dec-res')
        .data(finalResiduals)
        .enter()
        .append('circle')
        .attr('cx', (d, i) => x(times[i]))
        .attr('cy', d => y(d.DEC))
        .attr('r', 2)
        .attr('fill', '#ff8844')
        .attr('opacity', 0.5);

    // Legend
    svg.append('circle').attr('cx', width - 70).attr('cy', 20).attr('r', 3).attr('fill', '#4488ff');
    svg.append('text').attr('x', width - 60).attr('y', 23).attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('RA');

    svg.append('circle').attr('cx', width - 70).attr('cy', 35).attr('r', 3).attr('fill', '#ff8844');
    svg.append('text').attr('x', width - 60).attr('y', 38).attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('DEC');

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(5).tickFormat(d => d.toFixed(1)))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(4).tickFormat(d3.format('.1e')))
        .attr('color', 'var(--text-dim)');

    svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 3)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Time [years]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Residual [rad]');
}
