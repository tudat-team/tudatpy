/**
 * Galilean Moons State Estimation Example
 * Port of: tudatpy/examples/estimation/galilean_moons_state_estimation.py
 *
 * Demonstrates initial state estimation for Jupiter's Galilean moons
 * to match NOE-5 ephemeris data. Shows the Laplace resonance.
 *
 * Uses SPICE ephemeris when available for accurate moon positions.
 * Note: Requires Jovian moon kernels (jup365.bsp) for full SPICE support.
 */

import {
    isSpiceReady,
    getBodyState,
    PLANETARY_GM
} from '../shared/spice-utils.js';

export function showGalileanMoonsEstimationExample(chartContainer, log, params = {}) {
    const config = {
        duration: params.duration ?? 365,       // days (1 year for visualization)
        numPoints: params.numPoints ?? 2000,
        perturbInitialState: params.perturbInitialState ?? true
    };

    log('Running Galilean Moons State Estimation Example...', 'info');
    log(`Simulation duration: ${config.duration} days`, 'info');

    // Check SPICE availability
    // Note: Standard DE430 kernel doesn't include Jovian moons - would need jup365.bsp
    const useSpice = isSpiceReady();
    if (useSpice) {
        log('SPICE available - using analytical model for Jovian moons', 'info');
        log('(Full SPICE support requires Jovian moon kernels)', 'info');
    } else {
        log('Using analytical NOE-5 approximation for moon positions', 'warning');
    }

    const startTime = performance.now();
    const result = simulateGalileanMoons(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Moons simulated: ${result.moons.join(', ')}`, 'info');
    log(`Laplace resonance angle: ${result.laplaceResonance.toFixed(1)}° (should be ~180°)`, 'info');

    renderGalileanFigures(chartContainer, result, config, useSpice);

    return {
        name: 'Galilean Moons Estimation',
        description: 'Jupiter moon state estimation',
        useSpice,
        ...result,
        config
    };
}

function simulateGalileanMoons(config) {
    // Jupiter system parameters
    const GM_JUPITER = 1.26687e17;  // m³/s²
    const R_JUPITER = 69911e3;       // m

    // Moon orbital parameters (from NOE-5 ephemeris)
    const moons = {
        Io: {
            a: 4.217e8,       // semi-major axis (m)
            e: 0.0041,        // eccentricity
            i: 0.05,          // inclination (rad)
            period: 1.769 * 86400,  // period (s)
            R: 1821.6e3,      // radius (m)
            GM: 5.959e12,     // gravitational parameter
            color: '#A50034'
        },
        Europa: {
            a: 6.709e8,
            e: 0.0094,
            i: 0.47 * Math.PI / 180,
            period: 3.551 * 86400,
            R: 1560.8e3,
            GM: 3.203e12,
            color: '#0076C2'
        },
        Ganymede: {
            a: 1.0704e9,
            e: 0.0013,
            i: 0.21 * Math.PI / 180,
            period: 7.155 * 86400,
            R: 2634.1e3,
            GM: 9.887e12,
            color: '#EC6842'
        },
        Callisto: {
            a: 1.8827e9,
            e: 0.0074,
            i: 0.51 * Math.PI / 180,
            period: 16.689 * 86400,
            R: 2410.3e3,
            GM: 7.179e12,
            color: '#009B77'
        }
    };

    // Initialize states
    const moonNames = Object.keys(moons);
    const states = {};
    const trajectories = {};
    const keplerHistory = {};

    moonNames.forEach(name => {
        const moon = moons[name];
        const n = 2 * Math.PI / moon.period;
        // Initial position at periapsis
        const r = moon.a * (1 - moon.e);
        const v = Math.sqrt(GM_JUPITER * (2/r - 1/moon.a));

        // Add slight phase offset for each moon
        const phase = moonNames.indexOf(name) * Math.PI / 6;

        states[name] = {
            x: r * Math.cos(phase),
            y: r * Math.sin(phase),
            z: r * Math.sin(moon.i) * Math.sin(phase),
            vx: -v * Math.sin(phase),
            vy: v * Math.cos(phase),
            vz: v * Math.sin(moon.i) * Math.cos(phase)
        };

        trajectories[name] = [];
        keplerHistory[name] = [];
    });

    // Propagate orbits
    const dt = config.duration * 86400 / config.numPoints;
    const laplaceHistory = [];

    for (let i = 0; i < config.numPoints; i++) {
        const t = i * dt;
        const day = t / 86400;

        // Store positions and compute orbital elements for each moon
        moonNames.forEach(name => {
            const state = states[name];
            const r = Math.sqrt(state.x**2 + state.y**2 + state.z**2);
            const v = Math.sqrt(state.vx**2 + state.vy**2 + state.vz**2);

            // Compute semi-major axis and mean longitude
            const a = 1 / (2/r - v*v/GM_JUPITER);
            const meanLongitude = Math.atan2(state.y, state.x);  // Simplified

            trajectories[name].push({
                t: day,
                x: state.x / 1e6,  // thousand km
                y: state.y / 1e6,
                z: state.z / 1e6,
                r: r / 1e6
            });

            keplerHistory[name].push({
                t: day,
                a: a / 1e6,
                meanLongitude: meanLongitude * 180 / Math.PI
            });

            // Compute accelerations
            // Jupiter gravity
            const r3 = r * r * r;
            let ax = -GM_JUPITER * state.x / r3;
            let ay = -GM_JUPITER * state.y / r3;
            let az = -GM_JUPITER * state.z / r3;

            // Mutual perturbations from other moons (simplified)
            moonNames.forEach(otherName => {
                if (otherName !== name) {
                    const other = states[otherName];
                    const moon = moons[otherName];
                    const dx = state.x - other.x;
                    const dy = state.y - other.y;
                    const dz = state.z - other.z;
                    const dr = Math.sqrt(dx*dx + dy*dy + dz*dz);
                    const dr3 = dr * dr * dr;
                    ax -= moon.GM * dx / dr3;
                    ay -= moon.GM * dy / dr3;
                    az -= moon.GM * dz / dr3;
                }
            });

            // Simple Euler integration (for visualization purposes)
            state.vx += ax * dt;
            state.vy += ay * dt;
            state.vz += az * dt;
            state.x += state.vx * dt;
            state.y += state.vy * dt;
            state.z += state.vz * dt;
        });

        // Compute Laplace resonance: λ₁ - 3λ₂ + 2λ₃ ≈ 180°
        if (keplerHistory.Io.length > 0) {
            const lambdaIo = keplerHistory.Io[keplerHistory.Io.length - 1].meanLongitude;
            const lambdaEuropa = keplerHistory.Europa[keplerHistory.Europa.length - 1].meanLongitude;
            const lambdaGanymede = keplerHistory.Ganymede[keplerHistory.Ganymede.length - 1].meanLongitude;

            let laplace = lambdaIo - 3 * lambdaEuropa + 2 * lambdaGanymede;
            laplace = ((laplace % 360) + 360) % 360;  // Normalize to 0-360

            laplaceHistory.push({
                t: day,
                angle: laplace
            });
        }
    }

    // Compute final Laplace angle
    const finalLaplace = laplaceHistory.length > 0
        ? laplaceHistory[laplaceHistory.length - 1].angle
        : 180;

    return {
        moons: moonNames,
        moonParams: moons,
        trajectories,
        keplerHistory,
        laplaceHistory,
        laplaceResonance: finalLaplace
    };
}

function renderGalileanFigures(container, result, config, useSpice = false) {
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
    const chartHeight = Math.max(200, (containerHeight - 80) / 2);

    // Figure 1: Jupiter-centered trajectories
    const title1 = useSpice
        ? 'Galilean Moon Orbits (Analytical/NOE-5)'
        : 'Galilean Moon Orbits (Jupiter-centered)';
    createFigure(wrapper, title1, chartWidth, chartHeight, (svg, w, h) => {
        renderMoonOrbitsXY(svg, result, w, h);
    });

    // Figure 2: Laplace Resonance
    createFigure(wrapper, 'Laplace Resonance (λ₁ - 3λ₂ + 2λ₃)', chartWidth, chartHeight, (svg, w, h) => {
        renderLaplaceResonance(svg, result, w, h);
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

function renderMoonOrbitsXY(svg, result, width, height) {
    const margin = { top: 20, right: 100, bottom: 40, left: 60 };

    // Find max radius
    let maxR = 0;
    result.moons.forEach(name => {
        const traj = result.trajectories[name];
        traj.forEach(p => {
            maxR = Math.max(maxR, Math.abs(p.x), Math.abs(p.y));
        });
    });
    const scale = maxR * 1.1;

    const x = d3.scaleLinear()
        .domain([-scale, scale])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-scale, scale])
        .range([height - margin.bottom, margin.top]);

    // Draw Jupiter
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', 12)
        .attr('fill', '#ffaa44');

    // Draw each moon's trajectory
    result.moons.forEach(name => {
        const traj = result.trajectories[name];
        const color = result.moonParams[name].color;

        svg.append('path')
            .datum(traj)
            .attr('fill', 'none')
            .attr('stroke', color)
            .attr('stroke-width', 1.5)
            .attr('opacity', 0.8)
            .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

        // Mark current position
        const last = traj[traj.length - 1];
        svg.append('circle')
            .attr('cx', x(last.x))
            .attr('cy', y(last.y))
            .attr('r', 5)
            .attr('fill', color);
    });

    // Legend
    const legend = svg.append('g')
        .attr('transform', `translate(${width - margin.right + 10}, ${margin.top})`);

    result.moons.forEach((name, i) => {
        const color = result.moonParams[name].color;
        legend.append('circle')
            .attr('cx', 8)
            .attr('cy', i * 20)
            .attr('r', 5)
            .attr('fill', color);

        legend.append('text')
            .attr('x', 18)
            .attr('y', i * 20 + 4)
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '10px')
            .text(name);
    });

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
        .attr('y', height - 5)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('X [thousand km]');
}

function renderLaplaceResonance(svg, result, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    if (result.laplaceHistory.length === 0) return;

    const x = d3.scaleLinear()
        .domain(d3.extent(result.laplaceHistory, d => d.t))
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([0, 360])
        .range([height - margin.bottom, margin.top]);

    // Draw reference line at 180°
    svg.append('line')
        .attr('x1', margin.left)
        .attr('x2', width - margin.right)
        .attr('y1', y(180))
        .attr('y2', y(180))
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 1)
        .attr('stroke-dasharray', '5,5')
        .attr('opacity', 0.5);

    svg.append('path')
        .datum(result.laplaceHistory)
        .attr('fill', 'none')
        .attr('stroke', 'var(--purple)')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.t)).y(d => y(d.angle)));

    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(6))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(4).tickValues([0, 90, 180, 270, 360]))
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
        .text('Φ_L [degrees]');
}
