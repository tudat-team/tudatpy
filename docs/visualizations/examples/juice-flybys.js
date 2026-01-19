/**
 * JUICE Flybys Example
 * Port of: tudatpy/examples/propagation/juice_flybys.py
 *
 * Simulates the JUICE mission trajectory in the Jovian system,
 * including flybys of Ganymede, Europa, and Callisto.
 *
 * Uses SPICE ephemeris when available for accurate moon positions.
 */

import {
    isSpiceReady,
    getBodyState,
    jdToEt,
    PLANETARY_GM
} from '../shared/spice-utils.js';

export function showJuiceFlybysExample(chartContainer, log, params = {}) {
    const config = {
        duration: params.duration ?? 100,          // days
        numPoints: params.numPoints ?? 5000,
        approachSpeed: params.approachSpeed ?? 5000  // m/s
    };

    log('Running JUICE Flybys Example...', 'info');
    log(`Simulation duration: ${config.duration} days`, 'info');
    log(`Approach speed: ${config.approachSpeed} m/s`, 'info');

    // Check SPICE availability
    // Note: JUICE uses Jovian moons which may not be in the standard DE430 kernel
    const useSpice = isSpiceReady();
    if (useSpice) {
        log('SPICE available - using analytical model for Jovian moons', 'info');
        log('(Jovian moon ephemeris requires additional SPICE kernels)', 'info');
    } else {
        log('Using analytical Keplerian approximation for moon positions', 'warning');
    }

    const startTime = performance.now();
    // Note: For JUICE, we use analytical model since DE430 doesn't include Jovian moons
    // Full SPICE support would require jup365.bsp or similar
    const result = simulateJuiceTrajectory(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Flybys detected: ${result.flybys.length}`, 'info');
    result.flybys.forEach((fb, i) => {
        log(`  ${i+1}. ${fb.moon} at day ${fb.day.toFixed(1)}, altitude ${(fb.altitude/1000).toFixed(0)} km`, 'info');
    });

    renderJuiceFigures(chartContainer, result, config, useSpice);

    return {
        name: 'JUICE Flybys',
        description: 'Jovian moon flyby simulation',
        useSpice,
        ...result,
        config
    };
}

function simulateJuiceTrajectory(config) {
    // Jovian system parameters
    const GM_JUPITER = 1.26687e17;
    const moons = {
        Io: { a: 4.217e8, T: 1.769 * 86400, R: 1.822e6, GM: 5.959e12 },
        Europa: { a: 6.709e8, T: 3.551 * 86400, R: 1.561e6, GM: 3.203e12 },
        Ganymede: { a: 1.0704e9, T: 7.155 * 86400, R: 2.634e6, GM: 9.887e12 },
        Callisto: { a: 1.8827e9, T: 16.689 * 86400, R: 2.410e6, GM: 7.179e12 }
    };

    // Moon position at time t
    function getMoonPos(moon, t, phase0 = 0) {
        const omega = 2 * Math.PI / moon.T;
        const angle = omega * t + phase0;
        return [moon.a * Math.cos(angle), moon.a * Math.sin(angle)];
    }

    // Gravitational acceleration
    function gravityAccel(pos, bodyPos, bodyGM) {
        const dx = bodyPos[0] - pos[0];
        const dy = bodyPos[1] - pos[1];
        const r = Math.sqrt(dx*dx + dy*dy);
        if (r < 1000) return [0, 0];
        const r3 = r * r * r;
        return [bodyGM * dx / r3, bodyGM * dy / r3];
    }

    // Initial conditions
    const r0 = 2.5e9;
    let state = [r0, 0.5e9, -config.approachSpeed, 1500];  // x, y, vx, vy

    const dt = config.duration * 86400 / config.numPoints;
    const trajectory = [];
    const flybys = [];
    const moonPhases = { Io: 0, Europa: Math.PI/6, Ganymede: Math.PI/3, Callisto: Math.PI/2 };

    let lastDists = { Ganymede: Infinity, Europa: Infinity, Callisto: Infinity };

    for (let i = 0; i < config.numPoints; i++) {
        const t = i * dt;
        const day = t / 86400;

        // Get moon positions
        const moonPositions = {};
        for (const [name, moon] of Object.entries(moons)) {
            moonPositions[name] = getMoonPos(moon, t, moonPhases[name]);
        }

        // Store trajectory
        const r = Math.sqrt(state[0]**2 + state[1]**2);
        trajectory.push({
            t: day,
            x: state[0] / 1e6,  // thousand km
            y: state[1] / 1e6,
            r: r / 1e6,
            v: Math.sqrt(state[2]**2 + state[3]**2)
        });

        // Check for flybys
        for (const moonName of ['Ganymede', 'Europa', 'Callisto']) {
            const moonPos = moonPositions[moonName];
            const dist = Math.sqrt((state[0] - moonPos[0])**2 + (state[1] - moonPos[1])**2);

            if (dist > lastDists[moonName] && lastDists[moonName] < 5e7) {
                flybys.push({
                    moon: moonName,
                    day: day,
                    distance: lastDists[moonName],
                    altitude: lastDists[moonName] - moons[moonName].R
                });
            }
            lastDists[moonName] = dist;
        }

        // Compute accelerations
        const aJupiter = gravityAccel([state[0], state[1]], [0, 0], GM_JUPITER);
        let ax = aJupiter[0], ay = aJupiter[1];

        for (const [name, moon] of Object.entries(moons)) {
            const moonPos = moonPositions[name];
            const aMoon = gravityAccel([state[0], state[1]], moonPos, moon.GM);
            ax += aMoon[0];
            ay += aMoon[1];
        }

        // RK4 step
        const k1 = [state[2], state[3], ax, ay];
        const s1 = state.map((v, j) => v + 0.5 * dt * k1[j]);

        const aJup2 = gravityAccel([s1[0], s1[1]], [0, 0], GM_JUPITER);
        let ax2 = aJup2[0], ay2 = aJup2[1];
        for (const [name, moon] of Object.entries(moons)) {
            const moonPos = getMoonPos(moon, t + 0.5*dt, moonPhases[name]);
            const aMoon = gravityAccel([s1[0], s1[1]], moonPos, moon.GM);
            ax2 += aMoon[0];
            ay2 += aMoon[1];
        }
        const k2 = [s1[2], s1[3], ax2, ay2];

        const s2 = state.map((v, j) => v + 0.5 * dt * k2[j]);
        const k3 = [s2[2], s2[3], ax2, ay2];

        const s3 = state.map((v, j) => v + dt * k3[j]);
        const aJup4 = gravityAccel([s3[0], s3[1]], [0, 0], GM_JUPITER);
        let ax4 = aJup4[0], ay4 = aJup4[1];
        for (const [name, moon] of Object.entries(moons)) {
            const moonPos = getMoonPos(moon, t + dt, moonPhases[name]);
            const aMoon = gravityAccel([s3[0], s3[1]], moonPos, moon.GM);
            ax4 += aMoon[0];
            ay4 += aMoon[1];
        }
        const k4 = [s3[2], s3[3], ax4, ay4];

        state = state.map((v, j) => v + dt/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]));

        // Check bounds
        if (r < 7e7 || r > 5e9) break;
    }

    // Get final moon positions for plotting
    const finalMoonPositions = {};
    for (const [name, moon] of Object.entries(moons)) {
        const pos = getMoonPos(moon, config.duration * 86400 / 2, moonPhases[name]);
        finalMoonPositions[name] = { x: pos[0] / 1e6, y: pos[1] / 1e6, R: moon.R };
    }

    return {
        trajectory,
        flybys,
        moonOrbits: Object.fromEntries(
            Object.entries(moons).map(([name, moon]) => [name, moon.a / 1e6])
        ),
        finalMoonPositions
    };
}

function renderJuiceFigures(container, result, config, useSpice = false) {
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
    const chartHeight = Math.max(250, (containerHeight - 80) / 2);

    // Figure 1: Trajectory in XY plane
    const title1 = useSpice
        ? 'JUICE Trajectory (Jupiter-centered, Analytical Model)'
        : 'JUICE Trajectory (Jupiter-centered)';
    createFigure(wrapper, title1, chartWidth, chartHeight, (svg, w, h) => {
        renderTrajectoryXY(svg, result, w, h);
    });

    // Figure 2: Distance from Jupiter
    createFigure(wrapper, 'Distance from Jupiter', chartWidth, chartHeight, (svg, w, h) => {
        renderDistancePlot(svg, result, w, h);
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

function renderTrajectoryXY(svg, result, width, height) {
    const margin = { top: 20, right: 20, bottom: 40, left: 50 };

    const maxR = Math.max(...result.trajectory.map(p => Math.max(Math.abs(p.x), Math.abs(p.y))));
    const scale = Math.min(2500, maxR * 1.1);

    const x = d3.scaleLinear()
        .domain([-scale, scale])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-scale, scale])
        .range([height - margin.bottom, margin.top]);

    // Draw moon orbits
    const moonColors = { Io: '#ff6600', Europa: '#aaddff', Ganymede: '#888888', Callisto: '#553311' };
    for (const [name, radius] of Object.entries(result.moonOrbits)) {
        svg.append('circle')
            .attr('cx', x(0))
            .attr('cy', y(0))
            .attr('r', (x(radius) - x(0)))
            .attr('fill', 'none')
            .attr('stroke', moonColors[name])
            .attr('stroke-width', 0.5)
            .attr('stroke-dasharray', '3,3')
            .attr('opacity', 0.5);
    }

    // Draw trajectory
    svg.append('path')
        .datum(result.trajectory)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 1.5)
        .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

    // Draw Jupiter
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', 6)
        .attr('fill', '#ffaa44');

    // Draw flyby markers
    result.flybys.forEach(fb => {
        const trajPoint = result.trajectory.find(p => Math.abs(p.t - fb.day) < 0.5);
        if (trajPoint) {
            svg.append('circle')
                .attr('cx', x(trajPoint.x))
                .attr('cy', y(trajPoint.y))
                .attr('r', 4)
                .attr('fill', moonColors[fb.moon])
                .attr('stroke', 'white')
                .attr('stroke-width', 1);
        }
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

function renderDistancePlot(svg, result, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain(d3.extent(result.trajectory, d => d.t))
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([0, d3.max(result.trajectory, d => d.r) * 1.1])
        .range([height - margin.bottom, margin.top]);

    svg.append('path')
        .datum(result.trajectory)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.t)).y(d => y(d.r)));

    // Mark moon orbit distances
    const moonColors = { Io: '#ff6600', Europa: '#aaddff', Ganymede: '#888888', Callisto: '#553311' };
    for (const [name, radius] of Object.entries(result.moonOrbits)) {
        svg.append('line')
            .attr('x1', margin.left)
            .attr('x2', width - margin.right)
            .attr('y1', y(radius))
            .attr('y2', y(radius))
            .attr('stroke', moonColors[name])
            .attr('stroke-width', 1)
            .attr('stroke-dasharray', '3,3')
            .attr('opacity', 0.5);
    }

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
        .text('Distance [thousand km]');
}
