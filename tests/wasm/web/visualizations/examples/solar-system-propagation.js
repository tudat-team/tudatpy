/**
 * Solar System Propagation Example
 * Port of: tudatpy/examples/propagation/solar_system_propagation.py
 *
 * Demonstrates multi-body solar system propagation over 5 years.
 * Shows planetary orbits and their evolution.
 *
 * Uses SPICE ephemeris when available for accurate planetary positions.
 */

import {
    isSpiceReady,
    getBodyState,
    jdToEt,
    PLANETARY_SMA,
    PLANETARY_PERIOD
} from '../shared/spice-utils.js';

export function showSolarSystemExample(chartContainer, log, params = {}) {
    const config = {
        duration: params.duration ?? 5 * 365.25 * 86400,  // 5 years in seconds
        numPoints: params.numPoints ?? 2000,
        centralBody: params.centralBody ?? 'Sun',
        startEpoch: params.startEpoch ?? 0  // J2000 epoch (ET = 0)
    };

    log('Running Solar System Propagation Example...', 'info');
    log(`Duration: ${(config.duration / (365.25 * 86400)).toFixed(1)} years`, 'info');
    log(`Central body: ${config.centralBody}`, 'info');

    // Check SPICE availability
    const useSpice = isSpiceReady();
    if (useSpice) {
        log('Using SPICE ephemeris for planetary positions', 'success');
    } else {
        log('SPICE not available - using analytical Keplerian approximation', 'warning');
    }

    const startTime = performance.now();
    const result = useSpice
        ? simulateSolarSystemSpice(config, log)
        : simulateSolarSystemAnalytical(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    for (const planet of Object.keys(result.trajectories)) {
        log(`${planet}: ${result.trajectories[planet].length} points`, 'info');
    }

    renderSolarSystemFigures(chartContainer, result, config);

    return {
        name: 'Solar System',
        description: 'Multi-body solar system propagation',
        useSpice,
        ...result,
        config
    };
}

/**
 * Simulate solar system using SPICE ephemeris
 */
function simulateSolarSystemSpice(config, log) {
    const AU = 149597870.7;  // km
    const planetNames = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn'];
    const planetColors = {
        Mercury: '#b5b5b5',
        Venus: '#e6c87a',
        Earth: '#4a9fff',
        Mars: '#ff6b4a',
        Jupiter: '#e8a87c',
        Saturn: '#e8d4a8'
    };

    const trajectories = {};
    const dt = config.duration / config.numPoints;

    for (const name of planetNames) {
        const trajectory = [];

        for (let i = 0; i < config.numPoints; i++) {
            const t = i * dt;
            const et = config.startEpoch + t;

            // Get state from SPICE
            const state = getBodyState(name, 'Sun', et, 'ECLIPJ2000');

            if (state) {
                // Convert from meters to AU
                const x = state.x / (AU * 1000);
                const y = state.y / (AU * 1000);
                const z = state.z / (AU * 1000);
                const r = Math.sqrt(x * x + y * y + z * z);

                trajectory.push({
                    t: t / (365.25 * 86400),  // years
                    x, y, z,  // AU
                    r
                });
            }
        }

        if (trajectory.length > 0) {
            trajectories[name] = trajectory;
        } else {
            log(`Warning: No SPICE data for ${name}`, 'warning');
        }
    }

    // Build planets object for compatibility
    const planets = {};
    for (const name of planetNames) {
        planets[name] = {
            color: planetColors[name],
            a: (PLANETARY_SMA[name] || 1e11) / (AU * 1000),  // Convert to AU
            period: (PLANETARY_PERIOD[name] || 365.25 * 86400) / 86400  // Convert to days
        };
    }

    return { trajectories, planets, source: 'SPICE' };
}

/**
 * Simulate solar system using analytical Keplerian elements (fallback)
 */
function simulateSolarSystemAnalytical(config) {
    // Simplified Keplerian elements for inner planets (J2000)
    const planets = {
        Mercury: { a: 0.387, e: 0.206, i: 7.0, omega: 29.12, Omega: 48.33, period: 87.97, color: '#b5b5b5' },
        Venus: { a: 0.723, e: 0.007, i: 3.39, omega: 54.85, Omega: 76.68, period: 224.7, color: '#e6c87a' },
        Earth: { a: 1.000, e: 0.017, i: 0.0, omega: 102.9, Omega: 0, period: 365.25, color: '#4a9fff' },
        Mars: { a: 1.524, e: 0.093, i: 1.85, omega: 286.5, Omega: 49.56, period: 687.0, color: '#ff6b4a' }
    };

    // Outer planets
    const outerPlanets = {
        Jupiter: { a: 5.203, e: 0.049, i: 1.30, omega: 274.2, Omega: 100.5, period: 4332.6, color: '#e8a87c' },
        Saturn: { a: 9.537, e: 0.054, i: 2.48, omega: 338.7, Omega: 113.7, period: 10759, color: '#e8d4a8' }
    };

    const trajectories = {};
    const dt = config.duration / config.numPoints;

    // Generate trajectories for each planet
    for (const [name, params] of Object.entries({...planets, ...outerPlanets})) {
        const trajectory = [];
        const n = 2 * Math.PI / (params.period * 86400);  // Mean motion rad/s
        const incRad = params.i * Math.PI / 180;
        const omegaRad = params.omega * Math.PI / 180;
        const OmegaRad = params.Omega * Math.PI / 180;

        for (let i = 0; i < config.numPoints; i++) {
            const t = i * dt;

            // Mean anomaly
            const M = n * t;

            // Solve Kepler's equation (Newton-Raphson)
            let E = M;
            for (let j = 0; j < 10; j++) {
                E = E - (E - params.e * Math.sin(E) - M) / (1 - params.e * Math.cos(E));
            }

            // True anomaly
            const nu = 2 * Math.atan2(
                Math.sqrt(1 + params.e) * Math.sin(E / 2),
                Math.sqrt(1 - params.e) * Math.cos(E / 2)
            );

            // Orbital radius
            const r = params.a * (1 - params.e * Math.cos(E));

            // Position in orbital plane
            const xOrb = r * Math.cos(nu);
            const yOrb = r * Math.sin(nu);

            // Transform to ecliptic coordinates
            const cosO = Math.cos(OmegaRad);
            const sinO = Math.sin(OmegaRad);
            const cosI = Math.cos(incRad);
            const sinI = Math.sin(incRad);
            const cosW = Math.cos(omegaRad);
            const sinW = Math.sin(omegaRad);

            const x = (cosO * cosW - sinO * sinW * cosI) * xOrb + (-cosO * sinW - sinO * cosW * cosI) * yOrb;
            const y = (sinO * cosW + cosO * sinW * cosI) * xOrb + (-sinO * sinW + cosO * cosW * cosI) * yOrb;
            const z = (sinW * sinI) * xOrb + (cosW * sinI) * yOrb;

            trajectory.push({
                t: t / (365.25 * 86400),  // years
                x, y, z,  // AU
                r
            });
        }

        trajectories[name] = trajectory;
    }

    return { trajectories, planets: {...planets, ...outerPlanets}, source: 'analytical' };
}

function renderSolarSystemFigures(container, result, config) {
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

    // Figure 1: Inner solar system (top-down view)
    const title1 = result.source === 'SPICE'
        ? 'Inner Solar System (SPICE Ephemeris)'
        : 'Inner Solar System (XY Ecliptic Plane)';
    createFigure(wrapper, title1, chartWidth, chartHeight, (svg, w, h) => {
        renderOrbitalView(svg, result, ['Mercury', 'Venus', 'Earth', 'Mars'], w, h, 2.0);
    });

    // Figure 2: Outer solar system
    const title2 = result.source === 'SPICE'
        ? 'Outer Solar System (SPICE Ephemeris)'
        : 'Outer Solar System (XY Ecliptic Plane)';
    createFigure(wrapper, title2, chartWidth, chartHeight, (svg, w, h) => {
        renderOrbitalView(svg, result, ['Earth', 'Mars', 'Jupiter', 'Saturn'], w, h, 12.0);
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

function renderOrbitalView(svg, result, planetNames, width, height, maxAU) {
    const margin = { top: 20, right: 20, bottom: 40, left: 50 };

    const x = d3.scaleLinear()
        .domain([-maxAU, maxAU])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-maxAU, maxAU])
        .range([height - margin.bottom, margin.top]);

    // Draw Sun at center
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', 8)
        .attr('fill', '#ffcc00');

    // Draw orbits for selected planets
    for (const name of planetNames) {
        const trajectory = result.trajectories[name];
        const planet = result.planets[name];
        if (!trajectory) continue;

        // Draw orbit path
        svg.append('path')
            .datum(trajectory)
            .attr('fill', 'none')
            .attr('stroke', planet.color)
            .attr('stroke-width', 1.5)
            .attr('stroke-opacity', 0.7)
            .attr('d', d3.line()
                .x(d => x(d.x))
                .y(d => y(d.y)));

        // Draw planet at final position
        const lastPos = trajectory[trajectory.length - 1];
        svg.append('circle')
            .attr('cx', x(lastPos.x))
            .attr('cy', y(lastPos.y))
            .attr('r', 5)
            .attr('fill', planet.color);

        // Label
        svg.append('text')
            .attr('x', x(lastPos.x) + 8)
            .attr('y', y(lastPos.y) + 4)
            .attr('fill', planet.color)
            .attr('font-size', '10px')
            .text(name);
    }

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
        .text('X [AU]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Y [AU]');
}
