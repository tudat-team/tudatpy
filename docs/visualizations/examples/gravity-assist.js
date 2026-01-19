/**
 * Gravity Assist (Flyby) Example
 * Demonstrates planetary flyby mechanics.
 *
 * Shows how spacecraft can gain or lose velocity using planetary flybys.
 * Visualizes the hyperbolic trajectory and velocity change.
 *
 * Uses SPICE ephemeris when available for accurate planetary parameters.
 */

import {
    isSpiceReady,
    getGM,
    getRadius,
    PLANETARY_GM,
    PLANETARY_RADIUS
} from '../shared/spice-utils.js';

export function showGravityAssistExample(chartContainer, log, params = {}) {
    const config = {
        planet: params.planet ?? 'Jupiter',
        vInfIn: params.vInfIn ?? 10,           // km/s incoming excess velocity
        flybyAltitude: params.flybyAltitude ?? 500000,  // km above surface
        numPoints: params.numPoints ?? 300
    };

    // Check SPICE availability
    const useSpice = isSpiceReady();

    // Planet parameters - try SPICE first, fall back to hardcoded
    let planets;
    if (useSpice) {
        planets = {
            Venus: { mu: (getGM('Venus') || 324859e9) / 1e9, radius: (getRadius('Venus') || 6052000) / 1000, vOrbit: 35.02, color: '#e6c87a' },
            Earth: { mu: (getGM('Earth') || 398600e9) / 1e9, radius: (getRadius('Earth') || 6378000) / 1000, vOrbit: 29.78, color: '#4a9fff' },
            Mars: { mu: (getGM('Mars') || 42828e9) / 1e9, radius: (getRadius('Mars') || 3396000) / 1000, vOrbit: 24.07, color: '#ff6b4a' },
            Jupiter: { mu: (getGM('Jupiter') || 126686534e9) / 1e9, radius: (getRadius('Jupiter') || 71492000) / 1000, vOrbit: 13.07, color: '#e8a87c' },
            Saturn: { mu: (getGM('Saturn') || 37931187e9) / 1e9, radius: (getRadius('Saturn') || 60268000) / 1000, vOrbit: 9.69, color: '#e8d4a8' }
        };
    } else {
        planets = {
            Venus: { mu: 324859, radius: 6052, vOrbit: 35.02, color: '#e6c87a' },
            Earth: { mu: 398600, radius: 6378, vOrbit: 29.78, color: '#4a9fff' },
            Mars: { mu: 42828, radius: 3396, vOrbit: 24.07, color: '#ff6b4a' },
            Jupiter: { mu: 126686534, radius: 71492, vOrbit: 13.07, color: '#e8a87c' },
            Saturn: { mu: 37931187, radius: 60268, vOrbit: 9.69, color: '#e8d4a8' }
        };
    }

    const planet = planets[config.planet] || planets.Jupiter;

    log('Running Gravity Assist Example...', 'info');
    if (useSpice) {
        log('Using SPICE for planetary parameters', 'success');
    } else {
        log('SPICE not available - using hardcoded planetary parameters', 'warning');
    }
    log(`Flyby planet: ${config.planet}`, 'info');
    log(`V∞ incoming: ${config.vInfIn} km/s`, 'info');
    log(`Flyby altitude: ${config.flybyAltitude} km`, 'info');

    const startTime = performance.now();
    const result = computeGravityAssist(config, planet);
    const elapsed = performance.now() - startTime;
    log(`Computation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Turn angle: ${result.turnAngle.toFixed(1)}°`, 'info');
    log(`V∞ outgoing: ${result.vInfOut.toFixed(2)} km/s`, 'info');
    log(`Velocity change (heliocentric): ${result.deltaVHelio.toFixed(2)} km/s`, 'info');
    log(`Periapsis velocity: ${result.vPeriapsis.toFixed(2)} km/s`, 'info');

    renderGravityAssistFigures(chartContainer, result, config, planet, useSpice);

    return {
        name: 'Gravity Assist',
        description: `${config.planet} flyby trajectory`,
        useSpice,
        ...result,
        config
    };
}

function computeGravityAssist(config, planet) {
    const rPeriapsis = planet.radius + config.flybyAltitude;
    const vInf = config.vInfIn;

    // Hyperbolic excess velocity
    const vInfSq = vInf * vInf;

    // Eccentricity of hyperbolic trajectory
    const e = 1 + rPeriapsis * vInfSq / planet.mu;

    // Semi-major axis (negative for hyperbola)
    const a = -planet.mu / vInfSq;

    // Turn angle (total deflection)
    const delta = 2 * Math.asin(1 / e);
    const turnAngle = delta * 180 / Math.PI;

    // Periapsis velocity
    const vPeriapsis = Math.sqrt(vInfSq + 2 * planet.mu / rPeriapsis);

    // Impact parameter
    const b = rPeriapsis * Math.sqrt(1 + 2 * planet.mu / (rPeriapsis * vInfSq));

    // Velocity change in heliocentric frame
    // For maximum effect (retrograde flyby), the turn is in the direction of planet motion
    const deltaVHelio = 2 * vInf * Math.sin(delta / 2);

    // Generate hyperbolic trajectory
    const trajectory = [];
    const asymptoteAngle = Math.acos(-1 / e);

    // True anomaly range for hyperbola: -asymptoteAngle to +asymptoteAngle
    for (let i = 0; i <= config.numPoints; i++) {
        const theta = -asymptoteAngle * 0.95 + (i / config.numPoints) * asymptoteAngle * 1.9;
        const r = a * (1 - e * e) / (1 + e * Math.cos(theta));

        if (r > 0 && r < rPeriapsis * 10) {
            trajectory.push({
                x: r * Math.cos(theta),
                y: r * Math.sin(theta),
                r,
                theta: theta * 180 / Math.PI
            });
        }
    }

    // Velocity profile along trajectory
    const velocityProfile = [];
    for (let i = 0; i <= config.numPoints; i++) {
        const theta = -asymptoteAngle * 0.95 + (i / config.numPoints) * asymptoteAngle * 1.9;
        const r = a * (1 - e * e) / (1 + e * Math.cos(theta));

        if (r > 0 && r < rPeriapsis * 10) {
            const v = Math.sqrt(planet.mu * (2/r - 1/a));
            velocityProfile.push({
                theta: theta * 180 / Math.PI,
                v,
                r: r - planet.radius  // altitude
            });
        }
    }

    return {
        e, a, turnAngle,
        vPeriapsis,
        vInfOut: vInf,  // Magnitude unchanged
        deltaVHelio,
        b,
        rPeriapsis,
        trajectory,
        velocityProfile,
        asymptoteAngle: asymptoteAngle * 180 / Math.PI
    };
}

function renderGravityAssistFigures(container, result, config, planet, useSpice = false) {
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

    // Figure 1: Flyby trajectory
    const title1 = useSpice
        ? `${config.planet} Flyby Trajectory (SPICE Parameters)`
        : `${config.planet} Flyby Trajectory (Planet-centered)`;
    createFigure(wrapper, title1, chartWidth, chartHeight, (svg, w, h) => {
        renderFlybyTrajectory(svg, result, planet, w, h);
    });

    // Figure 2: Velocity profile
    createFigure(wrapper, 'Velocity Profile Along Trajectory', chartWidth, Math.max(180, chartHeight * 0.6), (svg, w, h) => {
        renderVelocityProfile(svg, result.velocityProfile, w, h);
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

function renderFlybyTrajectory(svg, result, planet, width, height) {
    const margin = { top: 20, right: 20, bottom: 40, left: 50 };

    const maxR = result.rPeriapsis * 5;

    const x = d3.scaleLinear()
        .domain([-maxR * 0.5, maxR * 1.5])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-maxR, maxR])
        .range([height - margin.bottom, margin.top]);

    // Draw planet
    const planetScreenRadius = Math.abs(x(planet.radius) - x(0));
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', Math.max(5, planetScreenRadius))
        .attr('fill', planet.color);

    // Draw hyperbolic trajectory
    svg.append('path')
        .datum(result.trajectory)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

    // Draw asymptotes (dashed)
    const asymAngle = result.asymptoteAngle * Math.PI / 180;

    svg.append('line')
        .attr('x1', x(-maxR * 0.5 * Math.cos(asymAngle)))
        .attr('y1', y(-maxR * 0.5 * Math.sin(asymAngle)))
        .attr('x2', x(0))
        .attr('y2', y(0))
        .attr('stroke', 'var(--text-dim)')
        .attr('stroke-width', 1)
        .attr('stroke-dasharray', '4,4');

    svg.append('line')
        .attr('x1', x(0))
        .attr('y1', y(0))
        .attr('x2', x(maxR * Math.cos(asymAngle)))
        .attr('y2', y(maxR * Math.sin(asymAngle)))
        .attr('stroke', 'var(--text-dim)')
        .attr('stroke-width', 1)
        .attr('stroke-dasharray', '4,4');

    // Periapsis marker
    svg.append('circle')
        .attr('cx', x(result.rPeriapsis))
        .attr('cy', y(0))
        .attr('r', 6)
        .attr('fill', 'var(--red)');

    svg.append('text')
        .attr('x', x(result.rPeriapsis) + 10)
        .attr('y', y(0) - 10)
        .attr('fill', 'var(--red)')
        .attr('font-size', '10px')
        .text('Periapsis');

    // Velocity vectors (incoming and outgoing)
    const arrowSize = maxR * 0.15;

    // Incoming V∞
    svg.append('line')
        .attr('x1', x(-maxR * 0.3))
        .attr('y1', y(result.b))
        .attr('x2', x(-maxR * 0.3 + arrowSize))
        .attr('y2', y(result.b))
        .attr('stroke', 'var(--green)')
        .attr('stroke-width', 3)
        .attr('marker-end', 'url(#arrowGreen)');

    // Outgoing V∞
    const outAngle = result.turnAngle * Math.PI / 180;
    svg.append('line')
        .attr('x1', x(maxR * 0.5))
        .attr('y1', y(-result.b * 0.5))
        .attr('x2', x(maxR * 0.5 + arrowSize * Math.cos(outAngle)))
        .attr('y2', y(-result.b * 0.5 + arrowSize * Math.sin(outAngle)))
        .attr('stroke', 'var(--orange)')
        .attr('stroke-width', 3);

    // Arrow marker
    svg.append('defs').append('marker')
        .attr('id', 'arrowGreen')
        .attr('viewBox', '0 -5 10 10')
        .attr('refX', 8)
        .attr('refY', 0)
        .attr('markerWidth', 6)
        .attr('markerHeight', 6)
        .attr('orient', 'auto')
        .append('path')
        .attr('d', 'M0,-5L10,0L0,5')
        .attr('fill', 'var(--green)');

    // Labels
    svg.append('text')
        .attr('x', x(-maxR * 0.25))
        .attr('y', y(result.b) - 10)
        .attr('fill', 'var(--green)')
        .attr('font-size', '10px')
        .text('V∞ in');

    // Turn angle arc
    svg.append('text')
        .attr('x', width - margin.right - 80)
        .attr('y', margin.top + 20)
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '11px')
        .text(`Turn: ${result.turnAngle.toFixed(1)}°`);

    svg.append('text')
        .attr('x', width - margin.right - 80)
        .attr('y', margin.top + 38)
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '11px')
        .text(`ΔV: ${result.deltaVHelio.toFixed(2)} km/s`);

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

function renderVelocityProfile(svg, data, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain(d3.extent(data, d => d.theta))
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([0, d3.max(data, d => d.v) * 1.1])
        .range([height - margin.bottom, margin.top]);

    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.theta)).y(d => y(d.v)));

    // Mark periapsis (theta = 0)
    svg.append('line')
        .attr('x1', x(0))
        .attr('x2', x(0))
        .attr('y1', margin.top)
        .attr('y2', height - margin.bottom)
        .attr('stroke', 'var(--red)')
        .attr('stroke-width', 1)
        .attr('stroke-dasharray', '4,4');

    svg.append('text')
        .attr('x', x(0) + 5)
        .attr('y', margin.top + 15)
        .attr('fill', 'var(--red)')
        .attr('font-size', '9px')
        .text('Periapsis');

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
        .text('True Anomaly [deg]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Velocity [km/s]');
}
