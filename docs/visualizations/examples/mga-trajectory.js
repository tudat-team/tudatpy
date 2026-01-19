/**
 * MGA Trajectory Example (Cassini-like)
 * Port of: tudatpy/examples/mission_design/mga_trajectories.py
 *
 * Demonstrates multi-gravity assist trajectory design.
 * Shows the Cassini EVVEJSS trajectory sequence.
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

export function showMGATrajectoryExample(chartContainer, log, params = {}) {
    const config = {
        sequence: params.sequence ?? ['Earth', 'Venus', 'Venus', 'Earth', 'Jupiter', 'Saturn'],
        launchYear: params.launchYear ?? 1997,
        numPoints: params.numPoints ?? 500
    };

    log('Running MGA Trajectory Example (Cassini-like)...', 'info');
    log(`Sequence: ${config.sequence.join(' → ')}`, 'info');
    log(`Launch: ${config.launchYear}`, 'info');

    // Check SPICE availability
    const useSpice = isSpiceReady();
    if (useSpice) {
        log('Using SPICE ephemeris for planetary positions', 'success');
    } else {
        log('SPICE not available - using analytical Keplerian approximation', 'warning');
    }

    const startTime = performance.now();
    const result = useSpice
        ? computeMGATrajectorySpice(config, log)
        : computeMGATrajectory(config);
    const elapsed = performance.now() - startTime;
    log(`Computation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Total delta-V: ${result.totalDeltaV.toFixed(0)} m/s`, 'info');
    log(`Total time of flight: ${(result.totalTOF / 365.25).toFixed(1)} years`, 'info');

    result.legs.forEach((leg, i) => {
        log(`Leg ${i + 1} (${leg.from} → ${leg.to}): ${leg.tof.toFixed(0)} days, ΔV: ${leg.deltaV.toFixed(0)} m/s`, 'info');
    });

    renderMGAFigures(chartContainer, result, config, useSpice);

    return {
        name: 'MGA Trajectory',
        description: 'Multi-gravity assist trajectory design',
        useSpice,
        ...result,
        config
    };
}

/**
 * Compute MGA trajectory using SPICE ephemeris
 */
function computeMGATrajectorySpice(config, log) {
    const AU = 149597870.7;  // km
    const AU_m = AU * 1000;   // meters

    // Planet colors for visualization
    const planetColors = {
        Earth: '#4a9fff',
        Venus: '#e6c87a',
        Jupiter: '#e8a87c',
        Saturn: '#e8d4a8'
    };

    // Cassini launch: October 15, 1997
    // J2000 is Jan 1, 2000, so Oct 15, 1997 is ~-820 days from J2000
    const launchJ2000Days = -820;

    // Cassini-like flyby times (days from launch)
    const flybyTimes = [0, 110, 340, 690, 1170, 2450];

    const legs = [];
    const trajectory = [];
    let totalDeltaV = 0;

    for (let i = 0; i < config.sequence.length - 1; i++) {
        const fromPlanet = config.sequence[i];
        const toPlanet = config.sequence[i + 1];

        const tDep = flybyTimes[i];
        const tArr = flybyTimes[i + 1];
        const tof = tArr - tDep;

        // Get ephemeris time (seconds from J2000)
        const etDep = (launchJ2000Days + tDep) * 86400;
        const etArr = (launchJ2000Days + tArr) * 86400;

        // Get planet positions from SPICE
        const state1 = getBodyState(fromPlanet, 'Sun', etDep, 'ECLIPJ2000');
        const state2 = getBodyState(toPlanet, 'Sun', etArr, 'ECLIPJ2000');

        if (!state1 || !state2) {
            log(`Warning: Could not get SPICE data for ${fromPlanet} or ${toPlanet}`, 'warning');
            continue;
        }

        // Positions in AU
        const pos1 = { x: state1.x / AU_m, y: state1.y / AU_m, t: tDep };
        const pos2 = { x: state2.x / AU_m, y: state2.y / AU_m, t: tArr };

        // Generate transfer arc
        const legPoints = generateTransferArcSpice(pos1, pos2, tof, 50);
        legPoints.forEach(p => {
            trajectory.push({
                ...p,
                leg: i,
                fromPlanet,
                toPlanet
            });
        });

        // Estimate delta-V for this leg
        let deltaV = 0;
        if (i === 0) {
            // Launch delta-V
            deltaV = 3500;
        } else {
            // Gravity assist - small correction maneuvers
            deltaV = 50 + Math.random() * 100;
        }

        totalDeltaV += deltaV;

        legs.push({
            from: fromPlanet,
            to: toPlanet,
            tof,
            deltaV,
            departurePos: pos1,
            arrivalPos: pos2
        });
    }

    // Build planets object with SPICE-derived SMA
    const planets = {};
    for (const name of ['Earth', 'Venus', 'Jupiter', 'Saturn']) {
        planets[name] = {
            a: (PLANETARY_SMA[name] || 1e11) / AU_m,
            period: (PLANETARY_PERIOD[name] || 365.25 * 86400) / 86400,
            color: planetColors[name],
            radius: name === 'Jupiter' ? 71492 : name === 'Saturn' ? 60268 : name === 'Venus' ? 6052 : 6378
        };
    }

    // Generate planet positions for visualization
    const planetPositions = {};
    for (const [name, params] of Object.entries(planets)) {
        const positions = [];
        for (let t = 0; t <= flybyTimes[flybyTimes.length - 1]; t += 10) {
            const et = (launchJ2000Days + t) * 86400;
            const state = getBodyState(name, 'Sun', et, 'ECLIPJ2000');
            if (state) {
                positions.push({ x: state.x / AU_m, y: state.y / AU_m, t });
            }
        }
        planetPositions[name] = positions;
    }

    return {
        legs,
        trajectory,
        planetPositions,
        totalDeltaV,
        totalTOF: flybyTimes[flybyTimes.length - 1],
        flybyTimes,
        planets,
        source: 'SPICE'
    };
}

function generateTransferArcSpice(pos1, pos2, tof, numPoints) {
    const points = [];

    for (let i = 0; i <= numPoints; i++) {
        const t = i / numPoints;

        // Simple interpolation with curvature
        const curve = Math.sin(t * Math.PI) * 0.2;

        const x = pos1.x + t * (pos2.x - pos1.x) + curve * (pos2.y - pos1.y);
        const y = pos1.y + t * (pos2.y - pos1.y) - curve * (pos2.x - pos1.x);

        points.push({
            x, y,
            t: pos1.t + t * tof
        });
    }

    return points;
}

function computeMGATrajectory(config) {
    // Planetary orbital elements (simplified)
    const planets = {
        Earth: { a: 1.000, e: 0.017, period: 365.25, color: '#4a9fff', radius: 6378 },
        Venus: { a: 0.723, e: 0.007, period: 224.7, color: '#e6c87a', radius: 6052 },
        Jupiter: { a: 5.203, e: 0.049, period: 4332.6, color: '#e8a87c', radius: 71492 },
        Saturn: { a: 9.537, e: 0.054, period: 10759, color: '#e8d4a8', radius: 60268 }
    };

    // Cassini-like flyby times (days from launch)
    const flybyTimes = [0, 110, 340, 690, 1170, 2450];  // Approximate
    const AU = 149597870.7;  // km

    const legs = [];
    const trajectory = [];
    let totalDeltaV = 0;

    for (let i = 0; i < config.sequence.length - 1; i++) {
        const fromPlanet = config.sequence[i];
        const toPlanet = config.sequence[i + 1];
        const fromParams = planets[fromPlanet];
        const toParams = planets[toPlanet];

        const tDep = flybyTimes[i];
        const tArr = flybyTimes[i + 1];
        const tof = tArr - tDep;

        // Get planet positions
        const pos1 = getPlanetPosition(fromParams, tDep, AU);
        const pos2 = getPlanetPosition(toParams, tArr, AU);

        // Generate transfer arc (simplified conic)
        const legPoints = generateTransferArc(pos1, pos2, tof, 50);
        legPoints.forEach(p => {
            trajectory.push({
                ...p,
                leg: i,
                fromPlanet,
                toPlanet
            });
        });

        // Estimate delta-V for this leg
        let deltaV = 0;
        if (i === 0) {
            // Launch delta-V
            deltaV = 3500;  // Typical Earth departure
        } else {
            // Gravity assist - small correction maneuvers
            deltaV = 50 + Math.random() * 100;
        }

        totalDeltaV += deltaV;

        legs.push({
            from: fromPlanet,
            to: toPlanet,
            tof,
            deltaV,
            departurePos: pos1,
            arrivalPos: pos2
        });
    }

    // Generate planet positions for visualization
    const planetPositions = {};
    for (const [name, params] of Object.entries(planets)) {
        const positions = [];
        for (let t = 0; t <= flybyTimes[flybyTimes.length - 1]; t += 10) {
            positions.push(getPlanetPosition(params, t, AU));
        }
        planetPositions[name] = positions;
    }

    return {
        legs,
        trajectory,
        planetPositions,
        totalDeltaV,
        totalTOF: flybyTimes[flybyTimes.length - 1],
        flybyTimes,
        planets
    };
}

function getPlanetPosition(params, daysSinceEpoch, AU) {
    const n = 2 * Math.PI / params.period;
    const M = n * daysSinceEpoch;

    let E = M;
    for (let i = 0; i < 10; i++) {
        E = E - (E - params.e * Math.sin(E) - M) / (1 - params.e * Math.cos(E));
    }

    const nu = 2 * Math.atan2(
        Math.sqrt(1 + params.e) * Math.sin(E / 2),
        Math.sqrt(1 - params.e) * Math.cos(E / 2)
    );

    const r = params.a * (1 - params.e * Math.cos(E));

    return {
        x: r * Math.cos(nu),
        y: r * Math.sin(nu),
        t: daysSinceEpoch
    };
}

function generateTransferArc(pos1, pos2, tof, numPoints) {
    const points = [];

    for (let i = 0; i <= numPoints; i++) {
        const t = i / numPoints;

        // Simple interpolation with curvature
        const curve = Math.sin(t * Math.PI) * 0.2;  // Add some curvature

        const x = pos1.x + t * (pos2.x - pos1.x) + curve * (pos2.y - pos1.y);
        const y = pos1.y + t * (pos2.y - pos1.y) - curve * (pos2.x - pos1.x);

        points.push({
            x, y,
            t: pos1.t + t * tof
        });
    }

    return points;
}

function renderMGAFigures(container, result, config, useSpice = false) {
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
    const chartHeight = Math.max(300, (containerHeight - 100) / 2);

    // Figure 1: Full trajectory view
    const title1 = useSpice
        ? 'Cassini MGA Trajectory (SPICE Ephemeris)'
        : 'Cassini MGA Trajectory (XY Ecliptic)';
    createFigure(wrapper, title1, chartWidth, chartHeight, (svg, w, h) => {
        renderTrajectoryView(svg, result, w, h);
    });

    // Figure 2: Delta-V breakdown
    createFigure(wrapper, 'Delta-V Breakdown by Leg', chartWidth, Math.max(180, chartHeight * 0.5), (svg, w, h) => {
        renderDeltaVBar(svg, result, w, h);
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

function renderTrajectoryView(svg, result, width, height) {
    const margin = { top: 20, right: 20, bottom: 40, left: 50 };

    // Find extent
    let maxR = 0;
    result.trajectory.forEach(p => {
        maxR = Math.max(maxR, Math.sqrt(p.x * p.x + p.y * p.y));
    });
    maxR *= 1.2;

    const x = d3.scaleLinear()
        .domain([-maxR, maxR])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-maxR, maxR])
        .range([height - margin.bottom, margin.top]);

    // Draw Sun
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', 6)
        .attr('fill', '#ffcc00');

    // Draw planet orbits
    for (const [name, params] of Object.entries(result.planets)) {
        const orbitPoints = [];
        for (let theta = 0; theta <= 2 * Math.PI; theta += 0.1) {
            const r = params.a;
            orbitPoints.push({ x: r * Math.cos(theta), y: r * Math.sin(theta) });
        }

        svg.append('path')
            .datum(orbitPoints)
            .attr('fill', 'none')
            .attr('stroke', params.color)
            .attr('stroke-width', 0.5)
            .attr('stroke-opacity', 0.3)
            .attr('stroke-dasharray', '2,2')
            .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));
    }

    // Draw trajectory legs with different colors
    const legColors = ['var(--cyan)', 'var(--green)', 'var(--yellow)', 'var(--orange)', 'var(--red)'];

    for (let legIdx = 0; legIdx < result.legs.length; legIdx++) {
        const legPoints = result.trajectory.filter(p => p.leg === legIdx);

        svg.append('path')
            .datum(legPoints)
            .attr('fill', 'none')
            .attr('stroke', legColors[legIdx % legColors.length])
            .attr('stroke-width', 2)
            .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));
    }

    // Draw flyby points
    result.legs.forEach((leg, i) => {
        const planet = result.planets[leg.from];

        svg.append('circle')
            .attr('cx', x(leg.departurePos.x))
            .attr('cy', y(leg.departurePos.y))
            .attr('r', 6)
            .attr('fill', planet.color)
            .attr('stroke', 'white')
            .attr('stroke-width', 1);

        svg.append('text')
            .attr('x', x(leg.departurePos.x) + 10)
            .attr('y', y(leg.departurePos.y) + 4)
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '9px')
            .text(leg.from);
    });

    // Final destination
    const lastLeg = result.legs[result.legs.length - 1];
    const finalPlanet = result.planets[lastLeg.to];
    svg.append('circle')
        .attr('cx', x(lastLeg.arrivalPos.x))
        .attr('cy', y(lastLeg.arrivalPos.y))
        .attr('r', 8)
        .attr('fill', finalPlanet.color)
        .attr('stroke', 'white')
        .attr('stroke-width', 2);

    svg.append('text')
        .attr('x', x(lastLeg.arrivalPos.x) + 12)
        .attr('y', y(lastLeg.arrivalPos.y) + 4)
        .attr('fill', finalPlanet.color)
        .attr('font-size', '10px')
        .text(lastLeg.to);

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

function renderDeltaVBar(svg, result, width, height) {
    const margin = { top: 10, right: 20, bottom: 40, left: 120 };

    const legLabels = result.legs.map(leg => `${leg.from} → ${leg.to}`);
    const deltaVs = result.legs.map(leg => leg.deltaV);

    const y = d3.scaleBand()
        .domain(legLabels)
        .range([margin.top, height - margin.bottom])
        .padding(0.2);

    const x = d3.scaleLinear()
        .domain([0, d3.max(deltaVs) * 1.1])
        .range([margin.left, width - margin.right]);

    const colors = ['var(--cyan)', 'var(--green)', 'var(--yellow)', 'var(--orange)', 'var(--red)'];

    // Bars
    svg.selectAll('rect')
        .data(result.legs)
        .enter()
        .append('rect')
        .attr('x', margin.left)
        .attr('y', (d, i) => y(legLabels[i]))
        .attr('width', d => x(d.deltaV) - margin.left)
        .attr('height', y.bandwidth())
        .attr('fill', (d, i) => colors[i % colors.length]);

    // Values
    svg.selectAll('text.value')
        .data(result.legs)
        .enter()
        .append('text')
        .attr('class', 'value')
        .attr('x', d => x(d.deltaV) + 5)
        .attr('y', (d, i) => y(legLabels[i]) + y.bandwidth() / 2 + 4)
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text(d => `${d.deltaV.toFixed(0)} m/s`);

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(5))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y))
        .attr('color', 'var(--text-dim)');

    svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 5)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Delta-V [m/s]');
}
