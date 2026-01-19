/**
 * Low-Thrust Earth-Mars Porkchop Plot Example
 * Port of: tudatpy/examples/mission_design/low_thrust_earth_mars_transfer_window.py
 *
 * Demonstrates porkchop plot generation for low-thrust transfers using
 * hodographic shaping. Shows ΔV requirements across departure/arrival windows.
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

export function showLowThrustPorkchopExample(chartContainer, log, params = {}) {
    const config = {
        departureDays: params.departureDays ?? 60,    // Grid size for departure
        arrivalDays: params.arrivalDays ?? 60,        // Grid size for arrival
        departureStart: params.departureStart ?? 0,   // Days from epoch
        arrivalStart: params.arrivalStart ?? 200,     // Min TOF
        arrivalEnd: params.arrivalEnd ?? 600          // Max TOF
    };

    log('Running Low-Thrust Porkchop Example...', 'info');
    log('Computing ΔV grid for Earth-Mars low-thrust transfers...', 'info');
    log(`Grid: ${config.departureDays} × ${config.arrivalDays} points`, 'info');

    // Check SPICE availability
    const useSpice = isSpiceReady();
    if (useSpice) {
        log('Using SPICE ephemeris for planetary positions', 'success');
    } else {
        log('SPICE not available - using analytical Keplerian approximation', 'warning');
    }

    const startTime = performance.now();
    const result = useSpice
        ? computeLowThrustPorkchopSpice(config, log)
        : computeLowThrustPorkchop(config);
    const elapsed = performance.now() - startTime;
    log(`Computation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Min ΔV: ${result.minDeltaV.toFixed(2)} km/s`, 'info');
    log(`Optimal departure: day ${result.optimalDeparture.toFixed(0)}`, 'info');
    log(`Optimal TOF: ${result.optimalTOF.toFixed(0)} days`, 'info');

    renderLowThrustFigures(chartContainer, result, config, useSpice);

    return {
        name: 'Low-Thrust Porkchop',
        description: 'Low-thrust Earth-Mars transfer windows',
        useSpice,
        ...result,
        config
    };
}

/**
 * Compute low-thrust porkchop plot using SPICE ephemeris
 */
function computeLowThrustPorkchopSpice(config, log) {
    const AU = 1.496e11;
    const AU_m = AU;  // meters
    const mu_sun = 1.327e20;

    // Grid setup
    const numDeparture = config.departureDays;
    const numArrival = config.arrivalDays;

    const departureWindow = 500;  // days
    const departureDates = [];
    for (let i = 0; i < numDeparture; i++) {
        departureDates.push(config.departureStart + (i / (numDeparture - 1)) * departureWindow);
    }

    const arrivalDates = [];
    for (let i = 0; i < numArrival; i++) {
        const tof = config.arrivalStart + (i / (numArrival - 1)) * (config.arrivalEnd - config.arrivalStart);
        arrivalDates.push(tof);
    }

    // Compute ΔV grid using SPICE
    const grid = [];
    let minDeltaV = Infinity;
    let optimalDeparture = 0;
    let optimalTOF = 0;

    for (let i = 0; i < numDeparture; i++) {
        const row = [];
        for (let j = 0; j < numArrival; j++) {
            const departureDay = departureDates[i];
            const tof = arrivalDates[j];
            const arrivalDay = departureDay + tof;

            // Convert to ephemeris time (days from J2000 to seconds)
            const etDep = departureDay * 86400;
            const etArr = arrivalDay * 86400;

            // Get planet positions from SPICE
            const earthState = getBodyState('Earth', 'Sun', etDep, 'ECLIPJ2000');
            const marsState = getBodyState('Mars', 'Sun', etArr, 'ECLIPJ2000');

            let deltaV;
            if (earthState && marsState) {
                const earthPos = {
                    x: earthState.x,
                    y: earthState.y,
                    r: Math.sqrt(earthState.x**2 + earthState.y**2 + earthState.z**2),
                    v: Math.sqrt(earthState.vx**2 + earthState.vy**2 + earthState.vz**2)
                };
                const marsPos = {
                    x: marsState.x,
                    y: marsState.y,
                    r: Math.sqrt(marsState.x**2 + marsState.y**2 + marsState.z**2),
                    v: Math.sqrt(marsState.vx**2 + marsState.vy**2 + marsState.vz**2)
                };
                deltaV = computeHodographicDeltaVSpice(earthPos, marsPos, tof, mu_sun);
            } else {
                deltaV = 100;  // Fallback
            }

            row.push({
                departure: departureDay,
                tof: tof,
                deltaV: deltaV
            });

            if (deltaV < minDeltaV) {
                minDeltaV = deltaV;
                optimalDeparture = departureDay;
                optimalTOF = tof;
            }
        }
        grid.push(row);
    }

    // Generate optimal trajectory using SPICE positions
    const etOptDep = optimalDeparture * 86400;
    const etOptArr = (optimalDeparture + optimalTOF) * 86400;
    const optEarthState = getBodyState('Earth', 'Sun', etOptDep, 'ECLIPJ2000');
    const optMarsState = getBodyState('Mars', 'Sun', etOptArr, 'ECLIPJ2000');

    let optimalTrajectory = [];
    if (optEarthState && optMarsState) {
        const earthPos = { x: optEarthState.x, y: optEarthState.y, r: Math.sqrt(optEarthState.x**2 + optEarthState.y**2) };
        const marsPos = { x: optMarsState.x, y: optMarsState.y, r: Math.sqrt(optMarsState.x**2 + optMarsState.y**2) };
        optimalTrajectory = generateLowThrustTrajectorySpice(earthPos, marsPos, optimalTOF, mu_sun, AU);
    }

    return {
        grid,
        departureDates,
        arrivalDates,
        minDeltaV,
        optimalDeparture,
        optimalTOF,
        optimalTrajectory,
        source: 'SPICE'
    };
}

function computeHodographicDeltaVSpice(pos1, pos2, tofDays, mu) {
    const tof = tofDays * 86400;

    const r1 = pos1.r;
    const r2 = pos2.r;

    // Transfer angle
    const cosTheta = (pos1.x * pos2.x + pos1.y * pos2.y) / (r1 * r2);
    const theta = Math.acos(Math.max(-1, Math.min(1, cosTheta)));

    // Characteristic velocity
    const a_transfer = (r1 + r2) / 2;
    const v_inf = Math.sqrt(mu / a_transfer);

    // Time factor
    const tau = tof / (2 * Math.PI * Math.sqrt(a_transfer**3 / mu));

    // Low-thrust approximation
    let deltaV = Math.abs(pos2.v - pos1.v) * 0.8;

    // Shaping contributions
    const radialComp = Math.abs(r2 - r1) / (r1 + r2) * v_inf * 0.5;
    const normalComp = Math.abs(Math.sin(theta)) * v_inf * 0.3;
    deltaV += radialComp + normalComp;

    // TOF penalties
    if (tofDays < 150) deltaV *= 1 + (150 - tofDays) / 50;
    if (tofDays > 500) deltaV *= 1 + (tofDays - 500) / 200;

    return deltaV / 1000;  // km/s
}

function generateLowThrustTrajectorySpice(pos1, pos2, tofDays, mu, AU) {
    const numPoints = 100;
    const trajectory = [];

    const r1 = pos1.r;
    const r2 = pos2.r;
    const angle1 = Math.atan2(pos1.y, pos1.x);
    const angle2 = Math.atan2(pos2.y, pos2.x);

    let totalAngle = angle2 - angle1;
    if (totalAngle < 0) totalAngle += 2 * Math.PI;
    totalAngle += 2 * Math.PI * 2;  // 2 revolutions for low-thrust

    for (let i = 0; i <= numPoints; i++) {
        const frac = i / numPoints;
        const angle = angle1 + frac * totalAngle;
        const r = r1 + frac * (r2 - r1);
        const oscillation = 0.05 * Math.sin(frac * Math.PI * 4) * (r2 - r1);

        trajectory.push({
            x: (r + oscillation) * Math.cos(angle) / AU,
            y: (r + oscillation) * Math.sin(angle) / AU
        });
    }

    return trajectory;
}

function computeLowThrustPorkchop(config) {
    // Planetary data (AU and days)
    const AU = 1.496e11;
    const mu_sun = 1.327e20;

    const earthOrbit = { a: 1.0 * AU, e: 0.0167, period: 365.25 };
    const marsOrbit = { a: 1.524 * AU, e: 0.0934, period: 687 };

    // Grid setup
    const numDeparture = config.departureDays;
    const numArrival = config.arrivalDays;

    const departureWindow = 500;  // days
    const departureDates = [];
    for (let i = 0; i < numDeparture; i++) {
        departureDates.push(config.departureStart + (i / (numDeparture - 1)) * departureWindow);
    }

    const arrivalDates = [];
    for (let i = 0; i < numArrival; i++) {
        const tof = config.arrivalStart + (i / (numArrival - 1)) * (config.arrivalEnd - config.arrivalStart);
        arrivalDates.push(tof);
    }

    // Compute ΔV grid
    const grid = [];
    let minDeltaV = Infinity;
    let optimalDeparture = 0;
    let optimalTOF = 0;

    for (let i = 0; i < numDeparture; i++) {
        const row = [];
        for (let j = 0; j < numArrival; j++) {
            const departureDay = departureDates[i];
            const tof = arrivalDates[j];
            const arrivalDay = departureDay + tof;

            // Get planet positions
            const earthPos = getPlanetPosition(departureDay, earthOrbit);
            const marsPos = getPlanetPosition(arrivalDay, marsOrbit);

            // Compute low-thrust ΔV using hodographic shaping approximation
            const deltaV = computeHodographicDeltaV(earthPos, marsPos, tof, mu_sun);

            row.push({
                departure: departureDay,
                tof: tof,
                deltaV: deltaV
            });

            if (deltaV < minDeltaV) {
                minDeltaV = deltaV;
                optimalDeparture = departureDay;
                optimalTOF = tof;
            }
        }
        grid.push(row);
    }

    // Generate optimal trajectory
    const optimalEarthPos = getPlanetPosition(optimalDeparture, earthOrbit);
    const optimalMarsPos = getPlanetPosition(optimalDeparture + optimalTOF, marsOrbit);
    const optimalTrajectory = generateLowThrustTrajectory(
        optimalEarthPos, optimalMarsPos, optimalTOF, mu_sun
    );

    return {
        grid,
        departureDates,
        arrivalDates,
        minDeltaV,
        optimalDeparture,
        optimalTOF,
        optimalTrajectory
    };
}

function getPlanetPosition(daysSinceEpoch, orbit) {
    const AU = 1.496e11;
    const meanAnomaly = 2 * Math.PI * daysSinceEpoch / orbit.period;

    // Solve Kepler's equation (simplified for low eccentricity)
    let E = meanAnomaly;
    for (let i = 0; i < 5; i++) {
        E = meanAnomaly + orbit.e * Math.sin(E);
    }

    const trueAnomaly = 2 * Math.atan2(
        Math.sqrt(1 + orbit.e) * Math.sin(E / 2),
        Math.sqrt(1 - orbit.e) * Math.cos(E / 2)
    );

    const r = orbit.a * (1 - orbit.e * Math.cos(E));

    return {
        x: r * Math.cos(trueAnomaly),
        y: r * Math.sin(trueAnomaly),
        r: r,
        v: Math.sqrt(1.327e20 * (2 / r - 1 / orbit.a))
    };
}

function computeHodographicDeltaV(pos1, pos2, tofDays, mu) {
    // Simplified hodographic shaping ΔV estimate
    const tof = tofDays * 86400;  // Convert to seconds
    const AU = 1.496e11;

    const r1 = pos1.r;
    const r2 = pos2.r;

    // Transfer angle
    const cosTheta = (pos1.x * pos2.x + pos1.y * pos2.y) / (r1 * r2);
    const theta = Math.acos(Math.max(-1, Math.min(1, cosTheta)));

    // Characteristic velocity (based on energy considerations)
    const a_transfer = (r1 + r2) / 2;
    const v_inf = Math.sqrt(mu / a_transfer);

    // Low-thrust specific: ΔV depends on transfer time and geometry
    const tau = tof / (2 * Math.PI * Math.sqrt(a_transfer * a_transfer * a_transfer / mu));

    // Hodographic shaping factor (simplified)
    const numRevs = 2;  // Multiple revolutions typical for low-thrust
    const shapeParam = Math.abs(theta + 2 * Math.PI * numRevs) / tau;

    // ΔV estimation for low-thrust
    let deltaV = Math.abs(pos2.v - pos1.v) * 0.8;  // Base ΔV

    // Add contributions from shaping
    const radialComp = Math.abs(r2 - r1) / (r1 + r2) * v_inf * 0.5;
    const normalComp = Math.abs(Math.sin(theta)) * v_inf * 0.3;

    deltaV += radialComp + normalComp;

    // Penalize very short or very long transfers
    if (tofDays < 150) {
        deltaV *= 1 + (150 - tofDays) / 50;
    }
    if (tofDays > 500) {
        deltaV *= 1 + (tofDays - 500) / 200;
    }

    // Convert to km/s
    return deltaV / 1000;
}

function generateLowThrustTrajectory(pos1, pos2, tofDays, mu) {
    const AU = 1.496e11;
    const numPoints = 100;
    const trajectory = [];

    const r1 = pos1.r;
    const r2 = pos2.r;
    const angle1 = Math.atan2(pos1.y, pos1.x);
    const angle2 = Math.atan2(pos2.y, pos2.x);

    // Handle angle wrapping for multiple revolutions
    let totalAngle = angle2 - angle1;
    if (totalAngle < 0) totalAngle += 2 * Math.PI;
    totalAngle += 2 * Math.PI * 2;  // 2 revolutions

    for (let i = 0; i <= numPoints; i++) {
        const frac = i / numPoints;

        // Spiral outward with hodographic shaping
        const angle = angle1 + frac * totalAngle;
        const r = r1 + frac * (r2 - r1);

        // Add some variation for realistic spiral
        const oscillation = 0.05 * Math.sin(frac * Math.PI * 4) * (r2 - r1);

        trajectory.push({
            x: (r + oscillation) * Math.cos(angle) / AU,
            y: (r + oscillation) * Math.sin(angle) / AU
        });
    }

    return trajectory;
}

function renderLowThrustFigures(container, result, config, useSpice = false) {
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

    // Figure 1: Porkchop plot (heatmap)
    const title1 = useSpice
        ? 'Low-Thrust Porkchop Plot (SPICE Ephemeris)'
        : 'Low-Thrust Porkchop Plot (ΔV)';
    createFigure(wrapper, title1, chartWidth, chartHeight, (svg, w, h) => {
        renderPorkchopHeatmap(svg, result, w, h);
    });

    // Figure 2: Optimal trajectory
    const title2 = useSpice
        ? 'Optimal Low-Thrust Trajectory (SPICE)'
        : 'Optimal Low-Thrust Trajectory';
    createFigure(wrapper, title2, chartWidth, chartHeight, (svg, w, h) => {
        renderOptimalTrajectory(svg, result, w, h);
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

function renderPorkchopHeatmap(svg, result, width, height) {
    const margin = { top: 10, right: 80, bottom: 40, left: 60 };

    const numX = result.grid.length;
    const numY = result.grid[0].length;

    const x = d3.scaleLinear()
        .domain([result.departureDates[0], result.departureDates[numX - 1]])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([result.arrivalDates[0], result.arrivalDates[numY - 1]])
        .range([height - margin.bottom, margin.top]);

    // Find ΔV range
    let minDV = Infinity, maxDV = 0;
    result.grid.forEach(row => {
        row.forEach(cell => {
            if (cell.deltaV < 50) {  // Filter outliers
                minDV = Math.min(minDV, cell.deltaV);
                maxDV = Math.max(maxDV, cell.deltaV);
            }
        });
    });

    // Color scale
    const colorScale = d3.scaleSequential(d3.interpolateViridis)
        .domain([maxDV, minDV]);  // Reversed so low ΔV is bright

    // Draw heatmap cells
    const cellWidth = (width - margin.left - margin.right) / numX;
    const cellHeight = (height - margin.top - margin.bottom) / numY;

    for (let i = 0; i < numX; i++) {
        for (let j = 0; j < numY; j++) {
            const cell = result.grid[i][j];
            const dv = Math.min(cell.deltaV, maxDV);  // Clamp outliers

            svg.append('rect')
                .attr('x', x(cell.departure) - cellWidth / 2)
                .attr('y', y(cell.tof) - cellHeight / 2)
                .attr('width', cellWidth)
                .attr('height', cellHeight)
                .attr('fill', colorScale(dv));
        }
    }

    // Mark optimal point
    svg.append('circle')
        .attr('cx', x(result.optimalDeparture))
        .attr('cy', y(result.optimalTOF))
        .attr('r', 6)
        .attr('fill', 'none')
        .attr('stroke', 'white')
        .attr('stroke-width', 2);

    // Color bar
    const colorBarWidth = 15;
    const colorBarHeight = height - margin.top - margin.bottom;
    const colorBarX = width - margin.right + 20;

    const colorBarScale = d3.scaleLinear()
        .domain([minDV, maxDV])
        .range([height - margin.bottom, margin.top]);

    const defs = svg.append('defs');
    const gradient = defs.append('linearGradient')
        .attr('id', 'colorbar-gradient')
        .attr('x1', '0%').attr('y1', '100%')
        .attr('x2', '0%').attr('y2', '0%');

    for (let i = 0; i <= 10; i++) {
        const t = i / 10;
        gradient.append('stop')
            .attr('offset', `${t * 100}%`)
            .attr('stop-color', colorScale(minDV + t * (maxDV - minDV)));
    }

    svg.append('rect')
        .attr('x', colorBarX)
        .attr('y', margin.top)
        .attr('width', colorBarWidth)
        .attr('height', colorBarHeight)
        .attr('fill', 'url(#colorbar-gradient)');

    svg.append('g')
        .attr('transform', `translate(${colorBarX + colorBarWidth}, 0)`)
        .call(d3.axisRight(colorBarScale).ticks(5).tickFormat(d => d.toFixed(1)))
        .attr('color', 'var(--text-dim)');

    svg.append('text')
        .attr('x', colorBarX + colorBarWidth / 2)
        .attr('y', margin.top - 5)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '8px')
        .text('ΔV km/s');

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(6))
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
        .text('Departure [days from epoch]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Time of Flight [days]');
}

function renderOptimalTrajectory(svg, result, width, height) {
    const margin = { top: 20, right: 80, bottom: 40, left: 60 };

    // Scale to fit trajectory
    const scale = 2.0;  // AU

    const x = d3.scaleLinear()
        .domain([-scale, scale])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-scale, scale])
        .range([height - margin.bottom, margin.top]);

    // Draw planet orbits
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', x(1.0) - x(0))
        .attr('fill', 'none')
        .attr('stroke', '#4488ff')
        .attr('stroke-width', 1)
        .attr('stroke-dasharray', '3,3')
        .attr('opacity', 0.5);

    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', x(1.524) - x(0))
        .attr('fill', 'none')
        .attr('stroke', '#ff6644')
        .attr('stroke-width', 1)
        .attr('stroke-dasharray', '3,3')
        .attr('opacity', 0.5);

    // Draw Sun
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', 5)
        .attr('fill', '#ffcc00');

    // Draw trajectory
    svg.append('path')
        .datum(result.optimalTrajectory)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

    // Mark start and end
    const start = result.optimalTrajectory[0];
    const end = result.optimalTrajectory[result.optimalTrajectory.length - 1];

    svg.append('circle')
        .attr('cx', x(start.x))
        .attr('cy', y(start.y))
        .attr('r', 5)
        .attr('fill', '#4488ff');

    svg.append('circle')
        .attr('cx', x(end.x))
        .attr('cy', y(end.y))
        .attr('r', 5)
        .attr('fill', '#ff6644');

    // Legend
    const legend = svg.append('g')
        .attr('transform', `translate(${width - margin.right + 10}, ${margin.top})`);

    const legendItems = [
        { name: 'Earth', color: '#4488ff' },
        { name: 'Mars', color: '#ff6644' },
        { name: 'Sun', color: '#ffcc00' }
    ];

    legendItems.forEach((item, i) => {
        legend.append('circle')
            .attr('cx', 5)
            .attr('cy', i * 16)
            .attr('r', 4)
            .attr('fill', item.color);

        legend.append('text')
            .attr('x', 14)
            .attr('y', i * 16 + 4)
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '9px')
            .text(item.name);
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
        .text('X [AU]');
}
