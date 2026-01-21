/**
 * Earth-Mars Transfer Window Example
 * Port of: tudatpy/examples/mission_design/earth_mars_transfer_window.py
 *
 * Demonstrates porkchop plot generation for Earth-Mars transfers.
 * Shows delta-V contours as function of departure and arrival dates.
 *
 * Uses binary SPICE kernels (via CALCEPH) when available for accurate planetary positions.
 * Falls back to precomputed JSON ephemeris or analytical approximations.
 */

import {
    isHighAccuracyEphemerisAvailable,
    isCalcephAvailable,
    loadPlanetaryEphemeris,
    getBodyState,
    jdToEt,
    PLANETARY_GM,
    loadPrecomputedEphemerisFromUrl
} from '../shared/spice-utils.js';

export async function showEarthMarsTransferExample(chartContainer, log, params = {}) {
    const config = {
        departureWindowDays: params.departureWindowDays ?? 400,  // days from J2000
        arrivalWindowDays: params.arrivalWindowDays ?? 400,      // days
        // Start from day 10 to stay within safe interpolation bounds for precomputed ephemeris
        // (8th order Lagrange interpolation requires ~5 data points margin on each end)
        departureStart: params.departureStart ?? 10,             // days from reference
        arrivalStart: params.arrivalStart ?? 210,                // minimum flight time
        resolution: params.resolution ?? 30,                     // grid points per axis
        // Binary kernel URL - can be overridden to use different kernels
        spkUrl: params.spkUrl ?? './data/de432s.bsp'            // JPL DE432s (10MB, covers 1950-2050)
    };

    log('Running Earth-Mars Transfer Window Example...', 'info');
    log(`Departure window: ${config.departureWindowDays} days`, 'info');
    log(`Arrival window: ${config.arrivalWindowDays} days`, 'info');
    log(`Resolution: ${config.resolution}x${config.resolution}`, 'info');

    // Check if high-accuracy ephemeris is already loaded
    let useHighAccuracy = isHighAccuracyEphemerisAvailable('Earth', 'Sun', 'J2000') &&
                          isHighAccuracyEphemerisAvailable('Mars', 'Sun', 'J2000');
    let ephemerisSource = 'cached';

    if (!useHighAccuracy) {
        // Priority 1: Try loading binary SPICE kernel via CALCEPH
        if (isCalcephAvailable()) {
            log('Loading binary SPICE kernel (CALCEPH)...', 'info');
            const result = await loadPlanetaryEphemeris(config.spkUrl, ['Earth', 'Mars'], 'Sun', 'J2000');

            if (result.success && result.loaded.length === 2) {
                useHighAccuracy = true;
                ephemerisSource = 'binary SPK (CALCEPH)';
                log(`Loaded binary kernel: ${result.loaded.join(', ')}`, 'success');
            } else if (result.loaded.length > 0) {
                log(`Partial load from binary kernel: ${result.loaded.join(', ')}`, 'warning');
                if (result.failed.length > 0) {
                    log(`Failed to load: ${result.failed.join(', ')}`, 'warning');
                }
            } else {
                log('Binary kernel not available, trying precomputed JSON...', 'info');
            }
        }

        // Priority 2: Fall back to precomputed JSON ephemeris
        if (!useHighAccuracy) {
            log('Loading precomputed JSON ephemeris...', 'info');
            const earthLoaded = await loadPrecomputedEphemerisFromUrl('./data/ephemeris/earth_sun_j2000.json');
            const marsLoaded = await loadPrecomputedEphemerisFromUrl('./data/ephemeris/mars_sun_j2000.json');
            useHighAccuracy = earthLoaded && marsLoaded;
            if (useHighAccuracy) {
                ephemerisSource = 'precomputed JSON';
            }
        }
    }

    if (useHighAccuracy) {
        log(`Using high-accuracy ephemeris (${ephemerisSource})`, 'success');
    } else {
        log('Using analytical Keplerian approximation (no ephemeris data available)', 'warning');
    }

    const startTime = performance.now();
    const result = useHighAccuracy
        ? computeTransferWindowSpice(config, log)
        : computeTransferWindowAnalytical(config);
    const elapsed = performance.now() - startTime;

    // Update source with more specific info
    if (useHighAccuracy) {
        result.source = ephemerisSource;
    }

    log(`Computation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Minimum delta-V: ${result.minDeltaV.toFixed(0)} m/s`, 'info');
    log(`Optimal departure: Day ${result.optimalDeparture.toFixed(0)}`, 'info');
    log(`Optimal arrival: Day ${result.optimalArrival.toFixed(0)}`, 'info');
    log(`Flight time: ${(result.optimalArrival - result.optimalDeparture).toFixed(0)} days`, 'info');

    renderPorkchopPlot(chartContainer, result, config);

    return {
        name: 'Earth-Mars Transfer',
        description: 'Interplanetary transfer window analysis',
        useHighAccuracy,
        ephemerisSource,
        ...result,
        config
    };
}

/**
 * Compute transfer window using SPICE ephemeris
 */
function computeTransferWindowSpice(config, log) {
    const AU = 149597870.7;  // km
    const muSun = PLANETARY_GM.Sun / 1e9;  // km³/s² (convert from m³/s²)

    const departureDays = [];
    const arrivalDays = [];
    const deltaVGrid = [];

    const dDep = config.departureWindowDays / (config.resolution - 1);
    const dArr = config.arrivalWindowDays / (config.resolution - 1);

    let minDeltaV = Infinity;
    let optimalDeparture = 0;
    let optimalArrival = 0;

    for (let i = 0; i < config.resolution; i++) {
        const tDep = config.departureStart + i * dDep;
        departureDays.push(tDep);

        const row = [];
        for (let j = 0; j < config.resolution; j++) {
            const tArr = config.arrivalStart + j * dArr;
            if (i === 0) arrivalDays.push(tArr);

            const tof = (tArr - tDep) * 86400;  // seconds

            if (tof < 100 * 86400 || tof > 400 * 86400) {
                // Invalid flight time
                row.push(NaN);
                continue;
            }

            // Convert days to ephemeris time (seconds since J2000)
            const etDep = tDep * 86400;
            const etArr = tArr * 86400;

            // Get planet positions from SPICE (using J2000 frame to match precomputed ephemeris)
            const earthState = getBodyState('Earth', 'Sun', etDep, 'J2000');
            const marsState = getBodyState('Mars', 'Sun', etArr, 'J2000');

            if (!earthState || !marsState) {
                row.push(NaN);
                continue;
            }

            // Convert from meters to km
            const earthPos = {
                x: earthState.x / 1000,
                y: earthState.y / 1000,
                z: earthState.z / 1000,
                vx: earthState.vx / 1000,
                vy: earthState.vy / 1000,
                vz: earthState.vz / 1000
            };
            const marsPos = {
                x: marsState.x / 1000,
                y: marsState.y / 1000,
                z: marsState.z / 1000,
                vx: marsState.vx / 1000,
                vy: marsState.vy / 1000,
                vz: marsState.vz / 1000
            };

            // Solve Lambert problem
            const deltaV = solveLambertApprox(earthPos, marsPos, tof, muSun, true);

            row.push(deltaV);

            if (deltaV < minDeltaV && !isNaN(deltaV)) {
                minDeltaV = deltaV;
                optimalDeparture = tDep;
                optimalArrival = tArr;
            }
        }
        deltaVGrid.push(row);
    }

    return { departureDays, arrivalDays, deltaVGrid, minDeltaV, optimalDeparture, optimalArrival, source: 'SPICE' };
}

/**
 * Compute transfer window using analytical Keplerian elements (fallback)
 */
function computeTransferWindowAnalytical(config) {
    // Simplified orbital elements
    const AU = 149597870.7;  // km
    const muSun = 1.32712440018e11;  // km³/s²

    const earthParams = { a: 1.000, e: 0.017, period: 365.25 };
    const marsParams = { a: 1.524, e: 0.093, period: 687.0 };

    const departureDays = [];
    const arrivalDays = [];
    const deltaVGrid = [];

    const dDep = config.departureWindowDays / (config.resolution - 1);
    const dArr = config.arrivalWindowDays / (config.resolution - 1);

    let minDeltaV = Infinity;
    let optimalDeparture = 0;
    let optimalArrival = 0;

    for (let i = 0; i < config.resolution; i++) {
        const tDep = config.departureStart + i * dDep;
        departureDays.push(tDep);

        const row = [];
        for (let j = 0; j < config.resolution; j++) {
            const tArr = config.arrivalStart + j * dArr;
            if (i === 0) arrivalDays.push(tArr);

            const tof = (tArr - tDep) * 86400;  // seconds

            if (tof < 100 * 86400 || tof > 400 * 86400) {
                // Invalid flight time
                row.push(NaN);
                continue;
            }

            // Get planet positions (simplified circular orbits)
            const earthPos = getPlanetPositionAnalytical(earthParams, tDep, AU);
            const marsPos = getPlanetPositionAnalytical(marsParams, tArr, AU);

            // Solve Lambert problem (simplified)
            const deltaV = solveLambertApprox(earthPos, marsPos, tof, muSun, false);

            row.push(deltaV);

            if (deltaV < minDeltaV && !isNaN(deltaV)) {
                minDeltaV = deltaV;
                optimalDeparture = tDep;
                optimalArrival = tArr;
            }
        }
        deltaVGrid.push(row);
    }

    return { departureDays, arrivalDays, deltaVGrid, minDeltaV, optimalDeparture, optimalArrival, source: 'analytical' };
}

function getPlanetPositionAnalytical(params, daysSinceEpoch, AU) {
    const n = 2 * Math.PI / params.period;  // rad/day
    const M = n * daysSinceEpoch;

    // Solve Kepler's equation
    let E = M;
    for (let i = 0; i < 10; i++) {
        E = E - (E - params.e * Math.sin(E) - M) / (1 - params.e * Math.cos(E));
    }

    const nu = 2 * Math.atan2(
        Math.sqrt(1 + params.e) * Math.sin(E / 2),
        Math.sqrt(1 - params.e) * Math.cos(E / 2)
    );

    const r = params.a * (1 - params.e * Math.cos(E)) * AU;

    return {
        x: r * Math.cos(nu),
        y: r * Math.sin(nu),
        z: 0,
        r: r
    };
}

function solveLambertApprox(r1, r2, tof, mu, hasVelocity = false) {
    // Simplified Lambert solver using Hohmann-like approximation
    const r1Mag = Math.sqrt(r1.x * r1.x + r1.y * r1.y + (r1.z || 0) * (r1.z || 0));
    const r2Mag = Math.sqrt(r2.x * r2.x + r2.y * r2.y + (r2.z || 0) * (r2.z || 0));

    // Transfer angle
    const dot = r1.x * r2.x + r1.y * r2.y + (r1.z || 0) * (r2.z || 0);
    const cosTheta = dot / (r1Mag * r2Mag);
    const theta = Math.acos(Math.max(-1, Math.min(1, cosTheta)));

    // Semi-major axis estimate
    const c = Math.sqrt(r1Mag * r1Mag + r2Mag * r2Mag - 2 * r1Mag * r2Mag * cosTheta);
    const s = (r1Mag + r2Mag + c) / 2;
    const aMin = s / 2;

    // Parabolic time of flight
    const tofPara = Math.sqrt(2) / 3 * Math.sqrt(s * s * s / mu) * (1 - Math.pow((s - c) / s, 1.5));

    // Estimate semi-major axis from time of flight ratio
    const tofRatio = tof / tofPara;
    let a = aMin * Math.pow(tofRatio, 0.5);

    // Vis-viva for departure and arrival velocities
    const v1Circ = Math.sqrt(mu / r1Mag);
    const v2Circ = Math.sqrt(mu / r2Mag);

    const v1Trans = Math.sqrt(mu * (2 / r1Mag - 1 / a));
    const v2Trans = Math.sqrt(mu * (2 / r2Mag - 1 / a));

    // Delta-V computation
    let dv1, dv2;

    if (hasVelocity && r1.vx !== undefined) {
        // Use actual planet velocities from SPICE
        const v1Planet = Math.sqrt(r1.vx * r1.vx + r1.vy * r1.vy + (r1.vz || 0) * (r1.vz || 0));
        const v2Planet = Math.sqrt(r2.vx * r2.vx + r2.vy * r2.vy + (r2.vz || 0) * (r2.vz || 0));
        dv1 = Math.abs(v1Trans - v1Planet);
        dv2 = Math.abs(v2Trans - v2Planet);
    } else {
        // Use circular velocity approximation
        dv1 = Math.abs(v1Trans - v1Circ);
        dv2 = Math.abs(v2Trans - v2Circ);
    }

    // Add penalty for very short/long transfers
    let penalty = 0;
    if (tof < 150 * 86400) penalty = (150 * 86400 - tof) / 86400 * 100;
    if (tof > 350 * 86400) penalty = (tof - 350 * 86400) / 86400 * 50;

    return (dv1 + dv2) * 1000 + penalty;  // Convert to m/s
}

function renderPorkchopPlot(container, result, config) {
    container.innerHTML = '';

    const containerRect = container.getBoundingClientRect();
    const containerWidth = containerRect.width || 600;
    const containerHeight = containerRect.height || 500;

    const wrapper = document.createElement('div');
    wrapper.style.cssText = `
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        padding: 15px;
        width: 100%;
        height: 100%;
        box-sizing: border-box;
    `;
    container.appendChild(wrapper);

    const title = document.createElement('div');
    title.style.cssText = 'font-family: "Orbitron", sans-serif; font-size: 14px; color: var(--cyan); margin-bottom: 10px; text-align: center;';
    // Show ephemeris source in title
    const sourceLabel = result.source === 'analytical' ? 'Analytical' :
                        result.source.includes('CALCEPH') ? 'Binary SPK' :
                        result.source.includes('JSON') ? 'Precomputed' : 'SPICE';
    title.textContent = `Earth-Mars Transfer (${sourceLabel})`;
    wrapper.appendChild(title);

    const chartDiv = document.createElement('div');
    wrapper.appendChild(chartDiv);

    const size = Math.min(containerWidth - 80, containerHeight - 120);
    const width = size;
    const height = size;
    const margin = { top: 20, right: 60, bottom: 50, left: 60 };

    const svg = d3.select(chartDiv)
        .append('svg')
        .attr('width', width)
        .attr('height', height);

    // Scales
    const x = d3.scaleLinear()
        .domain([result.departureDays[0], result.departureDays[result.departureDays.length - 1]])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([result.arrivalDays[0], result.arrivalDays[result.arrivalDays.length - 1]])
        .range([height - margin.bottom, margin.top]);

    // Color scale for delta-V
    const validValues = result.deltaVGrid.flat().filter(v => !isNaN(v));
    const minV = Math.min(...validValues);
    const maxV = Math.min(Math.max(...validValues), minV * 3);  // Cap at 3x minimum

    const color = d3.scaleSequential(d3.interpolateViridis)
        .domain([maxV, minV]);  // Reversed so lower is brighter

    // Draw heatmap cells
    const cellWidth = (width - margin.left - margin.right) / config.resolution;
    const cellHeight = (height - margin.top - margin.bottom) / config.resolution;

    for (let i = 0; i < config.resolution; i++) {
        for (let j = 0; j < config.resolution; j++) {
            const val = result.deltaVGrid[i][j];
            if (isNaN(val)) continue;

            svg.append('rect')
                .attr('x', margin.left + i * cellWidth)
                .attr('y', margin.top + (config.resolution - 1 - j) * cellHeight)
                .attr('width', cellWidth)
                .attr('height', cellHeight)
                .attr('fill', color(Math.min(val, maxV)));
        }
    }

    // Optimal point marker
    svg.append('circle')
        .attr('cx', x(result.optimalDeparture))
        .attr('cy', y(result.optimalArrival))
        .attr('r', 8)
        .attr('fill', 'none')
        .attr('stroke', 'var(--red)')
        .attr('stroke-width', 2);

    svg.append('text')
        .attr('x', x(result.optimalDeparture) + 12)
        .attr('y', y(result.optimalArrival) + 4)
        .attr('fill', 'var(--red)')
        .attr('font-size', '10px')
        .text(`${(result.minDeltaV/1000).toFixed(1)} km/s`);

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(6))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(6))
        .attr('color', 'var(--text-dim)');

    svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 10)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '11px')
        .text('Departure Day (from J2000)');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 15)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '11px')
        .text('Arrival Day (from J2000)');

    // Color bar
    const colorBarWidth = 15;
    const colorBarHeight = height - margin.top - margin.bottom;
    const colorBarX = width - margin.right + 10;

    const colorBarData = d3.range(0, 1.01, 0.02);
    colorBarData.forEach((t, i) => {
        svg.append('rect')
            .attr('x', colorBarX)
            .attr('y', margin.top + i * (colorBarHeight / colorBarData.length))
            .attr('width', colorBarWidth)
            .attr('height', colorBarHeight / colorBarData.length + 1)
            .attr('fill', color(minV + t * (maxV - minV)));
    });

    svg.append('text')
        .attr('x', colorBarX + colorBarWidth / 2)
        .attr('y', margin.top - 5)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-dim)')
        .attr('font-size', '9px')
        .text(`${(minV/1000).toFixed(1)}`);

    svg.append('text')
        .attr('x', colorBarX + colorBarWidth / 2)
        .attr('y', height - margin.bottom + 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-dim)')
        .attr('font-size', '9px')
        .text(`${(maxV/1000).toFixed(1)}`);

    svg.append('text')
        .attr('transform', `rotate(90)`)
        .attr('x', height / 2)
        .attr('y', -width + margin.right - 25)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('ΔV [km/s]');
}
