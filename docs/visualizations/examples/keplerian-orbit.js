/**
 * Keplerian Satellite Orbit Example
 * Port of: tudatpy/examples/propagation/keplerian_satellite_orbit.py
 *
 * Demonstrates basic propagation of a satellite under point-mass gravity.
 * Shows the classic two-body problem - ONE 3D trajectory plot matching the Python example.
 */

/**
 * Run the Keplerian orbit example and display results
 * @param {HTMLElement} chartContainer - Container for D3 charts
 * @param {Function} log - Logging function
 * @param {Object} params - Optional parameters to override defaults
 */
export function showKeplerianOrbitExample(chartContainer, log, params = {}) {
    // Default parameters matching the Python example (Delfi-C3)
    const config = {
        semiMajorAxis: params.semiMajorAxis ?? 6992.76221,  // km
        eccentricity: params.eccentricity ?? 0.00403294322,
        inclination: params.inclination ?? 98.0,  // degrees (approx from 1.71 rad)
        raan: params.raan ?? 21.94,  // degrees (approx from 0.383 rad)
        argPeriapsis: params.argPeriapsis ?? 75.2,  // degrees (approx from 1.312 rad)
        trueAnomaly: params.trueAnomaly ?? 175.9,  // degrees (approx from 3.07 rad)
        duration: params.duration ?? 86400,  // 1 day in seconds
        numPoints: params.numPoints ?? 1000
    };

    log('Running Keplerian Orbit Example...', 'info');
    log(`Semi-major axis: ${config.semiMajorAxis.toFixed(2)} km`, 'info');
    log(`Eccentricity: ${config.eccentricity.toFixed(6)}`, 'info');
    log(`Inclination: ${config.inclination.toFixed(2)} deg`, 'info');
    log(`Duration: ${(config.duration / 3600).toFixed(1)} hours`, 'info');

    // Check if WASM module is available
    if (typeof Module === 'undefined' || !Module.propagateKeplerOrbit) {
        log('WASM module not loaded', 'error');
        return { error: 'WASM module not available' };
    }

    // Propagate the orbit using Tudat WASM
    const startTime = performance.now();
    const result = Module.propagateKeplerOrbit(
        config.semiMajorAxis,
        config.eccentricity,
        config.inclination,
        config.raan,
        config.argPeriapsis,
        config.trueAnomaly,
        config.duration,
        config.numPoints
    );
    const elapsed = performance.now() - startTime;
    log(`Propagation completed in ${elapsed.toFixed(1)} ms`, 'success');

    // Parse results: [t, x, y, z, t, x, y, z, ...]
    const trajectory = [];
    for (let i = 0; i < result.length; i += 4) {
        trajectory.push({
            t: result[i],
            x: result[i + 1] / 1000,  // Convert to km
            y: result[i + 2] / 1000,
            z: result[i + 3] / 1000
        });
    }

    // Compute orbital parameters for display
    const r0 = Math.sqrt(trajectory[0].x**2 + trajectory[0].y**2 + trajectory[0].z**2);
    const rFinal = Math.sqrt(
        trajectory[trajectory.length-1].x**2 +
        trajectory[trajectory.length-1].y**2 +
        trajectory[trajectory.length-1].z**2
    );

    log(`Initial radius: ${r0.toFixed(2)} km`, 'info');
    log(`Final radius: ${rFinal.toFixed(2)} km`, 'info');

    // Render single 3D trajectory chart (matching Python example)
    render3DTrajectory(chartContainer, trajectory, config);

    return {
        name: 'Keplerian Orbit',
        description: 'Two-body satellite orbit propagation',
        trajectory,
        config
    };
}

/**
 * Render a single 3D trajectory visualization using isometric projection
 * Matches the Python example's single 3D plot - fully responsive
 */
function render3DTrajectory(container, trajectory, config) {
    // Clear existing content
    container.innerHTML = '';

    // Get actual container dimensions
    const containerRect = container.getBoundingClientRect();
    const containerWidth = containerRect.width || 600;
    const containerHeight = containerRect.height || 500;

    // Use the smaller dimension to keep aspect ratio square-ish
    const size = Math.min(containerWidth - 40, containerHeight - 60);
    const width = size;
    const height = size;
    const margin = { top: 30, right: 30, bottom: 50, left: 50 };

    // Create wrapper that fills container
    const wrapper = document.createElement('div');
    wrapper.style.cssText = `
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        width: 100%;
        height: 100%;
        overflow: hidden;
    `;
    container.appendChild(wrapper);

    // Title
    const title = document.createElement('div');
    title.style.cssText = 'font-family: "Orbitron", sans-serif; font-size: 14px; color: var(--cyan); margin-bottom: 10px; text-align: center;';
    title.textContent = 'Delfi-C3 trajectory around Earth';
    wrapper.appendChild(title);

    // Chart container - fills remaining space
    const chartDiv = document.createElement('div');
    chartDiv.style.cssText = `width: ${width}px; height: ${height}px;`;
    wrapper.appendChild(chartDiv);

    const svg = d3.select(chartDiv)
        .append('svg')
        .attr('width', width)
        .attr('height', height);

    // Isometric projection angles
    const angleX = -25 * Math.PI / 180;
    const angleZ = 35 * Math.PI / 180;

    function project(x, y, z) {
        const x1 = x * Math.cos(angleZ) - y * Math.sin(angleZ);
        const y1 = x * Math.sin(angleZ) + y * Math.cos(angleZ);
        const z1 = z;
        const x2 = x1;
        const y2 = y1 * Math.cos(angleX) - z1 * Math.sin(angleX);
        const z2 = y1 * Math.sin(angleX) + z1 * Math.cos(angleX);
        return { x: x2, y: z2, depth: y2 };
    }

    const projectedPoints = trajectory.map(p => project(p.x, p.y, p.z));

    const xExtent = d3.extent(projectedPoints, d => d.x);
    const yExtent = d3.extent(projectedPoints, d => d.y);
    const range = Math.max(xExtent[1] - xExtent[0], yExtent[1] - yExtent[0]) * 0.6;
    const xMid = (xExtent[0] + xExtent[1]) / 2;
    const yMid = (yExtent[0] + yExtent[1]) / 2;

    const x = d3.scaleLinear()
        .domain([xMid - range, xMid + range])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([yMid - range, yMid + range])
        .range([height - margin.bottom, margin.top]);

    // Earth
    const earthProjected = project(0, 0, 0);
    const earthRadius = 6378.137;
    const earthScreenRadius = Math.abs(x(earthRadius) - x(0)) * 0.8;

    const earthGradient = svg.append('defs')
        .append('radialGradient')
        .attr('id', 'earthGradient')
        .attr('cx', '35%')
        .attr('cy', '35%');
    earthGradient.append('stop').attr('offset', '0%').attr('stop-color', '#4a9fff');
    earthGradient.append('stop').attr('offset', '100%').attr('stop-color', '#0066cc');

    svg.append('circle')
        .attr('cx', x(earthProjected.x))
        .attr('cy', y(earthProjected.y))
        .attr('r', earthScreenRadius)
        .attr('fill', 'url(#earthGradient)')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 1);

    // Orbit path
    const line = d3.line()
        .x(d => x(d.x))
        .y(d => y(d.y));

    svg.append('path')
        .datum(projectedPoints)
        .attr('fill', 'none')
        .attr('stroke', 'var(--orange)')
        .attr('stroke-width', 2)
        .attr('stroke-dasharray', '8,4')
        .attr('d', line);

    // Satellite position
    const lastPoint = projectedPoints[projectedPoints.length - 1];
    svg.append('circle')
        .attr('cx', x(lastPoint.x))
        .attr('cy', y(lastPoint.y))
        .attr('r', 5)
        .attr('fill', 'var(--cyan)')
        .attr('stroke', 'white')
        .attr('stroke-width', 1);

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(5).tickFormat(d => `${(d/1000).toFixed(0)}k`))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(5).tickFormat(d => `${(d/1000).toFixed(0)}k`))
        .attr('color', 'var(--text-dim)');

    // Axis labels
    svg.append('text')
        .attr('x', width / 2)
        .attr('y', height - 8)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '11px')
        .text('x [km]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 14)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '11px')
        .text('z [km]');

    // Legend
    const legend = svg.append('g')
        .attr('transform', `translate(${width - margin.right - 80}, ${margin.top})`);

    legend.append('line')
        .attr('x1', 0).attr('y1', 0)
        .attr('x2', 20).attr('y2', 0)
        .attr('stroke', 'var(--orange)')
        .attr('stroke-width', 2)
        .attr('stroke-dasharray', '6,3');

    legend.append('text')
        .attr('x', 25).attr('y', 4)
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Delfi-C3');

    legend.append('circle')
        .attr('cx', 10).attr('cy', 18)
        .attr('r', 5)
        .attr('fill', 'url(#earthGradient)');

    legend.append('text')
        .attr('x', 25).attr('y', 22)
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Earth');
}
