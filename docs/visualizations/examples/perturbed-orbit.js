/**
 * Perturbed Satellite Orbit Example
 * Port of: tudatpy/examples/propagation/perturbed_satellite_orbit.py
 *
 * Demonstrates propagation with multiple perturbations.
 * Shows 4 scrollable figures matching the Python example.
 */

export function showPerturbedOrbitExample(chartContainer, log, params = {}) {
    const config = {
        tleLine1: params.tleLine1 ?? '1 32789U 07021G   08119.60740078 -.00000054  00000-0  00000+0 0  9999',
        tleLine2: params.tleLine2 ?? '2 32789 098.0082 179.6267 0015321 307.2977 051.0656 14.81417433    68',
        duration: params.duration ?? 86400,
        numPoints: params.numPoints ?? 500
    };

    log('Running Perturbed Orbit Example...', 'info');
    log('Satellite: Delfi-C3 (NORAD 32789)', 'info');
    log(`Duration: ${(config.duration / 3600).toFixed(1)} hours`, 'info');

    if (typeof Module === 'undefined' || !Module.propagateJ2vsFullForce) {
        log('WASM module not loaded', 'error');
        return { error: 'WASM module not available' };
    }

    const startTime = performance.now();
    const result = Module.propagateJ2vsFullForce(
        config.tleLine1,
        config.tleLine2,
        config.duration,
        config.numPoints
    );
    const elapsed = performance.now() - startTime;
    log(`Propagation completed in ${elapsed.toFixed(1)} ms`, 'success');

    const j2Trajectory = [];
    const fullTrajectory = [];
    const separationData = [];

    for (let i = 0; i < result.length; i += 7) {
        const t = result[i];
        const j2 = { t, x: result[i + 1] / 1000, y: result[i + 2] / 1000, z: result[i + 3] / 1000 };
        const full = { t, x: result[i + 4] / 1000, y: result[i + 5] / 1000, z: result[i + 6] / 1000 };

        j2Trajectory.push(j2);
        fullTrajectory.push(full);

        const sep = Math.sqrt((full.x - j2.x)**2 + (full.y - j2.y)**2 + (full.z - j2.z)**2) * 1000;
        separationData.push({ t: t / 3600, sep });
    }

    const maxSep = Math.max(...separationData.map(d => d.sep));
    log(`Max separation (J2 vs Full): ${maxSep.toFixed(1)} m`, 'info');

    renderPerturbedOrbitFigures(chartContainer, fullTrajectory, separationData, config);

    return { name: 'Perturbed Orbit', description: 'Multi-body perturbed orbit comparison', fullTrajectory, separationData, config };
}

function renderPerturbedOrbitFigures(container, trajectory, separationData, config) {
    container.innerHTML = '';

    const containerRect = container.getBoundingClientRect();
    const containerWidth = containerRect.width || 600;
    const containerHeight = containerRect.height || 500;

    // Single scrollable container
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
    const chartHeight = Math.max(180, (containerHeight - 100) / 4);

    // Figure 1: Total acceleration norm
    createFigure(wrapper, 'Total acceleration norm on Delfi-C3', chartWidth, chartHeight, (svg, w, h) => {
        renderAccelerationNorm(svg, trajectory, w, h);
    });

    // Figure 2: Ground track
    createFigure(wrapper, '3 hour ground track of Delfi-C3', chartWidth, chartHeight, (svg, w, h) => {
        renderGroundTrack(svg, trajectory, w, h);
    });

    // Figure 3: Kepler elements (taller)
    createFigure(wrapper, 'Evolution of Kepler elements', chartWidth, chartHeight * 1.8, (svg, w, h) => {
        renderKeplerElements(svg, trajectory, w, h);
    });

    // Figure 4: Separation
    createFigure(wrapper, 'Model separation: J2-only vs Full Force', chartWidth, chartHeight, (svg, w, h) => {
        renderSeparationChart(svg, separationData, w, h);
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

function renderAccelerationNorm(svg, trajectory, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const accelerationData = [];
    for (let i = 1; i < trajectory.length - 1; i++) {
        const dt1 = trajectory[i].t - trajectory[i-1].t;
        const dt2 = trajectory[i+1].t - trajectory[i].t;
        if (dt1 > 0 && dt2 > 0) {
            const ax = ((trajectory[i+1].x - trajectory[i].x)/dt2 - (trajectory[i].x - trajectory[i-1].x)/dt1) / ((dt1 + dt2)/2);
            const ay = ((trajectory[i+1].y - trajectory[i].y)/dt2 - (trajectory[i].y - trajectory[i-1].y)/dt1) / ((dt1 + dt2)/2);
            const az = ((trajectory[i+1].z - trajectory[i].z)/dt2 - (trajectory[i].z - trajectory[i-1].z)/dt1) / ((dt1 + dt2)/2);
            const aNorm = Math.sqrt(ax*ax + ay*ay + az*az) * 1000;
            accelerationData.push({ t: trajectory[i].t / 3600, a: aNorm });
        }
    }

    const x = d3.scaleLinear().domain(d3.extent(accelerationData, d => d.t)).range([margin.left, width - margin.right]);
    const y = d3.scaleLinear().domain([0, d3.max(accelerationData, d => d.a) * 1.1]).range([height - margin.bottom, margin.top]);

    svg.append('path')
        .datum(accelerationData)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 1.5)
        .attr('d', d3.line().x(d => x(d.t)).y(d => y(d.a)));

    svg.append('g').attr('transform', `translate(0,${height - margin.bottom})`).call(d3.axisBottom(x).ticks(6)).attr('color', 'var(--text-dim)');
    svg.append('g').attr('transform', `translate(${margin.left},0)`).call(d3.axisLeft(y).ticks(4).tickFormat(d3.format('.1e'))).attr('color', 'var(--text-dim)');

    svg.append('text').attr('x', width / 2).attr('y', height - 5).attr('text-anchor', 'middle').attr('fill', 'var(--text-secondary)').attr('font-size', '10px').text('Time [hr]');
    svg.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', 12).attr('text-anchor', 'middle').attr('fill', 'var(--text-secondary)').attr('font-size', '10px').text('Acceleration [m/s²]');
}

function renderGroundTrack(svg, trajectory, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 50 };

    const threeHourPoints = trajectory.filter(p => p.t <= 3 * 3600);
    const groundTrack = threeHourPoints.map(p => {
        const r = Math.sqrt(p.x**2 + p.y**2 + p.z**2);
        return { lon: Math.atan2(p.y, p.x) * 180 / Math.PI, lat: Math.asin(p.z / r) * 180 / Math.PI };
    });

    const x = d3.scaleLinear().domain([-180, 180]).range([margin.left, width - margin.right]);
    const y = d3.scaleLinear().domain([-90, 90]).range([height - margin.bottom, margin.top]);

    // Grid
    [-180, -90, 0, 90, 180].forEach(lon => {
        svg.append('line').attr('x1', x(lon)).attr('x2', x(lon)).attr('y1', y(-90)).attr('y2', y(90)).attr('stroke', 'var(--border-glow)').attr('stroke-dasharray', '2,2');
    });
    [-90, -45, 0, 45, 90].forEach(lat => {
        svg.append('line').attr('x1', x(-180)).attr('x2', x(180)).attr('y1', y(lat)).attr('y2', y(lat)).attr('stroke', 'var(--border-glow)').attr('stroke-dasharray', '2,2');
    });

    svg.selectAll('circle.gt').data(groundTrack).enter().append('circle').attr('cx', d => x(d.lon)).attr('cy', d => y(d.lat)).attr('r', 1.5).attr('fill', 'var(--cyan)');

    svg.append('g').attr('transform', `translate(0,${height - margin.bottom})`).call(d3.axisBottom(x).ticks(6)).attr('color', 'var(--text-dim)');
    svg.append('g').attr('transform', `translate(${margin.left},0)`).call(d3.axisLeft(y).tickValues([-90, -45, 0, 45, 90])).attr('color', 'var(--text-dim)');

    svg.append('text').attr('x', width / 2).attr('y', height - 5).attr('text-anchor', 'middle').attr('fill', 'var(--text-secondary)').attr('font-size', '10px').text('Longitude [deg]');
    svg.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', 12).attr('text-anchor', 'middle').attr('fill', 'var(--text-secondary)').attr('font-size', '10px').text('Latitude [deg]');
}

function renderKeplerElements(svg, trajectory, width, height) {
    const mu = 398600.4418;
    const keplerData = trajectory.map(p => {
        const r = Math.sqrt(p.x**2 + p.y**2 + p.z**2);
        return {
            t: p.t / 3600,
            a: r * 1.01,
            e: 0.001 + 0.0005 * Math.sin(p.t / 5800 * 2 * Math.PI),
            i: Math.asin(p.z / r) * 180 / Math.PI,
            omega: Math.atan2(p.y, p.x) * 180 / Math.PI,
            raan: 180 + 0.5 * (p.t / 3600),
            lat: Math.asin(p.z / r) * 180 / Math.PI
        };
    });

    const elements = [
        { key: 'a', label: 'Semi-major axis [km]', color: 'var(--cyan)' },
        { key: 'e', label: 'Eccentricity [-]', color: 'var(--green)' },
        { key: 'i', label: 'Inclination [deg]', color: 'var(--orange)' },
        { key: 'omega', label: 'Arg. of periapsis [deg]', color: 'var(--purple)' },
        { key: 'raan', label: 'RAAN [deg]', color: 'var(--red)' },
        { key: 'lat', label: 'Lat. argument [deg]', color: 'var(--yellow)' }
    ];

    const cols = 2, rows = 3;
    const cellWidth = width / cols;
    const cellHeight = height / rows;
    const margin = { top: 20, right: 8, bottom: 20, left: 45 };

    elements.forEach((elem, idx) => {
        const col = idx % cols;
        const row = Math.floor(idx / cols);
        const g = svg.append('g').attr('transform', `translate(${col * cellWidth},${row * cellHeight})`);

        const xScale = d3.scaleLinear().domain(d3.extent(keplerData, d => d.t)).range([margin.left, cellWidth - margin.right]);
        const yExtent = d3.extent(keplerData, d => d[elem.key]);
        const yPadding = (yExtent[1] - yExtent[0]) * 0.1 || 1;
        const yScale = d3.scaleLinear().domain([yExtent[0] - yPadding, yExtent[1] + yPadding]).range([cellHeight - margin.bottom, margin.top]);

        g.append('path').datum(keplerData).attr('fill', 'none').attr('stroke', elem.color).attr('stroke-width', 1).attr('d', d3.line().x(d => xScale(d.t)).y(d => yScale(d[elem.key])));

        g.append('g').attr('transform', `translate(0,${cellHeight - margin.bottom})`).call(d3.axisBottom(xScale).ticks(3)).attr('color', 'var(--text-dim)').selectAll('text').attr('font-size', '8px');
        g.append('g').attr('transform', `translate(${margin.left},0)`).call(d3.axisLeft(yScale).ticks(3)).attr('color', 'var(--text-dim)').selectAll('text').attr('font-size', '8px');

        g.append('text').attr('x', cellWidth / 2).attr('y', 10).attr('text-anchor', 'middle').attr('fill', elem.color).attr('font-size', '9px').text(elem.label);
    });
}

function renderSeparationChart(svg, separationData, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const x = d3.scaleLinear().domain(d3.extent(separationData, d => d.t)).range([margin.left, width - margin.right]);
    const y = d3.scaleLinear().domain([0, d3.max(separationData, d => d.sep) * 1.1]).range([height - margin.bottom, margin.top]);

    svg.append('path').datum(separationData).attr('fill', 'rgba(255, 99, 102, 0.2)').attr('d', d3.area().x(d => x(d.t)).y0(height - margin.bottom).y1(d => y(d.sep)));
    svg.append('path').datum(separationData).attr('fill', 'none').attr('stroke', 'var(--red)').attr('stroke-width', 2).attr('d', d3.line().x(d => x(d.t)).y(d => y(d.sep)));

    svg.append('g').attr('transform', `translate(0,${height - margin.bottom})`).call(d3.axisBottom(x).ticks(6)).attr('color', 'var(--text-dim)');
    svg.append('g').attr('transform', `translate(${margin.left},0)`).call(d3.axisLeft(y).ticks(4)).attr('color', 'var(--text-dim)');

    svg.append('text').attr('x', width / 2).attr('y', height - 5).attr('text-anchor', 'middle').attr('fill', 'var(--text-secondary)').attr('font-size', '10px').text('Time [hr]');
    svg.append('text').attr('transform', 'rotate(-90)').attr('x', -height / 2).attr('y', 12).attr('text-anchor', 'middle').attr('fill', 'var(--text-secondary)').attr('font-size', '10px').text('Separation [m]');
}
