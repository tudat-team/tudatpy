/**
 * CR3BP Manifolds Example
 * Port of: tudatpy/examples/propagation/impact_manifolds_lpo_cr3bp.py
 *
 * Demonstrates stable and unstable manifolds of libration point orbits
 * in the Circular Restricted Three-Body Problem.
 */

export function showCR3BPManifoldsExample(chartContainer, log, params = {}) {
    const config = {
        mu: params.mu ?? 0.01215,           // Earth-Moon mass ratio
        haloAmplitude: params.haloAmplitude ?? 0.02,  // Halo orbit amplitude
        manifoldTime: params.manifoldTime ?? 3.0,     // Propagation time (normalized)
        numTrajectories: params.numTrajectories ?? 20
    };

    log('Running CR3BP Manifolds Example...', 'info');
    log(`Mass ratio μ = ${config.mu.toFixed(5)} (Earth-Moon)`, 'info');
    log(`Halo amplitude: ${config.haloAmplitude}`, 'info');
    log(`Number of manifold trajectories: ${config.numTrajectories}`, 'info');

    const startTime = performance.now();
    const result = computeManifolds(config);
    const elapsed = performance.now() - startTime;
    log(`Computation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`L1 location: x = ${result.L1.toFixed(4)}`, 'info');
    log(`Halo orbit period: ${result.haloPeriod.toFixed(4)} (normalized)`, 'info');

    renderManifoldsFigures(chartContainer, result, config);

    return {
        name: 'CR3BP Manifolds',
        description: 'Stable/unstable manifolds of halo orbits',
        ...result,
        config
    };
}

function computeManifolds(config) {
    const { mu, haloAmplitude, manifoldTime, numTrajectories } = config;

    // Find L1 Lagrange point (between primaries)
    // L1 is at x where: x - (1-μ)/(x+μ)² + μ/(x-1+μ)² = 0
    // Approximate L1 for Earth-Moon
    const L1 = 1 - mu - Math.pow(mu / 3, 1/3);

    // Approximate halo orbit around L1 (Richardson third-order)
    const gamma = Math.pow(mu / 3, 1/3);
    const c2 = (mu + (1 - mu) * gamma**3 / (1 - gamma)**3) / gamma**3;
    const lambda = Math.sqrt((2 - c2 + Math.sqrt(9 * c2**2 - 8 * c2)) / 2);

    // Halo orbit (simplified circular approximation in xy plane)
    const Ax = haloAmplitude;
    const Az = haloAmplitude * 0.5;
    const period = 2 * Math.PI / lambda;

    // Generate halo orbit points
    const haloOrbit = [];
    const haloPoints = 100;
    for (let i = 0; i <= haloPoints; i++) {
        const tau = i / haloPoints * period;
        const x = L1 - Ax * Math.cos(lambda * tau);
        const y = Ax * Math.sin(lambda * tau) * 2;
        const z = Az * Math.cos(lambda * tau);
        haloOrbit.push({ x, y, z, tau });
    }

    // CR3BP equations of motion
    function cr3bpDerivatives(state) {
        const [x, y, z, vx, vy, vz] = state;

        const r1 = Math.sqrt((x + mu)**2 + y**2 + z**2);
        const r2 = Math.sqrt((x - 1 + mu)**2 + y**2 + z**2);

        const ax = 2 * vy + x - (1 - mu) * (x + mu) / r1**3 - mu * (x - 1 + mu) / r2**3;
        const ay = -2 * vx + y - (1 - mu) * y / r1**3 - mu * y / r2**3;
        const az = -(1 - mu) * z / r1**3 - mu * z / r2**3;

        return [vx, vy, vz, ax, ay, az];
    }

    // RK4 integrator
    function rk4Step(state, dt) {
        const k1 = cr3bpDerivatives(state);
        const s1 = state.map((v, i) => v + 0.5 * dt * k1[i]);
        const k2 = cr3bpDerivatives(s1);
        const s2 = state.map((v, i) => v + 0.5 * dt * k2[i]);
        const k3 = cr3bpDerivatives(s2);
        const s3 = state.map((v, i) => v + dt * k3[i]);
        const k4 = cr3bpDerivatives(s3);
        return state.map((v, i) => v + dt / 6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
    }

    // Generate manifold trajectories
    const stableManifolds = [];
    const unstableManifolds = [];
    const dt = 0.01;
    const steps = Math.floor(manifoldTime / dt);

    // Sample points along halo orbit
    for (let i = 0; i < numTrajectories; i++) {
        const tau = i / numTrajectories * period;
        const haloIdx = Math.floor(i / numTrajectories * haloPoints);
        const haloPoint = haloOrbit[haloIdx];

        // Approximate eigenvector direction (outward from L1)
        const dx = haloPoint.x - L1;
        const dy = haloPoint.y;
        const norm = Math.sqrt(dx**2 + dy**2) || 0.001;
        const eps = 0.001;  // Perturbation magnitude

        // Unstable manifold (forward propagation)
        let stateU = [
            haloPoint.x + eps * dx / norm,
            haloPoint.y + eps * dy / norm,
            haloPoint.z,
            -lambda * Ax * Math.sin(lambda * tau) + eps * 0.1,
            lambda * Ax * 2 * Math.cos(lambda * tau),
            -lambda * Az * Math.sin(lambda * tau)
        ];

        const unstableTraj = [{ x: stateU[0], y: stateU[1], z: stateU[2] }];
        for (let j = 0; j < steps; j++) {
            stateU = rk4Step(stateU, dt);
            if (j % 5 === 0) {
                unstableTraj.push({ x: stateU[0], y: stateU[1], z: stateU[2] });
            }
            // Stop if too far
            if (Math.abs(stateU[0]) > 2 || Math.abs(stateU[1]) > 2) break;
        }
        unstableManifolds.push(unstableTraj);

        // Stable manifold (backward propagation)
        let stateS = [
            haloPoint.x - eps * dx / norm,
            haloPoint.y - eps * dy / norm,
            haloPoint.z,
            -lambda * Ax * Math.sin(lambda * tau) - eps * 0.1,
            lambda * Ax * 2 * Math.cos(lambda * tau),
            -lambda * Az * Math.sin(lambda * tau)
        ];

        const stableTraj = [{ x: stateS[0], y: stateS[1], z: stateS[2] }];
        for (let j = 0; j < steps; j++) {
            stateS = rk4Step(stateS, -dt);  // Backward
            if (j % 5 === 0) {
                stableTraj.push({ x: stateS[0], y: stateS[1], z: stateS[2] });
            }
            if (Math.abs(stateS[0]) > 2 || Math.abs(stateS[1]) > 2) break;
        }
        stableManifolds.push(stableTraj);
    }

    return {
        L1,
        haloOrbit,
        haloPeriod: period,
        stableManifolds,
        unstableManifolds,
        mu,
        primaryPositions: { m1: [-mu, 0], m2: [1 - mu, 0] }
    };
}

function renderManifoldsFigures(container, result, config) {
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
    const chartHeight = Math.max(300, containerHeight - 100);

    // Main figure: XY view of manifolds
    createFigure(wrapper, 'CR3BP Manifolds (XY Projection)', chartWidth, chartHeight, (svg, w, h) => {
        renderManifoldsXY(svg, result, w, h);
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

function renderManifoldsXY(svg, result, width, height) {
    const margin = { top: 20, right: 20, bottom: 40, left: 50 };

    const x = d3.scaleLinear()
        .domain([0.6, 1.2])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-0.4, 0.4])
        .range([height - margin.bottom, margin.top]);

    // Draw stable manifolds (green)
    result.stableManifolds.forEach(traj => {
        if (traj.length > 1) {
            svg.append('path')
                .datum(traj)
                .attr('fill', 'none')
                .attr('stroke', 'var(--green)')
                .attr('stroke-width', 0.5)
                .attr('opacity', 0.6)
                .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));
        }
    });

    // Draw unstable manifolds (red)
    result.unstableManifolds.forEach(traj => {
        if (traj.length > 1) {
            svg.append('path')
                .datum(traj)
                .attr('fill', 'none')
                .attr('stroke', 'var(--red)')
                .attr('stroke-width', 0.5)
                .attr('opacity', 0.6)
                .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));
        }
    });

    // Draw halo orbit (yellow)
    svg.append('path')
        .datum(result.haloOrbit)
        .attr('fill', 'none')
        .attr('stroke', 'var(--yellow)')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

    // Draw L1 point
    svg.append('circle')
        .attr('cx', x(result.L1))
        .attr('cy', y(0))
        .attr('r', 4)
        .attr('fill', 'var(--cyan)');

    svg.append('text')
        .attr('x', x(result.L1) + 8)
        .attr('y', y(0) + 4)
        .attr('fill', 'var(--cyan)')
        .attr('font-size', '10px')
        .text('L1');

    // Draw primaries
    svg.append('circle')
        .attr('cx', x(result.primaryPositions.m1[0]))
        .attr('cy', y(0))
        .attr('r', 8)
        .attr('fill', 'var(--blue)');

    svg.append('circle')
        .attr('cx', x(result.primaryPositions.m2[0]))
        .attr('cy', y(0))
        .attr('r', 4)
        .attr('fill', 'var(--text-dim)');

    // Legend
    const legend = svg.append('g')
        .attr('transform', `translate(${width - margin.right - 100}, ${margin.top})`);

    [
        { color: 'var(--green)', label: 'Stable' },
        { color: 'var(--red)', label: 'Unstable' },
        { color: 'var(--yellow)', label: 'Halo Orbit' }
    ].forEach((item, i) => {
        legend.append('line')
            .attr('x1', 0).attr('y1', i * 16)
            .attr('x2', 20).attr('y2', i * 16)
            .attr('stroke', item.color)
            .attr('stroke-width', 2);

        legend.append('text')
            .attr('x', 25).attr('y', i * 16 + 4)
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '10px')
            .text(item.label);
    });

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
        .attr('y', height - 5)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('x (normalized)');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('y (normalized)');
}
