/**
 * Linear Sensitivity Analysis Example
 * Port of: tudatpy/examples/propagation/linear_sensitivity_analysis.py
 *
 * Demonstrates variational equations and state transition matrix.
 * Shows how small perturbations in initial conditions propagate.
 */

export function showLinearSensitivityExample(chartContainer, log, params = {}) {
    const config = {
        initialAltitude: params.initialAltitude ?? 400,  // km
        duration: params.duration ?? 86400,               // 1 day
        numPoints: params.numPoints ?? 500,
        perturbationSize: params.perturbationSize ?? 100  // meters
    };

    log('Running Linear Sensitivity Analysis...', 'info');
    log(`Initial altitude: ${config.initialAltitude} km`, 'info');
    log(`Perturbation: ${config.perturbationSize} m`, 'info');
    log(`Duration: ${(config.duration / 3600).toFixed(1)} hours`, 'info');

    const startTime = performance.now();
    const result = computeSensitivity(config);
    const elapsed = performance.now() - startTime;
    log(`Analysis completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Max position sensitivity: ${result.maxPosSensitivity.toFixed(1)}`, 'info');
    log(`Max velocity sensitivity: ${result.maxVelSensitivity.toFixed(4)}`, 'info');

    renderSensitivityFigures(chartContainer, result, config);

    return {
        name: 'Linear Sensitivity',
        description: 'State transition matrix analysis',
        ...result,
        config
    };
}

function computeSensitivity(config) {
    const RE = 6378.137;  // Earth radius km
    const mu = 3.986004418e5;  // km³/s²

    // Initial state (circular orbit)
    const r0 = RE + config.initialAltitude;
    const v0 = Math.sqrt(mu / r0);
    const period = 2 * Math.PI * Math.sqrt(r0 * r0 * r0 / mu);

    // Reference trajectory
    const refTrajectory = [];
    const perturbedTrajectories = [];

    // Perturbation directions (position x, y, z and velocity x, y, z)
    const perturbations = [
        { name: 'δx', dx: 1, dy: 0, dz: 0, dvx: 0, dvy: 0, dvz: 0 },
        { name: 'δy', dx: 0, dy: 1, dz: 0, dvx: 0, dvy: 0, dvz: 0 },
        { name: 'δz', dx: 0, dy: 0, dz: 1, dvx: 0, dvy: 0, dvz: 0 },
        { name: 'δvx', dx: 0, dy: 0, dz: 0, dvx: 1, dvy: 0, dvz: 0 },
        { name: 'δvy', dx: 0, dy: 0, dz: 0, dvx: 0, dvy: 1, dvz: 0 },
        { name: 'δvz', dx: 0, dy: 0, dz: 0, dvx: 0, dvy: 0, dvz: 1 }
    ];

    const dt = config.duration / config.numPoints;
    const pertSize = config.perturbationSize / 1000;  // Convert to km
    const velPertSize = 0.001;  // 1 mm/s for velocity perturbations

    // State derivative function
    const derivative = (state) => {
        const r = Math.sqrt(state[0]**2 + state[1]**2 + state[2]**2);
        const a = -mu / (r * r * r);
        return [
            state[3], state[4], state[5],
            a * state[0], a * state[1], a * state[2]
        ];
    };

    // RK4 integrator
    const rk4Step = (state, h) => {
        const k1 = derivative(state);
        const s2 = state.map((v, i) => v + h/2 * k1[i]);
        const k2 = derivative(s2);
        const s3 = state.map((v, i) => v + h/2 * k2[i]);
        const k3 = derivative(s3);
        const s4 = state.map((v, i) => v + h * k3[i]);
        const k4 = derivative(s4);
        return state.map((v, i) => v + h/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
    };

    // Initial reference state
    let refState = [r0, 0, 0, 0, v0, 0];

    // Initialize perturbed states
    const perturbedStates = perturbations.map(p => [
        r0 + p.dx * pertSize,
        p.dy * pertSize,
        p.dz * pertSize,
        p.dvx * velPertSize,
        v0 + p.dvy * velPertSize,
        p.dvz * velPertSize
    ]);

    // State transition matrix elements over time
    const stmHistory = [];
    let maxPosSensitivity = 0;
    let maxVelSensitivity = 0;

    const stepSize = 10;  // Integration step in seconds

    for (let i = 0; i < config.numPoints; i++) {
        const t = i * dt;

        // Store reference trajectory
        refTrajectory.push({
            t: t / 3600,  // hours
            x: refState[0],
            y: refState[1],
            z: refState[2]
        });

        // Compute STM elements from perturbed trajectories
        const stmRow = { t: t / 3600 };

        perturbations.forEach((p, j) => {
            const delta = [
                perturbedStates[j][0] - refState[0],
                perturbedStates[j][1] - refState[1],
                perturbedStates[j][2] - refState[2],
                perturbedStates[j][3] - refState[3],
                perturbedStates[j][4] - refState[4],
                perturbedStates[j][5] - refState[5]
            ];

            // Normalize by perturbation size
            const scale = j < 3 ? pertSize : velPertSize;
            const normalized = delta.map(d => d / scale);

            stmRow[`pos_${p.name}`] = Math.sqrt(normalized[0]**2 + normalized[1]**2 + normalized[2]**2);
            stmRow[`vel_${p.name}`] = Math.sqrt(normalized[3]**2 + normalized[4]**2 + normalized[5]**2);

            if (j < 3) {
                maxPosSensitivity = Math.max(maxPosSensitivity, stmRow[`pos_${p.name}`]);
            } else {
                maxVelSensitivity = Math.max(maxVelSensitivity, stmRow[`vel_${p.name}`]);
            }
        });

        stmHistory.push(stmRow);

        // Integrate reference state
        const numSteps = Math.ceil(dt / stepSize);
        const h = dt / numSteps;
        for (let s = 0; s < numSteps; s++) {
            refState = rk4Step(refState, h);
        }

        // Integrate perturbed states
        for (let j = 0; j < perturbedStates.length; j++) {
            for (let s = 0; s < numSteps; s++) {
                perturbedStates[j] = rk4Step(perturbedStates[j], h);
            }
        }
    }

    return {
        refTrajectory,
        stmHistory,
        perturbations,
        maxPosSensitivity,
        maxVelSensitivity,
        period: period / 3600  // hours
    };
}

function renderSensitivityFigures(container, result, config) {
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
    const chartHeight = Math.max(180, (containerHeight - 100) / 3);

    // Figure 1: Position sensitivity to position perturbations
    createFigure(wrapper, 'Position Sensitivity to Initial Position', chartWidth, chartHeight, (svg, w, h) => {
        renderMultiLine(svg, result.stmHistory, 't', ['pos_δx', 'pos_δy', 'pos_δz'],
            ['var(--cyan)', 'var(--green)', 'var(--orange)'],
            ['δx₀', 'δy₀', 'δz₀'],
            'Time [hours]', 'δr/δr₀', w, h);
    });

    // Figure 2: Position sensitivity to velocity perturbations
    createFigure(wrapper, 'Position Sensitivity to Initial Velocity', chartWidth, chartHeight, (svg, w, h) => {
        renderMultiLine(svg, result.stmHistory, 't', ['pos_δvx', 'pos_δvy', 'pos_δvz'],
            ['var(--purple)', 'var(--red)', 'var(--yellow)'],
            ['δvx₀', 'δvy₀', 'δvz₀'],
            'Time [hours]', 'δr/δv₀ [s]', w, h);
    });

    // Figure 3: Velocity sensitivity
    createFigure(wrapper, 'Velocity Sensitivity to Initial Velocity', chartWidth, chartHeight, (svg, w, h) => {
        renderMultiLine(svg, result.stmHistory, 't', ['vel_δvx', 'vel_δvy', 'vel_δvz'],
            ['var(--purple)', 'var(--red)', 'var(--yellow)'],
            ['δvx₀', 'δvy₀', 'δvz₀'],
            'Time [hours]', 'δv/δv₀', w, h);
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

function renderMultiLine(svg, data, xKey, yKeys, colors, labels, xLabel, yLabel, width, height) {
    const margin = { top: 10, right: 100, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain(d3.extent(data, d => d[xKey]))
        .range([margin.left, width - margin.right]);

    // Find y extent across all series
    let yMin = Infinity, yMax = -Infinity;
    yKeys.forEach(key => {
        data.forEach(d => {
            if (d[key] !== undefined) {
                yMin = Math.min(yMin, d[key]);
                yMax = Math.max(yMax, d[key]);
            }
        });
    });
    const yPadding = (yMax - yMin) * 0.1 || 1;

    const y = d3.scaleLinear()
        .domain([Math.max(0, yMin - yPadding), yMax + yPadding])
        .range([height - margin.bottom, margin.top]);

    // Draw lines
    yKeys.forEach((key, i) => {
        const lineData = data.filter(d => d[key] !== undefined);
        svg.append('path')
            .datum(lineData)
            .attr('fill', 'none')
            .attr('stroke', colors[i])
            .attr('stroke-width', 1.5)
            .attr('d', d3.line().x(d => x(d[xKey])).y(d => y(d[key])));
    });

    // Legend
    const legend = svg.append('g')
        .attr('transform', `translate(${width - margin.right + 10}, ${margin.top})`);

    labels.forEach((label, i) => {
        legend.append('line')
            .attr('x1', 0).attr('y1', i * 18)
            .attr('x2', 20).attr('y2', i * 18)
            .attr('stroke', colors[i])
            .attr('stroke-width', 2);

        legend.append('text')
            .attr('x', 25).attr('y', i * 18 + 4)
            .attr('fill', 'var(--text-secondary)')
            .attr('font-size', '10px')
            .text(label);
    });

    // Axes
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
        .text(xLabel);

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text(yLabel);
}
