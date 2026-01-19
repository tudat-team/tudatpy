/**
 * Earth-Moon Low Thrust Transfer Example
 * Port of: tudatpy/examples/propagation/thrust_between_Earth_Moon.py
 *
 * Demonstrates constant thrust propagation in the Earth-Moon system.
 * Shows vehicle trajectory, mass depletion, and altitude evolution.
 */

export function showEarthMoonThrustExample(chartContainer, log, params = {}) {
    const config = {
        thrustMagnitude: params.thrustMagnitude ?? 10,      // N
        specificImpulse: params.specificImpulse ?? 5000,    // s
        initialMass: params.initialMass ?? 5000,            // kg
        duration: params.duration ?? 30,                     // days
        numPoints: params.numPoints ?? 3000
    };

    log('Running Earth-Moon Thrust Transfer Example...', 'info');
    log(`Thrust: ${config.thrustMagnitude} N, Isp: ${config.specificImpulse} s`, 'info');
    log(`Initial mass: ${config.initialMass} kg`, 'info');

    const startTime = performance.now();
    const result = simulateEarthMoonTransfer(config);
    const elapsed = performance.now() - startTime;
    log(`Simulation completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Final altitude: ${(result.finalAltitude / 1000).toFixed(0)} km`, 'info');
    log(`Final mass: ${result.finalMass.toFixed(1)} kg`, 'info');
    log(`Propellant used: ${(config.initialMass - result.finalMass).toFixed(1)} kg`, 'info');
    log(`Termination: ${result.termination}`, 'info');

    renderEarthMoonFigures(chartContainer, result, config);

    return {
        name: 'Earth-Moon Thrust',
        description: 'Low-thrust transfer in Earth-Moon system',
        ...result,
        config
    };
}

function simulateEarthMoonTransfer(config) {
    // Constants
    const G = 6.674e-11;
    const M_EARTH = 5.972e24;
    const M_MOON = 7.342e22;
    const M_SUN = 1.989e30;
    const R_EARTH = 6.371e6;
    const MOON_DISTANCE = 3.844e8;  // Average distance to Moon
    const MOON_PERIOD = 27.3 * 86400;  // Moon orbital period
    const g0 = 9.80665;

    const mu_earth = G * M_EARTH;
    const mu_moon = G * M_MOON;
    const mu_sun = G * M_SUN;

    // Initial state: 8000 km from Earth center, 7.5 km/s velocity
    let x = 8.0e6, y = 0, z = 0;
    let vx = 0, vy = 7.5e3, vz = 0;
    let mass = config.initialMass;

    const dt = config.duration * 86400 / config.numPoints;
    const trajectory = [];
    const moonTrajectory = [];

    // Mass flow rate
    const massFlowRate = config.thrustMagnitude / (config.specificImpulse * g0);

    let termination = 'time';
    const maxAltitude = 100e6;  // 100,000 km
    const minMass = 4000;       // kg

    for (let i = 0; i < config.numPoints; i++) {
        const t = i * dt;
        const day = t / 86400;

        // Current position
        const r_earth = Math.sqrt(x*x + y*y + z*z);
        const altitude = r_earth - R_EARTH;

        // Moon position (simplified circular orbit)
        const moonAngle = 2 * Math.PI * t / MOON_PERIOD;
        const moonX = MOON_DISTANCE * Math.cos(moonAngle);
        const moonY = MOON_DISTANCE * Math.sin(moonAngle);
        const moonZ = 0;

        // Distance to Moon
        const dx_moon = x - moonX;
        const dy_moon = y - moonY;
        const dz_moon = z - moonZ;
        const r_moon = Math.sqrt(dx_moon*dx_moon + dy_moon*dy_moon + dz_moon*dz_moon);

        // Store trajectory
        trajectory.push({
            t: day,
            x: x / 1e6,  // thousand km
            y: y / 1e6,
            z: z / 1e6,
            altitude: altitude / 1000,  // km
            mass: mass,
            r: r_earth / 1e6,
            v: Math.sqrt(vx*vx + vy*vy + vz*vz) / 1000  // km/s
        });

        if (i % 50 === 0) {
            moonTrajectory.push({
                x: moonX / 1e6,
                y: moonY / 1e6,
                z: moonZ / 1e6
            });
        }

        // Check termination conditions
        if (altitude > maxAltitude) {
            termination = 'altitude';
            break;
        }
        if (mass < minMass) {
            termination = 'mass';
            break;
        }

        // Gravitational acceleration from Earth
        const r3_earth = r_earth * r_earth * r_earth;
        let ax = -mu_earth * x / r3_earth;
        let ay = -mu_earth * y / r3_earth;
        let az = -mu_earth * z / r3_earth;

        // Gravitational acceleration from Moon
        const r3_moon = r_moon * r_moon * r_moon;
        ax -= mu_moon * dx_moon / r3_moon;
        ay -= mu_moon * dy_moon / r3_moon;
        az -= mu_moon * dz_moon / r3_moon;

        // Thrust acceleration (aligned with velocity)
        const v = Math.sqrt(vx*vx + vy*vy + vz*vz);
        if (v > 0 && mass > minMass) {
            const thrust_accel = config.thrustMagnitude / mass;
            ax += thrust_accel * vx / v;
            ay += thrust_accel * vy / v;
            az += thrust_accel * vz / v;
        }

        // RK4 integration
        const k1_x = vx, k1_y = vy, k1_z = vz;
        const k1_vx = ax, k1_vy = ay, k1_vz = az;

        const x2 = x + 0.5 * dt * k1_x;
        const y2 = y + 0.5 * dt * k1_y;
        const z2 = z + 0.5 * dt * k1_z;
        const vx2 = vx + 0.5 * dt * k1_vx;
        const vy2 = vy + 0.5 * dt * k1_vy;
        const vz2 = vz + 0.5 * dt * k1_vz;

        const r2 = Math.sqrt(x2*x2 + y2*y2 + z2*z2);
        const r2_3 = r2 * r2 * r2;
        let ax2 = -mu_earth * x2 / r2_3;
        let ay2 = -mu_earth * y2 / r2_3;
        let az2 = -mu_earth * z2 / r2_3;

        const v2 = Math.sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
        if (v2 > 0) {
            const thrust_accel2 = config.thrustMagnitude / mass;
            ax2 += thrust_accel2 * vx2 / v2;
            ay2 += thrust_accel2 * vy2 / v2;
            az2 += thrust_accel2 * vz2 / v2;
        }

        const k2_x = vx2, k2_y = vy2, k2_z = vz2;
        const k2_vx = ax2, k2_vy = ay2, k2_vz = az2;

        // Update state (simplified RK2)
        x += dt * k2_x;
        y += dt * k2_y;
        z += dt * k2_z;
        vx += dt * k2_vx;
        vy += dt * k2_vy;
        vz += dt * k2_vz;

        // Update mass
        mass -= massFlowRate * dt;
        if (mass < minMass) mass = minMass;
    }

    const finalPoint = trajectory[trajectory.length - 1];

    return {
        trajectory,
        moonTrajectory,
        finalAltitude: finalPoint.altitude * 1000,  // m
        finalMass: finalPoint.mass,
        termination
    };
}

function renderEarthMoonFigures(container, result, config) {
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

    // Figure 1: XY Trajectory
    createFigure(wrapper, 'Earth-Moon Transfer Trajectory', chartWidth, chartHeight, (svg, w, h) => {
        renderTrajectoryXY(svg, result, w, h);
    });

    // Figure 2: Altitude vs Time
    createFigure(wrapper, 'Altitude vs Time', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.trajectory, 't', 'altitude',
            'Time [days]', 'Altitude [km]', 'var(--cyan)', w, h);
    });

    // Figure 3: Mass vs Time
    createFigure(wrapper, 'Vehicle Mass vs Time', chartWidth, chartHeight, (svg, w, h) => {
        renderLineChart(svg, result.trajectory, 't', 'mass',
            'Time [days]', 'Mass [kg]', 'var(--orange)', w, h);
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
    const margin = { top: 20, right: 20, bottom: 40, left: 60 };

    // Find bounds including Moon trajectory
    let maxR = Math.max(...result.trajectory.map(p => Math.max(Math.abs(p.x), Math.abs(p.y))));
    if (result.moonTrajectory.length > 0) {
        maxR = Math.max(maxR, ...result.moonTrajectory.map(p => Math.max(Math.abs(p.x), Math.abs(p.y))));
    }
    const scale = maxR * 1.1;

    const x = d3.scaleLinear()
        .domain([-scale, scale])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-scale, scale])
        .range([height - margin.bottom, margin.top]);

    // Draw Moon orbit (dashed)
    if (result.moonTrajectory.length > 0) {
        svg.append('circle')
            .attr('cx', x(0))
            .attr('cy', y(0))
            .attr('r', x(384.4) - x(0))  // Moon distance in thousand km
            .attr('fill', 'none')
            .attr('stroke', '#666')
            .attr('stroke-width', 1)
            .attr('stroke-dasharray', '5,5')
            .attr('opacity', 0.5);
    }

    // Draw Earth
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', 8)
        .attr('fill', '#4488ff');

    // Draw vehicle trajectory
    svg.append('path')
        .datum(result.trajectory)
        .attr('fill', 'none')
        .attr('stroke', 'var(--green)')
        .attr('stroke-width', 1.5)
        .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

    // Draw Moon trajectory
    if (result.moonTrajectory.length > 1) {
        svg.append('path')
            .datum(result.moonTrajectory)
            .attr('fill', 'none')
            .attr('stroke', '#888')
            .attr('stroke-width', 1)
            .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));
    }

    // Mark start and end
    const start = result.trajectory[0];
    const end = result.trajectory[result.trajectory.length - 1];

    svg.append('circle')
        .attr('cx', x(start.x))
        .attr('cy', y(start.y))
        .attr('r', 4)
        .attr('fill', 'var(--cyan)');

    svg.append('circle')
        .attr('cx', x(end.x))
        .attr('cy', y(end.y))
        .attr('r', 4)
        .attr('fill', 'var(--orange)');

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

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Y [thousand km]');
}

function renderLineChart(svg, data, xKey, yKey, xLabel, yLabel, color, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain(d3.extent(data, d => d[xKey]))
        .range([margin.left, width - margin.right]);

    const yExtent = d3.extent(data, d => d[yKey]);
    const yPadding = (yExtent[1] - yExtent[0]) * 0.1 || 1;

    const y = d3.scaleLinear()
        .domain([yExtent[0] - yPadding, yExtent[1] + yPadding])
        .range([height - margin.bottom, margin.top]);

    svg.append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', color)
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d[xKey])).y(d => y(d[yKey])));

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
