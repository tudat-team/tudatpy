/**
 * Cassini MGA Trajectory Optimization Example
 * Port of: tudatpy/examples/mission_design/cassini1_mga_optimization.py
 *
 * Demonstrates trajectory optimization for a simplified Cassini mission
 * using multiple gravity assists (Earth-Venus-Venus-Earth-Jupiter-Saturn).
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

export function showCassiniMGAExample(chartContainer, log, params = {}) {
    const config = {
        populationSize: params.populationSize ?? 30,
        generations: params.generations ?? 200,
        showEvolution: params.showEvolution ?? true
    };

    log('Running Cassini MGA Optimization Example...', 'info');
    log('Route: Earth → Venus → Venus → Earth → Jupiter → Saturn', 'info');
    log(`Optimization: ${config.populationSize} individuals, ${config.generations} generations`, 'info');

    // Check SPICE availability
    const useSpice = isSpiceReady();
    if (useSpice) {
        log('Using SPICE ephemeris for planetary positions', 'success');
    } else {
        log('SPICE not available - using analytical Keplerian approximation', 'warning');
    }

    const startTime = performance.now();
    const result = useSpice
        ? optimizeCassiniTrajectorySpice(config, log)
        : optimizeCassiniTrajectory(config);
    const elapsed = performance.now() - startTime;
    log(`Optimization completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Best ΔV: ${(result.bestDeltaV / 1000).toFixed(2)} km/s`, 'info');
    log(`Total flight time: ${result.totalFlightTime.toFixed(0)} days`, 'info');

    renderCassiniFigures(chartContainer, result, config, useSpice);

    return {
        name: 'Cassini MGA',
        description: 'Multi-gravity assist trajectory',
        useSpice,
        ...result,
        config
    };
}

function optimizeCassiniTrajectory(config) {
    // Planetary data (simplified circular orbits)
    const AU = 1.496e11;  // Astronomical Unit in meters
    const planets = {
        Earth: { a: 1.0 * AU, period: 365.25 },
        Venus: { a: 0.723 * AU, period: 224.7 },
        Jupiter: { a: 5.203 * AU, period: 4332.6 },
        Saturn: { a: 9.537 * AU, period: 10759.2 }
    };

    // Transfer sequence
    const sequence = ['Earth', 'Venus', 'Venus', 'Earth', 'Jupiter', 'Saturn'];

    // Bounds for design variables (days)
    const departureDateBounds = [0, 1000];  // Departure window
    const legTimeBounds = [
        [30, 400],    // Earth-Venus
        [100, 470],   // Venus-Venus
        [30, 400],    // Venus-Earth
        [400, 2000],  // Earth-Jupiter
        [1000, 6000]  // Jupiter-Saturn
    ];

    // Simple evolutionary optimization (Differential Evolution style)
    const dim = 6;  // 1 departure + 5 leg times
    const pop = [];
    const fitnessHistory = [];

    // Initialize population
    for (let i = 0; i < config.populationSize; i++) {
        const individual = [
            departureDateBounds[0] + Math.random() * (departureDateBounds[1] - departureDateBounds[0])
        ];
        for (let j = 0; j < 5; j++) {
            individual.push(
                legTimeBounds[j][0] + Math.random() * (legTimeBounds[j][1] - legTimeBounds[j][0])
            );
        }
        pop.push({
            x: individual,
            fitness: evaluateTrajectory(individual, sequence, planets)
        });
    }

    // Evolution
    for (let gen = 0; gen < config.generations; gen++) {
        // Sort by fitness
        pop.sort((a, b) => a.fitness - b.fitness);

        // Record best
        fitnessHistory.push({
            generation: gen,
            bestFitness: pop[0].fitness,
            avgFitness: pop.reduce((s, p) => s + p.fitness, 0) / pop.length
        });

        // Generate new population using differential evolution
        const newPop = [pop[0]];  // Keep best

        for (let i = 1; i < config.populationSize; i++) {
            // Select 3 random individuals
            const a = pop[Math.floor(Math.random() * config.populationSize)];
            const b = pop[Math.floor(Math.random() * config.populationSize)];
            const c = pop[Math.floor(Math.random() * config.populationSize)];

            // Mutation
            const F = 0.5;  // Scaling factor
            const CR = 0.7; // Crossover rate
            const trial = [];

            for (let j = 0; j < dim; j++) {
                if (Math.random() < CR) {
                    let val = a.x[j] + F * (b.x[j] - c.x[j]);
                    // Bounds check
                    if (j === 0) {
                        val = Math.max(departureDateBounds[0], Math.min(departureDateBounds[1], val));
                    } else {
                        val = Math.max(legTimeBounds[j-1][0], Math.min(legTimeBounds[j-1][1], val));
                    }
                    trial.push(val);
                } else {
                    trial.push(pop[i].x[j]);
                }
            }

            const trialFitness = evaluateTrajectory(trial, sequence, planets);
            if (trialFitness < pop[i].fitness) {
                newPop.push({ x: trial, fitness: trialFitness });
            } else {
                newPop.push(pop[i]);
            }
        }

        pop.length = 0;
        pop.push(...newPop);
    }

    // Get best solution
    pop.sort((a, b) => a.fitness - b.fitness);
    const best = pop[0];

    // Generate trajectory for visualization
    const trajectory = generateTrajectory(best.x, sequence, planets);

    // Calculate total flight time
    const totalFlightTime = best.x.slice(1).reduce((a, b) => a + b, 0);

    return {
        bestDeltaV: best.fitness,
        bestParameters: best.x,
        totalFlightTime,
        trajectory,
        flybyPoints: trajectory.flybys,
        fitnessHistory,
        sequence
    };
}

/**
 * Optimize Cassini trajectory using SPICE ephemeris
 */
function optimizeCassiniTrajectorySpice(config, log) {
    const AU = 1.496e11;  // Astronomical Unit in meters

    // Transfer sequence
    const sequence = ['Earth', 'Venus', 'Venus', 'Earth', 'Jupiter', 'Saturn'];

    // Cassini launch window: October 1997
    // J2000 epoch is Jan 1, 2000, so Oct 15, 1997 is roughly -820 days before J2000
    const launchJ2000Days = -820;  // Days from J2000 (Oct 1997)

    // Bounds for design variables (days)
    const departureDateBounds = [launchJ2000Days - 100, launchJ2000Days + 200];
    const legTimeBounds = [
        [30, 400],    // Earth-Venus
        [100, 470],   // Venus-Venus
        [30, 400],    // Venus-Earth
        [400, 2000],  // Earth-Jupiter
        [1000, 6000]  // Jupiter-Saturn
    ];

    // Planet colors for visualization
    const planetColors = {
        Earth: '#4488ff',
        Venus: '#ffaa44',
        Jupiter: '#ff6644',
        Saturn: '#ffdd88'
    };

    // Simple evolutionary optimization (Differential Evolution style)
    const dim = 6;  // 1 departure + 5 leg times
    const pop = [];
    const fitnessHistory = [];

    // Initialize population
    for (let i = 0; i < config.populationSize; i++) {
        const individual = [
            departureDateBounds[0] + Math.random() * (departureDateBounds[1] - departureDateBounds[0])
        ];
        for (let j = 0; j < 5; j++) {
            individual.push(
                legTimeBounds[j][0] + Math.random() * (legTimeBounds[j][1] - legTimeBounds[j][0])
            );
        }
        pop.push({
            x: individual,
            fitness: evaluateTrajectorySpice(individual, sequence, AU)
        });
    }

    // Evolution
    for (let gen = 0; gen < config.generations; gen++) {
        // Sort by fitness
        pop.sort((a, b) => a.fitness - b.fitness);

        // Record best
        fitnessHistory.push({
            generation: gen,
            bestFitness: pop[0].fitness,
            avgFitness: pop.reduce((s, p) => s + p.fitness, 0) / pop.length
        });

        // Generate new population using differential evolution
        const newPop = [pop[0]];  // Keep best

        for (let i = 1; i < config.populationSize; i++) {
            // Select 3 random individuals
            const a = pop[Math.floor(Math.random() * config.populationSize)];
            const b = pop[Math.floor(Math.random() * config.populationSize)];
            const c = pop[Math.floor(Math.random() * config.populationSize)];

            // Mutation
            const F = 0.5;  // Scaling factor
            const CR = 0.7; // Crossover rate
            const trial = [];

            for (let j = 0; j < dim; j++) {
                if (Math.random() < CR) {
                    let val = a.x[j] + F * (b.x[j] - c.x[j]);
                    // Bounds check
                    if (j === 0) {
                        val = Math.max(departureDateBounds[0], Math.min(departureDateBounds[1], val));
                    } else {
                        val = Math.max(legTimeBounds[j-1][0], Math.min(legTimeBounds[j-1][1], val));
                    }
                    trial.push(val);
                } else {
                    trial.push(pop[i].x[j]);
                }
            }

            const trialFitness = evaluateTrajectorySpice(trial, sequence, AU);
            if (trialFitness < pop[i].fitness) {
                newPop.push({ x: trial, fitness: trialFitness });
            } else {
                newPop.push(pop[i]);
            }
        }

        pop.length = 0;
        pop.push(...newPop);
    }

    // Get best solution
    pop.sort((a, b) => a.fitness - b.fitness);
    const best = pop[0];

    // Generate trajectory for visualization using SPICE
    const trajectory = generateTrajectorySpice(best.x, sequence, AU);

    // Calculate total flight time
    const totalFlightTime = best.x.slice(1).reduce((a, b) => a + b, 0);

    // Build planets object for compatibility
    const planets = {};
    for (const name of ['Earth', 'Venus', 'Jupiter', 'Saturn']) {
        planets[name] = {
            a: (PLANETARY_SMA[name] || 1e11) / AU,
            period: (PLANETARY_PERIOD[name] || 365.25 * 86400) / 86400,
            color: planetColors[name]
        };
    }

    return {
        bestDeltaV: best.fitness,
        bestParameters: best.x,
        totalFlightTime,
        trajectory,
        flybyPoints: trajectory.flybys,
        fitnessHistory,
        sequence,
        source: 'SPICE'
    };
}

/**
 * Evaluate trajectory using SPICE ephemeris
 */
function evaluateTrajectorySpice(params, sequence, AU) {
    let totalDeltaV = 0;
    const departureDay = params[0];

    // Calculate node times (days from J2000)
    const nodeTimes = [departureDay];
    for (let i = 1; i < params.length; i++) {
        nodeTimes.push(nodeTimes[i-1] + params[i]);
    }

    const mu = 1.327e20;  // Sun's gravitational parameter (m³/s²)

    // Evaluate each leg
    for (let leg = 0; leg < sequence.length - 1; leg++) {
        const from = sequence[leg];
        const to = sequence[leg + 1];
        const tof = params[leg + 1];  // Time of flight in days

        const t1 = nodeTimes[leg];
        const t2 = nodeTimes[leg + 1];

        // Convert to ephemeris time (seconds from J2000)
        const et1 = t1 * 86400;
        const et2 = t2 * 86400;

        // Get planet positions from SPICE
        const state1 = getBodyState(from, 'Sun', et1, 'ECLIPJ2000');
        const state2 = getBodyState(to, 'Sun', et2, 'ECLIPJ2000');

        if (!state1 || !state2) {
            return 1e10;  // Penalty if SPICE query fails
        }

        // Positions in meters
        const pos1 = { x: state1.x, y: state1.y, z: state1.z };
        const pos2 = { x: state2.x, y: state2.y, z: state2.z };

        const r1 = Math.sqrt(pos1.x*pos1.x + pos1.y*pos1.y + pos1.z*pos1.z);
        const r2 = Math.sqrt(pos2.x*pos2.x + pos2.y*pos2.y + pos2.z*pos2.z);

        // Transfer angle
        const dotProd = pos1.x*pos2.x + pos1.y*pos2.y + pos1.z*pos2.z;
        const cosTheta = dotProd / (r1 * r2);
        const theta = Math.acos(Math.max(-1, Math.min(1, cosTheta)));

        // Approximate transfer geometry (chord)
        const c = Math.sqrt(r1*r1 + r2*r2 - 2*r1*r2*cosTheta);
        const s = (r1 + r2 + c) / 2;

        // Circular velocities
        const v_circular_from = Math.sqrt(mu / r1);
        const v_circular_to = Math.sqrt(mu / r2);

        // Transfer velocities (simplified Hohmann-like)
        const a_transfer = (r1 + r2) / 2;
        const v_transfer_1 = Math.sqrt(mu * (2/r1 - 1/a_transfer));
        const v_transfer_2 = Math.sqrt(mu * (2/r2 - 1/a_transfer));

        // ΔV at departure
        const dv1 = Math.abs(v_transfer_1 - v_circular_from);

        // ΔV at arrival (or gravity assist)
        let dv2;
        if (leg === sequence.length - 2) {
            // Final insertion
            dv2 = Math.abs(v_transfer_2 - v_circular_to) * 0.5;
        } else {
            // Gravity assist - reduce ΔV cost
            dv2 = Math.abs(v_transfer_2 - v_circular_to) * 0.1;
        }

        totalDeltaV += dv1 + dv2;
    }

    // Add penalty for unrealistic time of flight
    const totalTof = params.slice(1).reduce((a, b) => a + b, 0);
    if (totalTof < 1000 || totalTof > 8000) {
        totalDeltaV += 10000;
    }

    return totalDeltaV;
}

/**
 * Generate trajectory points using SPICE ephemeris
 */
function generateTrajectorySpice(params, sequence, AU) {
    const trajectory = [];
    const flybys = [];

    // Calculate node times (days from J2000)
    const nodeTimes = [params[0]];
    for (let i = 1; i < params.length; i++) {
        nodeTimes.push(nodeTimes[i-1] + params[i]);
    }

    // Generate points for each leg
    for (let leg = 0; leg < sequence.length - 1; leg++) {
        const from = sequence[leg];
        const to = sequence[leg + 1];
        const tof = params[leg + 1];

        const t1 = nodeTimes[leg];
        const t2 = nodeTimes[leg + 1];

        // Get planet positions from SPICE
        const et1 = t1 * 86400;
        const et2 = t2 * 86400;

        const state1 = getBodyState(from, 'Sun', et1, 'ECLIPJ2000');
        const state2 = getBodyState(to, 'Sun', et2, 'ECLIPJ2000');

        if (!state1 || !state2) continue;

        const pos1 = { x: state1.x / AU, y: state1.y / AU };
        const pos2 = { x: state2.x / AU, y: state2.y / AU };

        // Record flyby
        flybys.push({
            planet: from,
            x: pos1.x,
            y: pos1.y,
            day: t1
        });

        // Interpolate trajectory (curved arc)
        const numPoints = Math.min(100, Math.ceil(tof / 5));
        for (let i = 0; i <= numPoints; i++) {
            const frac = i / numPoints;
            const angle1 = Math.atan2(pos1.y, pos1.x);
            const angle2 = Math.atan2(pos2.y, pos2.x);
            let dAngle = angle2 - angle1;
            if (dAngle > Math.PI) dAngle -= 2 * Math.PI;
            if (dAngle < -Math.PI) dAngle += 2 * Math.PI;

            const r1 = Math.sqrt(pos1.x*pos1.x + pos1.y*pos1.y);
            const r2 = Math.sqrt(pos2.x*pos2.x + pos2.y*pos2.y);

            const angle = angle1 + frac * dAngle;
            const r = r1 + frac * (r2 - r1);

            trajectory.push({
                x: r * Math.cos(angle),
                y: r * Math.sin(angle),
                leg: leg
            });
        }
    }

    // Add final destination
    const lastT = nodeTimes[nodeTimes.length - 1];
    const lastState = getBodyState(sequence[sequence.length - 1], 'Sun', lastT * 86400, 'ECLIPJ2000');
    if (lastState) {
        flybys.push({
            planet: sequence[sequence.length - 1],
            x: lastState.x / AU,
            y: lastState.y / AU,
            day: lastT
        });
    }

    return { points: trajectory, flybys };
}

function evaluateTrajectory(params, sequence, planets) {
    // Simplified trajectory evaluation using patched conics
    let totalDeltaV = 0;
    const departureDate = params[0];

    // Calculate node times
    const nodeTimes = [departureDate];
    for (let i = 1; i < params.length; i++) {
        nodeTimes.push(nodeTimes[i-1] + params[i]);
    }

    // Evaluate each leg
    for (let leg = 0; leg < sequence.length - 1; leg++) {
        const from = sequence[leg];
        const to = sequence[leg + 1];
        const tof = params[leg + 1];  // Time of flight in days

        // Get planet positions at departure and arrival
        const t1 = nodeTimes[leg];
        const t2 = nodeTimes[leg + 1];

        const pos1 = getPlanetPosition(from, t1, planets);
        const pos2 = getPlanetPosition(to, t2, planets);

        // Calculate transfer orbit (simplified Lambert problem)
        const r1 = Math.sqrt(pos1.x*pos1.x + pos1.y*pos1.y);
        const r2 = Math.sqrt(pos2.x*pos2.x + pos2.y*pos2.y);

        // Transfer angle
        const cosTheta = (pos1.x*pos2.x + pos1.y*pos2.y) / (r1 * r2);
        const theta = Math.acos(Math.max(-1, Math.min(1, cosTheta)));

        // Semi-major axis of transfer orbit (vis-viva approximation)
        const mu = 1.327e20;  // Sun's gravitational parameter (m³/s²)
        const tofSec = tof * 86400;

        // Approximate transfer geometry
        const c = Math.sqrt(r1*r1 + r2*r2 - 2*r1*r2*cosTheta);  // Chord
        const s = (r1 + r2 + c) / 2;  // Semi-perimeter

        // Minimum energy transfer
        const a_min = s / 2;
        const tof_min = Math.PI * Math.sqrt(a_min*a_min*a_min / mu);

        // Estimate ΔV based on transfer geometry
        const v_circular_from = Math.sqrt(mu / r1);
        const v_circular_to = Math.sqrt(mu / r2);

        // Transfer velocities (simplified Hohmann-like)
        const a_transfer = (r1 + r2) / 2;
        const v_transfer_1 = Math.sqrt(mu * (2/r1 - 1/a_transfer));
        const v_transfer_2 = Math.sqrt(mu * (2/r2 - 1/a_transfer));

        // ΔV at departure (excess velocity relative to planet)
        const dv1 = Math.abs(v_transfer_1 - v_circular_from);

        // ΔV at arrival (or gravity assist)
        let dv2;
        if (leg === sequence.length - 2) {
            // Final insertion
            dv2 = Math.abs(v_transfer_2 - v_circular_to) * 0.5;  // Capture ΔV
        } else {
            // Gravity assist - reduce ΔV cost
            dv2 = Math.abs(v_transfer_2 - v_circular_to) * 0.1;  // GA efficiency
        }

        totalDeltaV += dv1 + dv2;
    }

    // Add penalty for unrealistic time of flight
    const totalTof = params.slice(1).reduce((a, b) => a + b, 0);
    if (totalTof < 1000 || totalTof > 8000) {
        totalDeltaV += 10000;  // Penalty
    }

    return totalDeltaV;
}

function getPlanetPosition(planet, daysSinceEpoch, planets) {
    const data = planets[planet];
    const angle = 2 * Math.PI * daysSinceEpoch / data.period;
    return {
        x: data.a * Math.cos(angle),
        y: data.a * Math.sin(angle)
    };
}

function generateTrajectory(params, sequence, planets) {
    const AU = 1.496e11;
    const trajectory = [];
    const flybys = [];

    // Calculate node times
    const nodeTimes = [params[0]];
    for (let i = 1; i < params.length; i++) {
        nodeTimes.push(nodeTimes[i-1] + params[i]);
    }

    // Generate points for each leg
    for (let leg = 0; leg < sequence.length - 1; leg++) {
        const from = sequence[leg];
        const to = sequence[leg + 1];
        const tof = params[leg + 1];

        const t1 = nodeTimes[leg];
        const t2 = nodeTimes[leg + 1];

        const pos1 = getPlanetPosition(from, t1, planets);
        const pos2 = getPlanetPosition(to, t2, planets);

        // Record flyby
        flybys.push({
            planet: from,
            x: pos1.x / AU,
            y: pos1.y / AU,
            day: t1
        });

        // Interpolate trajectory (simplified arc)
        const numPoints = Math.min(100, Math.ceil(tof / 5));
        for (let i = 0; i <= numPoints; i++) {
            const frac = i / numPoints;
            // Curved path (not straight line)
            const angle1 = Math.atan2(pos1.y, pos1.x);
            const angle2 = Math.atan2(pos2.y, pos2.x);
            let dAngle = angle2 - angle1;
            if (dAngle > Math.PI) dAngle -= 2 * Math.PI;
            if (dAngle < -Math.PI) dAngle += 2 * Math.PI;

            const r1 = Math.sqrt(pos1.x*pos1.x + pos1.y*pos1.y);
            const r2 = Math.sqrt(pos2.x*pos2.x + pos2.y*pos2.y);

            const angle = angle1 + frac * dAngle;
            const r = r1 + frac * (r2 - r1);  // Linear interpolation of radius

            trajectory.push({
                x: r * Math.cos(angle) / AU,
                y: r * Math.sin(angle) / AU,
                leg: leg
            });
        }
    }

    // Add final destination
    const lastPos = getPlanetPosition(sequence[sequence.length - 1], nodeTimes[nodeTimes.length - 1], planets);
    flybys.push({
        planet: sequence[sequence.length - 1],
        x: lastPos.x / AU,
        y: lastPos.y / AU,
        day: nodeTimes[nodeTimes.length - 1]
    });

    return { points: trajectory, flybys };
}

function renderCassiniFigures(container, result, config, useSpice = false) {
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

    // Figure 1: Trajectory
    const title1 = useSpice
        ? 'Cassini Transfer Trajectory (SPICE Ephemeris)'
        : 'Cassini Transfer Trajectory (Sun-centered)';
    createFigure(wrapper, title1, chartWidth, chartHeight, (svg, w, h) => {
        renderCassiniTrajectory(svg, result, w, h);
    });

    // Figure 2: Optimization convergence
    createFigure(wrapper, 'Optimization Convergence', chartWidth, chartHeight, (svg, w, h) => {
        renderConvergence(svg, result, w, h);
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

function renderCassiniTrajectory(svg, result, width, height) {
    const margin = { top: 20, right: 100, bottom: 40, left: 60 };

    // Scale to fit Saturn's orbit
    const scale = 12;  // AU

    const x = d3.scaleLinear()
        .domain([-scale, scale])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-scale, scale])
        .range([height - margin.bottom, margin.top]);

    // Draw planet orbits
    const orbitColors = {
        Earth: '#4488ff',
        Venus: '#ffaa44',
        Jupiter: '#ff6644',
        Saturn: '#ffdd88'
    };

    const orbitRadii = {
        Earth: 1.0,
        Venus: 0.723,
        Jupiter: 5.203,
        Saturn: 9.537
    };

    Object.entries(orbitRadii).forEach(([planet, r]) => {
        svg.append('circle')
            .attr('cx', x(0))
            .attr('cy', y(0))
            .attr('r', x(r) - x(0))
            .attr('fill', 'none')
            .attr('stroke', orbitColors[planet])
            .attr('stroke-width', 1)
            .attr('stroke-dasharray', '3,3')
            .attr('opacity', 0.4);
    });

    // Draw Sun
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', 6)
        .attr('fill', '#ffcc00');

    // Draw trajectory
    svg.append('path')
        .datum(result.trajectory.points)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 1.5)
        .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

    // Draw flyby markers
    const flybyColors = {
        Earth: '#4488ff',
        Venus: '#ffaa44',
        Jupiter: '#ff6644',
        Saturn: '#aaaaaa'
    };

    result.flybyPoints.forEach((fb, i) => {
        svg.append('circle')
            .attr('cx', x(fb.x))
            .attr('cy', y(fb.y))
            .attr('r', 5)
            .attr('fill', flybyColors[fb.planet] || 'white')
            .attr('stroke', 'white')
            .attr('stroke-width', 1);
    });

    // Legend
    const legend = svg.append('g')
        .attr('transform', `translate(${width - margin.right + 10}, ${margin.top})`);

    const legendItems = [
        { name: 'Earth', color: '#4488ff' },
        { name: 'Venus', color: '#ffaa44' },
        { name: 'Jupiter', color: '#ff6644' },
        { name: 'Saturn', color: '#aaaaaa' }
    ];

    legendItems.forEach((item, i) => {
        legend.append('circle')
            .attr('cx', 8)
            .attr('cy', i * 18)
            .attr('r', 4)
            .attr('fill', item.color);

        legend.append('text')
            .attr('x', 18)
            .attr('y', i * 18 + 4)
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

function renderConvergence(svg, result, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain([0, result.fitnessHistory.length])
        .range([margin.left, width - margin.right]);

    const yExtent = d3.extent(result.fitnessHistory, d => d.bestFitness / 1000);
    const y = d3.scaleLinear()
        .domain([yExtent[0] * 0.9, Math.min(yExtent[1] * 1.1, 25)])
        .range([height - margin.bottom, margin.top]);

    // Best fitness
    svg.append('path')
        .datum(result.fitnessHistory)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 2)
        .attr('d', d3.line()
            .x(d => x(d.generation))
            .y(d => y(d.bestFitness / 1000)));

    // Mark champion
    const champion = result.fitnessHistory.reduce((best, d) =>
        d.bestFitness < best.bestFitness ? d : best
    );

    svg.append('circle')
        .attr('cx', x(champion.generation))
        .attr('cy', y(champion.bestFitness / 1000))
        .attr('r', 5)
        .attr('fill', 'var(--orange)');

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
        .text('Generation');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('ΔV [km/s]');
}
