/**
 * Hodographic Shaping MGA Optimization Example
 * Port of: tudatpy/examples/mission_design/hodographic_shaping_mga_optimization.py
 *
 * Demonstrates low-thrust MGA trajectory optimization using hodographic
 * shaping for Earth-Mars-Earth-Jupiter transfer.
 *
 * Uses SPICE ephemeris when available for accurate planetary positions.
 */

import {
    isSpiceReady,
    getBodyState,
    PLANETARY_SMA,
    PLANETARY_PERIOD
} from '../shared/spice-utils.js';

export function showHodographicShapingMGAExample(chartContainer, log, params = {}) {
    const config = {
        populationSize: params.populationSize ?? 100,
        numGenerations: params.numGenerations ?? 30,
        transferSequence: params.transferSequence ?? ['Earth', 'Mars', 'Earth', 'Jupiter']
    };

    log('Running Hodographic Shaping MGA Optimization...', 'info');
    log(`Transfer: ${config.transferSequence.join(' → ')}`, 'info');
    log(`Population: ${config.populationSize}, Generations: ${config.numGenerations}`, 'info');

    // Check SPICE availability
    const useSpice = isSpiceReady();
    if (useSpice) {
        log('Using SPICE ephemeris for planetary positions', 'success');
    } else {
        log('SPICE not available - using analytical Keplerian approximation', 'warning');
    }

    const startTime = performance.now();
    const result = useSpice
        ? optimizeHodographicMGASpice(config, log)
        : optimizeHodographicMGA(config);
    const elapsed = performance.now() - startTime;
    log(`Optimization completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Best ΔV: ${(result.bestDeltaV / 1000).toFixed(2)} km/s`, 'info');
    log(`Total time of flight: ${(result.totalTOF / 86400).toFixed(0)} days`, 'info');
    log(`Departure: MJD ${(result.departureEpoch / 86400).toFixed(0)}`, 'info');
    result.legTOFs.forEach((tof, i) => {
        log(`  Leg ${i + 1} TOF: ${(tof / 86400).toFixed(0)} days`, 'info');
    });

    renderHodographicMGAFigures(chartContainer, result, config, useSpice);

    return {
        name: 'Hodographic Shaping MGA',
        description: 'Low-thrust multi-gravity assist optimization',
        useSpice,
        ...result,
        config
    };
}

/**
 * Optimize hodographic MGA trajectory using SPICE ephemeris
 */
function optimizeHodographicMGASpice(config, log) {
    const AU = 1.496e11;
    const JULIAN_DAY = 86400;
    const GM_SUN = 1.327e20;

    // Get planet position from SPICE
    function getPlanetPositionSpice(name, epoch) {
        const et = epoch;  // Already in seconds from J2000
        const state = getBodyState(name, 'Sun', et, 'ECLIPJ2000');
        if (state) {
            return {
                x: state.x,
                y: state.y,
                vx: state.vx,
                vy: state.vy
            };
        }
        // Fallback to analytical if SPICE fails
        const fallbackA = PLANETARY_SMA[name] || AU;
        const fallbackT = PLANETARY_PERIOD[name] || 365.25 * JULIAN_DAY;
        const n = 2 * Math.PI / fallbackT;
        const angle = n * epoch;
        return {
            x: fallbackA * Math.cos(angle),
            y: fallbackA * Math.sin(angle),
            vx: -fallbackA * n * Math.sin(angle),
            vy: fallbackA * n * Math.cos(angle)
        };
    }

    // Hodographic shaping approximation for low-thrust leg
    function computeLegDeltaV(r1, r2, tof, numRevs = 0) {
        const a_transfer = (r1 + r2) / 2;
        const v1 = Math.sqrt(GM_SUN / r1);
        const v2 = Math.sqrt(GM_SUN / r2);
        const v_transfer = Math.sqrt(GM_SUN * (2 / r1 - 1 / a_transfer));
        const minTOF = Math.PI * Math.sqrt(a_transfer**3 / GM_SUN);
        const tofRatio = tof / minTOF;

        let deltaV = Math.abs(v_transfer - v1) + Math.abs(v2 - v_transfer);
        deltaV *= (1 + 0.1 * numRevs);

        if (tofRatio < 1) deltaV *= 2 * (1 - tofRatio) + 1;
        else if (tofRatio > 3) deltaV *= 1 + 0.1 * (tofRatio - 3);

        return deltaV;
    }

    // Optimization bounds
    const bounds = {
        departureEpoch: [9000 * JULIAN_DAY, 9200 * JULIAN_DAY],
        legTOF: [200 * JULIAN_DAY, 1200 * JULIAN_DAY],
        numRevolutions: [0, 2]
    };

    // Fitness function using SPICE positions
    function fitness(individual) {
        const departureEpoch = individual.departureEpoch;
        const legTOFs = individual.legTOFs;
        const numRevs = individual.numRevs;

        let totalDeltaV = 0;
        let epoch = departureEpoch;
        const nodePositions = [];
        const nodeEpochs = [epoch];

        for (let i = 0; i < config.transferSequence.length - 1; i++) {
            const fromPlanet = config.transferSequence[i];
            const toPlanet = config.transferSequence[i + 1];

            const pos1 = getPlanetPositionSpice(fromPlanet, epoch);
            epoch += legTOFs[i];
            const pos2 = getPlanetPositionSpice(toPlanet, epoch);

            const r1 = Math.sqrt(pos1.x**2 + pos1.y**2);
            const r2 = Math.sqrt(pos2.x**2 + pos2.y**2);

            const legDV = computeLegDeltaV(r1, r2, legTOFs[i], numRevs[i]);
            totalDeltaV += legDV;

            nodePositions.push({ planet: fromPlanet, ...pos1, epoch: nodeEpochs[i] });
            nodeEpochs.push(epoch);
        }

        const lastPlanet = config.transferSequence[config.transferSequence.length - 1];
        const lastPos = getPlanetPositionSpice(lastPlanet, epoch);
        nodePositions.push({ planet: lastPlanet, ...lastPos, epoch });

        return { totalDeltaV, nodePositions, nodeEpochs };
    }

    // Generate random individual
    function randomIndividual() {
        const numLegs = config.transferSequence.length - 1;
        return {
            departureEpoch: bounds.departureEpoch[0] + Math.random() * (bounds.departureEpoch[1] - bounds.departureEpoch[0]),
            legTOFs: Array(numLegs).fill(0).map(() =>
                bounds.legTOF[0] + Math.random() * (bounds.legTOF[1] - bounds.legTOF[0])
            ),
            swingbyAlts: Array(numLegs - 1).fill(1e6),
            numRevs: Array(numLegs).fill(0).map(() =>
                Math.floor(Math.random() * (bounds.numRevolutions[1] + 1))
            )
        };
    }

    // Mutation
    function mutate(individual, rate = 0.1) {
        const mutated = JSON.parse(JSON.stringify(individual));
        const numLegs = config.transferSequence.length - 1;

        if (Math.random() < rate) {
            mutated.departureEpoch += (Math.random() - 0.5) * 20 * JULIAN_DAY;
            mutated.departureEpoch = Math.max(bounds.departureEpoch[0],
                Math.min(bounds.departureEpoch[1], mutated.departureEpoch));
        }

        for (let i = 0; i < numLegs; i++) {
            if (Math.random() < rate) {
                mutated.legTOFs[i] += (Math.random() - 0.5) * 100 * JULIAN_DAY;
                mutated.legTOFs[i] = Math.max(bounds.legTOF[0],
                    Math.min(bounds.legTOF[1], mutated.legTOFs[i]));
            }
            if (Math.random() < rate) {
                mutated.numRevs[i] = Math.floor(Math.random() * (bounds.numRevolutions[1] + 1));
            }
        }

        return mutated;
    }

    // Crossover
    function crossover(parent1, parent2) {
        const child = JSON.parse(JSON.stringify(parent1));
        const numLegs = config.transferSequence.length - 1;

        if (Math.random() < 0.5) child.departureEpoch = parent2.departureEpoch;
        for (let i = 0; i < numLegs; i++) {
            if (Math.random() < 0.5) child.legTOFs[i] = parent2.legTOFs[i];
            if (Math.random() < 0.5) child.numRevs[i] = parent2.numRevs[i];
        }

        return child;
    }

    // Run genetic algorithm
    let population = Array(config.populationSize).fill(null).map(() => randomIndividual());
    const fitnessHistory = [];
    let bestIndividual = null;
    let bestFitness = Infinity;

    for (let gen = 0; gen < config.numGenerations; gen++) {
        const evaluated = population.map(ind => ({
            individual: ind,
            result: fitness(ind)
        }));

        evaluated.sort((a, b) => a.result.totalDeltaV - b.result.totalDeltaV);

        if (evaluated[0].result.totalDeltaV < bestFitness) {
            bestFitness = evaluated[0].result.totalDeltaV;
            bestIndividual = evaluated[0].individual;
        }

        fitnessHistory.push({
            generation: gen,
            best: evaluated[0].result.totalDeltaV,
            mean: evaluated.reduce((s, e) => s + e.result.totalDeltaV, 0) / evaluated.length
        });

        const elite = evaluated.slice(0, Math.floor(config.populationSize * 0.1));
        const newPop = elite.map(e => e.individual);

        while (newPop.length < config.populationSize) {
            const parent1 = elite[Math.floor(Math.random() * elite.length)].individual;
            const parent2 = elite[Math.floor(Math.random() * elite.length)].individual;
            let child = crossover(parent1, parent2);
            child = mutate(child, 0.2);
            newPop.push(child);
        }

        population = newPop;
    }

    const bestResult = fitness(bestIndividual);
    const trajectory = generateTrajectorySpice(bestIndividual, bestResult, AU);

    return {
        bestDeltaV: bestResult.totalDeltaV,
        departureEpoch: bestIndividual.departureEpoch,
        legTOFs: bestIndividual.legTOFs,
        totalTOF: bestIndividual.legTOFs.reduce((a, b) => a + b, 0),
        numRevolutions: bestIndividual.numRevs,
        nodePositions: bestResult.nodePositions,
        fitnessHistory,
        trajectory,
        planetOrbits: Object.fromEntries(
            config.transferSequence.map(name => [name, (PLANETARY_SMA[name] || AU) / AU])
        ),
        source: 'SPICE'
    };
}

function generateTrajectorySpice(individual, result, AU) {
    const trajectory = [];
    const numPoints = 200;

    let epoch = individual.departureEpoch;

    for (let legIdx = 0; legIdx < individual.legTOFs.length; legIdx++) {
        const startNode = result.nodePositions[legIdx];
        const endNode = result.nodePositions[legIdx + 1];
        const tof = individual.legTOFs[legIdx];

        for (let i = 0; i <= numPoints / individual.legTOFs.length; i++) {
            const t = i / (numPoints / individual.legTOFs.length);
            const currentEpoch = epoch + t * tof;

            const r1 = Math.sqrt(startNode.x**2 + startNode.y**2);
            const r2 = Math.sqrt(endNode.x**2 + endNode.y**2);
            const theta1 = Math.atan2(startNode.y, startNode.x);
            const theta2 = Math.atan2(endNode.y, endNode.x);

            let dTheta = theta2 - theta1;
            if (dTheta > Math.PI) dTheta -= 2 * Math.PI;
            if (dTheta < -Math.PI) dTheta += 2 * Math.PI;
            dTheta += individual.numRevs[legIdx] * 2 * Math.PI;

            const theta = theta1 + dTheta * t;
            const r = r1 + (r2 - r1) * t + 0.1 * AU * Math.sin(Math.PI * t);

            trajectory.push({
                x: r * Math.cos(theta) / AU,
                y: r * Math.sin(theta) / AU,
                leg: legIdx,
                t: (currentEpoch - individual.departureEpoch) / 86400
            });
        }

        epoch += tof;
    }

    return trajectory;
}

function optimizeHodographicMGA(config) {
    const AU = 1.496e11;
    const JULIAN_DAY = 86400;
    const GM_SUN = 1.327e20;

    // Planetary parameters (simplified circular coplanar orbits)
    const planets = {
        Earth: { a: 1.0 * AU, T: 365.25 * JULIAN_DAY, R: 6.371e6, GM: 3.986e14 },
        Mars: { a: 1.524 * AU, T: 687 * JULIAN_DAY, R: 3.390e6, GM: 4.283e13 },
        Jupiter: { a: 5.203 * AU, T: 4333 * JULIAN_DAY, R: 6.995e7, GM: 1.267e17 },
        Venus: { a: 0.723 * AU, T: 224.7 * JULIAN_DAY, R: 6.052e6, GM: 3.249e14 }
    };

    // Get planet position at epoch
    function getPlanetPosition(name, epoch, phase0 = 0) {
        const planet = planets[name];
        const n = 2 * Math.PI / planet.T;
        const angle = n * epoch + phase0;
        const x = planet.a * Math.cos(angle);
        const y = planet.a * Math.sin(angle);
        const vx = -planet.a * n * Math.sin(angle);
        const vy = planet.a * n * Math.cos(angle);
        return { x, y, vx, vy };
    }

    // Hodographic shaping approximation for low-thrust leg
    function computeLegDeltaV(r1, r2, tof, numRevs = 0) {
        // Simplified hodographic shaping cost estimation
        const a_transfer = (r1 + r2) / 2;
        const e_transfer = Math.abs(r2 - r1) / (r2 + r1);

        // Required velocity change (simplified)
        const v1 = Math.sqrt(GM_SUN / r1);
        const v2 = Math.sqrt(GM_SUN / r2);
        const v_transfer = Math.sqrt(GM_SUN * (2 / r1 - 1 / a_transfer));

        // Time factor (longer TOF = lower thrust = lower instantaneous DV but may increase total)
        const minTOF = Math.PI * Math.sqrt(a_transfer ** 3 / GM_SUN);
        const tofRatio = tof / minTOF;

        // Low-thrust approximation
        let deltaV = Math.abs(v_transfer - v1) + Math.abs(v2 - v_transfer);

        // Adjust for revolutions
        deltaV *= (1 + 0.1 * numRevs);

        // Adjust for TOF
        if (tofRatio < 1) {
            deltaV *= 2 * (1 - tofRatio) + 1;  // Penalize too short TOF
        } else if (tofRatio > 3) {
            deltaV *= 1 + 0.1 * (tofRatio - 3);  // Slight penalty for very long TOF
        }

        return deltaV;
    }

    // Compute gravity assist delta-V capability
    function swingbyDeltaV(planet, vInf, periapsisAlt) {
        const rp = planets[planet].R + periapsisAlt;
        const GM = planets[planet].GM;
        const delta = 2 * Math.asin(1 / (1 + rp * vInf ** 2 / GM));
        return 2 * vInf * Math.sin(delta / 2);
    }

    // Optimization bounds
    const bounds = {
        departureEpoch: [9000 * JULIAN_DAY, 9200 * JULIAN_DAY],
        legTOF: [200 * JULIAN_DAY, 1200 * JULIAN_DAY],
        swingbyAltitude: [200e3, 2e8],  // 200 km to 200,000 km
        numRevolutions: [0, 2]
    };

    // Initial planet phases
    const phases = {
        Earth: 0,
        Mars: Math.PI / 4,
        Jupiter: Math.PI / 2
    };

    // Fitness function
    function fitness(individual) {
        const departureEpoch = individual.departureEpoch;
        const legTOFs = individual.legTOFs;
        const swingbyAlts = individual.swingbyAlts;
        const numRevs = individual.numRevs;

        let totalDeltaV = 0;
        let epoch = departureEpoch;
        const nodePositions = [];
        const nodeEpochs = [epoch];

        // Compute each leg
        for (let i = 0; i < config.transferSequence.length - 1; i++) {
            const fromPlanet = config.transferSequence[i];
            const toPlanet = config.transferSequence[i + 1];

            const pos1 = getPlanetPosition(fromPlanet, epoch, phases[fromPlanet] || 0);
            epoch += legTOFs[i];
            const pos2 = getPlanetPosition(toPlanet, epoch, phases[toPlanet] || 0);

            const r1 = Math.sqrt(pos1.x ** 2 + pos1.y ** 2);
            const r2 = Math.sqrt(pos2.x ** 2 + pos2.y ** 2);

            // Leg delta-V
            const legDV = computeLegDeltaV(r1, r2, legTOFs[i], numRevs[i]);
            totalDeltaV += legDV;

            nodePositions.push({ planet: fromPlanet, ...pos1, epoch: nodeEpochs[i] });
            nodeEpochs.push(epoch);
        }

        // Add final node
        const lastPlanet = config.transferSequence[config.transferSequence.length - 1];
        const lastPos = getPlanetPosition(lastPlanet, epoch, phases[lastPlanet] || 0);
        nodePositions.push({ planet: lastPlanet, ...lastPos, epoch });

        return { totalDeltaV, nodePositions, nodeEpochs };
    }

    // Generate random individual
    function randomIndividual() {
        const numLegs = config.transferSequence.length - 1;
        const numSwingbys = numLegs - 1;

        return {
            departureEpoch: bounds.departureEpoch[0] + Math.random() * (bounds.departureEpoch[1] - bounds.departureEpoch[0]),
            legTOFs: Array(numLegs).fill(0).map(() =>
                bounds.legTOF[0] + Math.random() * (bounds.legTOF[1] - bounds.legTOF[0])
            ),
            swingbyAlts: Array(numSwingbys).fill(0).map(() =>
                Math.exp(Math.log(bounds.swingbyAltitude[0]) +
                    Math.random() * (Math.log(bounds.swingbyAltitude[1]) - Math.log(bounds.swingbyAltitude[0])))
            ),
            numRevs: Array(numLegs).fill(0).map(() =>
                Math.floor(bounds.numRevolutions[0] + Math.random() * (bounds.numRevolutions[1] - bounds.numRevolutions[0] + 1))
            )
        };
    }

    // Mutation
    function mutate(individual, rate = 0.1) {
        const mutated = JSON.parse(JSON.stringify(individual));
        const numLegs = config.transferSequence.length - 1;

        if (Math.random() < rate) {
            mutated.departureEpoch += (Math.random() - 0.5) * 20 * JULIAN_DAY;
            mutated.departureEpoch = Math.max(bounds.departureEpoch[0],
                Math.min(bounds.departureEpoch[1], mutated.departureEpoch));
        }

        for (let i = 0; i < numLegs; i++) {
            if (Math.random() < rate) {
                mutated.legTOFs[i] += (Math.random() - 0.5) * 100 * JULIAN_DAY;
                mutated.legTOFs[i] = Math.max(bounds.legTOF[0],
                    Math.min(bounds.legTOF[1], mutated.legTOFs[i]));
            }
            if (Math.random() < rate) {
                mutated.numRevs[i] = Math.floor(Math.random() * (bounds.numRevolutions[1] + 1));
            }
        }

        return mutated;
    }

    // Crossover
    function crossover(parent1, parent2) {
        const child = JSON.parse(JSON.stringify(parent1));
        const numLegs = config.transferSequence.length - 1;

        if (Math.random() < 0.5) child.departureEpoch = parent2.departureEpoch;

        for (let i = 0; i < numLegs; i++) {
            if (Math.random() < 0.5) child.legTOFs[i] = parent2.legTOFs[i];
            if (Math.random() < 0.5) child.numRevs[i] = parent2.numRevs[i];
        }

        return child;
    }

    // Run genetic algorithm
    let population = Array(config.populationSize).fill(null).map(() => randomIndividual());
    const fitnessHistory = [];
    let bestIndividual = null;
    let bestFitness = Infinity;

    for (let gen = 0; gen < config.numGenerations; gen++) {
        // Evaluate fitness
        const evaluated = population.map(ind => ({
            individual: ind,
            result: fitness(ind)
        }));

        // Sort by fitness
        evaluated.sort((a, b) => a.result.totalDeltaV - b.result.totalDeltaV);

        // Track best
        if (evaluated[0].result.totalDeltaV < bestFitness) {
            bestFitness = evaluated[0].result.totalDeltaV;
            bestIndividual = evaluated[0].individual;
        }

        fitnessHistory.push({
            generation: gen,
            best: evaluated[0].result.totalDeltaV,
            mean: evaluated.reduce((s, e) => s + e.result.totalDeltaV, 0) / evaluated.length
        });

        // Selection and reproduction
        const elite = evaluated.slice(0, Math.floor(config.populationSize * 0.1));
        const newPop = elite.map(e => e.individual);

        while (newPop.length < config.populationSize) {
            const parent1 = elite[Math.floor(Math.random() * elite.length)].individual;
            const parent2 = elite[Math.floor(Math.random() * elite.length)].individual;
            let child = crossover(parent1, parent2);
            child = mutate(child, 0.2);
            newPop.push(child);
        }

        population = newPop;
    }

    // Final evaluation of best
    const bestResult = fitness(bestIndividual);

    // Generate trajectory for visualization
    const trajectory = generateTrajectory(bestIndividual, bestResult, planets, phases, GM_SUN, AU);

    return {
        bestDeltaV: bestResult.totalDeltaV,
        departureEpoch: bestIndividual.departureEpoch,
        legTOFs: bestIndividual.legTOFs,
        totalTOF: bestIndividual.legTOFs.reduce((a, b) => a + b, 0),
        numRevolutions: bestIndividual.numRevs,
        nodePositions: bestResult.nodePositions,
        fitnessHistory,
        trajectory,
        planetOrbits: Object.fromEntries(
            Object.entries(planets).map(([name, p]) => [name, p.a / AU])
        )
    };
}

function generateTrajectory(individual, result, planets, phases, GM_SUN, AU) {
    const trajectory = [];
    const numPoints = 200;

    let epoch = individual.departureEpoch;

    for (let legIdx = 0; legIdx < individual.legTOFs.length; legIdx++) {
        const startNode = result.nodePositions[legIdx];
        const endNode = result.nodePositions[legIdx + 1];
        const tof = individual.legTOFs[legIdx];

        for (let i = 0; i <= numPoints / individual.legTOFs.length; i++) {
            const t = i / (numPoints / individual.legTOFs.length);
            const currentEpoch = epoch + t * tof;

            // Simple interpolation with curve
            const r1 = Math.sqrt(startNode.x ** 2 + startNode.y ** 2);
            const r2 = Math.sqrt(endNode.x ** 2 + endNode.y ** 2);
            const theta1 = Math.atan2(startNode.y, startNode.x);
            const theta2 = Math.atan2(endNode.y, endNode.x);

            // Ensure proper angle interpolation
            let dTheta = theta2 - theta1;
            if (dTheta > Math.PI) dTheta -= 2 * Math.PI;
            if (dTheta < -Math.PI) dTheta += 2 * Math.PI;

            // Add revolutions
            dTheta += individual.numRevs[legIdx] * 2 * Math.PI;

            const theta = theta1 + dTheta * t;
            const r = r1 + (r2 - r1) * t + 0.1 * AU * Math.sin(Math.PI * t);  // Slight curve

            trajectory.push({
                x: r * Math.cos(theta) / AU,
                y: r * Math.sin(theta) / AU,
                leg: legIdx,
                t: (currentEpoch - individual.departureEpoch) / 86400  // days
            });
        }

        epoch += tof;
    }

    return trajectory;
}

function renderHodographicMGAFigures(container, result, config, useSpice = false) {
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
    const chartHeight = Math.max(220, (containerHeight - 80) / 2);

    // Figure 1: Trajectory
    const title1 = useSpice
        ? 'Low-Thrust MGA Trajectory (SPICE Ephemeris)'
        : 'Low-Thrust MGA Trajectory';
    createFigure(wrapper, title1, chartWidth, chartHeight, (svg, w, h) => {
        renderTrajectory(svg, result, config, w, h);
    });

    // Figure 2: Fitness evolution
    createFigure(wrapper, 'Optimization Convergence', chartWidth, chartHeight, (svg, w, h) => {
        renderFitnessEvolution(svg, result, w, h);
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

function renderTrajectory(svg, result, config, width, height) {
    const margin = { top: 20, right: 120, bottom: 35, left: 55 };

    const maxR = Math.max(
        d3.max(result.trajectory, d => Math.sqrt(d.x ** 2 + d.y ** 2)),
        ...Object.values(result.planetOrbits)
    ) * 1.1;

    const x = d3.scaleLinear()
        .domain([-maxR, maxR])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([-maxR, maxR])
        .range([height - margin.bottom, margin.top]);

    // Planet orbit colors
    const planetColors = {
        Earth: '#4488ff',
        Mars: '#ff6644',
        Jupiter: '#ffaa44',
        Venus: '#ffdd88'
    };

    // Draw planet orbits
    Object.entries(result.planetOrbits).forEach(([name, radius]) => {
        if (config.transferSequence.includes(name)) {
            const points = [];
            for (let i = 0; i <= 64; i++) {
                const angle = (i / 64) * 2 * Math.PI;
                points.push({ x: radius * Math.cos(angle), y: radius * Math.sin(angle) });
            }

            svg.append('path')
                .datum(points)
                .attr('fill', 'none')
                .attr('stroke', planetColors[name])
                .attr('stroke-width', 0.5)
                .attr('stroke-dasharray', '3,3')
                .attr('opacity', 0.5)
                .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));
        }
    });

    // Draw Sun
    svg.append('circle')
        .attr('cx', x(0))
        .attr('cy', y(0))
        .attr('r', 6)
        .attr('fill', '#ffcc00');

    // Draw trajectory colored by leg
    const legColors = ['var(--cyan)', '#00ff88', '#ff88ff'];
    for (let leg = 0; leg < config.transferSequence.length - 1; leg++) {
        const legData = result.trajectory.filter(d => d.leg === leg);
        svg.append('path')
            .datum(legData)
            .attr('fill', 'none')
            .attr('stroke', legColors[leg % legColors.length])
            .attr('stroke-width', 2)
            .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));
    }

    // Draw node markers
    result.nodePositions.forEach((node, i) => {
        const nx = node.x / 1.496e11;
        const ny = node.y / 1.496e11;

        svg.append('circle')
            .attr('cx', x(nx))
            .attr('cy', y(ny))
            .attr('r', 5)
            .attr('fill', planetColors[node.planet] || 'white')
            .attr('stroke', 'white')
            .attr('stroke-width', 1);
    });

    // Legend
    let legendY = 20;
    config.transferSequence.forEach((planet, i) => {
        if (i < config.transferSequence.length - 1) {
            svg.append('line')
                .attr('x1', width - 110)
                .attr('x2', width - 90)
                .attr('y1', legendY)
                .attr('y2', legendY)
                .attr('stroke', legColors[i % legColors.length])
                .attr('stroke-width', 2);
            svg.append('text')
                .attr('x', width - 85)
                .attr('y', legendY + 4)
                .attr('fill', 'var(--text-secondary)')
                .attr('font-size', '9px')
                .text(`${planet}→${config.transferSequence[i + 1]}`);
            legendY += 15;
        }
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
        .attr('x', width / 2 - 30)
        .attr('y', height - 3)
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

function renderFitnessEvolution(svg, result, width, height) {
    const margin = { top: 15, right: 80, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain([0, result.fitnessHistory.length - 1])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([
            d3.min(result.fitnessHistory, d => d.best) * 0.9,
            d3.max(result.fitnessHistory, d => d.mean) * 1.1
        ])
        .range([height - margin.bottom, margin.top]);

    // Best fitness line
    svg.append('path')
        .datum(result.fitnessHistory)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.generation)).y(d => y(d.best)));

    // Mean fitness line
    svg.append('path')
        .datum(result.fitnessHistory)
        .attr('fill', 'none')
        .attr('stroke', '#888888')
        .attr('stroke-width', 1)
        .attr('stroke-dasharray', '4,2')
        .attr('d', d3.line().x(d => x(d.generation)).y(d => y(d.mean)));

    // Legend
    svg.append('line').attr('x1', width - 70).attr('x2', width - 50).attr('y1', 20).attr('y2', 20)
        .attr('stroke', 'var(--cyan)').attr('stroke-width', 2);
    svg.append('text').attr('x', width - 45).attr('y', 23)
        .attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('Best');

    svg.append('line').attr('x1', width - 70).attr('x2', width - 50).attr('y1', 35).attr('y2', 35)
        .attr('stroke', '#888888').attr('stroke-width', 1).attr('stroke-dasharray', '4,2');
    svg.append('text').attr('x', width - 45).attr('y', 38)
        .attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('Mean');

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(6).tickFormat(d3.format('d')))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(5).tickFormat(d => (d / 1000).toFixed(0)))
        .attr('color', 'var(--text-dim)');

    svg.append('text')
        .attr('x', width / 2 - 10)
        .attr('y', height - 3)
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
