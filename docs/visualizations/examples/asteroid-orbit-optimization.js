/**
 * Asteroid Orbit Optimization Example
 * Port of: tudatpy/examples/pygmo/asteroid_orbit_optimization
 *
 * Demonstrates multi-objective optimization for asteroid mission design.
 * Shows Pareto front evolution and design space exploration.
 */

export function showAsteroidOrbitOptimizationExample(chartContainer, log, params = {}) {
    const config = {
        numGenerations: params.numGenerations ?? 30,
        populationSize: params.populationSize ?? 50,
        targetAsteroid: params.targetAsteroid ?? 'Eros'
    };

    log('Running Asteroid Orbit Optimization...', 'info');
    log(`Target: ${config.targetAsteroid}`, 'info');
    log(`Population: ${config.populationSize}, Generations: ${config.numGenerations}`, 'info');

    const startTime = performance.now();
    const result = runAsteroidOptimization(config);
    const elapsed = performance.now() - startTime;
    log(`Optimization completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Pareto front size: ${result.paretoFront.length} solutions`, 'info');
    log(`Min ΔV: ${result.minDeltaV.toFixed(0)} m/s`, 'info');
    log(`Min TOF: ${result.minTOF.toFixed(0)} days`, 'info');

    renderAsteroidOptFigures(chartContainer, result, config);

    return {
        name: 'Asteroid Orbit Optimization',
        description: 'Multi-objective mission design',
        ...result,
        config
    };
}

function runAsteroidOptimization(config) {
    // Design space bounds
    const bounds = {
        departureDate: { min: 0, max: 365 * 2 },      // days from epoch
        tof: { min: 200, max: 600 },                   // time of flight days
        launchVInf: { min: 0, max: 5 }                 // km/s
    };

    // Generate initial population
    let population = [];
    for (let i = 0; i < config.populationSize; i++) {
        const individual = {
            departureDate: bounds.departureDate.min + Math.random() * (bounds.departureDate.max - bounds.departureDate.min),
            tof: bounds.tof.min + Math.random() * (bounds.tof.max - bounds.tof.min),
            launchVInf: bounds.launchVInf.min + Math.random() * (bounds.launchVInf.max - bounds.launchVInf.min)
        };
        const objectives = evaluateAsteroidMission(individual);
        population.push({ ...individual, ...objectives });
    }

    const paretoHistory = [];
    const convergenceHistory = [];

    // NSGA-II style evolution
    for (let gen = 0; gen < config.numGenerations; gen++) {
        // Non-dominated sorting
        const fronts = nonDominatedSort(population);
        const paretoFront = fronts[0];

        paretoHistory.push({
            generation: gen,
            frontSize: paretoFront.length,
            front: paretoFront.map(p => ({ deltaV: p.deltaV, tof: p.tof }))
        });

        convergenceHistory.push({
            generation: gen,
            avgDeltaV: population.reduce((s, p) => s + p.deltaV, 0) / population.length,
            minDeltaV: Math.min(...population.map(p => p.deltaV)),
            avgTOF: population.reduce((s, p) => s + p.tof, 0) / population.length
        });

        // Generate offspring
        const offspring = [];
        while (offspring.length < config.populationSize) {
            // Tournament selection
            const parent1 = tournamentSelect(population, fronts);
            const parent2 = tournamentSelect(population, fronts);

            // Crossover
            const child = {
                departureDate: Math.random() < 0.5 ? parent1.departureDate : parent2.departureDate,
                tof: Math.random() < 0.5 ? parent1.tof : parent2.tof,
                launchVInf: Math.random() < 0.5 ? parent1.launchVInf : parent2.launchVInf
            };

            // Mutation
            if (Math.random() < 0.2) {
                child.departureDate += (Math.random() - 0.5) * 60;
                child.departureDate = Math.max(bounds.departureDate.min, Math.min(bounds.departureDate.max, child.departureDate));
            }
            if (Math.random() < 0.2) {
                child.tof += (Math.random() - 0.5) * 100;
                child.tof = Math.max(bounds.tof.min, Math.min(bounds.tof.max, child.tof));
            }
            if (Math.random() < 0.2) {
                child.launchVInf += (Math.random() - 0.5) * 1;
                child.launchVInf = Math.max(bounds.launchVInf.min, Math.min(bounds.launchVInf.max, child.launchVInf));
            }

            const objectives = evaluateAsteroidMission(child);
            offspring.push({ ...child, ...objectives });
        }

        // Combine and select
        const combined = [...population, ...offspring];
        const combinedFronts = nonDominatedSort(combined);

        // Select best individuals
        population = [];
        for (const front of combinedFronts) {
            if (population.length + front.length <= config.populationSize) {
                population.push(...front);
            } else {
                // Crowding distance selection
                const remaining = config.populationSize - population.length;
                const sorted = front.sort((a, b) => a.deltaV - b.deltaV);
                population.push(...sorted.slice(0, remaining));
                break;
            }
        }
    }

    // Final Pareto front
    const finalFronts = nonDominatedSort(population);
    const paretoFront = finalFronts[0];

    return {
        paretoFront,
        paretoHistory,
        convergenceHistory,
        population,
        minDeltaV: Math.min(...paretoFront.map(p => p.deltaV)),
        minTOF: Math.min(...paretoFront.map(p => p.tof)),
        bounds
    };
}

function evaluateAsteroidMission(individual) {
    // Simplified mission model
    // In reality this would use Lambert solver and actual asteroid ephemeris

    const { departureDate, tof, launchVInf } = individual;

    // Departure C3 (characteristic energy)
    const c3 = launchVInf * launchVInf;

    // Launch delta-V (Earth escape)
    const vEarth = 29.78;  // km/s
    const deltaVLaunch = Math.sqrt(c3 + 2 * 398600 / 6678) - Math.sqrt(398600 / 6678);

    // Transfer geometry penalty (simplified)
    const phaseAnglePenalty = 0.5 * Math.abs(Math.sin(departureDate * 2 * Math.PI / 365));

    // Arrival delta-V (simplified model based on TOF)
    const optimalTOF = 350;  // Approximate optimal for Eros-like target
    const tofPenalty = Math.abs(tof - optimalTOF) / 200;
    const deltaVArrival = 2.0 + tofPenalty * 3 + phaseAnglePenalty * 1.5;

    const totalDeltaV = (deltaVLaunch + deltaVArrival) * 1000;  // Convert to m/s

    return {
        deltaV: totalDeltaV,
        c3,
        deltaVLaunch: deltaVLaunch * 1000,
        deltaVArrival: deltaVArrival * 1000
    };
}

function nonDominatedSort(population) {
    const fronts = [[]];
    const dominationCount = new Map();
    const dominatedBy = new Map();

    for (const p of population) {
        dominationCount.set(p, 0);
        dominatedBy.set(p, []);
    }

    for (const p of population) {
        for (const q of population) {
            if (p === q) continue;

            if (dominates(p, q)) {
                dominatedBy.get(p).push(q);
            } else if (dominates(q, p)) {
                dominationCount.set(p, dominationCount.get(p) + 1);
            }
        }

        if (dominationCount.get(p) === 0) {
            fronts[0].push(p);
        }
    }

    let i = 0;
    while (fronts[i].length > 0) {
        const nextFront = [];
        for (const p of fronts[i]) {
            for (const q of dominatedBy.get(p)) {
                dominationCount.set(q, dominationCount.get(q) - 1);
                if (dominationCount.get(q) === 0) {
                    nextFront.push(q);
                }
            }
        }
        i++;
        fronts.push(nextFront);
    }

    return fronts.filter(f => f.length > 0);
}

function dominates(p, q) {
    // p dominates q if p is at least as good in all objectives and strictly better in at least one
    const pBetterDV = p.deltaV <= q.deltaV;
    const pBetterTOF = p.tof <= q.tof;
    const pStrictlyBetter = p.deltaV < q.deltaV || p.tof < q.tof;
    return pBetterDV && pBetterTOF && pStrictlyBetter;
}

function tournamentSelect(population, fronts) {
    const a = population[Math.floor(Math.random() * population.length)];
    const b = population[Math.floor(Math.random() * population.length)];

    // Find front indices
    let aFront = 0, bFront = 0;
    for (let i = 0; i < fronts.length; i++) {
        if (fronts[i].includes(a)) aFront = i;
        if (fronts[i].includes(b)) bFront = i;
    }

    return aFront < bFront ? a : (bFront < aFront ? b : (Math.random() < 0.5 ? a : b));
}

function renderAsteroidOptFigures(container, result, config) {
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
    const chartHeight = Math.max(220, (containerHeight - 100) / 2);

    // Figure 1: Pareto front
    createFigure(wrapper, 'Pareto Front (ΔV vs Time of Flight)', chartWidth, chartHeight, (svg, w, h) => {
        renderParetoFront(svg, result, w, h);
    });

    // Figure 2: Convergence
    createFigure(wrapper, 'Optimization Convergence', chartWidth, Math.max(180, chartHeight * 0.7), (svg, w, h) => {
        renderOptConvergence(svg, result.convergenceHistory, w, h);
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

function renderParetoFront(svg, result, width, height) {
    const margin = { top: 20, right: 20, bottom: 40, left: 60 };

    const x = d3.scaleLinear()
        .domain([result.bounds.tof.min, result.bounds.tof.max])
        .range([margin.left, width - margin.right]);

    const allDeltaV = result.population.map(p => p.deltaV);
    const y = d3.scaleLinear()
        .domain([Math.min(...allDeltaV) * 0.9, Math.max(...allDeltaV) * 1.1])
        .range([height - margin.bottom, margin.top]);

    // Draw all population points (faded)
    svg.selectAll('circle.pop')
        .data(result.population)
        .enter()
        .append('circle')
        .attr('class', 'pop')
        .attr('cx', d => x(d.tof))
        .attr('cy', d => y(d.deltaV))
        .attr('r', 3)
        .attr('fill', 'var(--text-dim)')
        .attr('opacity', 0.3);

    // Draw Pareto front points
    svg.selectAll('circle.pareto')
        .data(result.paretoFront)
        .enter()
        .append('circle')
        .attr('class', 'pareto')
        .attr('cx', d => x(d.tof))
        .attr('cy', d => y(d.deltaV))
        .attr('r', 5)
        .attr('fill', 'var(--cyan)')
        .attr('stroke', 'white')
        .attr('stroke-width', 1);

    // Connect Pareto front with line
    const sortedPareto = [...result.paretoFront].sort((a, b) => a.tof - b.tof);
    svg.append('path')
        .datum(sortedPareto)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 1.5)
        .attr('stroke-dasharray', '4,2')
        .attr('d', d3.line().x(d => x(d.tof)).y(d => y(d.deltaV)));

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
        .text('Time of Flight [days]');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('Total ΔV [m/s]');
}

function renderOptConvergence(svg, history, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain([0, history.length - 1])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([0, Math.max(...history.map(d => d.avgDeltaV)) * 1.1])
        .range([height - margin.bottom, margin.top]);

    // Average delta-V
    svg.append('path')
        .datum(history)
        .attr('fill', 'none')
        .attr('stroke', 'var(--orange)')
        .attr('stroke-width', 1.5)
        .attr('d', d3.line().x(d => x(d.generation)).y(d => y(d.avgDeltaV)));

    // Min delta-V
    svg.append('path')
        .datum(history)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 2)
        .attr('d', d3.line().x(d => x(d.generation)).y(d => y(d.minDeltaV)));

    // Legend
    const legend = svg.append('g')
        .attr('transform', `translate(${width - margin.right - 80}, ${margin.top + 10})`);

    legend.append('line').attr('x1', 0).attr('y1', 0).attr('x2', 20).attr('y2', 0).attr('stroke', 'var(--cyan)').attr('stroke-width', 2);
    legend.append('text').attr('x', 25).attr('y', 4).attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('Min ΔV');

    legend.append('line').attr('x1', 0).attr('y1', 16).attr('x2', 20).attr('y2', 16).attr('stroke', 'var(--orange)').attr('stroke-width', 1.5);
    legend.append('text').attr('x', 25).attr('y', 20).attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('Avg ΔV');

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
        .text('ΔV [m/s]');
}
