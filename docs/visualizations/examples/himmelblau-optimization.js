/**
 * Himmelblau Optimization Example
 * Port of: tudatpy/examples/pygmo/himmelblau_minimization.py
 *
 * Demonstrates optimization algorithm convergence on Himmelblau's function.
 * Shows the 2D function landscape and optimization path.
 */

export function showHimmelblauOptimizationExample(chartContainer, log, params = {}) {
    const config = {
        numGenerations: params.numGenerations ?? 50,
        populationSize: params.populationSize ?? 20,
        algorithm: params.algorithm ?? 'de'  // differential evolution
    };

    log('Running Himmelblau Optimization Example...', 'info');
    log(`Algorithm: Differential Evolution`, 'info');
    log(`Population size: ${config.populationSize}`, 'info');
    log(`Generations: ${config.numGenerations}`, 'info');

    const startTime = performance.now();
    const result = runOptimization(config);
    const elapsed = performance.now() - startTime;
    log(`Optimization completed in ${elapsed.toFixed(1)} ms`, 'success');

    log(`Best solution: (${result.bestX.toFixed(4)}, ${result.bestY.toFixed(4)})`, 'info');
    log(`Best fitness: ${result.bestFitness.toFixed(6)}`, 'info');
    log(`Converged to minimum #${result.minimumIndex + 1}`, 'info');

    renderHimmelblauFigures(chartContainer, result, config);

    return {
        name: 'Himmelblau Optimization',
        description: 'Multi-modal function optimization',
        ...result,
        config
    };
}

// Himmelblau's function: f(x,y) = (x² + y - 11)² + (x + y² - 7)²
function himmelblau(x, y) {
    return Math.pow(x*x + y - 11, 2) + Math.pow(x + y*y - 7, 2);
}

// Known minima of Himmelblau's function
const MINIMA = [
    { x: 3.0, y: 2.0 },
    { x: -2.805118, y: 3.131312 },
    { x: -3.779310, y: -3.283186 },
    { x: 3.584428, y: -1.848126 }
];

function runOptimization(config) {
    const bounds = { xMin: -5, xMax: 5, yMin: -5, yMax: 5 };

    // Initialize population randomly
    let population = [];
    for (let i = 0; i < config.populationSize; i++) {
        const x = bounds.xMin + Math.random() * (bounds.xMax - bounds.xMin);
        const y = bounds.yMin + Math.random() * (bounds.yMax - bounds.yMin);
        population.push({ x, y, fitness: himmelblau(x, y) });
    }

    // Track convergence history
    const convergenceHistory = [];
    const bestPath = [];

    // Differential Evolution parameters
    const F = 0.8;   // Mutation factor
    const CR = 0.9;  // Crossover probability

    for (let gen = 0; gen < config.numGenerations; gen++) {
        // Find best in current population
        population.sort((a, b) => a.fitness - b.fitness);
        const best = population[0];

        convergenceHistory.push({
            generation: gen,
            bestFitness: best.fitness,
            avgFitness: population.reduce((s, p) => s + p.fitness, 0) / population.length
        });

        bestPath.push({ x: best.x, y: best.y, gen });

        // Evolution step
        const newPopulation = [];
        for (let i = 0; i < config.populationSize; i++) {
            // Select three random individuals (different from i)
            const candidates = population.filter((_, idx) => idx !== i);
            const [a, b, c] = shuffleArray(candidates).slice(0, 3);

            // Mutation
            let mutantX = a.x + F * (b.x - c.x);
            let mutantY = a.y + F * (b.y - c.y);

            // Bound constraints
            mutantX = Math.max(bounds.xMin, Math.min(bounds.xMax, mutantX));
            mutantY = Math.max(bounds.yMin, Math.min(bounds.yMax, mutantY));

            // Crossover
            const trial = { x: population[i].x, y: population[i].y };
            if (Math.random() < CR) trial.x = mutantX;
            if (Math.random() < CR) trial.y = mutantY;
            trial.fitness = himmelblau(trial.x, trial.y);

            // Selection
            if (trial.fitness <= population[i].fitness) {
                newPopulation.push(trial);
            } else {
                newPopulation.push({ ...population[i] });
            }
        }

        population = newPopulation;
    }

    // Final best
    population.sort((a, b) => a.fitness - b.fitness);
    const best = population[0];

    // Find which minimum we converged to
    let minimumIndex = 0;
    let minDist = Infinity;
    MINIMA.forEach((m, i) => {
        const dist = Math.sqrt((best.x - m.x)**2 + (best.y - m.y)**2);
        if (dist < minDist) {
            minDist = dist;
            minimumIndex = i;
        }
    });

    // Generate function landscape for visualization
    const landscape = generateLandscape(bounds, 50);

    return {
        bestX: best.x,
        bestY: best.y,
        bestFitness: best.fitness,
        minimumIndex,
        convergenceHistory,
        bestPath,
        finalPopulation: population,
        landscape,
        bounds
    };
}

function shuffleArray(array) {
    const arr = [...array];
    for (let i = arr.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [arr[i], arr[j]] = [arr[j], arr[i]];
    }
    return arr;
}

function generateLandscape(bounds, resolution) {
    const data = [];
    const dx = (bounds.xMax - bounds.xMin) / resolution;
    const dy = (bounds.yMax - bounds.yMin) / resolution;

    for (let i = 0; i <= resolution; i++) {
        for (let j = 0; j <= resolution; j++) {
            const x = bounds.xMin + i * dx;
            const y = bounds.yMin + j * dy;
            const z = himmelblau(x, y);
            data.push({ x, y, z, i, j });
        }
    }
    return { data, resolution };
}

function renderHimmelblauFigures(container, result, config) {
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
    const chartHeight = Math.max(250, (containerHeight - 100) / 2);

    // Figure 1: Function landscape with optimization path
    createFigure(wrapper, 'Himmelblau Function Landscape', chartWidth, chartHeight, (svg, w, h) => {
        renderLandscape(svg, result, w, h);
    });

    // Figure 2: Convergence history
    createFigure(wrapper, 'Optimization Convergence', chartWidth, Math.max(180, chartHeight * 0.6), (svg, w, h) => {
        renderConvergence(svg, result.convergenceHistory, w, h);
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

function renderLandscape(svg, result, width, height) {
    const margin = { top: 20, right: 20, bottom: 40, left: 50 };
    const { bounds, landscape, bestPath, finalPopulation } = result;

    const x = d3.scaleLinear()
        .domain([bounds.xMin, bounds.xMax])
        .range([margin.left, width - margin.right]);

    const y = d3.scaleLinear()
        .domain([bounds.yMin, bounds.yMax])
        .range([height - margin.bottom, margin.top]);

    // Color scale for function values (log scale for better visualization)
    const zValues = landscape.data.map(d => d.z);
    const zMax = Math.min(200, Math.max(...zValues));  // Cap for better color distribution

    const color = d3.scaleSequential(d3.interpolateViridis)
        .domain([Math.log(zMax + 1), 0]);

    // Draw heatmap
    const cellWidth = (width - margin.left - margin.right) / (landscape.resolution + 1);
    const cellHeight = (height - margin.top - margin.bottom) / (landscape.resolution + 1);

    landscape.data.forEach(d => {
        const logZ = Math.log(Math.min(d.z, zMax) + 1);
        svg.append('rect')
            .attr('x', margin.left + d.i * cellWidth)
            .attr('y', margin.top + (landscape.resolution - d.j) * cellHeight)
            .attr('width', cellWidth + 1)
            .attr('height', cellHeight + 1)
            .attr('fill', color(logZ));
    });

    // Draw known minima
    MINIMA.forEach((m, i) => {
        svg.append('circle')
            .attr('cx', x(m.x))
            .attr('cy', y(m.y))
            .attr('r', 8)
            .attr('fill', 'none')
            .attr('stroke', 'white')
            .attr('stroke-width', 2);

        svg.append('text')
            .attr('x', x(m.x))
            .attr('y', y(m.y) - 12)
            .attr('text-anchor', 'middle')
            .attr('fill', 'white')
            .attr('font-size', '10px')
            .text(`M${i + 1}`);
    });

    // Draw optimization path
    svg.append('path')
        .datum(bestPath)
        .attr('fill', 'none')
        .attr('stroke', 'var(--red)')
        .attr('stroke-width', 2)
        .attr('stroke-dasharray', '4,2')
        .attr('d', d3.line().x(d => x(d.x)).y(d => y(d.y)));

    // Draw final best point
    svg.append('circle')
        .attr('cx', x(result.bestX))
        .attr('cy', y(result.bestY))
        .attr('r', 6)
        .attr('fill', 'var(--red)')
        .attr('stroke', 'white')
        .attr('stroke-width', 2);

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
        .text('x');

    svg.append('text')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', 12)
        .attr('text-anchor', 'middle')
        .attr('fill', 'var(--text-secondary)')
        .attr('font-size', '10px')
        .text('y');
}

function renderConvergence(svg, history, width, height) {
    const margin = { top: 10, right: 15, bottom: 35, left: 60 };

    const x = d3.scaleLinear()
        .domain([0, history.length - 1])
        .range([margin.left, width - margin.right]);

    const yMax = Math.max(...history.map(d => d.avgFitness));
    const y = d3.scaleLog()
        .domain([0.0001, yMax * 1.2])
        .range([height - margin.bottom, margin.top]);

    // Average fitness line
    svg.append('path')
        .datum(history)
        .attr('fill', 'none')
        .attr('stroke', 'var(--orange)')
        .attr('stroke-width', 1.5)
        .attr('d', d3.line()
            .x(d => x(d.generation))
            .y(d => y(Math.max(0.0001, d.avgFitness))));

    // Best fitness line
    svg.append('path')
        .datum(history)
        .attr('fill', 'none')
        .attr('stroke', 'var(--cyan)')
        .attr('stroke-width', 2)
        .attr('d', d3.line()
            .x(d => x(d.generation))
            .y(d => y(Math.max(0.0001, d.bestFitness))));

    // Legend
    const legend = svg.append('g')
        .attr('transform', `translate(${width - margin.right - 80}, ${margin.top + 10})`);

    legend.append('line').attr('x1', 0).attr('y1', 0).attr('x2', 20).attr('y2', 0).attr('stroke', 'var(--cyan)').attr('stroke-width', 2);
    legend.append('text').attr('x', 25).attr('y', 4).attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('Best');

    legend.append('line').attr('x1', 0).attr('y1', 16).attr('x2', 20).attr('y2', 16).attr('stroke', 'var(--orange)').attr('stroke-width', 1.5);
    legend.append('text').attr('x', 25).attr('y', 20).attr('fill', 'var(--text-secondary)').attr('font-size', '9px').text('Average');

    // Axes
    svg.append('g')
        .attr('transform', `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(6))
        .attr('color', 'var(--text-dim)');

    svg.append('g')
        .attr('transform', `translate(${margin.left},0)`)
        .call(d3.axisLeft(y).ticks(4, '.0e'))
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
        .text('Fitness (log scale)');
}
