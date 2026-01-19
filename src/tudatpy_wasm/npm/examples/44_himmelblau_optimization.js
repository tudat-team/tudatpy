/**
 * Example 44: Himmelblau Function Optimization
 *
 * Ported from: examples/tudatpy/pygmo/himmelblau_minimization.py
 *
 * This example demonstrates optimization fundamentals using the
 * Himmelblau function, a classic test function for optimization
 * algorithms with multiple global minima.
 *
 * Key concepts:
 * - User-defined optimization problems
 * - Fitness function evaluation
 * - Bounds and constraints
 * - Multiple optima detection
 * - Evolutionary algorithms
 *
 * Run with: node 44_himmelblau_optimization.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Himmelblau Function Optimization ===\n');

    const tudat = await createTudatModule();

    /**
     * Himmelblau function:
     * f(x, y) = (x² + y - 11)² + (x + y² - 7)²
     *
     * This function has four identical global minima at:
     * (3.0, 2.0)
     * (-2.805118, 3.131312)
     * (-3.779310, -3.283186)
     * (3.584428, -1.848126)
     *
     * All minima have f(x,y) = 0
     */
    function himmelblau(x, y) {
        const term1 = Math.pow(x * x + y - 11, 2);
        const term2 = Math.pow(x + y * y - 7, 2);
        return term1 + term2;
    }

    // Problem bounds
    const xMin = -5, xMax = 5;
    const yMin = -5, yMax = 5;

    console.log('Himmelblau Function: f(x,y) = (x² + y - 11)² + (x + y² - 7)²');
    console.log(`Search bounds: x ∈ [${xMin}, ${xMax}], y ∈ [${yMin}, ${yMax}]\n`);

    // Known optima
    const knownOptima = [
        { x: 3.0, y: 2.0 },
        { x: -2.805118, y: 3.131312 },
        { x: -3.779310, y: -3.283186 },
        { x: 3.584428, y: -1.848126 }
    ];

    console.log('Known Global Minima (f = 0):');
    knownOptima.forEach((opt, i) => {
        const f = himmelblau(opt.x, opt.y);
        console.log(`  ${i + 1}. (${opt.x.toFixed(6)}, ${opt.y.toFixed(6)}) -> f = ${f.toExponential(2)}`);
    });

    // ========================================
    // Simple Genetic Algorithm Implementation
    // ========================================
    console.log('\n--- Running Genetic Algorithm ---\n');

    class GeneticAlgorithm {
        constructor(fitnessFunc, bounds, popSize = 50, mutationRate = 0.1) {
            this.fitnessFunc = fitnessFunc;
            this.bounds = bounds;  // [[xMin, xMax], [yMin, yMax]]
            this.popSize = popSize;
            this.mutationRate = mutationRate;
            this.dimension = bounds.length;
        }

        randomIndividual() {
            return this.bounds.map(b => b[0] + Math.random() * (b[1] - b[0]));
        }

        initializePopulation() {
            return Array(this.popSize).fill(null).map(() => this.randomIndividual());
        }

        evaluate(population) {
            return population.map(ind => ({
                genes: ind,
                fitness: this.fitnessFunc(ind[0], ind[1])
            }));
        }

        select(evaluated) {
            // Tournament selection
            const tournamentSize = 3;
            const selected = [];

            for (let i = 0; i < this.popSize; i++) {
                const tournament = [];
                for (let j = 0; j < tournamentSize; j++) {
                    const idx = Math.floor(Math.random() * evaluated.length);
                    tournament.push(evaluated[idx]);
                }
                tournament.sort((a, b) => a.fitness - b.fitness);
                selected.push([...tournament[0].genes]);
            }

            return selected;
        }

        crossover(parent1, parent2) {
            // Arithmetic crossover
            const alpha = Math.random();
            return parent1.map((g, i) => alpha * g + (1 - alpha) * parent2[i]);
        }

        mutate(individual) {
            return individual.map((g, i) => {
                if (Math.random() < this.mutationRate) {
                    const range = this.bounds[i][1] - this.bounds[i][0];
                    const mutation = (Math.random() - 0.5) * range * 0.2;
                    return Math.max(this.bounds[i][0], Math.min(this.bounds[i][1], g + mutation));
                }
                return g;
            });
        }

        evolve(population) {
            const evaluated = this.evaluate(population);
            const selected = this.select(evaluated);

            const newPop = [];
            for (let i = 0; i < this.popSize; i += 2) {
                const p1 = selected[i];
                const p2 = selected[(i + 1) % this.popSize];

                let child1 = this.crossover(p1, p2);
                let child2 = this.crossover(p2, p1);

                child1 = this.mutate(child1);
                child2 = this.mutate(child2);

                newPop.push(child1, child2);
            }

            return newPop.slice(0, this.popSize);
        }

        run(generations) {
            let population = this.initializePopulation();
            let bestEver = { fitness: Infinity };
            const history = [];

            for (let gen = 0; gen < generations; gen++) {
                const evaluated = this.evaluate(population);
                const best = evaluated.reduce((a, b) => a.fitness < b.fitness ? a : b);

                if (best.fitness < bestEver.fitness) {
                    bestEver = { ...best };
                }

                history.push({
                    generation: gen,
                    bestFitness: best.fitness,
                    avgFitness: evaluated.reduce((s, e) => s + e.fitness, 0) / evaluated.length
                });

                population = this.evolve(population);
            }

            return { bestEver, history };
        }
    }

    // Run GA multiple times to find different optima
    const bounds = [[xMin, xMax], [yMin, yMax]];
    const ga = new GeneticAlgorithm(himmelblau, bounds, 50, 0.15);

    console.log('GA Parameters:');
    console.log('  Population size: 50');
    console.log('  Generations: 100');
    console.log('  Mutation rate: 15%');
    console.log('  Selection: Tournament (size 3)');
    console.log('  Crossover: Arithmetic\n');

    // Run multiple times to find different optima
    const numRuns = 8;
    const foundOptima = [];

    console.log('Running optimization...\n');

    for (let run = 0; run < numRuns; run++) {
        const result = ga.run(100);
        const x = result.bestEver.genes[0];
        const y = result.bestEver.genes[1];
        const f = result.bestEver.fitness;

        // Check if this is a new optimum (not already found)
        const isNew = !foundOptima.some(opt =>
            Math.sqrt((opt.x - x)**2 + (opt.y - y)**2) < 0.5
        );

        if (isNew && f < 0.01) {
            foundOptima.push({ x, y, f, run: run + 1 });
        }
    }

    console.log('Found Optima:');
    foundOptima.forEach((opt, i) => {
        // Find closest known optimum
        let minDist = Infinity;
        let closestIdx = 0;
        knownOptima.forEach((known, j) => {
            const dist = Math.sqrt((known.x - opt.x)**2 + (known.y - opt.y)**2);
            if (dist < minDist) {
                minDist = dist;
                closestIdx = j + 1;
            }
        });

        console.log(`  ${i + 1}. (${opt.x.toFixed(4)}, ${opt.y.toFixed(4)}) -> f = ${opt.f.toExponential(2)} (near optimum ${closestIdx})`);
    });

    console.log(`\nFound ${foundOptima.length} of 4 global optima`);

    // ========================================
    // Differential Evolution for comparison
    // ========================================
    console.log('\n--- Running Differential Evolution ---\n');

    class DifferentialEvolution {
        constructor(fitnessFunc, bounds, popSize = 30, F = 0.8, CR = 0.9) {
            this.fitnessFunc = fitnessFunc;
            this.bounds = bounds;
            this.popSize = popSize;
            this.F = F;   // Mutation factor
            this.CR = CR; // Crossover rate
            this.dimension = bounds.length;
        }

        randomIndividual() {
            return this.bounds.map(b => b[0] + Math.random() * (b[1] - b[0]));
        }

        run(generations) {
            // Initialize population
            let population = Array(this.popSize).fill(null).map(() => this.randomIndividual());
            let fitness = population.map(ind => this.fitnessFunc(ind[0], ind[1]));

            let bestIdx = fitness.indexOf(Math.min(...fitness));
            let bestEver = { genes: [...population[bestIdx]], fitness: fitness[bestIdx] };

            for (let gen = 0; gen < generations; gen++) {
                for (let i = 0; i < this.popSize; i++) {
                    // Select three random individuals (not i)
                    const candidates = [...Array(this.popSize).keys()].filter(j => j !== i);
                    const [a, b, c] = candidates.sort(() => Math.random() - 0.5).slice(0, 3);

                    // Mutation
                    const mutant = population[a].map((g, d) => {
                        let val = g + this.F * (population[b][d] - population[c][d]);
                        return Math.max(this.bounds[d][0], Math.min(this.bounds[d][1], val));
                    });

                    // Crossover
                    const jRand = Math.floor(Math.random() * this.dimension);
                    const trial = population[i].map((g, d) =>
                        (Math.random() < this.CR || d === jRand) ? mutant[d] : g
                    );

                    // Selection
                    const trialFitness = this.fitnessFunc(trial[0], trial[1]);
                    if (trialFitness < fitness[i]) {
                        population[i] = trial;
                        fitness[i] = trialFitness;

                        if (trialFitness < bestEver.fitness) {
                            bestEver = { genes: [...trial], fitness: trialFitness };
                        }
                    }
                }
            }

            return bestEver;
        }
    }

    const de = new DifferentialEvolution(himmelblau, bounds, 30, 0.8, 0.9);

    console.log('DE Parameters:');
    console.log('  Population size: 30');
    console.log('  Generations: 100');
    console.log('  Mutation factor F: 0.8');
    console.log('  Crossover rate CR: 0.9\n');

    const foundDE = [];
    for (let run = 0; run < numRuns; run++) {
        const result = de.run(100);
        const x = result.genes[0];
        const y = result.genes[1];
        const f = result.fitness;

        const isNew = !foundDE.some(opt =>
            Math.sqrt((opt.x - x)**2 + (opt.y - y)**2) < 0.5
        );

        if (isNew && f < 0.01) {
            foundDE.push({ x, y, f });
        }
    }

    console.log('DE Found Optima:');
    foundDE.forEach((opt, i) => {
        console.log(`  ${i + 1}. (${opt.x.toFixed(4)}, ${opt.y.toFixed(4)}) -> f = ${opt.f.toExponential(2)}`);
    });
    console.log(`\nFound ${foundDE.length} of 4 global optima`);

    // ========================================
    // Function Landscape Analysis
    // ========================================
    console.log('\n--- Function Landscape Analysis ---\n');

    // Sample the function on a grid
    const gridSize = 21;
    const xVals = Array(gridSize).fill(0).map((_, i) => xMin + i * (xMax - xMin) / (gridSize - 1));
    const yVals = Array(gridSize).fill(0).map((_, i) => yMin + i * (yMax - yMin) / (gridSize - 1));

    let minVal = Infinity, maxVal = -Infinity;
    for (const x of xVals) {
        for (const y of yVals) {
            const f = himmelblau(x, y);
            minVal = Math.min(minVal, f);
            maxVal = Math.max(maxVal, f);
        }
    }

    console.log('Function Statistics:');
    console.log(`  Minimum value: ${minVal.toFixed(6)}`);
    console.log(`  Maximum value (in domain): ${maxVal.toFixed(1)}`);
    console.log(`  Dynamic range: ${(maxVal / (minVal + 1e-10)).toExponential(2)}`);

    // Gradient analysis at a non-optimal point
    const testX = 0, testY = 0;
    const eps = 1e-6;
    const gradX = (himmelblau(testX + eps, testY) - himmelblau(testX - eps, testY)) / (2 * eps);
    const gradY = (himmelblau(testX, testY + eps) - himmelblau(testX, testY - eps)) / (2 * eps);

    console.log(`\nGradient at (0, 0):`);
    console.log(`  ∂f/∂x = ${gradX.toFixed(2)}`);
    console.log(`  ∂f/∂y = ${gradY.toFixed(2)}`);
    console.log(`  |∇f| = ${Math.sqrt(gradX*gradX + gradY*gradY).toFixed(2)}`);

    // Summary
    console.log('\n--- Summary ---\n');
    console.log('The Himmelblau function is useful for testing because:');
    console.log('  1. Multiple global optima test multimodal search');
    console.log('  2. Symmetric structure helps verify algorithm correctness');
    console.log('  3. Known analytical optima enable accuracy assessment');
    console.log('  4. Smooth landscape suitable for gradient-based methods');

    console.log('\n=== Himmelblau Optimization Complete ===');
}

main().catch(console.error);
