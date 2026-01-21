#!/usr/bin/env node
// Quick test for orbit determination convergence

const createModule = require('../../build-wasm/tests/wasm/web/tudat_wasm_test.js');

async function test() {
    const Module = await createModule();

    const truthState = {
        semiMajorAxis: 6793,
        eccentricity: 0.0001,
        inclination: 51.6,
        raan: 45.0,
        argPeriapsis: 90.0,
        trueAnomaly: 0.0
    };

    const initialGuess = {
        semiMajorAxis: 6793 + 5,
        eccentricity: 0.0001,
        inclination: 51.6 + 0.1,
        raan: 45.0 + 0.2,
        argPeriapsis: 90.0,
        trueAnomaly: 0.5
    };

    const period = 5400;
    const duration = period * 2;
    const numObservations = 50;
    const numOrbitSamples = 500;
    const noiseStdDev = 100;
    const maxIterations = 10;

    console.log('Running orbit determination...');
    console.log(`Truth SMA: ${truthState.semiMajorAxis} km`);
    console.log(`Initial guess SMA: ${initialGuess.semiMajorAxis} km`);
    console.log(`Observations: ${numObservations}, Noise: ${noiseStdDev} m`);
    console.log('');

    const result = Module.runOrbitDetermination(
        JSON.stringify(initialGuess),
        JSON.stringify(truthState),
        duration,
        numObservations,
        numOrbitSamples,
        noiseStdDev,
        maxIterations,
        'fullforce'
    );

    const numIterations = Math.round(result[0]);
    const numObs = Math.round(result[1]);
    const numSamples = Math.round(result[2]);

    console.log(`Converged in ${numIterations} iterations (samples: ${numSamples})`);

    const iterationDataSize = 1 + 6 + numObs * 3;
    let offset = 3;

    for (let iter = 0; iter < numIterations; iter++) {
        const rms = result[offset];
        console.log(`Iteration ${iter}: RMS = ${rms.toFixed(2)} m`);
        offset += iterationDataSize;
    }

    console.log('');
    if (result[offset - iterationDataSize] < noiseStdDev * 2) {
        console.log('SUCCESS: RMS converged to near noise level');
    } else {
        console.log('FAILURE: RMS did not converge');
        process.exit(1);
    }
}

test().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
