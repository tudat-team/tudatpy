const createModule = require('../../build-wasm/tests/wasm/web/tudat_wasm_test.js');

async function test() {
    const Module = await createModule();
    const truthState = { semiMajorAxis: 6793, eccentricity: 0.0001, inclination: 51.6, raan: 45.0, argPeriapsis: 90.0, trueAnomaly: 0.0 };
    const initialGuess = { semiMajorAxis: 6798, eccentricity: 0.0001, inclination: 51.7, raan: 45.2, argPeriapsis: 90.0, trueAnomaly: 0.5 };

    console.log('Testing with 500 observations...');
    const result = Module.runOrbitDetermination(JSON.stringify(initialGuess), JSON.stringify(truthState), 10800, 500, 500, 100, 10, 'fullforce');

    const numIterations = Math.round(result[0]);
    const numObs = Math.round(result[1]);
    const numSamples = Math.round(result[2]);
    console.log('Converged in ' + numIterations + ' iterations, numObs=' + numObs + ', numSamples=' + numSamples);

    let offset = 3;
    const iterSize = 1 + 6 + numObs * 3;
    for (let iter = 0; iter < numIterations; iter++) {
        console.log('Iteration ' + iter + ': RMS = ' + result[offset].toFixed(2) + ' m');
        offset += iterSize;
    }
}
test();
