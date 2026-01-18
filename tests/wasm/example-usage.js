/**
 * Example: Using Tudat WASM from JavaScript/Node.js
 *
 * This example shows how to:
 * 1. Load the Tudat WASM module
 * 2. Mount data files to the virtual filesystem
 * 3. Call Tudat functions from JavaScript
 *
 * Build Tudat with Emscripten first:
 *   mkdir build-wasm && cd build-wasm
 *   cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake_modules/toolchain-emscripten.cmake
 *   make -j$(nproc)
 *
 * Run this example:
 *   node tests/wasm/example-usage.js
 */

const fs = require('fs');
const path = require('path');

// Path to the built WASM module (adjust based on your build)
const WASM_MODULE_PATH = path.join(__dirname, '../../build-wasm/tests/wasm/tudat_wasm_test.js');

async function main() {
    console.log('='.repeat(50));
    console.log('  Tudat WASM JavaScript Usage Example');
    console.log('='.repeat(50));

    // Check if the WASM module exists
    if (!fs.existsSync(WASM_MODULE_PATH)) {
        console.error(`\nError: WASM module not found at ${WASM_MODULE_PATH}`);
        console.error('\nPlease build Tudat with Emscripten first:');
        console.error('  mkdir build-wasm && cd build-wasm');
        console.error('  cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake_modules/toolchain-emscripten.cmake');
        console.error('  make -j$(nproc)');
        process.exit(1);
    }

    console.log('\n[1] Loading WASM module...');

    // For a modularized build (MODULARIZE=1), you would load it like this:
    // const createTudatModule = require(WASM_MODULE_PATH);
    // const Module = await createTudatModule({
    //     // Optional: Pre-configure the module
    //     print: (text) => console.log('[Tudat]', text),
    //     printErr: (text) => console.error('[Tudat Error]', text),
    //
    //     // Pre-run callback to set up the virtual filesystem
    //     preRun: [function(Module) {
    //         // Create directories for data files
    //         Module.FS.mkdir('/tudat_data');
    //         Module.FS.mkdir('/tudat_data/ephemeris');
    //         Module.FS.mkdir('/tudat_data/spice_kernels');
    //
    //         // Mount local files to the virtual filesystem
    //         // Example: Mount a SPICE kernel
    //         // const kernelData = fs.readFileSync('/path/to/kernel.bsp');
    //         // Module.FS.writeFile('/tudat_data/spice_kernels/kernel.bsp', kernelData);
    //     }]
    // });

    // For the current non-modularized test build, we can run it directly:
    console.log('[2] Running the test module...');
    console.log('-'.repeat(50));

    // Execute the WASM test
    require(WASM_MODULE_PATH);
}

/**
 * Example: How to mount data files to the virtual filesystem
 *
 * When using Tudat WASM in a browser or Node.js, you need to provide
 * data files (ephemeris, SPICE kernels, etc.) via Emscripten's virtual FS.
 *
 * Browser example with pre-loaded data:
 * ```javascript
 * const Module = await createTudatModule({
 *     preRun: [(Module) => {
 *         // Create the data directory structure
 *         Module.FS.mkdir('/tudat_data');
 *         Module.FS.mkdir('/tudat_data/spice_kernels');
 *
 *         // Write data files fetched from your server
 *         const kernelData = new Uint8Array(await fetch('/data/de430.bsp').then(r => r.arrayBuffer()));
 *         Module.FS.writeFile('/tudat_data/spice_kernels/de430.bsp', kernelData);
 *     }]
 * });
 * ```
 *
 * Node.js example:
 * ```javascript
 * const fs = require('fs');
 * const createTudatModule = require('./tudat_module.js');
 *
 * const Module = await createTudatModule({
 *     preRun: [(Module) => {
 *         Module.FS.mkdir('/tudat_data');
 *         Module.FS.mkdir('/tudat_data/spice_kernels');
 *
 *         // Read local file and write to virtual FS
 *         const kernelData = fs.readFileSync('./data/de430.bsp');
 *         Module.FS.writeFile('/tudat_data/spice_kernels/de430.bsp', kernelData);
 *     }]
 * });
 * ```
 */

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
