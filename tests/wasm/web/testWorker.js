// Web Worker for running WASM tests off the main thread
// This keeps the UI responsive while tests execute

let wasmModule = null;
let testCount = 0;
let passCount = 0;
let failCount = 0;

// Override console output capture
self.Module = {
    print: function(text) {
        processOutput(text);
    },
    printErr: function(text) {
        processOutput(text);
    },
    onRuntimeInitialized: function() {
        self.postMessage({ type: 'status', message: 'WASM runtime initialized' });
    }
};

function processOutput(text) {
    if (!text || typeof text !== 'string') return;

    // Send raw output to main thread for console display
    self.postMessage({ type: 'output', text: text });

    // Parse test results
    if (text.startsWith('[PASS]')) {
        testCount++;
        passCount++;
        const testName = text.substring(7).trim();
        self.postMessage({
            type: 'result',
            passed: true,
            name: testName,
            current: testCount,
            passCount: passCount,
            failCount: failCount
        });
    } else if (text.startsWith('[FAIL]')) {
        testCount++;
        failCount++;
        const testName = text.substring(7).trim();
        self.postMessage({
            type: 'result',
            passed: false,
            name: testName,
            current: testCount,
            passCount: passCount,
            failCount: failCount
        });
    } else if (text.includes('===') && text.includes('Tests')) {
        // Category header
        self.postMessage({ type: 'category', name: text });
    } else if (text.includes('Total:')) {
        // Summary line - extract total count
        const match = text.match(/Total:\s*(\d+)/);
        if (match) {
            self.postMessage({ type: 'total', count: parseInt(match[1]) });
        }
    }
}

self.onmessage = async function(e) {
    const { type, wasmUrl } = e.data;

    if (type === 'load') {
        try {
            self.postMessage({ type: 'status', message: 'Loading WASM module...' });

            // Import the WASM module
            importScripts(wasmUrl);

            // Wait for module to be ready
            await waitForModule();
            wasmModule = Module;

            self.postMessage({ type: 'loaded' });
        } catch (error) {
            self.postMessage({ type: 'error', message: error.message });
        }
    } else if (type === 'run') {
        if (!wasmModule) {
            self.postMessage({ type: 'error', message: 'WASM module not loaded' });
            return;
        }

        // Reset counters
        testCount = 0;
        passCount = 0;
        failCount = 0;

        self.postMessage({ type: 'started' });

        try {
            // Run the tests
            if (wasmModule.callMain) {
                wasmModule.callMain([]);
            } else if (wasmModule._main) {
                wasmModule._main();
            } else {
                self.postMessage({ type: 'error', message: 'No entry point found' });
                return;
            }

            // Small delay to ensure all output is flushed
            await new Promise(r => setTimeout(r, 100));

            self.postMessage({
                type: 'finished',
                total: testCount,
                passed: passCount,
                failed: failCount
            });
        } catch (error) {
            self.postMessage({ type: 'error', message: error.message });
        }
    }
};

function waitForModule() {
    return new Promise((resolve, reject) => {
        const timeout = setTimeout(() => reject(new Error('Module load timeout')), 30000);

        const check = () => {
            if (typeof Module !== 'undefined' && Module.calledRun) {
                clearTimeout(timeout);
                resolve();
            } else if (typeof Module !== 'undefined' && Module.onRuntimeInitialized) {
                const original = Module.onRuntimeInitialized;
                Module.onRuntimeInitialized = () => {
                    original();
                    clearTimeout(timeout);
                    resolve();
                };
            } else {
                setTimeout(check, 50);
            }
        };
        check();
    });
}
