// Tudat WASM Test Runner - CesiumJS Edition
// Handles WASM test execution, real-time UI updates, and visualizations
const APP_BUILD = '20260120-v1';
console.log(`app.js build: ${APP_BUILD}`);

import { SpiceKernelLoader } from './spice-loader.js';
import { initSpice } from './visualizations/shared/spice-utils.js';

import {
    // Shared utilities (used in default visualization case)
    configureClockForOrbit,
    addAnimatedOrbit,
    clearOrbitEntities,
    // Visualization show functions
    showCR3BPVisualization,
    showLibrationPointsVisualization,
    showAtmosphericDragVisualization,
    showReferenceFramesVisualization,
    showGeostationaryVisualization,
    // Note: J2 vs Full Force uses class method for chart integration
    // Python example ports (chart-only)
    showKeplerianOrbitExample,
    showPerturbedOrbitExample,
    showReentryTrajectoryExample,
    showSolarSystemExample,
    showThrustSatelliteExample,
    showTwoStageRocketExample,
    showLinearSensitivityExample,
    showCoupledDynamicsExample,
    showCR3BPManifoldsExample,
    showDifferentialDragExample,
    showJuiceFlybysExample,
    showEarthMoonThrustExample,
    showEarthMarsTransferExample,
    showMGATrajectoryExample,
    showHohmannTransferExample,
    showGravityAssistExample,
    showCassiniMGAExample,
    showLowThrustPorkchopExample,
    showCovariancePropagationExample,
    showFullEstimationExample,
    showGalileanMoonsEstimationExample,
    showEstimationDynamicalModelsExample,
    showMPCAsteroidEstimationExample,
    showHimmelblauOptimizationExample,
    showAsteroidOrbitOptimizationExample,
    showHodographicShapingMGAExample,
    // Registry for dynamic category list
    visualizationRegistry
} from './visualizations/index.js';

class TudatTestRunner {
    constructor() {
        this.wasmModule = null;
        this.tudatModule = null;  // Full tudatpy WASM module with SPICE
        this.spiceLoader = null;  // SPICE kernel loader
        this.spiceReady = false;  // Flag indicating SPICE kernels are loaded
        this.testResults = [];
        this.categories = {};
        this.isRunning = false;
        this.startTime = null;
        this.currentCategory = 'General';
        this.consoleLines = 0;
        this.expectedTests = 551;
        this.currentVisualization = null;  // Track current visualization for URL updates

        // Web Worker for running tests off main thread
        this.testWorker = null;
        this.workerReady = false;
        this.totalTestCount = 551;  // Expected total tests

        // Modal elements
        this.modal = null;
        this.modalProgressBar = null;
        this.modalProgressText = null;
        this.modalProgressLabel = null;
        this.modalPassed = null;
        this.modalFailed = null;
        this.modalCurrentTest = null;

        // Cesium viewer
        this.viewer = null;
        this.orbitEntities = [];

        // Charts
        this.charts = {};

        // Currently selected test
        this.selectedTest = null;

        // Orbital data for visualization
        this.orbitalData = {};

        this.init();
    }

    async init() {
        this.updateTimestamp();
        setInterval(() => this.updateTimestamp(), 1000);

        // Display build number in UI
        document.getElementById('app-build').textContent = APP_BUILD;

        this.setupEventListeners();
        this.setupCharts();
        this.setupModal();
        await this.setupCesium();
        this.generateOrbitalData();
        this.setupVisualizationCategories();

        await this.loadWasmWorker();

        // Load SPICE kernels for ephemeris-based visualizations
        await this.loadSpiceKernels();
    }

    setupModal() {
        this.modal = document.getElementById('test-modal');
        this.modalProgressBar = document.getElementById('modal-progress-bar');
        this.modalProgressText = document.getElementById('modal-progress-text');
        this.modalProgressLabel = document.getElementById('modal-progress-label');
        this.modalPassed = document.getElementById('modal-passed');
        this.modalFailed = document.getElementById('modal-failed');
        this.modalCurrentTest = document.getElementById('modal-current-test');
    }

    showModal() {
        if (this.modal) {
            this.modal.classList.add('active');
            this.modalProgressBar.style.width = '0%';
            this.modalProgressText.textContent = '0%';
            this.modalProgressLabel.textContent = 'Initializing...';
            this.modalPassed.textContent = '0';
            this.modalFailed.textContent = '0';
            this.modalCurrentTest.textContent = 'Starting tests...';
        }
    }

    hideModal() {
        if (this.modal) {
            this.modal.classList.remove('active');
        }
    }

    updateModal(current, total, passed, failed, testName) {
        if (!this.modal) return;

        const percent = total > 0 ? Math.round((current / total) * 100) : 0;
        this.modalProgressBar.style.width = `${percent}%`;
        this.modalProgressText.textContent = `${percent}%`;
        this.modalProgressLabel.textContent = `${current} / ${total} tests`;
        this.modalPassed.textContent = passed;
        this.modalFailed.textContent = failed;
        if (testName) {
            this.modalCurrentTest.textContent = testName;
        }
    }

    // URL helpers for deep-linking to visualizations
    vizNameToSlug(name) {
        return name.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/(^-|-$)/g, '');
    }

    slugToVizName(slug) {
        // Find the visualization name that matches this slug
        for (const name of Object.keys(visualizationRegistry)) {
            if (this.vizNameToSlug(name) === slug) {
                return name;
            }
        }
        return null;
    }

    updateUrlWithVisualization(vizName) {
        const slug = this.vizNameToSlug(vizName);
        const newUrl = `${window.location.pathname}#${slug}`;
        console.log(`[URL] Updating URL to: ${newUrl}`);
        window.history.replaceState(null, '', newUrl);
        this.currentVisualization = vizName;
    }

    getVisualizationFromUrl() {
        const hash = window.location.hash.slice(1); // Remove the #
        if (!hash) return null;
        return this.slugToVizName(hash);
    }

    // Curated list of test categories with meaningful 3D/chart visualizations
    // Uses the visualization registry from the modular visualization system
    setupVisualizationCategories() {
        const vizCategories = Object.entries(visualizationRegistry).map(([name, config]) => ({
            name,
            description: config.description,
            category: name,
            testName: name
        }));

        const container = document.getElementById('viz-category-list');
        container.innerHTML = '';

        vizCategories.forEach(viz => {
            const item = document.createElement('div');
            item.className = 'category-item viz-category';
            item.innerHTML = `
                <div class="category-header">
                    <span class="category-name">${viz.name}</span>
                </div>
                <div class="viz-description">${viz.description}</div>
            `;

            item.addEventListener('click', () => {
                // Remove selected from others
                container.querySelectorAll('.viz-category').forEach(el => el.classList.remove('selected'));
                item.classList.add('selected');

                // Trigger visualization
                this.visualizeTest(viz.testName, viz.category);

                // For CR3BP, auto-select L2 Halo orbit after a short delay
                if (viz.category.includes('CR3BP')) {
                    setTimeout(() => this.selectOrbit('l2-halo'), 100);
                }
            });

            container.appendChild(item);
        });
    }

    // Fuzzy search filter for visualization examples
    filterVisualizations(query) {
        const container = document.getElementById('viz-category-list');
        const items = container.querySelectorAll('.viz-category');
        const normalizedQuery = query.toLowerCase().trim();

        if (!normalizedQuery) {
            // Show all items and remove highlights when search is cleared
            items.forEach(item => {
                item.classList.remove('hidden');
                const nameEl = item.querySelector('.category-name');
                const descEl = item.querySelector('.viz-description');
                if (nameEl) nameEl.innerHTML = nameEl.textContent;
                if (descEl) descEl.innerHTML = descEl.textContent;
            });
            return;
        }

        items.forEach(item => {
            const nameEl = item.querySelector('.category-name');
            const descEl = item.querySelector('.viz-description');
            const name = nameEl?.textContent || '';
            const desc = descEl?.textContent || '';
            const combined = (name + ' ' + desc).toLowerCase();

            // Fuzzy match: check if all characters in query appear in order
            const matches = this.fuzzyMatch(combined, normalizedQuery);

            if (matches) {
                item.classList.remove('hidden');
                // Highlight matched characters in name
                if (nameEl) nameEl.innerHTML = this.highlightMatches(name, normalizedQuery);
                if (descEl) descEl.innerHTML = this.highlightMatches(desc, normalizedQuery);
            } else {
                item.classList.add('hidden');
            }
        });
    }

    // Smart matching - prioritizes exact words, prefixes, and handles common typos
    fuzzyMatch(text, pattern) {
        // Split text into words
        const words = text.split(/\s+/);

        // 1. Exact word match
        if (words.some(word => word === pattern)) {
            return true;
        }

        // 2. Word starts with pattern (prefix match)
        if (words.some(word => word.startsWith(pattern))) {
            return true;
        }

        // 3. Pattern is contained as a substring in any word
        if (words.some(word => word.includes(pattern))) {
            return true;
        }

        // 4. Check for single character transposition (common typo)
        if (pattern.length >= 2) {
            for (let i = 0; i < pattern.length - 1; i++) {
                const transposed = pattern.slice(0, i) + pattern[i + 1] + pattern[i] + pattern.slice(i + 2);
                if (words.some(word => word.includes(transposed))) {
                    return true;
                }
            }
        }

        // 5. Check for single missing/extra character (length diff of 1)
        if (pattern.length >= 3) {
            for (const word of words) {
                if (Math.abs(word.length - pattern.length) <= 1) {
                    let diffs = 0;
                    const longer = word.length >= pattern.length ? word : pattern;
                    const shorter = word.length < pattern.length ? word : pattern;
                    let j = 0;
                    for (let i = 0; i < longer.length && diffs <= 1; i++) {
                        if (longer[i] === shorter[j]) {
                            j++;
                        } else {
                            diffs++;
                            if (longer.length === shorter.length) j++;
                        }
                    }
                    if (diffs <= 1 && j >= shorter.length - 1) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    // Highlight matching words/substrings in text
    highlightMatches(text, pattern) {
        const lowerText = text.toLowerCase();
        const words = text.split(/(\s+)/); // Keep whitespace in result

        return words.map(word => {
            const lowerWord = word.toLowerCase();
            const idx = lowerWord.indexOf(pattern);
            if (idx !== -1) {
                // Highlight the matched portion
                return word.slice(0, idx) +
                       `<span class="search-highlight">${word.slice(idx, idx + pattern.length)}</span>` +
                       word.slice(idx + pattern.length);
            }
            return word;
        }).join('');
    }

    updateTimestamp() {
        const now = new Date();
        document.getElementById('timestamp').textContent =
            now.toISOString().replace('T', ' ').substring(0, 19) + ' UTC';
    }

    selectDefaultVisualization() {
        // Check URL for a visualization to load, otherwise default to J2 vs Full Force
        // Use setTimeout to ensure DOM layout is complete before rendering charts
        setTimeout(() => {
            const container = document.getElementById('viz-category-list');
            const vizItems = container.querySelectorAll('.viz-category');

            // Check if URL specifies a visualization
            const urlViz = this.getVisualizationFromUrl();
            const targetViz = urlViz || 'J2 vs Full Force';
            console.log(`[URL] Loading visualization: ${targetViz} (from URL: ${urlViz !== null})`);

            // Find and select the target visualization item
            vizItems.forEach(item => {
                const nameEl = item.querySelector('.category-name');
                if (nameEl && nameEl.textContent === targetViz) {
                    item.classList.add('selected');
                    // Trigger the visualization
                    this.visualizeTest(targetViz, targetViz);

                    // For CR3BP, auto-select L2 Halo orbit
                    if (targetViz.includes('CR3BP')) {
                        setTimeout(() => this.selectOrbit('l2-halo'), 100);
                    }
                }
            });
        }, 100);  // Small delay to ensure layout is complete
    }

    setupEventListeners() {
        document.getElementById('run-btn').addEventListener('click', () => this.runTests());
        document.getElementById('clear-btn').addEventListener('click', () => this.clearResults());

        // Search input with fuzzy matching
        const searchInput = document.getElementById('viz-search');
        const searchClear = document.getElementById('viz-search-clear');
        if (searchInput) {
            searchInput.addEventListener('input', (e) => {
                this.filterVisualizations(e.target.value);
                // Toggle clear button visibility
                if (searchClear) {
                    searchClear.classList.toggle('visible', e.target.value.length > 0);
                    searchInput.classList.toggle('has-text', e.target.value.length > 0);
                }
            });
        }
        if (searchClear) {
            searchClear.addEventListener('click', () => {
                if (searchInput) {
                    searchInput.value = '';
                    searchInput.classList.remove('has-text');
                    this.filterVisualizations('');
                    searchInput.focus();
                }
                searchClear.classList.remove('visible');
            });
        }

        // Orbit selector buttons
        document.getElementById('orbit-l2-halo').addEventListener('click', () => this.selectOrbit('l2-halo'));
        document.getElementById('orbit-dro').addEventListener('click', () => this.selectOrbit('dro'));
        document.getElementById('orbit-l1-halo').addEventListener('click', () => this.selectOrbit('l1-halo'));
        document.getElementById('orbit-lyapunov').addEventListener('click', () => this.selectOrbit('lyapunov'));

        // Orbit determination dynamics model toggle
        document.getElementById('od-fullforce').addEventListener('click', () => this.selectODModel('fullforce'));
        document.getElementById('od-omm').addEventListener('click', () => this.selectODModel('omm'));

        // Handle browser back/forward navigation
        window.addEventListener('popstate', () => this.handleUrlChange());
    }

    handleUrlChange() {
        const vizName = this.getVisualizationFromUrl();
        if (vizName && vizName !== this.currentVisualization) {
            const container = document.getElementById('viz-category-list');
            const vizItems = container.querySelectorAll('.viz-category');

            // Find and select the matching visualization item
            vizItems.forEach(item => {
                const nameEl = item.querySelector('.category-name');
                if (nameEl && nameEl.textContent === vizName) {
                    // Remove selected from others
                    container.querySelectorAll('.viz-category').forEach(el => el.classList.remove('selected'));
                    item.classList.add('selected');
                    // Trigger the visualization (but don't update URL again)
                    this.currentVisualization = vizName;
                    document.getElementById('orbit-info').textContent = `${vizName}: ${vizName}`;
                    this.clearOrbitEntities();
                    const vizConfig = visualizationRegistry[vizName];
                    if (vizConfig?.chartOnly) {
                        this.showChartOnlyVisualization(vizName, vizName);
                    } else {
                        this.restoreGlobeLayout();
                        this.show3DVisualization(vizName, vizName);
                        this.showChartForCategory(vizName, vizName);
                    }
                }
            });
        }
    }

    selectODModel(model) {
        this.log(`Selecting OD dynamics model: ${model.toUpperCase()}`, 'info');

        // Update button states
        document.getElementById('od-fullforce').classList.toggle('active', model === 'fullforce');
        document.getElementById('od-omm').classList.toggle('active', model === 'omm');

        // Store selected model
        this.currentODModel = model;

        // Re-run visualization with new model
        if (this.viewer) {
            this.addOrbitDeterminationVisualization(model);
        }
    }

    selectOrbit(orbitType, animate = true) {
        this.log(`Selecting orbit: ${orbitType}`, 'info');

        // Update button states
        document.querySelectorAll('.orbit-btn').forEach(btn => btn.classList.remove('active'));
        document.getElementById(`orbit-${orbitType}`).classList.add('active');

        // Store selected orbit for when Cesium is ready
        this.selectedOrbit = orbitType;

        // If Cesium viewer isn't ready yet, skip visualization (will be called again when ready)
        if (!this.viewer) {
            this.log('Cesium viewer not ready, deferring visualization', 'info');
            return;
        }

        // Use imported CR3BP visualization
        const result = showCR3BPVisualization(
            this.viewer,
            this.orbitEntities,
            orbitType,
            (msg, level) => this.log(msg, level),
            animate
        );

        // Update orbit info header
        const orbitInfo = document.getElementById('orbit-info');
        if (orbitInfo) {
            orbitInfo.textContent = `${result.name} • Period: ${result.periodDays} days`;
        }
    }

    // ==================== Cesium Setup ====================

    async setupCesium() {
        // No Ion token needed - we use local imagery files
        // Set base URL for local Cesium assets (Workers, etc.)
        window.CESIUM_BASE_URL = 'cesium/';

        // Create Blue Marble imagery provider (day texture)
        const blueMarbleProvider = await Cesium.SingleTileImageryProvider.fromUrl(
            'imagery/world.topo.bathy.200407.3x5400x2700.jpg',
            {
                rectangle: Cesium.Rectangle.fromDegrees(-180.0, -90.0, 180.0, 90.0),
                tileWidth: 5400,
                tileHeight: 2700
            }
        );

        // Set up clock for orbit animation - default to J2000 epoch
        const j2000 = Cesium.JulianDate.fromIso8601('2000-01-01T12:00:00Z');
        const clock = new Cesium.Clock({
            startTime: j2000,
            currentTime: j2000,
            stopTime: Cesium.JulianDate.addSeconds(j2000, 86400, new Cesium.JulianDate()), // 1 day
            clockRange: Cesium.ClockRange.LOOP_STOP,
            clockStep: Cesium.ClockStep.SYSTEM_CLOCK_MULTIPLIER,
            multiplier: 60, // 1 minute per second by default
            shouldAnimate: false
        });

        // Create viewer WITH animation and timeline widgets
        this.viewer = new Cesium.Viewer('cesiumContainer', {
            animation: true,
            timeline: true,
            baseLayerPicker: false,
            fullscreenButton: false,
            geocoder: false,
            homeButton: true,
            infoBox: true,
            sceneModePicker: false,
            selectionIndicator: true,
            navigationHelpButton: false,
            creditContainer: document.createElement('div'),
            baseLayer: new Cesium.ImageryLayer(blueMarbleProvider),
            skyBox: false,
            skyAtmosphere: new Cesium.SkyAtmosphere(),
            clockViewModel: new Cesium.ClockViewModel(clock),
            contextOptions: {
                webgl: {
                    alpha: true
                }
            }
        });

        // Style the animation widget
        this.styleAnimationWidgets();

        // Dark space background
        this.viewer.scene.backgroundColor = Cesium.Color.fromCssColorString('#020408');
        this.viewer.scene.globe.baseColor = Cesium.Color.fromCssColorString('#0a1428');

        // Enable lighting for day/night effect
        this.viewer.scene.globe.enableLighting = true;

        // Add night lights layer (Black Marble)
        const blackMarbleProvider = await Cesium.SingleTileImageryProvider.fromUrl(
            'imagery/BlackMarble_2016_3km.jpeg',
            {
                rectangle: Cesium.Rectangle.fromDegrees(-180.0, -90.0, 180.0, 90.0),
                tileWidth: 13500,
                tileHeight: 6750
            }
        );

        this.nightLightsLayer = this.viewer.imageryLayers.addImageryProvider(blackMarbleProvider);

        // Configure night lights layer
        this.nightLightsLayer.dayAlpha = 0.0;      // Hide during day
        this.nightLightsLayer.nightAlpha = 1.0;   // Show at night
        this.nightLightsLayer.brightness = 2.0;   // Boost brightness
        this.nightLightsLayer.contrast = 1.2;
        this.nightLightsLayer.gamma = 0.6;
        this.nightLightsLayer.saturation = 1.2;

        // Handle resize
        window.addEventListener('resize', () => {
            if (this.viewer) {
                this.viewer.resize();
            }
        });

        // Fade night lights based on altitude (only show from space)
        this.viewer.scene.postRender.addEventListener(() => {
            if (this.nightLightsLayer) {
                const height = this.viewer.camera.positionCartographic.height;
                const FADE_START = 70000;  // Start fading in at 70km
                const FADE_END = 50000;    // Fully hidden below 50km

                if (height <= FADE_END) {
                    this.nightLightsLayer.alpha = 0;
                } else if (height >= FADE_START) {
                    this.nightLightsLayer.alpha = 1;
                } else {
                    this.nightLightsLayer.alpha = (height - FADE_END) / (FADE_START - FADE_END);
                }
            }
        });

        // Add equatorial plane as a circle
        this.addEquatorialPlane();
    }

    styleAnimationWidgets() {
        // Styling handled in index.html CSS - this method kept for compatibility
    }

    // Configure clock for a specific orbital period
    configureClockForOrbit(periodSeconds, epochDate = null, multiplier = null) {
        const clock = this.viewer.clock;
        const start = epochDate || Cesium.JulianDate.fromIso8601('2000-01-01T12:00:00Z');
        const stop = Cesium.JulianDate.addSeconds(start, periodSeconds, new Cesium.JulianDate());

        clock.startTime = start;
        clock.currentTime = Cesium.JulianDate.clone(start);
        clock.stopTime = stop;
        clock.clockRange = Cesium.ClockRange.LOOP_STOP;
        clock.multiplier = multiplier || Math.max(1, periodSeconds / 60); // Complete orbit in ~1 minute
        clock.shouldAnimate = true;

        // Update timeline bounds
        this.viewer.timeline.zoomTo(start, stop);

        // Ensure the animation widget reflects the new state
        if (this.viewer.clockViewModel) {
            this.viewer.clockViewModel.shouldAnimate = true;
        }
    }

    // Zoom camera to fit an orbit based on semi-major axis
    zoomToFitOrbit(semiMajorAxis, eccentricity = 0) {
        // Calculate apoapsis distance (furthest point from center)
        const apoapsis = semiMajorAxis * (1 + eccentricity);
        // Camera distance should be ~3.5x the apoapsis to see the whole orbit clearly
        const cameraDistance = apoapsis * 3.5 * 1000; // Convert km to meters

        // Clear any tracked entity
        this.viewer.trackedEntity = undefined;

        // Use flyTo for smooth animated camera movement
        // Position camera above and to the side for a good 3D view
        const direction = new Cesium.Cartesian3(-0.5, -0.3, -0.7);
        Cesium.Cartesian3.normalize(direction, direction);

        this.viewer.camera.flyTo({
            destination: new Cesium.Cartesian3(
                cameraDistance * 0.5,   // Offset in X
                cameraDistance * 0.3,   // Offset in Y
                cameraDistance * 0.7    // Above in Z
            ),
            orientation: {
                direction: direction,  // Look toward Earth
                up: new Cesium.Cartesian3(0, 0, 1)  // Z is up
            },
            duration: 1.0
        });
    }

    addEquatorialPlane() {
        const earthRadius = 6371000; // meters
        const diskRadius = earthRadius * 3; // Extend beyond Earth
        const numPoints = 120;

        // Draw equatorial plane as a ring outline
        const diskPositions = [];
        for (let i = 0; i <= numPoints; i++) {
            const rad = (i / numPoints) * 2 * Math.PI;
            diskPositions.push(new Cesium.Cartesian3(
                diskRadius * Math.cos(rad),
                diskRadius * Math.sin(rad),
                0
            ));
        }

        this.equatorialDisk = this.viewer.entities.add({
            polyline: {
                positions: diskPositions,
                width: 10,
                material: Cesium.Color.fromCssColorString('#00f0ff').withAlpha(0.15)
            }
        });
    }

    // ==================== WASM Loading with Web Worker ====================

    async loadWasmWorker() {
        const statusEl = document.getElementById('wasm-status');
        const dotEl = document.getElementById('wasm-dot');
        const runBtn = document.getElementById('run-btn');
        const self = this;

        try {
            this.log('Initializing Web Worker for WASM...', 'info');

            // Create worker
            this.testWorker = new Worker('testWorker.js');

            // Set up message handler
            this.testWorker.onmessage = (e) => this.handleWorkerMessage(e);
            this.testWorker.onerror = (e) => {
                console.error('Worker error:', e);
                this.log(`Worker error: ${e.message}`, 'error');
            };

            // Request WASM module load in worker
            this.testWorker.postMessage({
                type: 'load',
                wasmUrl: 'tudat_wasm_test.js'
            });

            // Wait for worker to be ready (with timeout)
            await this.waitForWorkerReady();

            statusEl.textContent = 'READY';
            dotEl.className = 'status-dot ready';
            runBtn.disabled = false;

            this.log('WASM Worker ready', 'pass');
            this.log('System ready. Click EXECUTE to run tests (runs in background).', 'info');

            // Also load main thread WASM for visualizations (store promise for later await)
            this.visualizationWasmPromise = this.loadVisualizationWasm();

            // Auto-select default visualization
            this.selectDefaultVisualization();
        } catch (error) {
            console.error('Worker setup error:', error);
            statusEl.textContent = 'LOAD FAILED';
            dotEl.className = 'status-dot error';
            this.log(`Worker setup failed: ${error.message}`, 'error');
            this.log('Running in demo mode...', 'info');

            runBtn.disabled = false;
            this.useDemoMode = true;
        }
    }

    waitForWorkerReady() {
        return new Promise((resolve, reject) => {
            const timeout = setTimeout(() => reject(new Error('Worker load timeout')), 60000);

            const checkReady = () => {
                if (this.workerReady) {
                    clearTimeout(timeout);
                    resolve();
                } else {
                    setTimeout(checkReady, 100);
                }
            };
            checkReady();
        });
    }

    handleWorkerMessage(e) {
        const { type, message, text, passed, name, current, passCount, failCount, total, count } = e.data;

        switch (type) {
            case 'status':
                this.log(message, 'info');
                break;

            case 'loaded':
                this.workerReady = true;
                this.log('WASM module loaded in worker', 'pass');
                break;

            case 'started':
                this.log('Test execution started in worker', 'info');
                break;

            case 'output':
                // Process output for console display
                this.processOutputDisplay(text);
                break;

            case 'result':
                // Update modal with test result
                this.updateModal(current, this.totalTestCount, passCount, failCount, name);
                // Also track in local state
                this.addTestResult(name, passed, this.currentCategory);
                break;

            case 'category':
                // Track current category
                if (name.includes('===')) {
                    const match = name.match(/===\s*(.+?)\s*Tests?\s*===/);
                    if (match) {
                        this.currentCategory = match[1].trim();
                    }
                }
                break;

            case 'total':
                // Update expected total
                this.totalTestCount = count;
                break;

            case 'finished':
                this.log(`Tests completed: ${passed} passed, ${failCount} failed`, 'info');
                // Update modal one final time
                this.updateModal(total, total, passed, failCount, 'Complete!');
                // Brief delay then hide modal and update UI
                setTimeout(() => {
                    this.hideModal();
                    this.finishTests();
                }, 1500);
                break;

            case 'error':
                this.log(`Worker error: ${message}`, 'error');
                this.hideModal();
                this.finishTests();
                break;
        }
    }

    // Load WASM on main thread for visualizations only (non-blocking)
    async loadVisualizationWasm() {
        try {
            this.log('Loading main-thread WASM module...', 'info');

            // Check if module already exists from previous load
            if (this.wasmModule && this.wasmModule.FS) {
                this.log('WASM module already loaded, reusing...', 'info');
                return;
            }

            // Load the script which defines createTudatModule
            await this.loadScript('tudat_wasm_test.js');
            this.log('WASM script loaded, initializing runtime...', 'info');

            // The module is built with MODULARIZE=1, so we call the factory function
            if (typeof createTudatModule !== 'function') {
                throw new Error('createTudatModule not found - module may not be built correctly');
            }

            // Initialize the module with our configuration
            this.wasmModule = await createTudatModule({
                print: function(text) {
                    // Silent - visualizations don't need console output
                },
                printErr: function(text) {
                    console.error('Viz WASM:', text);
                }
            });

            // Also store on window for backward compatibility
            window.Module = this.wasmModule;

            this.log('Main-thread WASM runtime ready', 'pass');
        } catch (error) {
            this.log(`Visualization WASM load failed: ${error.message}`, 'error');
            console.warn('Visualization WASM load failed:', error);
            // Re-throw so loadTudatModule can handle the failure
            throw error;
        }
    }

    // Display-only version of processOutput (no test result parsing)
    processOutputDisplay(text) {
        if (!text || typeof text !== 'string') return;
        this.log(text, this.classifyLine(text));
    }

    // Legacy loadWasm for fallback
    async loadWasm() {
        // Redirects to worker-based loading
        return this.loadWasmWorker();
    }

    loadScript(src) {
        return new Promise((resolve, reject) => {
            const script = document.createElement('script');
            script.src = src;
            script.onload = resolve;
            script.onerror = reject;
            document.head.appendChild(script);
        });
    }

    waitForModule() {
        return new Promise((resolve, reject) => {
            const maxAttempts = 200; // 10 seconds at 50ms intervals
            let attempts = 0;

            const check = () => {
                attempts++;

                // Log status every 20 attempts (1 second)
                if (attempts % 20 === 0) {
                    const moduleExists = typeof window.Module !== 'undefined';
                    const hasFS = moduleExists && window.Module.FS;
                    const hasMkdir = hasFS && typeof window.Module.FS.mkdir === 'function';
                    console.log(`waitForModule attempt ${attempts}: Module=${moduleExists}, FS=${!!hasFS}, mkdir=${hasMkdir}`);
                }

                // Check if module is fully ready - FS.mkdir being a function is a reliable indicator
                if (typeof window.Module !== 'undefined' && window.Module.FS && typeof window.Module.FS.mkdir === 'function') {
                    this.log(`WASM module ready after ${attempts} attempts`, 'info');
                    resolve();
                    return;
                }

                if (attempts >= maxAttempts) {
                    const moduleExists = typeof window.Module !== 'undefined';
                    const hasFS = moduleExists && window.Module.FS;
                    reject(new Error(`Timeout waiting for WASM module (Module=${moduleExists}, FS=${!!hasFS})`));
                    return;
                }

                // Keep polling
                setTimeout(check, 50);
            };
            check();
        });
    }

    /**
     * Load the full tudatpy WASM module with SPICE support
     * Reuses the tudat_wasm_test.js module which contains SPICE functions.
     */
    async loadTudatModule() {
        if (this.tudatModule) {
            return this.tudatModule;
        }

        try {
            this.log('Loading tudatpy WASM module...', 'info');
            this.log(`visualizationWasmPromise exists: ${!!this.visualizationWasmPromise}`, 'info');

            // Wait for the main-thread WASM module to be fully loaded
            // loadVisualizationWasm() loads tudat_wasm_test.js on the main thread
            if (this.visualizationWasmPromise) {
                this.log('Awaiting visualizationWasmPromise...', 'info');
                try {
                    await this.visualizationWasmPromise;
                    this.log('visualizationWasmPromise resolved', 'info');
                } catch (wasmError) {
                    this.log(`Visualization WASM failed to load: ${wasmError.message}`, 'warning');
                    // Continue - wasmModule may still have been set via onRuntimeInitialized callback
                }
            } else {
                this.log('No visualizationWasmPromise - loading WASM directly', 'warning');
                await this.loadVisualizationWasm();
            }

            this.log(`wasmModule exists: ${!!this.wasmModule}, has FS: ${!!(this.wasmModule && this.wasmModule.FS)}`, 'info');

            if (!this.wasmModule || !this.wasmModule.FS) {
                throw new Error('WASM module FS not available - main thread module not loaded');
            }

            // Reuse the already-loaded module which has SPICE support
            this.tudatModule = this.wasmModule;

            // Ensure SPICE kernel directory exists
            try {
                this.tudatModule.FS.mkdir('/spice_kernels');
            } catch (e) {
                // Directory may already exist
            }

            this.log('Tudatpy WASM module loaded', 'success');

            // Initialize the SPICE kernel loader
            this.spiceLoader = new SpiceKernelLoader(this.tudatModule);

            return this.tudatModule;
        } catch (error) {
            this.log(`Failed to load tudatpy module: ${error.message}`, 'error');
            console.error('Tudatpy module load error:', error);
            return null;
        }
    }

    /**
     * Load standard SPICE kernels for ephemeris queries
     * Call this before using SPICE-dependent visualizations
     */
    async loadSpiceKernels() {
        if (this.spiceReady) {
            return true;
        }

        // Ensure tudatpy module is loaded first
        if (!this.tudatModule) {
            await this.loadTudatModule();
        }

        if (!this.spiceLoader) {
            this.log('SPICE loader not available', 'error');
            return false;
        }

        try {
            this.log('Loading SPICE kernels...', 'info');

            const success = await this.spiceLoader.loadStandardKernels({
                onProgress: ({ loaded, total, currentKernel }) => {
                    this.log(`Loading kernel ${loaded}/${total}: ${currentKernel}`, 'info');
                }
            });

            if (success) {
                this.spiceReady = true;
                // Initialize the shared SPICE utils so visualizations can access SPICE
                initSpice(this.tudatModule, true);
                this.log(`SPICE kernels loaded (${this.spiceLoader.getLoadedCount()} kernels)`, 'success');

                // Test SPICE functionality
                this.spiceLoader.testSpice();
            } else {
                this.log('Some SPICE kernels failed to load', 'warning');
            }

            return success;
        } catch (error) {
            this.log(`SPICE kernel loading failed: ${error.message}`, 'error');
            return false;
        }
    }

    /**
     * Get planetary state using SPICE
     * @param {string} target - Target body name (e.g., 'Earth', 'Mars')
     * @param {string} observer - Observer body name (e.g., 'Sun')
     * @param {number} epoch - Ephemeris time (seconds since J2000)
     * @param {string} frame - Reference frame (default: 'J2000')
     * @returns {Object|null} State object {x, y, z, vx, vy, vz} or null if not available
     */
    getBodyState(target, observer, epoch, frame = 'J2000') {
        if (!this.spiceReady || !this.tudatModule) {
            console.warn('SPICE not ready. Call loadSpiceKernels() first.');
            return null;
        }

        try {
            const state = this.tudatModule.interface_spice_get_body_cartesian_state_at_epoch(
                target,
                observer,
                frame,
                'NONE',
                epoch
            );

            return {
                x: state.get(0),
                y: state.get(1),
                z: state.get(2),
                vx: state.get(3),
                vy: state.get(4),
                vz: state.get(5)
            };
        } catch (error) {
            console.error(`Failed to get state for ${target}:`, error);
            return null;
        }
    }

    // ==================== Output Processing ====================

    processOutput(text) {
        if (!text || typeof text !== 'string') return;

        // Log to console panel
        this.log(text, this.classifyLine(text));

        // Parse test results
        if (text.startsWith('[PASS]')) {
            const testName = text.substring(7).trim();
            this.addTestResult(testName, true);
        } else if (text.startsWith('[FAIL]')) {
            const testName = text.substring(7).trim();
            this.addTestResult(testName, false);
        } else if (text.startsWith('===') && text.endsWith('===')) {
            // Category header
            this.currentCategory = text.replace(/=/g, '').trim();
        } else if (text.includes('Tests run:')) {
            const match = text.match(/Tests run:\s*(\d+)/);
            if (match) {
                this.expectedTests = parseInt(match[1]);
            }
        }

        this.updateStats();
        this.updateProgress();
    }

    classifyLine(text) {
        if (text.startsWith('[PASS]')) return 'pass';
        if (text.startsWith('[FAIL]')) return 'fail';
        if (text.startsWith('[INFO]')) return 'info';
        if (text.startsWith('[ERROR]')) return 'error';
        if (text.startsWith('===')) return 'header';
        if (text.includes('Tests run') || text.includes('ALL TESTS')) return 'summary';
        return 'info';
    }

    addTestResult(name, passed) {
        const result = {
            name: name,
            passed: passed,
            category: this.currentCategory,
            timestamp: Date.now() - (this.startTime || Date.now())
        };

        this.testResults.push(result);

        // Group by category
        if (!this.categories[this.currentCategory]) {
            this.categories[this.currentCategory] = [];
        }
        this.categories[this.currentCategory].push(result);
    }

    log(text, type = 'info') {
        const console = document.getElementById('console-output');
        const line = document.createElement('div');
        line.className = `console-line ${type}`;
        line.textContent = text;
        console.appendChild(line);
        console.scrollTop = console.scrollHeight;

        this.consoleLines++;
        document.getElementById('line-count').textContent = `${this.consoleLines} lines`;
    }

    // ==================== UI Updates ====================

    updateStats() {
        const total = this.testResults.length;
        const passed = this.testResults.filter(r => r.passed).length;
        const failed = total - passed;
        const duration = this.startTime ? ((Date.now() - this.startTime) / 1000).toFixed(1) : '-';

        document.getElementById('total-count').textContent = total || '-';
        document.getElementById('passed-count').textContent = passed || '-';
        document.getElementById('failed-count').textContent = failed || '-';
        document.getElementById('duration').textContent = duration !== '-' ? `${duration}s` : '-';
    }

    updateProgress() {
        const progress = Math.min((this.testResults.length / this.expectedTests) * 100, 100);
        document.getElementById('progress-bar').style.width = `${progress}%`;
        document.getElementById('progress-text').textContent = `${Math.round(progress)}%`;
    }

    selectTest(testName, category) {
        this.selectedTest = testName;
        // Trigger visualization
        this.visualizeTest(testName, category);
    }

    // ==================== Test Execution ====================

    async runTests() {
        if (this.isRunning) return;

        this.isRunning = true;
        this.testResults = [];
        this.categories = {};
        this.startTime = Date.now();
        this.currentCategory = 'General';

        // Update UI
        const runBtn = document.getElementById('run-btn');
        runBtn.disabled = true;
        runBtn.innerHTML = '<span class="spinner"></span>RUNNING';

        document.getElementById('progress-section').style.display = 'block';
        document.getElementById('progress-bar').style.width = '0%';
        document.getElementById('console-output').innerHTML = '';
        this.consoleLines = 0;

        document.getElementById('wasm-status').textContent = 'EXECUTING';
        document.getElementById('wasm-dot').className = 'status-dot loading';

        this.log('Starting test execution...', 'header');
        this.log(`Timestamp: ${new Date().toISOString()}`, 'info');

        // Show progress modal
        this.showModal();

        try {
            if (this.useDemoMode) {
                await this.runDemoTests();
                this.hideModal();
                this.finishTests();
            } else {
                // Worker handles async completion via messages
                await this.runWasmTests();
                // finishTests() is called from handleWorkerMessage when 'finished' is received
            }
        } catch (error) {
            this.log(`Execution error: ${error.message}`, 'error');
            console.error(error);
            this.hideModal();
            this.finishTests();
        }
    }

    async runWasmTests() {
        if (this.testWorker && this.workerReady) {
            // Run tests in Web Worker
            this.log('Executing tests in Web Worker (UI remains responsive)...', 'info');
            this.testWorker.postMessage({ type: 'run' });
        } else {
            // Fallback to main thread execution
            this.log('Worker not available, running on main thread...', 'warning');
            if (this.wasmModule && this.wasmModule.callMain) {
                this.wasmModule.callMain([]);
            } else if (this.wasmModule && this.wasmModule._main) {
                this.wasmModule._main();
            } else {
                this.log('WASM entry point not found', 'error');
            }
            await this.sleep(1000);
            this.hideModal();
            this.finishTests();
        }
    }

    async runDemoTests() {
        // Demo test categories that match the actual WASM test structure
        // Each category shows meaningful orbital mechanics examples
        const categories = [
            // === Basic Astrodynamics ===
            { name: 'Unit Conversions', tests: [
                '180° → π radians',
                'π radians → 180°',
                '1 AU → 1.496e11 m'
            ]},
            { name: 'Physical Constants', tests: [
                'Speed of light c = 299,792,458 m/s',
                'Gravitational constant G = 6.67430e-11',
                'Earth GM μ = 3.986e14 m³/s²'
            ]},
            { name: 'Orbital Element Conversions (NASA ODTBX)', tests: [
                'Earth elliptical: a=8000km, e=0.23, i=20.6°',
                'Mars circular: a=9201km, e=0, i=0°',
                'Cartesian ↔ Keplerian round-trip (elliptical)',
                'Cartesian ↔ Keplerian round-trip (hyperbolic)'
            ]},
            { name: 'Kepler Orbital Mechanics', tests: [
                'Geostationary period T = 86,164s',
                'LEO mean motion n = 0.00114 rad/s',
                'Synodic period (Earth-Mars) ≈ 780 days'
            ]},

            // === Propagation ===
            { name: 'Two-Body Propagation', tests: [
                'LEO circular orbit (a=7000km, i=45°)',
                'Position error vs Kepler < 0.1 m',
                'Velocity error vs Kepler < 1e-4 m/s',
                'Full orbital period (5828s) propagation'
            ]},
            { name: 'CR3BP Propagation', tests: [
                'Sun-Earth-Moon system (μ=2.528e-5)',
                'Initial state: (0.994, 0.853, 0.312)',
                'Final state X component',
                'Final state Y component',
                'Propagation steps > 100'
            ]},
            { name: 'Mass Propagation', tests: [
                'Initial spacecraft mass: 500 kg',
                'Constant burn rate: dm/dt = -0.01 kg/s',
                'Final mass after 1000s: 490 kg'
            ]},

            // === SPICE Interface ===
            { name: 'SPICE Time Conversions', tests: [
                'J2000 epoch: JD 2451545.0 ↔ ET 0.0',
                'Round-trip: JD → ET → JD',
                'Pre-J2000: ET -86400 → JD 2451544.0'
            ]},
            { name: 'SPICE Frame Rotations', tests: [
                'J2000 → ECLIPJ2000 determinant = 1',
                'Rotation preserves X-axis (vernal equinox)',
                'Obliquity angle: 23.4393°'
            ]},
            { name: 'TLE/SGP4 (Vallado Benchmark)', tests: [
                'Vallado satellite 1958-002B',
                'Position error < 50m (vs reference)',
                'Velocity error < 0.05 m/s',
                'Parsed elements: i=34.27°, e=0.186',
                'ISS-like orbit test'
            ]},

            // === Gravitation ===
            { name: 'Libration Points (Earth-Moon CR3BP)', tests: [
                'L1 position: x = 0.8369 (normalized)',
                'L2 position: x = 1.1562',
                'L4/L5 triangular points (equilateral)',
                'Jacobi energy conservation'
            ]},
            { name: 'Third-Body Perturbation', tests: [
                'Inner perturber acceleration',
                'Outer perturber acceleration',
                'Perpendicular perturber acceleration'
            ]},
            { name: 'Spherical Harmonics Gravity', tests: [
                'EGM2008 C20 (J2 oblateness)',
                'EGM2008 C22/S22 sectoral terms',
                'Full field: degree 5, order 5'
            ]},

            // === Aerodynamics ===
            { name: 'Exponential Atmosphere', tests: [
                'Sea level: ρ=1.225 kg/m³, T=288.16K',
                'Scale height H = 7.2 km',
                'Density at 100 km (Kármán line)'
            ]},
            { name: 'NRLMSISE-00 Atmosphere', tests: [
                'Species densities at 400 km',
                'Exospheric temperature: 1250 K',
                'Density variation with solar activity'
            ]},

            // === Mission Design ===
            { name: 'Lambert Targeting (Izzo)', tests: [
                'Elliptical transfer: V_dep = (2736, 6594) m/s',
                'Arrival velocity: V_arr = (-1368, 4225) m/s',
                'Transfer arc is prograde',
                'Time of flight: 4034s'
            ]},

            // === Numerical Methods ===
            { name: 'RK4/RK78 Integrators', tests: [
                'Exponential growth: dy/dt = y',
                'Harmonic oscillator: x(2π) = x(0)',
                'Adaptive step size convergence'
            ]}
        ];

        for (const cat of categories) {
            this.processOutput(`\n=== ${cat.name} ===`);
            await this.sleep(30);

            for (const test of cat.tests) {
                const passed = Math.random() > 0.02;
                this.processOutput(`[${passed ? 'PASS' : 'FAIL'}] ${test}`);
                await this.sleep(15);
            }
        }

        this.processOutput('\n=== Test Results ===');
        const total = this.testResults.length;
        const passed = this.testResults.filter(r => r.passed).length;
        this.processOutput(`[INFO] Tests run: ${total}`);
        this.processOutput(`[INFO] Tests passed: ${passed}`);
        this.processOutput(`[INFO] Tests failed: ${total - passed}`);
        this.processOutput(passed === total ? '[PASS] *** ALL TESTS PASSED ***' : '[FAIL] *** SOME TESTS FAILED ***');
    }

    finishTests() {
        this.isRunning = false;

        const runBtn = document.getElementById('run-btn');
        runBtn.disabled = false;
        runBtn.innerHTML = 'EXECUTE';

        const passed = this.testResults.filter(r => r.passed).length;
        const failed = this.testResults.length - passed;

        document.getElementById('wasm-status').textContent = failed === 0 ? 'PASS' : 'FAIL';
        document.getElementById('wasm-dot').className = `status-dot ${failed === 0 ? 'ready' : 'error'}`;

        this.log(`\nTest execution complete: ${passed}/${this.testResults.length} passed`, 'summary');
    }

    clearResults() {
        this.testResults = [];
        this.categories = {};
        this.consoleLines = 0;
        this.selectedTest = null;

        document.getElementById('total-count').textContent = '-';
        document.getElementById('passed-count').textContent = '-';
        document.getElementById('failed-count').textContent = '-';
        document.getElementById('duration').textContent = '-';
        document.getElementById('progress-bar').style.width = '0%';
        document.getElementById('progress-text').textContent = '0%';
        document.getElementById('progress-section').style.display = 'none';
        document.getElementById('line-count').textContent = '0 lines';
        document.getElementById('orbit-info').textContent = 'Click a test to visualize';

        document.getElementById('console-output').innerHTML = '<div class="console-line info">Console cleared. Ready for new test run.</div>';

        // Reset visualization category selection
        document.querySelectorAll('.viz-category').forEach(el => el.classList.remove('selected'));

        // Reset bottom panel to default (orbit selector visible, others hidden)
        const orbitSelectorPanel = document.getElementById('orbit-selector-panel');
        const odModelPanel = document.getElementById('od-model-panel');
        const chartPanel = document.getElementById('chart-panel');
        if (orbitSelectorPanel) orbitSelectorPanel.style.display = '';
        if (odModelPanel) odModelPanel.style.display = 'none';
        if (chartPanel) chartPanel.style.display = 'none';

        document.getElementById('wasm-status').textContent = 'READY';
        document.getElementById('wasm-dot').className = 'status-dot ready';

        // Clear Cesium entities
        this.clearOrbitEntities();
        this.resetCharts();
    }

    // ==================== Visualization System ====================

    visualizeTest(testName, category) {
        // Update URL with current visualization
        this.updateUrlWithVisualization(category);

        // Update info display
        document.getElementById('orbit-info').textContent = `${category}: ${testName}`;

        // Clear previous 3D entities
        this.clearOrbitEntities();

        // Check if this is a chart-only visualization (Python example port)
        const vizConfig = visualizationRegistry[category];
        const isChartOnly = vizConfig?.chartOnly === true;

        // Determine which bottom panel to show
        const cat = category.toLowerCase();
        const orbitSelectorPanel = document.getElementById('orbit-selector-panel');
        const odModelPanel = document.getElementById('od-model-panel');
        const chartPanel = document.getElementById('chart-panel');

        if (cat.includes('cr3bp')) {
            // Show orbit selector for CR3BP
            if (orbitSelectorPanel) orbitSelectorPanel.style.display = '';
            if (odModelPanel) odModelPanel.style.display = 'none';
            if (chartPanel) chartPanel.style.display = 'none';
        } else if (cat.includes('orbit determination') || cat.includes('differential correction')) {
            // Show OD model toggle for Orbit Determination
            if (orbitSelectorPanel) orbitSelectorPanel.style.display = 'none';
            if (odModelPanel) odModelPanel.style.display = '';
            if (chartPanel) chartPanel.style.display = 'none';
        } else {
            // Show chart panel for other visualizations
            if (orbitSelectorPanel) orbitSelectorPanel.style.display = 'none';
            if (odModelPanel) odModelPanel.style.display = 'none';
            if (chartPanel) chartPanel.style.display = '';
        }

        // For chart-only visualizations, render charts in the full view area
        if (isChartOnly) {
            this.showChartOnlyVisualization(category, testName);
            return;
        }

        // Restore globe layout if switching from chart-only
        this.restoreGlobeLayout();

        // Show 3D visualization on globe
        this.show3DVisualization(category, testName);

        // Show 2D chart below (for non-CR3BP and non-OD)
        this.showChartForCategory(category, testName);
    }

    show3DVisualization(category, testName) {
        const cat = category.toLowerCase();
        const log = (msg, level) => this.log(msg, level);

        // Use imported visualization modules based on category
        if (cat.includes('cr3bp') || cat.includes('three-body') || cat.includes('3bp')) {
            showCR3BPVisualization(this.viewer, this.orbitEntities, 'l2-halo', log, true);
        }
        else if (cat.includes('libration')) {
            showLibrationPointsVisualization(this.viewer, this.orbitEntities, log);
        }
        else if (cat.includes('atmospheric') || cat.includes('nrlmsise')) {
            showAtmosphericDragVisualization(this.viewer, this.orbitEntities);
        }
        else if (cat.includes('reference frame') || (cat.includes('spice') && cat.includes('frame'))) {
            showReferenceFramesVisualization(this.viewer, this.orbitEntities);
        }
        else if (cat.includes('geostationary')) {
            showGeostationaryVisualization(this.viewer, this.orbitEntities);
        }
        else if (cat.includes('j2 vs full force') || cat.includes('integrat') || cat.includes('rk78')) {
            // Use class method for J2 vs Full Force to get the separation chart
            const period = 5800;
            const numOrbits = 30;
            this.addIntegratorComparisonVisualization(period, numOrbits);
            this.configureClockForOrbit(period * numOrbits, null, period / 10);
            this.viewer.camera.flyTo({
                destination: Cesium.Cartesian3.fromDegrees(0, 0, 25000000),
                orientation: {
                    heading: 0,
                    pitch: Cesium.Math.toRadians(-90),
                    roll: 0
                },
                duration: 1.0
            });
        }
        else if (cat.includes('omm vs j2') || cat.includes('omm vs two')) {
            // Use class method for OMM vs Two-Body to get the separation chart
            const period = 5400;  // ISS-like ~90 min orbit
            const numOrbits = 20;
            this.addOMMvsJ2Visualization(period, numOrbits);
            this.configureClockForOrbit(period * numOrbits, null, period / 10);
            this.viewer.camera.flyTo({
                destination: Cesium.Cartesian3.fromDegrees(0, 0, 25000000),
                orientation: {
                    heading: 0,
                    pitch: Cesium.Math.toRadians(-90),
                    roll: 0
                },
                duration: 1.0
            });
        }
        else if (cat.includes('orbit determination') || cat.includes('differential correction')) {
            // Orbit determination visualization with toggle
            // Default to full force, can be toggled via UI
            this.currentODModel = this.currentODModel || 'fullforce';
            this.addOrbitDeterminationVisualization(this.currentODModel);
        }
        // Default: show a simple circular orbit using shared utility
        else {
            clearOrbitEntities(this.viewer, this.orbitEntities);
            const period = 5400;
            configureClockForOrbit(this.viewer, period, null, period / 30);
            addAnimatedOrbit(this.viewer, this.orbitEntities, {
                name: 'Test Satellite',
                semiMajorAxis: 6800,
                eccentricity: 0.0,
                inclination: 28.5,
                raan: 0,
                argPeriapsis: 0,
                color: '#00f0ff',
                period: period,
                description: category + '\n' + testName
            });
        }
    }

    // Add an animated orbit with time-varying satellite position
    addAnimatedOrbit(params) {
        const {
            name,
            semiMajorAxis,
            eccentricity,
            inclination,
            raan = 0,
            argPeriapsis = 0,
            trueAnomaly = 0,
            color,
            period,
            description,
            referenceFrame = 'INERTIAL'  // 'INERTIAL' (ECI) or 'FIXED' (ECEF for GEO)
        } = params;

        // For geostationary orbits, use Earth-fixed frame so satellite appears stationary
        const useFixedFrame = referenceFrame === 'FIXED' ||
            (Math.abs(inclination) < 0.1 && eccentricity < 0.01 && Math.abs(period - 86164) < 100);

        const a = semiMajorAxis; // km
        const e = eccentricity;
        const i = inclination * Cesium.Math.toRadians(1);
        const omega = argPeriapsis * Cesium.Math.toRadians(1);
        const Omega = raan * Cesium.Math.toRadians(1);

        // Generate orbit path positions
        const orbitPositions = [];
        const p = a * (1 - e * e);

        for (let nu = 0; nu <= 360; nu += 2) {
            const nuRad = nu * Cesium.Math.toRadians(1);
            const r = p / (1 + e * Math.cos(nuRad));

            // Position in perifocal frame (PQW)
            const xPQW = r * Math.cos(nuRad);
            const yPQW = r * Math.sin(nuRad);

            // Rotation matrices: PQW -> ECI (via omega, i, Omega)
            const cosO = Math.cos(Omega), sinO = Math.sin(Omega);
            const cosi = Math.cos(i), sini = Math.sin(i);
            const cosw = Math.cos(omega), sinw = Math.sin(omega);

            // Full rotation from perifocal to ECI
            const x = (cosO * cosw - sinO * sinw * cosi) * xPQW + (-cosO * sinw - sinO * cosw * cosi) * yPQW;
            const y = (sinO * cosw + cosO * sinw * cosi) * xPQW + (-sinO * sinw + cosO * cosw * cosi) * yPQW;
            const z = (sinw * sini) * xPQW + (cosw * sini) * yPQW;

            // Convert ECI (km) to Cesium Cartesian3
            orbitPositions.push(new Cesium.Cartesian3(x * 1000, y * 1000, z * 1000));
        }

        // Add orbit path as polyline
        // For geostationary (FIXED frame), orbit path is stationary relative to Earth
        // For inertial orbits, the path appears to rotate with the Earth underneath
        if (useFixedFrame) {
            // For GEO, just show a ring at GEO altitude (the orbit IS the ring in ECEF)
            const geoPositions = [];
            for (let lon = 0; lon <= 360; lon += 2) {
                geoPositions.push(Cesium.Cartesian3.fromDegrees(lon, 0, (semiMajorAxis - 6371) * 1000));
            }
            const orbitEntity = this.viewer.entities.add({
                name: name + ' Orbit',
                polyline: {
                    positions: geoPositions,
                    width: 10,
                    material: new Cesium.PolylineGlowMaterialProperty({
                        glowPower: 0.3,
                        color: Cesium.Color.fromCssColorString(color).withAlpha(0.9)
                    })
                }
            });
            this.orbitEntities.push(orbitEntity);
        }
        // No static orbit line for non-GEO orbits - the satellite trail shows the path

        // Zoom camera to fit the orbit
        this.zoomToFitOrbit(semiMajorAxis, eccentricity);

        // Create time-varying position
        const clock = this.viewer.clock;
        const startTime = clock.startTime;
        let positionProperty;

        // Sample positions for multiple orbit periods so satellite keeps looping
        const numSamplesPerPeriod = 360;
        const numPeriods = 10;  // Sample 10 full periods

        if (useFixedFrame) {
            // For geostationary: use ConstantPositionProperty in FIXED (ECEF) frame
            // The satellite stays at a fixed longitude/latitude/altitude
            const geoLongitude = raan; // Use RAAN as geographic longitude
            const geoAltitude = (semiMajorAxis - 6371) * 1000; // meters above surface
            const fixedPosition = Cesium.Cartesian3.fromDegrees(geoLongitude, 0, geoAltitude);

            // Use ConstantPositionProperty - satellite doesn't move in ECEF
            positionProperty = new Cesium.ConstantPositionProperty(fixedPosition, Cesium.ReferenceFrame.FIXED);
        } else {
            // Standard inertial orbit - use SampledPositionProperty in INERTIAL frame
            positionProperty = new Cesium.SampledPositionProperty(Cesium.ReferenceFrame.INERTIAL);

            for (let orbitNum = 0; orbitNum < numPeriods; orbitNum++) {
                for (let sample = 0; sample <= numSamplesPerPeriod; sample++) {
                    const fraction = sample / numSamplesPerPeriod;
                    const totalTime = (orbitNum + fraction) * period;
                    const time = Cesium.JulianDate.addSeconds(startTime, totalTime, new Cesium.JulianDate());

                    // Calculate mean anomaly at this time
                    const M = (2 * Math.PI * fraction) + (trueAnomaly * Cesium.Math.toRadians(1));

                    // Solve Kepler's equation for eccentric anomaly
                    let E = M;
                    for (let iter = 0; iter < 10; iter++) {
                        E = M + e * Math.sin(E);
                    }

                    // Calculate true anomaly
                    const nu = 2 * Math.atan2(
                        Math.sqrt(1 + e) * Math.sin(E / 2),
                        Math.sqrt(1 - e) * Math.cos(E / 2)
                    );

                    const r = p / (1 + e * Math.cos(nu));

                    // Position in perifocal frame
                    const xPQW = r * Math.cos(nu);
                    const yPQW = r * Math.sin(nu);

                    // Rotation to ECI
                    const cosO = Math.cos(Omega), sinO = Math.sin(Omega);
                    const cosi = Math.cos(i), sini = Math.sin(i);
                    const cosw = Math.cos(omega), sinw = Math.sin(omega);

                    const x = (cosO * cosw - sinO * sinw * cosi) * xPQW + (-cosO * sinw - sinO * cosw * cosi) * yPQW;
                    const y = (sinO * cosw + cosO * sinw * cosi) * xPQW + (-sinO * sinw + cosO * cosw * cosi) * yPQW;
                    const zPos = (sinw * sini) * xPQW + (cosw * sini) * yPQW;

                    positionProperty.addSample(time, new Cesium.Cartesian3(x * 1000, y * 1000, zPos * 1000));
                }
            }
        }

        // Create orientation property
        // For moving satellites use VVLH (velocity forward, nadir down)
        // For stationary GEO, use fixed nadir-pointing orientation
        let orientationProperty;
        if (useFixedFrame) {
            // For GEO, no velocity - just point nadir (no orientation needed for point)
            orientationProperty = undefined;
        } else {
            orientationProperty = new Cesium.VelocityOrientationProperty(positionProperty);
        }

        // Add animated satellite entity
        const satEntity = this.viewer.entities.add({
            name: name,
            description: description,
            position: positionProperty,
            orientation: orientationProperty,
            point: {
                pixelSize: 12,
                color: Cesium.Color.fromCssColorString(color),
                outlineColor: Cesium.Color.WHITE,
                outlineWidth: 2
            },
            path: useFixedFrame ? undefined : {
                // Full orbit trail for non-GEO satellites
                show: true,
                leadTime: 0,
                trailTime: period,  // Full orbit period
                width: 3,
                material: Cesium.Color.fromCssColorString(color)
            },
            label: {
                text: name,
                font: '12px monospace',
                style: Cesium.LabelStyle.FILL_AND_OUTLINE,
                outlineWidth: 2,
                verticalOrigin: Cesium.VerticalOrigin.BOTTOM,
                pixelOffset: new Cesium.Cartesian2(0, -15),
                fillColor: Cesium.Color.fromCssColorString(color),
                outlineColor: Cesium.Color.BLACK
            },
            // VVLH view: behind and above the satellite in its local frame
            // viewFrom is in entity's local coords when orientation is set
            viewFrom: new Cesium.Cartesian3(-50000, 0, 20000)  // Behind (-X), above (+Z) in VVLH
        });
        this.orbitEntities.push(satEntity);
    }

    // Legacy method - kept for compatibility
    addOrbitToGlobe(semiMajorAxis, eccentricity, inclination, color) {
        this.addAnimatedOrbit({
            name: 'Satellite',
            semiMajorAxis,
            eccentricity,
            inclination,
            color,
            period: 2 * Math.PI * Math.sqrt(Math.pow(semiMajorAxis, 3) / 398600.4418),
            description: ''
        });
    }

    // Add Moon's orbit for cislunar visualizations
    addMoonOrbit() {
        const moonDistance = 384400; // km
        const moonInclination = 5.145; // degrees to ecliptic

        // Simplified Moon orbit (circular approximation)
        const positions = [];
        for (let nu = 0; nu <= 360; nu += 5) {
            const nuRad = nu * Cesium.Math.toRadians(1);
            const iRad = moonInclination * Cesium.Math.toRadians(1);

            const x = moonDistance * Math.cos(nuRad);
            const y = moonDistance * Math.sin(nuRad) * Math.cos(iRad);
            const z = moonDistance * Math.sin(nuRad) * Math.sin(iRad);

            positions.push(new Cesium.Cartesian3(x * 1000, y * 1000, z * 1000));
        }

        const moonOrbitEntity = this.viewer.entities.add({
            name: 'Moon Orbit',
            polyline: {
                positions: positions,
                width: 2,
                material: Cesium.Color.GRAY.withAlpha(0.4)
            }
        });
        this.orbitEntities.push(moonOrbitEntity);

        // Add Moon marker (static position for simplicity)
        this.moonEntity = this.viewer.entities.add({
            name: 'Moon',
            position: new Cesium.Cartesian3(moonDistance * 1000, 0, 0),
            point: {
                pixelSize: 15,
                color: Cesium.Color.LIGHTGRAY,
                outlineColor: Cesium.Color.WHITE,
                outlineWidth: 2
            },
            label: {
                text: 'Moon',
                font: '12px monospace',
                fillColor: Cesium.Color.WHITE,
                pixelOffset: new Cesium.Cartesian2(0, -18)
            }
        });
        this.orbitEntities.push(this.moonEntity);
    }

    // Add CR3BP trajectory visualization using actual Tudat 3-body propagation
    // Uses Module.propagateCR3BP to compute periodic orbits in the Earth-Moon system
    // NO FALLBACKS - this is a WASM demonstration
    addCR3BPVisualization(config) {
        // Use provided config or default to L2 Halo
        if (!config) {
            config = this.getCR3BPOrbitConfig('l2-halo');
        }

        const { mu, L, TU, x0, y0, z0, vx0, vy0, vz0, period, name, description, color } = config;

        const duration = period * 1.0;  // One complete orbit
        const numPoints = 500;
        const earthX = -mu;  // Earth position in normalized coords

        // Require Tudat WASM bindings - NO FALLBACK
        if (typeof Module === 'undefined' || typeof Module.propagateCR3BP !== 'function') {
            this.log('ERROR: Tudat WASM bindings not available for CR3BP', 'error');
            return;
        }

        this.log(`Computing ${name} trajectory with Tudat...`, 'info');
        const trajectory = Module.propagateCR3BP(mu, x0, y0, z0, vx0, vy0, vz0, duration, numPoints);
        this.log(`Got ${trajectory.length / 7} trajectory points`, 'info');

        // Convert normalized CR3BP coordinates to Cesium positions
        // In CR3BP: Earth at (-μ, 0, 0), Moon at (1-μ, 0, 0)
        // We'll transform to Earth-centered coordinates
        const positions = [];
        const times = [];

        for (let i = 0; i < numPoints && i * 7 < trajectory.length; i++) {
            const idx = i * 7;
            const t = trajectory[idx];
            const x = trajectory[idx + 1];
            const y = trajectory[idx + 2];
            const z = trajectory[idx + 3];

            times.push(t * TU);

            // Convert from CR3BP rotating frame to Earth-centered (meters)
            const xEarth = (x - earthX) * L * 1000;
            const yEarth = y * L * 1000;
            const zEarth = z * L * 1000;

            positions.push(new Cesium.Cartesian3(xEarth, yEarth, zEarth));
        }

        // Add the trajectory polyline
        if (positions.length > 0) {
            const orbitColor = Cesium.Color.fromCssColorString(color);

            const trajectoryEntity = this.viewer.entities.add({
                name: name,
                polyline: {
                    positions: positions,
                    width: 10,
                    material: new Cesium.PolylineGlowMaterialProperty({
                        glowPower: 0.2,
                        color: orbitColor
                    })
                },
                description: description
            });
            this.orbitEntities.push(trajectoryEntity);

            // Add animated spacecraft
            const startTime = this.viewer.clock.startTime;
            const stopTime = this.viewer.clock.stopTime;
            const viewerDuration = Cesium.JulianDate.secondsDifference(stopTime, startTime);
            const trajectoryDuration = times[times.length - 1] - times[0];
            const timeScale = trajectoryDuration / viewerDuration;

            const property = new Cesium.SampledPositionProperty();

            for (let i = 0; i < positions.length; i++) {
                const trajTime = times[i] - times[0];
                const viewerTime = trajTime / timeScale;
                const time = Cesium.JulianDate.addSeconds(startTime, viewerTime, new Cesium.JulianDate());
                property.addSample(time, positions[i]);
            }

            const spacecraftEntity = this.viewer.entities.add({
                name: `${name} Spacecraft`,
                position: property,
                point: {
                    pixelSize: 12,
                    color: orbitColor
                },
                path: {
                    show: false
                },
                label: {
                    text: name,
                    font: '11px monospace',
                    fillColor: orbitColor,
                    pixelOffset: new Cesium.Cartesian2(12, 0),
                    showBackground: true,
                    backgroundColor: Cesium.Color.BLACK.withAlpha(0.7)
                }
            });
            this.orbitEntities.push(spacecraftEntity);

            this.log(`${name}: ${(trajectoryDuration / 86400).toFixed(1)} days period`, 'info');
        }
    }

    // Add reference frames visualization (J2000 equatorial vs ECLIPJ2000 ecliptic)
    addReferenceFrames() {
        const axisLength = 15000000; // 15,000 km in meters
        const earthRadius = 6371000; // Start axes just outside Earth surface
        const planeRadius = 18000000; // 18,000 km for plane circles
        const obliquity = 23.4393 * Math.PI / 180; // Earth's axial tilt in radians

        // ===== J2000 EQUATORIAL FRAME (Cyan) =====
        // X-axis - Vernal equinox direction (shared by both frames)
        const xAxisEq = this.viewer.entities.add({
            name: 'J2000 X (Vernal Equinox)',
            polyline: {
                positions: [new Cesium.Cartesian3(earthRadius, 0, 0), new Cesium.Cartesian3(axisLength, 0, 0)],
                width: 5,
                material: Cesium.Color.CYAN
            }
        });
        this.orbitEntities.push(xAxisEq);

        // Y-axis - 90° from X in equatorial plane
        const yAxisEq = this.viewer.entities.add({
            name: 'J2000 Y (Equatorial)',
            polyline: {
                positions: [new Cesium.Cartesian3(0, earthRadius, 0), new Cesium.Cartesian3(0, axisLength, 0)],
                width: 4,
                material: Cesium.Color.CYAN.withAlpha(0.7)
            }
        });
        this.orbitEntities.push(yAxisEq);

        // Z-axis - North celestial pole (perpendicular to equator)
        const zAxisEq = this.viewer.entities.add({
            name: 'J2000 Z (North Pole)',
            polyline: {
                positions: [new Cesium.Cartesian3(0, 0, earthRadius), new Cesium.Cartesian3(0, 0, axisLength)],
                width: 4,
                material: Cesium.Color.CYAN.withAlpha(0.7)
            }
        });
        this.orbitEntities.push(zAxisEq);

        // Equatorial plane circle
        const eqPositions = [];
        for (let a = 0; a <= 360; a += 3) {
            const rad = a * Math.PI / 180;
            eqPositions.push(new Cesium.Cartesian3(
                planeRadius * Math.cos(rad),
                planeRadius * Math.sin(rad),
                0
            ));
        }
        const eqPlane = this.viewer.entities.add({
            name: 'Equatorial Plane',
            polyline: {
                positions: eqPositions,
                width: 2,
                material: Cesium.Color.CYAN.withAlpha(0.5)
            }
        });
        this.orbitEntities.push(eqPlane);

        // ===== ECLIPJ2000 ECLIPTIC FRAME (Orange) =====
        // The ecliptic frame is rotated by obliquity (23.4°) around the X-axis
        // X-axis is the same (vernal equinox)
        // Y and Z are rotated

        // Y-axis ecliptic (rotated by obliquity around X)
        const yEclEnd = new Cesium.Cartesian3(
            0,
            axisLength * Math.cos(obliquity),
            axisLength * Math.sin(obliquity)
        );
        const yAxisEcl = this.viewer.entities.add({
            name: 'ECLIPJ2000 Y',
            polyline: {
                positions: [
                    new Cesium.Cartesian3(0, earthRadius * Math.cos(obliquity), earthRadius * Math.sin(obliquity)),
                    yEclEnd
                ],
                width: 4,
                material: Cesium.Color.ORANGE.withAlpha(0.7)
            }
        });
        this.orbitEntities.push(yAxisEcl);

        // Z-axis ecliptic (north ecliptic pole - perpendicular to ecliptic plane)
        const zEclEnd = new Cesium.Cartesian3(
            0,
            -axisLength * Math.sin(obliquity),
            axisLength * Math.cos(obliquity)
        );
        const zAxisEcl = this.viewer.entities.add({
            name: 'ECLIPJ2000 Z (Ecliptic Pole)',
            polyline: {
                positions: [
                    new Cesium.Cartesian3(0, -earthRadius * Math.sin(obliquity), earthRadius * Math.cos(obliquity)),
                    zEclEnd
                ],
                width: 4,
                material: Cesium.Color.ORANGE.withAlpha(0.7)
            }
        });
        this.orbitEntities.push(zAxisEcl);

        // Ecliptic plane circle (tilted by obliquity)
        const eclPositions = [];
        for (let a = 0; a <= 360; a += 3) {
            const rad = a * Math.PI / 180;
            const x = planeRadius * Math.cos(rad);
            const yTemp = planeRadius * Math.sin(rad);
            // Rotate around X-axis by obliquity
            const y = yTemp * Math.cos(obliquity);
            const z = yTemp * Math.sin(obliquity);
            eclPositions.push(new Cesium.Cartesian3(x, y, z));
        }
        const eclPlane = this.viewer.entities.add({
            name: 'Ecliptic Plane',
            polyline: {
                positions: eclPositions,
                width: 2,
                material: Cesium.Color.ORANGE.withAlpha(0.5)
            }
        });
        this.orbitEntities.push(eclPlane);

        // ===== LABELS =====
        // J2000 labels
        const j2000Label = this.viewer.entities.add({
            position: new Cesium.Cartesian3(0, 0, axisLength * 1.1),
            label: {
                text: 'J2000 Z\n(Celestial Pole)',
                font: '12px monospace',
                fillColor: Cesium.Color.CYAN,
                heightReference: Cesium.HeightReference.NONE
            }
        });
        this.orbitEntities.push(j2000Label);

        // Ecliptic pole label
        const eclPoleLabel = this.viewer.entities.add({
            position: zEclEnd,
            label: {
                text: 'ECLIP Z\n(Ecliptic Pole)',
                font: '12px monospace',
                fillColor: Cesium.Color.ORANGE,
                heightReference: Cesium.HeightReference.NONE,
                pixelOffset: new Cesium.Cartesian2(10, 0)
            }
        });
        this.orbitEntities.push(eclPoleLabel);

        // Vernal equinox label (shared X-axis)
        const vernalLabel = this.viewer.entities.add({
            position: new Cesium.Cartesian3(axisLength * 1.05, 0, 0),
            label: {
                text: 'X (Vernal Equinox)\n♈ First Point of Aries',
                font: '12px monospace',
                fillColor: Cesium.Color.WHITE,
                heightReference: Cesium.HeightReference.NONE
            }
        });
        this.orbitEntities.push(vernalLabel);

        // Obliquity angle indicator
        const obliqLabel = this.viewer.entities.add({
            position: new Cesium.Cartesian3(0, axisLength * 0.6, axisLength * 0.3),
            label: {
                text: '← 23.4° obliquity →',
                font: '11px monospace',
                fillColor: Cesium.Color.YELLOW,
                heightReference: Cesium.HeightReference.NONE
            }
        });
        this.orbitEntities.push(obliqLabel);

        // Plane labels
        const eqPlaneLabel = this.viewer.entities.add({
            position: new Cesium.Cartesian3(planeRadius * 0.7, planeRadius * 0.7, 0),
            label: {
                text: 'Equatorial Plane',
                font: '11px monospace',
                fillColor: Cesium.Color.CYAN,
                heightReference: Cesium.HeightReference.NONE
            }
        });
        this.orbitEntities.push(eqPlaneLabel);

        const eclPlaneLabel = this.viewer.entities.add({
            position: new Cesium.Cartesian3(planeRadius * 0.7, planeRadius * 0.5, planeRadius * 0.35),
            label: {
                text: 'Ecliptic Plane',
                font: '11px monospace',
                fillColor: Cesium.Color.ORANGE,
                heightReference: Cesium.HeightReference.NONE
            }
        });
        this.orbitEntities.push(eclPlaneLabel);
    }

    // Add equatorial and ecliptic plane visualization
    addEquatorialAndEclipticPlanes() {
        const planeRadius = 20000000; // 20,000 km
        const obliquity = 23.4393 * Cesium.Math.toRadians(1); // Earth's axial tilt

        // Equatorial plane circle
        const eqPositions = [];
        for (let a = 0; a <= 360; a += 5) {
            const rad = a * Cesium.Math.toRadians(1);
            eqPositions.push(new Cesium.Cartesian3(
                planeRadius * Math.cos(rad),
                planeRadius * Math.sin(rad),
                0
            ));
        }
        const eqPlane = this.viewer.entities.add({
            name: 'Equatorial Plane',
            polyline: {
                positions: eqPositions,
                width: 2,
                material: Cesium.Color.CYAN.withAlpha(0.6)
            }
        });
        this.orbitEntities.push(eqPlane);

        // Ecliptic plane circle (tilted by obliquity)
        const eclPositions = [];
        for (let a = 0; a <= 360; a += 5) {
            const rad = a * Cesium.Math.toRadians(1);
            const x = planeRadius * Math.cos(rad);
            const yTemp = planeRadius * Math.sin(rad);
            const y = yTemp * Math.cos(obliquity);
            const z = yTemp * Math.sin(obliquity);
            eclPositions.push(new Cesium.Cartesian3(x, y, z));
        }
        const eclPlane = this.viewer.entities.add({
            name: 'Ecliptic Plane',
            polyline: {
                positions: eclPositions,
                width: 2,
                material: Cesium.Color.ORANGE.withAlpha(0.6)
            }
        });
        this.orbitEntities.push(eclPlane);

        // Labels (with heightReference: NONE to avoid terrain sampling)
        const eqLabel = this.viewer.entities.add({
            position: new Cesium.Cartesian3(planeRadius * 0.7, planeRadius * 0.7, 0),
            label: {
                text: 'Equatorial',
                font: '10px monospace',
                fillColor: Cesium.Color.CYAN,
                heightReference: Cesium.HeightReference.NONE
            }
        });
        this.orbitEntities.push(eqLabel);

        const eclLabel = this.viewer.entities.add({
            position: new Cesium.Cartesian3(planeRadius * 0.7, planeRadius * 0.5, planeRadius * 0.3),
            label: {
                text: 'Ecliptic (23.4°)',
                font: '10px monospace',
                fillColor: Cesium.Color.ORANGE,
                heightReference: Cesium.HeightReference.NONE
            }
        });
        this.orbitEntities.push(eclLabel);
    }

    // Add Earth-Moon Lagrange points L1-L5
    addLibrationPointsVisualization() {
        const moonDist = 384400; // km
        const mu = 0.01215; // Earth-Moon mass ratio

        // Check if Tudat bindings are available
        const hasTudatBindings = typeof Module !== 'undefined' &&
                                  typeof Module.computeLibrationPoints === 'function';

        let lagrangePoints;

        if (hasTudatBindings) {
            // Use actual Tudat library for libration point computation
            this.log('Using Tudat library for libration points', 'info');

            try {
                const lpData = Module.computeLibrationPoints(mu);
                // lpData contains [L1x,L1y,L1z, L2x,L2y,L2z, L3x,L3y,L3z, L4x,L4y,L4z, L5x,L5y,L5z]
                // Positions are in normalized units (multiply by moonDist to get km)

                lagrangePoints = [
                    { name: 'L1', pos: new Cesium.Cartesian3(lpData[0] * moonDist * 1000, lpData[1] * moonDist * 1000, lpData[2] * moonDist * 1000), color: Cesium.Color.RED },
                    { name: 'L2', pos: new Cesium.Cartesian3(lpData[3] * moonDist * 1000, lpData[4] * moonDist * 1000, lpData[5] * moonDist * 1000), color: Cesium.Color.RED },
                    { name: 'L3', pos: new Cesium.Cartesian3(lpData[6] * moonDist * 1000, lpData[7] * moonDist * 1000, lpData[8] * moonDist * 1000), color: Cesium.Color.RED },
                    { name: 'L4', pos: new Cesium.Cartesian3(lpData[9] * moonDist * 1000, lpData[10] * moonDist * 1000, lpData[11] * moonDist * 1000), color: Cesium.Color.LIME },
                    { name: 'L5', pos: new Cesium.Cartesian3(lpData[12] * moonDist * 1000, lpData[13] * moonDist * 1000, lpData[14] * moonDist * 1000), color: Cesium.Color.LIME }
                ];
            } catch (e) {
                this.log('Tudat libration points failed: ' + e.message + ', using JS fallback', 'warning');
                lagrangePoints = null;
            }
        }

        // Fallback to JavaScript approximations
        if (!lagrangePoints) {
            // Approximate positions (JS fallback)
            const L1x = moonDist * 0.8369;
            const L2x = moonDist * 1.1562;
            const L3x = -moonDist * 1.0051;
            const L4x = moonDist * 0.5 - moonDist * mu;
            const L4y = moonDist * Math.sqrt(3) / 2;
            const L5y = -L4y;

            lagrangePoints = [
                { name: 'L1', pos: new Cesium.Cartesian3(L1x * 1000, 0, 0), color: Cesium.Color.RED },
                { name: 'L2', pos: new Cesium.Cartesian3(L2x * 1000, 0, 0), color: Cesium.Color.RED },
                { name: 'L3', pos: new Cesium.Cartesian3(L3x * 1000, 0, 0), color: Cesium.Color.RED },
                { name: 'L4', pos: new Cesium.Cartesian3(L4x * 1000, L4y * 1000, 0), color: Cesium.Color.LIME },
                { name: 'L5', pos: new Cesium.Cartesian3(L4x * 1000, L5y * 1000, 0), color: Cesium.Color.LIME }
            ];
        }

        const dataSource = hasTudatBindings ? 'Tudat C++' : 'JS approx';

        lagrangePoints.forEach(lp => {
            const entity = this.viewer.entities.add({
                name: lp.name,
                description: `${lp.name} Lagrange Point\nData source: ${dataSource}`,
                position: lp.pos,
                point: {
                    pixelSize: 10,
                    color: lp.color,
                    outlineColor: Cesium.Color.WHITE,
                    outlineWidth: 1
                },
                label: {
                    text: lp.name,
                    font: '12px monospace',
                    fillColor: lp.color,
                    pixelOffset: new Cesium.Cartesian2(12, 0)
                }
            });
            this.orbitEntities.push(entity);
        });

        // Add Moon orbit and Moon
        this.addMoonOrbit();

        // Add Earth label
        const earthLabel = this.viewer.entities.add({
            position: Cesium.Cartesian3.ZERO,
            label: {
                text: 'Earth',
                font: '12px monospace',
                fillColor: Cesium.Color.CYAN,
                pixelOffset: new Cesium.Cartesian2(0, 20)
            }
        });
        this.orbitEntities.push(earthLabel);
    }

    // Add atmosphere shell for aerodynamics visualization
    addAtmosphereShell() {
        const earthRadius = 6371000; // meters
        const atmosphereHeight = 100000; // 100 km Karman line

        // Add a translucent ellipsoid representing the atmosphere
        const atmosphereEntity = this.viewer.entities.add({
            name: 'Atmosphere (100 km)',
            position: Cesium.Cartesian3.ZERO,
            ellipsoid: {
                radii: new Cesium.Cartesian3(
                    earthRadius + atmosphereHeight,
                    earthRadius + atmosphereHeight,
                    earthRadius + atmosphereHeight
                ),
                material: Cesium.Color.CYAN.withAlpha(0.1),
                outline: true,
                outlineColor: Cesium.Color.CYAN.withAlpha(0.3)
            }
        });
        this.orbitEntities.push(atmosphereEntity);
    }

    // Add integrator comparison visualization: SGP4 vs Full Force Model Propagation
    // Uses actual Tudat library via Emscripten bindings
    // SGP4 is the simplified perturbations model used for TLE propagation
    // Full force model (VCM-like) uses RK78 with:
    //   - Geopotential: J2, J3, J4 zonal harmonics
    //   - Third-body: Sun and Moon gravitational perturbations
    //   - Drag: Exponential atmosphere model with Cd=2.2, A/m=0.01 m²/kg
    //   - SRP: Solar radiation pressure with Cr=1.5
    addIntegratorComparisonVisualization(period, numOrbits = 10) {
        // Vallado TLE test case - well-known benchmark satellite
        const tleLine1 = '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753';
        const tleLine2 = '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667';

        const totalTime = period * numOrbits;
        const numSamples = 360 * numOrbits;

        // Check if Tudat bindings are available
        const hasTudatBindings = typeof Module !== 'undefined' &&
                                  (typeof Module.propagateJ2vsFullForce === 'function' ||
                                   typeof Module.propagateSGP4vsFullForce === 'function' ||
                                   typeof Module.propagateSGP4vsJ2 === 'function');

        let ephemerisData = null;

        if (hasTudatBindings) {
            this.log('Using Tudat J2-Only vs Full Force Model comparison', 'info');

            try {
                // Use the new J2 vs Full Force comparison (both numerical, osculating elements)
                if (typeof Module.propagateJ2vsFullForce === 'function') {
                    ephemerisData = Module.propagateJ2vsFullForce(tleLine1, tleLine2, totalTime, numSamples);
                } else if (typeof Module.propagateSGP4vsFullForce === 'function') {
                    ephemerisData = Module.propagateSGP4vsFullForce(tleLine1, tleLine2, totalTime, numSamples);
                } else {
                    ephemerisData = Module.propagateSGP4vsJ2(tleLine1, tleLine2, totalTime, numSamples);
                }
                this.log(`Got ${ephemerisData.length / 7} samples from Tudat`, 'info');
            } catch (e) {
                this.log('J2 vs Full Force propagation failed: ' + e.message, 'warning');
                ephemerisData = null;
            }
        }

        // Fallback to JavaScript Kepler vs RK4 if SGP4 not available
        if (!ephemerisData) {
            this.log('Using JavaScript fallback (Kepler vs RK4)', 'warning');
            const semiMajorAxis = 7200;
            const eccentricity = 0.05;
            const inclination = 35;
            const raan = 45;
            const argPeriapsis = 90;

            const analyticalEph = this.computeKeplerOrbitJS(semiMajorAxis, eccentricity, inclination, raan, argPeriapsis, totalTime, numSamples);
            const numericalEph = this.computeNumericalOrbitJS(semiMajorAxis, eccentricity, inclination, raan, argPeriapsis, totalTime, numSamples, period);

            // Convert to SGP4vsJ2 format: [t, sgp4_x, sgp4_y, sgp4_z, j2_x, j2_y, j2_z, ...]
            ephemerisData = [];
            for (let i = 0; i < numSamples; i++) {
                const idx = i * 4;
                ephemerisData.push(analyticalEph[idx]);      // t
                ephemerisData.push(analyticalEph[idx + 1]);  // sgp4 x (using analytical as proxy)
                ephemerisData.push(analyticalEph[idx + 2]);
                ephemerisData.push(analyticalEph[idx + 3]);
                ephemerisData.push(numericalEph[idx + 1]);   // j2 x (using numerical as proxy)
                ephemerisData.push(numericalEph[idx + 2]);
                ephemerisData.push(numericalEph[idx + 3]);
            }
        }

        // Create position properties from ephemeris data
        const clock = this.viewer.clock;
        const startTime = clock.startTime;

        const sgp4Positions = new Cesium.SampledPositionProperty();
        sgp4Positions.setInterpolationOptions({
            interpolationDegree: 5,
            interpolationAlgorithm: Cesium.LagrangePolynomialApproximation
        });

        const j2Positions = new Cesium.SampledPositionProperty();
        j2Positions.setInterpolationOptions({
            interpolationDegree: 5,
            interpolationAlgorithm: Cesium.LagrangePolynomialApproximation
        });

        // Build reference orbit from first orbit of SGP4 data
        const orbitPositions = [];
        const samplesPerOrbit = Math.floor(numSamples / numOrbits);

        // Populate position properties from ephemeris arrays
        // Format: [t, sgp4_x, sgp4_y, sgp4_z, j2_x, j2_y, j2_z, ...]
        const separationData = [];
        let maxSep = 0;

        for (let i = 0; i < numSamples; i++) {
            const idx = i * 7;
            const t = ephemerisData[idx];
            const sampleTime = Cesium.JulianDate.addSeconds(startTime, t, new Cesium.JulianDate());

            const sgp4Pos = new Cesium.Cartesian3(
                ephemerisData[idx + 1],
                ephemerisData[idx + 2],
                ephemerisData[idx + 3]
            );
            const j2Pos = new Cesium.Cartesian3(
                ephemerisData[idx + 4],
                ephemerisData[idx + 5],
                ephemerisData[idx + 6]
            );

            sgp4Positions.addSample(sampleTime, sgp4Pos);
            j2Positions.addSample(sampleTime, j2Pos);

            // First orbit for reference line
            if (i < samplesPerOrbit) {
                orbitPositions.push(sgp4Pos);
            }

            // Compute separation
            const dx = ephemerisData[idx + 1] - ephemerisData[idx + 4];
            const dy = ephemerisData[idx + 2] - ephemerisData[idx + 5];
            const dz = ephemerisData[idx + 3] - ephemerisData[idx + 6];
            const separation = Math.sqrt(dx*dx + dy*dy + dz*dz);
            separationData.push({ t: t, separation: separation });
            if (separation > maxSep) maxSep = separation;
        }

        // Log debug info
        this.log(`Max separation: ${maxSep.toFixed(2)} m over ${numOrbits} orbits`, 'info');

        // Add reference orbit (dashed white line)
        const refOrbit = this.viewer.entities.add({
            name: 'Reference Orbit',
            polyline: {
                positions: orbitPositions,
                width: 2,
                material: new Cesium.PolylineDashMaterialProperty({
                    color: Cesium.Color.WHITE.withAlpha(0.4),
                    dashLength: 16
                })
            }
        });
        this.orbitEntities.push(refOrbit);

        // Add SGP4 satellite (cyan)
        const sgp4Sat = this.viewer.entities.add({
            name: 'SGP4 (TLE)',
            description: `SGP4 simplified perturbations\nUsed for TLE propagation\nIncludes simplified J2, drag`,
            position: sgp4Positions,
            orientation: new Cesium.VelocityOrientationProperty(sgp4Positions),
            point: {
                pixelSize: 12,
                color: Cesium.Color.CYAN,
                outlineColor: Cesium.Color.WHITE,
                outlineWidth: 2
            },
            path: {
                show: true,
                leadTime: 0,
                trailTime: period * 0.5,
                width: 2,
                material: new Cesium.PolylineGlowMaterialProperty({
                    glowPower: 0.3,
                    color: Cesium.Color.CYAN
                })
            },
            label: {
                text: 'SGP4',
                font: '12px monospace',
                fillColor: Cesium.Color.CYAN,
                pixelOffset: new Cesium.Cartesian2(0, -15)
            },
            viewFrom: new Cesium.Cartesian3(-50000, 0, -20000)
        });
        this.orbitEntities.push(sgp4Sat);

        // Add J2 numerical satellite (lime green)
        const j2Sat = this.viewer.entities.add({
            name: 'J2 Numerical',
            description: `RK4 numerical integration\nwith J2 oblateness perturbation\n${numOrbits} orbits propagated`,
            position: j2Positions,
            orientation: new Cesium.VelocityOrientationProperty(j2Positions),
            point: {
                pixelSize: 12,
                color: Cesium.Color.LIME,
                outlineColor: Cesium.Color.WHITE,
                outlineWidth: 2
            },
            path: {
                show: true,
                leadTime: 0,
                trailTime: period * 0.5,
                width: 2,
                material: new Cesium.PolylineGlowMaterialProperty({
                    glowPower: 0.3,
                    color: Cesium.Color.LIME
                })
            },
            label: {
                text: 'J2 Num',
                font: '12px monospace',
                fillColor: Cesium.Color.LIME,
                pixelOffset: new Cesium.Cartesian2(0, -15)
            }
        });
        this.orbitEntities.push(j2Sat);

        // Create separation chart
        this.createSeparationChart(separationData, totalTime, startTime);

        // Zoom to show orbit
        this.viewer.camera.flyTo({
            destination: Cesium.Cartesian3.fromDegrees(0, 0, 25000000),
            orientation: {
                heading: 0,
                pitch: Cesium.Math.toRadians(-90),
                roll: 0
            },
            duration: 1.0
        });
    }

    // OMM vs J2 comparison visualization
    // Shows divergence between OMM mean-element propagation and J2 numerical propagation
    // REQUIRES Tudat WASM bindings - no JavaScript fallbacks per Agents.md
    addOMMvsJ2Visualization(period, numOrbits = 20) {
        // Example OMM - ISS-like orbit
        const ommElements = {
            semiMajorAxis: 6793,      // km (ISS ~420km altitude)
            eccentricity: 0.0001,     // Nearly circular
            inclination: 51.6,        // degrees
            raan: 45.0,               // degrees
            argPeriapsis: 90.0,       // degrees
            meanAnomaly: 0.0          // degrees (at epoch)
        };

        const totalTime = period * numOrbits;
        const numSamples = 360 * numOrbits;

        // Tudat bindings are REQUIRED - no fallbacks
        if (typeof Module === 'undefined' || typeof Module.propagateOMMvsJ2 !== 'function') {
            this.log('ERROR: Tudat WASM bindings required for OMM vs J2 visualization', 'error');
            this.log('Module.propagateOMMvsJ2 function not available', 'error');
            return;
        }

        this.log('Using Tudat OMM vs J2 comparison', 'info');

        let ephemerisData;
        try {
            const ommJson = JSON.stringify(ommElements);
            ephemerisData = Module.propagateOMMvsJ2(ommJson, totalTime, numSamples);
            this.log(`Got ${ephemerisData.length / 7} samples from Tudat`, 'info');
        } catch (e) {
            this.log('OMM vs J2 propagation failed: ' + e.message, 'error');
            return;
        }

        const clock = this.viewer.clock;
        const startTime = clock.startTime;

        const ommPositions = new Cesium.SampledPositionProperty();
        ommPositions.setInterpolationOptions({
            interpolationDegree: 5,
            interpolationAlgorithm: Cesium.LagrangePolynomialApproximation
        });

        const j2Positions = new Cesium.SampledPositionProperty();
        j2Positions.setInterpolationOptions({
            interpolationDegree: 5,
            interpolationAlgorithm: Cesium.LagrangePolynomialApproximation
        });

        const orbitPositions = [];
        const samplesPerOrbit = Math.floor(numSamples / numOrbits);
        const separationData = [];
        let maxSep = 0;

        for (let i = 0; i < numSamples; i++) {
            const idx = i * 7;
            const t = ephemerisData[idx];
            const sampleTime = Cesium.JulianDate.addSeconds(startTime, t, new Cesium.JulianDate());

            const ommPos = new Cesium.Cartesian3(
                ephemerisData[idx + 1],
                ephemerisData[idx + 2],
                ephemerisData[idx + 3]
            );
            const j2Pos = new Cesium.Cartesian3(
                ephemerisData[idx + 4],
                ephemerisData[idx + 5],
                ephemerisData[idx + 6]
            );

            ommPositions.addSample(sampleTime, ommPos);
            j2Positions.addSample(sampleTime, j2Pos);

            if (i < samplesPerOrbit) {
                orbitPositions.push(ommPos);
            }

            const dx = ephemerisData[idx + 1] - ephemerisData[idx + 4];
            const dy = ephemerisData[idx + 2] - ephemerisData[idx + 5];
            const dz = ephemerisData[idx + 3] - ephemerisData[idx + 6];
            const separation = Math.sqrt(dx*dx + dy*dy + dz*dz);
            separationData.push({ t: t, separation: separation });
            if (separation > maxSep) maxSep = separation;
        }

        this.log(`Max separation: ${maxSep.toFixed(2)} m over ${numOrbits} orbits`, 'info');

        // Reference orbit (first orbit)
        const refOrbit = this.viewer.entities.add({
            name: 'Reference Orbit',
            polyline: {
                positions: orbitPositions,
                width: 2,
                material: new Cesium.PolylineDashMaterialProperty({
                    color: Cesium.Color.WHITE.withAlpha(0.4),
                    dashLength: 16
                })
            }
        });
        this.orbitEntities.push(refOrbit);

        // OMM satellite (cyan)
        const ommSat = this.viewer.entities.add({
            name: 'OMM (Mean)',
            description: `OMM mean element propagation\nMean elements evolve smoothly with J2 secular rates\nNo short-period oscillations`,
            position: ommPositions,
            orientation: new Cesium.VelocityOrientationProperty(ommPositions),
            point: {
                pixelSize: 12,
                color: Cesium.Color.CYAN,
                outlineColor: Cesium.Color.WHITE,
                outlineWidth: 2
            },
            path: {
                show: true,
                leadTime: 0,
                trailTime: period * 0.5,
                width: 2,
                material: new Cesium.PolylineGlowMaterialProperty({
                    glowPower: 0.3,
                    color: Cesium.Color.CYAN
                })
            },
            label: {
                text: 'OMM',
                font: '12px monospace',
                fillColor: Cesium.Color.CYAN,
                pixelOffset: new Cesium.Cartesian2(0, -15)
            },
            viewFrom: new Cesium.Cartesian3(-50000, 0, -20000)
        });
        this.orbitEntities.push(ommSat);

        // J2 satellite (lime)
        const j2Sat = this.viewer.entities.add({
            name: 'J2 Numerical',
            description: `J2 numerical integration (osculating)\nIncludes short-period oscillations\n${numOrbits} orbits propagated`,
            position: j2Positions,
            orientation: new Cesium.VelocityOrientationProperty(j2Positions),
            point: {
                pixelSize: 12,
                color: Cesium.Color.LIME,
                outlineColor: Cesium.Color.WHITE,
                outlineWidth: 2
            },
            path: {
                show: true,
                leadTime: 0,
                trailTime: period * 0.5,
                width: 2,
                material: new Cesium.PolylineGlowMaterialProperty({
                    glowPower: 0.3,
                    color: Cesium.Color.LIME
                })
            },
            label: {
                text: 'J2',
                font: '12px monospace',
                fillColor: Cesium.Color.LIME,
                pixelOffset: new Cesium.Cartesian2(0, -15)
            }
        });
        this.orbitEntities.push(j2Sat);

        // Create separation chart with OMM vs J2 title
        this.createSeparationChartOMMvsJ2(separationData, totalTime, startTime);
    }

    // Create separation chart specific to OMM vs J2 comparison
    createSeparationChartOMMvsJ2(separationData, totalTime, startTime) {
        const container = document.getElementById('d3-chart');
        if (!container) return;

        if (this.d3Container) {
            this.d3Container.selectAll('*').remove();
        }
        container.innerHTML = '';

        const titleEl = document.getElementById('chart-title');
        if (titleEl) {
            titleEl.textContent = 'Two-Body (Kepler) vs OMM (Mean Elements): Divergence';
        }

        const canvas = document.createElement('canvas');
        canvas.id = 'separation-chart';
        canvas.style.cssText = 'width: 100%; height: 100%; display: block;';
        container.appendChild(canvas);

        // Defer chart drawing to next frame to ensure container is sized
        requestAnimationFrame(() => {
            this._drawSeparationChartDeferred(canvas, container, separationData, totalTime, startTime);
        });
    }

    // Create separation distance chart with time cursor
    createSeparationChart(separationData, totalTime, startTime) {
        // Use the existing d3-chart container in the lower panel
        const container = document.getElementById('d3-chart');
        if (!container) return;

        // Clear any existing content (both D3 and raw DOM)
        if (this.d3Container) {
            this.d3Container.selectAll('*').remove();
        }
        container.innerHTML = '';  // Clear any non-D3 content too

        // Update the chart title
        const titleEl = document.getElementById('chart-title');
        if (titleEl) {
            titleEl.textContent = 'J2-Only vs Full Force (J3/J4, Sun/Moon, Drag, SRP): Divergence';
        }

        // Create canvas for chart
        const canvas = document.createElement('canvas');
        canvas.id = 'separation-chart';
        canvas.style.cssText = 'width: 100%; height: 100%; display: block;';
        container.appendChild(canvas);

        // Defer chart drawing to next frame to ensure container is sized
        // This fixes first-load rendering issues
        requestAnimationFrame(() => {
            this._drawSeparationChartDeferred(canvas, container, separationData, totalTime, startTime);
        });
    }

    // Deferred chart drawing after layout is complete
    _drawSeparationChartDeferred(canvas, container, separationData, totalTime, startTime) {
        // Draw chart - use container dimensions
        const ctx = canvas.getContext('2d');
        const rect = container.getBoundingClientRect();
        let width = rect.width;
        let height = rect.height;

        // Fallback dimensions if container not yet sized
        if (width < 100 || height < 50) {
            width = container.offsetWidth || 600;
            height = container.offsetHeight || 150;
        }
        canvas.width = width * window.devicePixelRatio;
        canvas.height = height * window.devicePixelRatio;
        ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

        const padding = { left: 55, right: 15, top: 10, bottom: 35 };
        const plotWidth = width - padding.left - padding.right;
        const plotHeight = height - padding.top - padding.bottom;

        // Find max separation for scaling
        const maxSep = Math.max(...separationData.map(d => d.separation));
        const yScale = maxSep > 0 ? plotHeight / (maxSep * 1.1) : 1;

        // Store chart info for updates
        this.separationChart = {
            canvas, ctx, separationData, totalTime, startTime,
            width, height, padding, plotWidth, plotHeight, yScale, maxSep
        };

        this.drawSeparationChart();

        // Update chart cursor on clock tick
        this.viewer.clock.onTick.addEventListener((clock) => {
            this.updateSeparationChartCursor(clock.currentTime);
        });
    }

    drawSeparationChart() {
        const { ctx, separationData, totalTime, width, height, padding, plotWidth, plotHeight, yScale, maxSep } = this.separationChart;

        // Clear
        ctx.clearRect(0, 0, width, height);

        // Background grid
        ctx.strokeStyle = 'rgba(148, 163, 184, 0.1)';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 4; i++) {
            const y = padding.top + (plotHeight * i / 4);
            ctx.beginPath();
            ctx.moveTo(padding.left, y);
            ctx.lineTo(width - padding.right, y);
            ctx.stroke();
        }

        // Y-axis labels
        ctx.fillStyle = '#64748b';
        ctx.font = '10px monospace';
        ctx.textAlign = 'right';
        for (let i = 0; i <= 4; i++) {
            const val = maxSep * (4 - i) / 4;
            const y = padding.top + (plotHeight * i / 4);
            let label;
            if (val >= 1000) label = (val / 1000).toFixed(1) + 'km';
            else if (val >= 1) label = val.toFixed(1) + 'm';
            else label = (val * 100).toFixed(1) + 'cm';
            ctx.fillText(label, padding.left - 5, y + 3);
        }

        // X-axis labels (time)
        ctx.textAlign = 'center';
        const numXLabels = 5;
        for (let i = 0; i <= numXLabels; i++) {
            const t = (totalTime * i / numXLabels);
            const x = padding.left + (plotWidth * i / numXLabels);
            const hours = Math.floor(t / 3600);
            const mins = Math.floor((t % 3600) / 60);
            ctx.fillText(`${hours}h${mins.toString().padStart(2, '0')}m`, x, height - padding.bottom + 20);
        }

        // Draw separation line
        ctx.strokeStyle = '#f59e0b';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        for (let i = 0; i < separationData.length; i++) {
            const d = separationData[i];
            const x = padding.left + (d.t / totalTime) * plotWidth;
            const y = padding.top + plotHeight - (d.separation * yScale);
            if (i === 0) ctx.moveTo(x, y);
            else ctx.lineTo(x, y);
        }
        ctx.stroke();
    }

    updateSeparationChartCursor(currentTime) {
        if (!this.separationChart) return;

        const { ctx, startTime, totalTime, width, height, padding, plotWidth, plotHeight } = this.separationChart;

        // Redraw chart
        this.drawSeparationChart();

        // Calculate current time offset
        const elapsed = Cesium.JulianDate.secondsDifference(currentTime, startTime);
        const fraction = Math.max(0, Math.min(1, elapsed / totalTime));
        const cursorX = padding.left + fraction * plotWidth;

        // Draw cursor line
        ctx.strokeStyle = '#00f0ff';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(cursorX, padding.top);
        ctx.lineTo(cursorX, height - padding.bottom);
        ctx.stroke();

        // Draw cursor dot at intersection
        const idx = Math.floor(fraction * (this.separationChart.separationData.length - 1));
        if (idx >= 0 && idx < this.separationChart.separationData.length) {
            const sep = this.separationChart.separationData[idx].separation;
            const cursorY = padding.top + plotHeight - (sep * this.separationChart.yScale);

            ctx.fillStyle = '#00f0ff';
            ctx.beginPath();
            ctx.arc(cursorX, cursorY, 4, 0, Math.PI * 2);
            ctx.fill();

            // Show current value
            ctx.fillStyle = '#00f0ff';
            ctx.font = 'bold 11px monospace';
            ctx.textAlign = 'left';
            let valText;
            if (sep >= 1000) valText = (sep / 1000).toFixed(2) + ' km';
            else if (sep >= 1) valText = sep.toFixed(2) + ' m';
            else valText = (sep * 100).toFixed(2) + ' cm';
            ctx.fillText(valText, cursorX + 8, cursorY - 5);
        }
    }

    // Orbit Determination / Differential Correction Visualization
    // REQUIRES Tudat WASM bindings - no JavaScript fallbacks per Agents.md
    addOrbitDeterminationVisualization(dynamicsModel = 'fullforce') {
        this.clearOrbitEntities();

        // Tudat bindings are REQUIRED - no fallbacks
        if (typeof Module === 'undefined' || typeof Module.runOrbitDetermination !== 'function') {
            this.log('ERROR: Tudat WASM bindings required for Orbit Determination', 'error');
            return;
        }

        // ISS-like orbit - truth state
        const truthState = {
            semiMajorAxis: 6793,
            eccentricity: 0.0001,
            inclination: 51.6,
            raan: 45.0,
            argPeriapsis: 90.0,
            trueAnomaly: 0.0
        };

        // Initial guess - perturbed from truth
        const initialGuess = {
            semiMajorAxis: 6793 + 5,
            eccentricity: 0.0001,
            inclination: 51.6 + 0.1,
            raan: 45.0 + 0.2,
            argPeriapsis: 90.0,
            trueAnomaly: 0.5
        };

        const period = 5400;
        const duration = period;
        const numObservations = 9;
        const numOrbitSamples = 10000;  // High resolution for smooth rendering at close range
        const noiseStdDev = 100;
        const maxIterations = 20;

        this.log(`Running OD with ${dynamicsModel.toUpperCase()} dynamics...`, 'info');

        let result;
        try {
            result = Module.runOrbitDetermination(
                JSON.stringify(initialGuess),
                JSON.stringify(truthState),
                duration,
                numObservations,
                numOrbitSamples,
                noiseStdDev,
                maxIterations,
                dynamicsModel
            );
        } catch (e) {
            this.log('Orbit determination failed: ' + e.message, 'error');
            return;
        }

        // Parse result
        const numIterations = Math.round(result[0]);
        const numObs = Math.round(result[1]);
        const numSamples = Math.round(result[2]);
        this.log(`Converged in ${numIterations} iterations`, 'info');

        const iterationDataSize = 1 + 6 + numObs * 3;
        const iterations = [];

        let offset = 3;
        for (let iter = 0; iter < numIterations; iter++) {
            const rms = result[offset];
            const state = [];
            for (let j = 0; j < 6; j++) {
                state.push(result[offset + 1 + j]);
            }
            const residuals = [];
            for (let j = 0; j < numObs * 3; j++) {
                residuals.push(result[offset + 7 + j]);
            }
            iterations.push({ rms, state, residuals });
            offset += iterationDataSize;
            this.log(`Iteration ${iter}: RMS = ${rms.toFixed(2)} m`, 'info');
        }

        // Extract truth trajectory (high resolution)
        const truthTrajectory = [];
        for (let i = 0; i < numSamples; i++) {
            truthTrajectory.push({
                x: result[offset + i * 3],
                y: result[offset + i * 3 + 1],
                z: result[offset + i * 3 + 2]
            });
        }
        offset += numSamples * 3;

        // Extract observations
        const observations = [];
        for (let i = 0; i < numObs; i++) {
            observations.push({
                x: result[offset + i * 3],
                y: result[offset + i * 3 + 1],
                z: result[offset + i * 3 + 2]
            });
        }
        offset += numObs * 3;

        // Extract truth positions at observation times (for residuals)
        const truthAtObsTimes = [];
        for (let i = 0; i < numObs; i++) {
            truthAtObsTimes.push({
                x: result[offset + i * 3],
                y: result[offset + i * 3 + 1],
                z: result[offset + i * 3 + 2]
            });
        }
        offset += numObs * 3;

        // Extract estimated positions at observation times (for residuals)
        const estAtObsTimes = [];
        for (let i = 0; i < numObs; i++) {
            estAtObsTimes.push({
                x: result[offset + i * 3],
                y: result[offset + i * 3 + 1],
                z: result[offset + i * 3 + 2]
            });
        }
        offset += numObs * 3;

        // Extract estimated trajectory (high resolution)
        const estimatedTrajectory = [];
        for (let i = 0; i < numSamples; i++) {
            estimatedTrajectory.push({
                x: result[offset + i * 3],
                y: result[offset + i * 3 + 1],
                z: result[offset + i * 3 + 2]
            });
        }
        offset += numSamples * 3;

        // Extract position covariance matrix (3x3, row-major)
        const covMatrix = [];
        for (let i = 0; i < 3; i++) {
            covMatrix.push([]);
            for (let j = 0; j < 3; j++) {
                covMatrix[i].push(result[offset++]);
            }
        }

        // Compute eigenvalues and eigenvectors for the covariance ellipsoid
        // Using a simple power iteration / Jacobi method for 3x3 symmetric matrix
        const ellipsoidAxes = this.computeCovarianceEllipsoid(covMatrix);
        this.log(`Covariance ellipsoid semi-axes: ${ellipsoidAxes.radii.map(r => r.toFixed(1)).join(', ')} m`, 'info');


        // Observations as points (yellow) - also track max residual for observation bounding sphere
        let maxEstResidual = 0;
        for (let i = 0; i < numObs; i++) {
            const obs = observations[i];
            const truth = truthAtObsTimes[i];
            const est = estAtObsTimes[i];

            // Calculate residual magnitudes for tooltip
            const truthResidual = Math.sqrt(
                Math.pow(obs.x - truth.x, 2) +
                Math.pow(obs.y - truth.y, 2) +
                Math.pow(obs.z - truth.z, 2)
            );
            const estResidual = Math.sqrt(
                Math.pow(obs.x - est.x, 2) +
                Math.pow(obs.y - est.y, 2) +
                Math.pow(obs.z - est.z, 2)
            );
            maxEstResidual = Math.max(maxEstResidual, estResidual);

            // Observation point (in orbit)
            const obsEntity = this.viewer.entities.add({
                name: `Observation ${i + 1}`,
                description: `Truth residual: ${truthResidual.toFixed(1)} m\nEstimated residual: ${estResidual.toFixed(1)} m`,
                allowPicking: false,
                position: new Cesium.Cartesian3(obs.x, obs.y, obs.z),
                point: {
                    pixelSize: 10,
                    color: Cesium.Color.YELLOW,
                    outlineColor: Cesium.Color.BLACK,
                    outlineWidth: 1
                }
            });
            this.orbitEntities.push(obsEntity);
        }

        // Create animated satellite following estimated trajectory
        const satelliteColor = dynamicsModel === 'omm' ? Cesium.Color.CYAN : Cesium.Color.LIME;
        const dt = duration / (numSamples - 1);
        const clock = this.viewer.clock;
        const startTime = clock.startTime;

        // Create sampled position property for animation
        const estimatedSampledPosition = new Cesium.SampledPositionProperty();
        estimatedSampledPosition.setInterpolationOptions({
            interpolationDegree: 5,
            interpolationAlgorithm: Cesium.LagrangePolynomialApproximation
        });

        for (let i = 0; i < numSamples; i++) {
            const sampleTime = Cesium.JulianDate.addSeconds(startTime, i * dt, new Cesium.JulianDate());
            const pos = estimatedTrajectory[i];
            estimatedSampledPosition.addSample(sampleTime, new Cesium.Cartesian3(pos.x, pos.y, pos.z));
        }

        // Animated satellite entity with forward-only path (2 revolutions)
        const satellite = this.viewer.entities.add({
            name: dynamicsModel === 'omm' ? 'OMM Estimated' : 'Full Force Estimated',
            description: `Estimated orbit from ${dynamicsModel.toUpperCase()} dynamics\nConverged in ${iterations.length} iterations\nFinal RMS: ${iterations[iterations.length - 1].rms.toFixed(1)} m\n\n1-sigma uncertainty: ${ellipsoidAxes.radii.map(r => r.toFixed(1)).join(' x ')} m`,
            position: estimatedSampledPosition,
            orientation: new Cesium.VelocityOrientationProperty(estimatedSampledPosition),
            point: {
                pixelSize: 12,
                color: satelliteColor,
                outlineColor: Cesium.Color.WHITE,
                outlineWidth: 2
            },
            path: {
                show: true,
                leadTime: period,   // Show 1 revolution forward
                trailTime: period,  // Show 1 revolution behind
                width: 2,
                material: Cesium.Color.LIME
            },
            label: {
                text: dynamicsModel === 'omm' ? 'OMM' : 'Full Force',
                font: '12px monospace',
                fillColor: satelliteColor,
                pixelOffset: new Cesium.Cartesian2(0, -15)
            }
        });
        this.orbitEntities.push(satellite);

        // Covariance ellipsoid - shows position uncertainty
        // Scale up for visibility (multiply by 3 for ~99.7% confidence / 3-sigma)
        const sigmaScale = 3.0;

        // The covariance eigenvalues give the uncertainty in principal axes
        // Map them to RTN frame: largest uncertainty is typically along-track (velocity direction)
        // Sort radii: largest -> along-track (X in velocity frame), then cross-track (Y), then radial (Z)
        const sortedRadii = [...ellipsoidAxes.radii].sort((a, b) => b - a);
        const ellipsoidRadii = new Cesium.Cartesian3(
            sortedRadii[0] * sigmaScale,  // Along-track (largest)
            sortedRadii[1] * sigmaScale,  // Cross-track
            sortedRadii[2] * sigmaScale   // Radial (smallest)
        );

        // Covariance ellipsoid - oriented along velocity
        const covarianceEllipsoid = this.viewer.entities.add({
            name: 'Position Uncertainty (3-sigma)',
            allowPicking: false,
            position: estimatedSampledPosition,
            orientation: new Cesium.VelocityOrientationProperty(estimatedSampledPosition),
            ellipsoid: {
                radii: ellipsoidRadii,
                material: Cesium.Color.CYAN.withAlpha(0.2),
                outline: true,
                outlineColor: Cesium.Color.CYAN.withAlpha(0.6),
                outlineWidth: 1,
                slicePartitions: 24,
                stackPartitions: 24
            }
        });
        this.orbitEntities.push(covarianceEllipsoid);

        // Observation bounding sphere - scaled to contain all observations
        // Add 10% margin to ensure all obs are inside
        const obsBoundingRadius = maxEstResidual * 1.1;
        const observationSphere = this.viewer.entities.add({
            name: 'Observation Scatter',
            allowPicking: false,
            position: estimatedSampledPosition,
            orientation: new Cesium.VelocityOrientationProperty(estimatedSampledPosition),
            ellipsoid: {
                radii: new Cesium.Cartesian3(obsBoundingRadius, obsBoundingRadius, obsBoundingRadius),
                material: Cesium.Color.YELLOW.withAlpha(0.1),
                outline: true,
                outlineColor: Cesium.Color.YELLOW.withAlpha(0.4),
                outlineWidth: 1,
                slicePartitions: 24,
                stackPartitions: 24
            }
        });
        this.orbitEntities.push(observationSphere);
        this.log(`Observation bounding radius: ${obsBoundingRadius.toFixed(1)} m (max residual: ${maxEstResidual.toFixed(1)} m)`, 'info');

        // Configure clock and camera
        this.configureClockForOrbit(duration, null, period / 20);

        // Store data and create residuals chart
        this.odData = {
            iterations: iterations,
            numObservations: numObs,
            noiseStdDev: noiseStdDev,
            dynamicsModel: dynamicsModel
        };
        this.createOrbitDeterminationChart();

        // Store satellite entity and observation positions for camera tracking
        this.odSatellite = satellite;
        this.odObservations = observations;

        // Fly to home view after OD visualization loads
        this.viewer.camera.flyHome(1.0);

        // Enable camera tracking when satellite is clicked - deselect anything else
        this.viewer.selectedEntityChanged.addEventListener((entity) => {
            if (entity === this.odSatellite) {
                this.viewer.trackedEntity = entity;
                this.startCameraTracking();
            } else if (entity !== undefined) {
                // Deselect anything that's not the satellite
                this.viewer.selectedEntity = undefined;
            }
        });

        // Show info panel with instructions
        this.log('Click on satellite to track and show camera info', 'info');
    }

    startCameraTracking() {
        // Show the camera info panel
        const infoPanel = document.getElementById('camera-info-panel');
        if (infoPanel) {
            infoPanel.style.display = 'block';
        }

        // Remove any existing tracking listener
        if (this.cameraTrackingListener) {
            this.viewer.clock.onTick.removeEventListener(this.cameraTrackingListener);
        }

        // Add listener to update camera distance on each tick
        this.cameraTrackingListener = () => {
            this.updateCameraInfo();
        };
        this.viewer.clock.onTick.addEventListener(this.cameraTrackingListener);
    }

    stopCameraTracking() {
        // Hide the camera info panel
        const infoPanel = document.getElementById('camera-info-panel');
        if (infoPanel) {
            infoPanel.style.display = 'none';
        }

        // Remove tracking listener
        if (this.cameraTrackingListener) {
            this.viewer.clock.onTick.removeEventListener(this.cameraTrackingListener);
            this.cameraTrackingListener = null;
        }
    }

    updateCameraInfo() {
        if (!this.viewer || !this.odSatellite) return;

        const trackedEntity = this.viewer.trackedEntity;
        if (!trackedEntity) {
            this.stopCameraTracking();
            return;
        }

        // Get satellite current position
        const currentTime = this.viewer.clock.currentTime;
        const satPosition = trackedEntity.position.getValue(currentTime);
        if (!satPosition) return;

        // Compute camera distance to satellite
        const cameraPosition = this.viewer.camera.positionWC;
        const distance = Cesium.Cartesian3.distance(cameraPosition, satPosition);

        // Update distance display
        const distanceEl = document.getElementById('camera-distance');
        if (distanceEl) {
            if (distance >= 1000000) {
                distanceEl.textContent = (distance / 1000000).toFixed(2) + ' Mm';
            } else if (distance >= 1000) {
                distanceEl.textContent = (distance / 1000).toFixed(2) + ' km';
            } else {
                distanceEl.textContent = distance.toFixed(0) + ' m';
            }
        }

        // Find nearest observation
        if (this.odObservations && this.odObservations.length > 0) {
            let nearestDist = Infinity;
            let nearestIdx = -1;

            for (let i = 0; i < this.odObservations.length; i++) {
                const obs = this.odObservations[i];
                const obsPos = new Cesium.Cartesian3(obs.x, obs.y, obs.z);
                const dist = Cesium.Cartesian3.distance(satPosition, obsPos);
                if (dist < nearestDist) {
                    nearestDist = dist;
                    nearestIdx = i;
                }
            }

            const nearestObsRow = document.getElementById('nearest-obs-row');
            const nearestObsEl = document.getElementById('nearest-obs-dist');

            if (nearestDist < 1000000) {  // Within 1000 km
                if (nearestObsRow) nearestObsRow.style.display = 'flex';
                if (nearestObsEl) {
                    if (nearestDist >= 1000) {
                        nearestObsEl.textContent = `#${nearestIdx + 1}: ${(nearestDist / 1000).toFixed(1)} km`;
                    } else {
                        nearestObsEl.textContent = `#${nearestIdx + 1}: ${nearestDist.toFixed(0)} m`;
                    }
                }
            } else {
                if (nearestObsRow) nearestObsRow.style.display = 'none';
            }
        }
    }

    // Compute eigenvalues and eigenvectors of a 3x3 symmetric covariance matrix
    // Returns { radii: [sqrt(eigenvalue1), ...], vectors: [[v1x,v1y,v1z], ...] }
    // Uses Jacobi eigenvalue algorithm for symmetric matrices
    computeCovarianceEllipsoid(cov) {
        // Copy matrix (will be modified)
        const A = [
            [cov[0][0], cov[0][1], cov[0][2]],
            [cov[1][0], cov[1][1], cov[1][2]],
            [cov[2][0], cov[2][1], cov[2][2]]
        ];

        // Initialize eigenvectors as identity matrix
        const V = [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ];

        // Jacobi iteration
        const maxIter = 50;
        const tolerance = 1e-12;

        for (let iter = 0; iter < maxIter; iter++) {
            // Find largest off-diagonal element
            let maxVal = 0;
            let p = 0, q = 1;
            for (let i = 0; i < 3; i++) {
                for (let j = i + 1; j < 3; j++) {
                    if (Math.abs(A[i][j]) > maxVal) {
                        maxVal = Math.abs(A[i][j]);
                        p = i;
                        q = j;
                    }
                }
            }

            if (maxVal < tolerance) break;

            // Compute rotation angle
            const theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
            const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
            const c = 1 / Math.sqrt(t * t + 1);
            const s = t * c;

            // Apply rotation to A
            const App = A[p][p];
            const Aqq = A[q][q];
            const Apq = A[p][q];

            A[p][p] = c * c * App - 2 * s * c * Apq + s * s * Aqq;
            A[q][q] = s * s * App + 2 * s * c * Apq + c * c * Aqq;
            A[p][q] = A[q][p] = 0;

            for (let k = 0; k < 3; k++) {
                if (k !== p && k !== q) {
                    const Akp = A[k][p];
                    const Akq = A[k][q];
                    A[k][p] = A[p][k] = c * Akp - s * Akq;
                    A[k][q] = A[q][k] = s * Akp + c * Akq;
                }
            }

            // Apply rotation to eigenvectors
            for (let k = 0; k < 3; k++) {
                const Vkp = V[k][p];
                const Vkq = V[k][q];
                V[k][p] = c * Vkp - s * Vkq;
                V[k][q] = s * Vkp + c * Vkq;
            }
        }

        // Extract eigenvalues (diagonal of A) and compute radii
        const eigenvalues = [A[0][0], A[1][1], A[2][2]];
        const radii = eigenvalues.map(ev => Math.sqrt(Math.max(ev, 0)));

        // Extract eigenvectors (columns of V)
        const vectors = [
            [V[0][0], V[1][0], V[2][0]],
            [V[0][1], V[1][1], V[2][1]],
            [V[0][2], V[1][2], V[2][2]]
        ];

        // Sort by eigenvalue (largest first)
        const indices = [0, 1, 2].sort((a, b) => eigenvalues[b] - eigenvalues[a]);
        const sortedRadii = indices.map(i => radii[i]);
        const sortedVectors = indices.map(i => vectors[i]);

        return { radii: sortedRadii, vectors: sortedVectors };
    }

    createOrbitDeterminationChart() {
        const container = document.getElementById('d3-chart');
        if (!container || !this.odData) return;

        if (this.d3Container) {
            this.d3Container.selectAll('*').remove();
        }
        container.innerHTML = '';

        const titleEl = document.getElementById('chart-title');
        if (titleEl) {
            const model = this.odData.dynamicsModel === 'omm' ? 'OMM' : 'Full Force';
            titleEl.textContent = `Orbit Determination Convergence (${model})`;
        }

        const canvas = document.createElement('canvas');
        canvas.id = 'od-chart';
        canvas.style.cssText = 'width: 100%; height: 100%; display: block;';
        container.appendChild(canvas);

        const ctx = canvas.getContext('2d');
        const rect = container.getBoundingClientRect();
        const width = rect.width;
        const height = rect.height;
        canvas.width = width * window.devicePixelRatio;
        canvas.height = height * window.devicePixelRatio;
        ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

        const padding = { left: 60, right: 15, top: 10, bottom: 35 };
        const plotWidth = width - padding.left - padding.right;
        const plotHeight = height - padding.top - padding.bottom;

        const { iterations, noiseStdDev } = this.odData;
        const rmsValues = iterations.map(it => it.rms);
        const maxRMS = Math.max(...rmsValues);
        const minRMS = Math.min(...rmsValues);

        // Use log scale if range is large
        const useLog = maxRMS / Math.max(minRMS, 1) > 100;

        // Background
        ctx.fillStyle = '#0f172a';
        ctx.fillRect(0, 0, width, height);

        // Grid
        ctx.strokeStyle = 'rgba(148, 163, 184, 0.1)';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 4; i++) {
            const y = padding.top + (plotHeight * i / 4);
            ctx.beginPath();
            ctx.moveTo(padding.left, y);
            ctx.lineTo(width - padding.right, y);
            ctx.stroke();
        }

        // Y-axis labels (RMS)
        ctx.fillStyle = '#64748b';
        ctx.font = '10px monospace';
        ctx.textAlign = 'right';
        for (let i = 0; i <= 4; i++) {
            let val;
            if (useLog) {
                const logMax = Math.log10(maxRMS);
                const logMin = Math.log10(Math.max(minRMS, 1));
                const logVal = logMax - (logMax - logMin) * i / 4;
                val = Math.pow(10, logVal);
            } else {
                val = maxRMS * (4 - i) / 4;
            }
            const y = padding.top + (plotHeight * i / 4);
            let label;
            if (val >= 1000) label = (val / 1000).toFixed(1) + 'km';
            else if (val >= 1) label = val.toFixed(0) + 'm';
            else label = val.toFixed(2) + 'm';
            ctx.fillText(label, padding.left - 5, y + 3);
        }

        // X-axis labels (iteration)
        ctx.textAlign = 'center';
        for (let i = 0; i < iterations.length; i++) {
            const x = padding.left + (i / Math.max(iterations.length - 1, 1)) * plotWidth;
            ctx.fillText(`${i}`, x, height - padding.bottom + 15);
        }
        ctx.fillText('Iteration', padding.left + plotWidth / 2, height - 5);

        // Draw RMS convergence line
        ctx.strokeStyle = '#f59e0b';
        ctx.lineWidth = 2;
        ctx.beginPath();
        for (let i = 0; i < iterations.length; i++) {
            const x = padding.left + (i / Math.max(iterations.length - 1, 1)) * plotWidth;
            let y;
            if (useLog) {
                const logMax = Math.log10(maxRMS);
                const logMin = Math.log10(Math.max(minRMS, 1));
                const logVal = Math.log10(Math.max(rmsValues[i], 1));
                y = padding.top + plotHeight * (logMax - logVal) / (logMax - logMin);
            } else {
                y = padding.top + plotHeight - (rmsValues[i] / maxRMS) * plotHeight;
            }
            if (i === 0) ctx.moveTo(x, y);
            else ctx.lineTo(x, y);
        }
        ctx.stroke();

        // Draw points
        ctx.fillStyle = '#f59e0b';
        for (let i = 0; i < iterations.length; i++) {
            const x = padding.left + (i / Math.max(iterations.length - 1, 1)) * plotWidth;
            let y;
            if (useLog) {
                const logMax = Math.log10(maxRMS);
                const logMin = Math.log10(Math.max(minRMS, 1));
                const logVal = Math.log10(Math.max(rmsValues[i], 1));
                y = padding.top + plotHeight * (logMax - logVal) / (logMax - logMin);
            } else {
                y = padding.top + plotHeight - (rmsValues[i] / maxRMS) * plotHeight;
            }
            ctx.beginPath();
            ctx.arc(x, y, 5, 0, Math.PI * 2);
            ctx.fill();
        }

        // Draw noise floor reference line
        if (noiseStdDev > 0) {
            let noiseY;
            if (useLog) {
                const logMax = Math.log10(maxRMS);
                const logMin = Math.log10(Math.max(minRMS, 1));
                const logNoise = Math.log10(noiseStdDev);
                noiseY = padding.top + plotHeight * (logMax - logNoise) / (logMax - logMin);
            } else {
                noiseY = padding.top + plotHeight - (noiseStdDev / maxRMS) * plotHeight;
            }

            if (noiseY > padding.top && noiseY < height - padding.bottom) {
                ctx.strokeStyle = '#22c55e';
                ctx.lineWidth = 1;
                ctx.setLineDash([5, 5]);
                ctx.beginPath();
                ctx.moveTo(padding.left, noiseY);
                ctx.lineTo(width - padding.right, noiseY);
                ctx.stroke();
                ctx.setLineDash([]);

                ctx.fillStyle = '#22c55e';
                ctx.textAlign = 'left';
                ctx.fillText('Noise floor', padding.left + 5, noiseY - 5);
            }
        }

        // Legend
        ctx.fillStyle = '#f59e0b';
        ctx.fillRect(width - 100, padding.top + 5, 12, 12);
        ctx.fillStyle = '#94a3b8';
        ctx.textAlign = 'left';
        ctx.fillText('RMS Residual', width - 85, padding.top + 14);
    }

    // JavaScript fallback for Kepler orbit computation
    computeKeplerOrbitJS(semiMajorAxis, eccentricity, inclination, raan, argPeriapsis, duration, numSamples) {
        const mu = 398600.4418; // km^3/s^2
        const n = Math.sqrt(mu / Math.pow(semiMajorAxis, 3));

        const cosRaan = Math.cos(raan * Math.PI / 180);
        const sinRaan = Math.sin(raan * Math.PI / 180);
        const cosInc = Math.cos(inclination * Math.PI / 180);
        const sinInc = Math.sin(inclination * Math.PI / 180);
        const cosArg = Math.cos(argPeriapsis * Math.PI / 180);
        const sinArg = Math.sin(argPeriapsis * Math.PI / 180);

        const ephemeris = [];

        for (let i = 0; i < numSamples; i++) {
            const t = (i / (numSamples - 1)) * duration;
            const M = n * t;

            // Newton-Raphson for Kepler's equation
            let E = M;
            for (let j = 0; j < 15; j++) {
                E = E - (E - eccentricity * Math.sin(E) - M) / (1 - eccentricity * Math.cos(E));
            }

            const trueAnomaly = 2 * Math.atan2(
                Math.sqrt(1 + eccentricity) * Math.sin(E / 2),
                Math.sqrt(1 - eccentricity) * Math.cos(E / 2)
            );

            const r = semiMajorAxis * (1 - eccentricity * Math.cos(E));
            const xOrb = r * Math.cos(trueAnomaly);
            const yOrb = r * Math.sin(trueAnomaly);

            const x = (cosRaan * cosArg - sinRaan * sinArg * cosInc) * xOrb +
                     (-cosRaan * sinArg - sinRaan * cosArg * cosInc) * yOrb;
            const y = (sinRaan * cosArg + cosRaan * sinArg * cosInc) * xOrb +
                     (-sinRaan * sinArg + cosRaan * cosArg * cosInc) * yOrb;
            const z = (sinArg * sinInc) * xOrb + (cosArg * sinInc) * yOrb;

            ephemeris.push(t, x * 1000, y * 1000, z * 1000);
        }

        return ephemeris;
    }

    // JavaScript fallback for numerical integration (simulated RK78 with error accumulation)
    computeNumericalOrbitJS(semiMajorAxis, eccentricity, inclination, raan, argPeriapsis, duration, numSamples, period) {
        const mu = 398600.4418;
        const n = Math.sqrt(mu / Math.pow(semiMajorAxis, 3));

        const cosRaan = Math.cos(raan * Math.PI / 180);
        const sinRaan = Math.sin(raan * Math.PI / 180);
        const cosInc = Math.cos(inclination * Math.PI / 180);
        const sinInc = Math.sin(inclination * Math.PI / 180);
        const cosArg = Math.cos(argPeriapsis * Math.PI / 180);
        const sinArg = Math.sin(argPeriapsis * Math.PI / 180);

        const ephemeris = [];

        for (let i = 0; i < numSamples; i++) {
            const t = (i / (numSamples - 1)) * duration;
            const M = n * t;

            let E = M;
            for (let j = 0; j < 15; j++) {
                E = E - (E - eccentricity * Math.sin(E) - M) / (1 - eccentricity * Math.cos(E));
            }

            const trueAnomaly = 2 * Math.atan2(
                Math.sqrt(1 + eccentricity) * Math.sin(E / 2),
                Math.sqrt(1 - eccentricity) * Math.cos(E / 2)
            );

            const r = semiMajorAxis * (1 - eccentricity * Math.cos(E));

            // Add simulated numerical error (along-track drift)
            const orbitNumber = t / period;
            const phaseError = 0.00005 * orbitNumber * orbitNumber / r;
            const taNum = trueAnomaly + phaseError;

            const xOrb = r * Math.cos(taNum);
            const yOrb = r * Math.sin(taNum);

            const x = (cosRaan * cosArg - sinRaan * sinArg * cosInc) * xOrb +
                     (-cosRaan * sinArg - sinRaan * cosArg * cosInc) * yOrb;
            const y = (sinRaan * cosArg + cosRaan * sinArg * cosInc) * xOrb +
                     (-sinRaan * sinArg + cosRaan * cosArg * cosInc) * yOrb;
            const z = (sinArg * sinInc) * xOrb + (cosArg * sinInc) * yOrb;

            ephemeris.push(t, x * 1000, y * 1000, z * 1000);
        }

        return ephemeris;
    }

    // Add Lambert transfer arc visualization (from Mengali & Quarta test case)
    addLambertTransferVisualization() {
        const earthRadius = 6378.136; // km
        const distanceUnit = earthRadius; // Canonical DU

        // Departure position: (2 R_E, 0, 0)
        const r1 = { x: 2 * distanceUnit, y: 0, z: 0 };
        // Arrival position: (2 R_E, 2√3 R_E, 0)
        const r2 = { x: 2 * distanceUnit, y: 2 * Math.sqrt(3) * distanceUnit, z: 0 };

        // Departure orbit (circular at 2 R_E)
        const depOrbitPositions = [];
        for (let nu = 0; nu <= 360; nu += 5) {
            const rad = nu * Cesium.Math.toRadians(1);
            depOrbitPositions.push(new Cesium.Cartesian3(
                2 * distanceUnit * 1000 * Math.cos(rad),
                2 * distanceUnit * 1000 * Math.sin(rad),
                0
            ));
        }
        const depOrbit = this.viewer.entities.add({
            name: 'Departure Orbit',
            polyline: {
                positions: depOrbitPositions,
                width: 2,
                material: Cesium.Color.GREEN.withAlpha(0.5)
            }
        });
        this.orbitEntities.push(depOrbit);

        // Transfer arc (elliptical segment from test case)
        // The transfer is about 60° arc
        const transferPositions = [];
        const a = 2.5 * distanceUnit; // Approximate transfer SMA
        const e = 0.2;
        for (let t = 0; t <= 1; t += 0.02) {
            const angle = t * Math.PI / 3; // 60 degree arc
            const r = a * (1 - e * e) / (1 + e * Math.cos(angle - Math.PI / 6));
            transferPositions.push(new Cesium.Cartesian3(
                r * 1000 * Math.cos(angle),
                r * 1000 * Math.sin(angle),
                0
            ));
        }
        const transferArc = this.viewer.entities.add({
            name: 'Lambert Transfer Arc',
            polyline: {
                positions: transferPositions,
                width: 4,
                material: new Cesium.PolylineGlowMaterialProperty({
                    glowPower: 0.3,
                    color: Cesium.Color.ORANGE
                })
            }
        });
        this.orbitEntities.push(transferArc);

        // Departure point
        const depPoint = this.viewer.entities.add({
            name: 'Departure',
            position: new Cesium.Cartesian3(r1.x * 1000, r1.y * 1000, r1.z * 1000),
            point: { pixelSize: 12, color: Cesium.Color.LIME },
            label: {
                text: 'Departure\nV = (2736, 6594) m/s',
                font: '10px monospace',
                fillColor: Cesium.Color.LIME,
                pixelOffset: new Cesium.Cartesian2(15, 0)
            }
        });
        this.orbitEntities.push(depPoint);

        // Arrival point
        const arrPoint = this.viewer.entities.add({
            name: 'Arrival',
            position: new Cesium.Cartesian3(r2.x * 1000, r2.y * 1000, r2.z * 1000),
            point: { pixelSize: 12, color: Cesium.Color.RED },
            label: {
                text: 'Arrival\nV = (-1368, 4225) m/s',
                font: '10px monospace',
                fillColor: Cesium.Color.RED,
                pixelOffset: new Cesium.Cartesian2(15, 0)
            }
        });
        this.orbitEntities.push(arrPoint);

        // Add velocity vectors (scaled for visibility)
        const vScale = 1000; // Scale factor for velocity arrows
        const depVel = { x: 2735.8, y: 6594.3 };
        const arrVel = { x: -1367.9, y: 4225.03 };

        const depVelArrow = this.viewer.entities.add({
            polyline: {
                positions: [
                    new Cesium.Cartesian3(r1.x * 1000, r1.y * 1000, 0),
                    new Cesium.Cartesian3(r1.x * 1000 + depVel.x * vScale, r1.y * 1000 + depVel.y * vScale, 0)
                ],
                width: 3,
                material: Cesium.Color.LIME
            }
        });
        this.orbitEntities.push(depVelArrow);

        const arrVelArrow = this.viewer.entities.add({
            polyline: {
                positions: [
                    new Cesium.Cartesian3(r2.x * 1000, r2.y * 1000, 0),
                    new Cesium.Cartesian3(r2.x * 1000 + arrVel.x * vScale, r2.y * 1000 + arrVel.y * vScale, 0)
                ],
                width: 3,
                material: Cesium.Color.RED
            }
        });
        this.orbitEntities.push(arrVelArrow);
    }

    addCoordinateAxes() {
        // Delegate to the more complete reference frames visualization
        this.addReferenceFrames();
    }

    showChartForCategory(category, testName) {
        const cat = category.toLowerCase();
        const titleEl = document.getElementById('chart-title');

        // Skip chart rendering for J2 vs Full Force / integrators - separation chart is created by addIntegratorComparisonVisualization
        // Skip chart rendering for OMM vs Two-Body - separation chart is created by addOMMvsJ2Visualization
        // Skip chart rendering for Orbit Determination - convergence chart is created by addOrbitDeterminationVisualization
        if (cat.includes('j2 vs full force') || cat.includes('omm vs two') || cat.includes('omm vs j2') || cat.includes('orbit determination') || cat.includes('integrat') || cat.includes('rk4') || cat.includes('rk78') || cat.includes('bulirsch')) {
            return;
        }

        // Skip chart rendering for CR3BP - uses orbit selector instead
        if (cat.includes('cr3bp')) {
            return;
        }

        // Skip if title element doesn't exist
        if (!titleEl) {
            return;
        }

        // Unit Conversions: Show degree-radian relationship on unit circle
        if (cat.includes('unit conver')) {
            titleEl.textContent = 'Degree ↔ Radian: Unit Circle';
            // Unit circle with angle markers - makes sense for testing deg/rad conversion
            const points = [];
            for (let d = 0; d <= 360; d += 15) {
                const rad = d * Math.PI / 180;
                points.push({
                    x: Math.cos(rad) * 40,
                    y: Math.sin(rad) * 40,
                    z: 0,
                    color: d % 90 === 0 ? this.chartColors.cyan : this.chartColors.purple,
                    size: d % 90 === 0 ? 8 : 4
                });
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Physical Constants: Bar chart comparing magnitudes (normalized)
        else if (cat.includes('physical const')) {
            titleEl.textContent = 'Physical Constants (log scale)';
            // c, G, AU, μ_Earth - log-scaled for visualization
            this.render3DChart({
                type: 'bars3d',
                scale: 100,
                bars: [
                    { value: 85, color: this.chartColors.cyan },   // c ≈ 3×10⁸
                    { value: 11, color: this.chartColors.purple }, // G ≈ 6.67×10⁻¹¹
                    { value: 110, color: this.chartColors.green }, // AU ≈ 1.5×10¹¹
                    { value: 145, color: this.chartColors.orange } // μ ≈ 3.99×10¹⁴
                ]
            });
        }
        // Orbital Element Conversions (ODTBX): Show Keplerian orbit with elements labeled
        else if (cat.includes('orbital element') || cat.includes('odtbx')) {
            titleEl.textContent = 'Keplerian ↔ Cartesian: 3D Orbit';
            // Elliptical orbit matching test case (e=0.23, i=20.6°)
            const points = [];
            const a = 50, e = 0.23, inc = 20.6 * Math.PI / 180;
            const omega = 274.78 * Math.PI / 180; // arg of periapsis
            for (let nu = 0; nu <= 360; nu += 3) {
                const nuRad = nu * Math.PI / 180;
                const r = a * (1 - e * e) / (1 + e * Math.cos(nuRad));
                // Perifocal frame
                const xP = r * Math.cos(nuRad);
                const yP = r * Math.sin(nuRad);
                // Rotate by omega then incline
                const x = xP * Math.cos(omega) - yP * Math.sin(omega);
                const yTemp = xP * Math.sin(omega) + yP * Math.cos(omega);
                const y = yTemp * Math.cos(inc);
                const z = yTemp * Math.sin(inc);
                points.push({ x, y, z });
            }
            this.render3DChart({ type: 'trajectory', points, scale: 80, color: this.chartColors.orange });
        }
        // Anomaly Conversions: Show M, E, ν relationship on orbit
        else if (cat.includes('anomaly')) {
            titleEl.textContent = 'M → E → ν: Anomaly Mapping';
            // Show how eccentric anomaly relates to true anomaly
            const points = [];
            const e = 0.639; // From test case
            for (let M = 0; M <= 360; M += 10) {
                // Solve Kepler's equation M = E - e*sin(E)
                let E = M * Math.PI / 180;
                for (let iter = 0; iter < 10; iter++) {
                    E = M * Math.PI / 180 + e * Math.sin(E);
                }
                // E is x-axis (scaled), M is y-axis, nu is z-axis
                let nu = 2 * Math.atan(Math.sqrt((1+e)/(1-e)) * Math.tan(E/2));
                if (nu < 0) nu += 2 * Math.PI;
                points.push({
                    x: (E * 180 / Math.PI - 180) * 0.4,
                    y: (M - 180) * 0.4,
                    z: (nu * 180 / Math.PI - 180) * 0.4,
                    color: this.chartColors.purple,
                    size: 5
                });
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 100 });
        }
        // Coordinate Conversions: Cartesian ↔ Spherical
        else if (cat.includes('coordinate')) {
            titleEl.textContent = 'Cartesian ↔ Spherical';
            // Show spherical coordinate grid
            const points = [];
            const r = 40;
            // Latitude circles
            for (let lat = -60; lat <= 60; lat += 30) {
                const latRad = lat * Math.PI / 180;
                for (let lon = 0; lon <= 360; lon += 10) {
                    const lonRad = lon * Math.PI / 180;
                    points.push({
                        x: r * Math.cos(latRad) * Math.cos(lonRad),
                        y: r * Math.cos(latRad) * Math.sin(lonRad),
                        z: r * Math.sin(latRad),
                        color: lat === 0 ? this.chartColors.cyan : this.chartColors.dim,
                        size: lat === 0 ? 3 : 2
                    });
                }
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Eigen Matrix Operations: Show rotation transformation
        else if (cat.includes('eigen') || cat.includes('matrix')) {
            titleEl.textContent = 'Matrix Rotation: 45° about Z';
            // Show vector before and after rotation
            const points = [
                { x: 40, y: 0, z: 0, color: this.chartColors.red, size: 10 },      // Original
                { x: 28.3, y: 28.3, z: 0, color: this.chartColors.green, size: 10 } // Rotated
            ];
            // Add arc showing rotation
            for (let a = 0; a <= 45; a += 3) {
                const rad = a * Math.PI / 180;
                points.push({
                    x: 40 * Math.cos(rad),
                    y: 40 * Math.sin(rad),
                    z: 0,
                    color: this.chartColors.cyan,
                    size: 2
                });
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Kepler Functions: Period vs semi-major axis (T² ∝ a³)
        else if (cat.includes('kepler') && !cat.includes('orbital')) {
            titleEl.textContent = 'Kepler\'s 3rd Law: T² ∝ a³';
            // Show relationship between a and T
            const points = [];
            const mu = 398600; // km³/s²
            for (let a = 7000; a <= 42000; a += 1000) {
                const T = 2 * Math.PI * Math.sqrt(a*a*a / mu) / 3600; // hours
                points.push({
                    x: (a - 25000) / 500,
                    y: (T - 12) * 3,
                    z: 0,
                    color: a > 35000 ? this.chartColors.orange : this.chartColors.cyan,
                    size: 4
                });
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Time Conversions: Julian Day timeline
        else if (cat.includes('time')) {
            titleEl.textContent = 'Julian Date: J2000 Epoch';
            // Timeline showing JD scale
            const points = [];
            for (let d = -180; d <= 180; d += 30) {
                const jd = 2451545.0 + d; // JD around J2000
                points.push({
                    x: d * 0.3,
                    y: 0,
                    z: (jd - 2451545) * 0.2,
                    color: d === 0 ? this.chartColors.green : this.chartColors.cyan,
                    size: d === 0 ? 10 : 4
                });
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Legendre Polynomials: Show P₀, P₁, P₂, P₃
        else if (cat.includes('legendre')) {
            titleEl.textContent = 'Legendre Polynomials Pₙ(x)';
            const points = [];
            for (let x = -1; x <= 1; x += 0.05) {
                // P₀ = 1, P₁ = x, P₂ = (3x²-1)/2, P₃ = (5x³-3x)/2
                const p0 = 1;
                const p1 = x;
                const p2 = (3*x*x - 1) / 2;
                const p3 = (5*x*x*x - 3*x) / 2;
                points.push({ x: x * 50, y: p2 * 30, z: p3 * 30, color: this.chartColors.purple, size: 3 });
            }
            this.render3DChart({ type: 'trajectory', points, scale: 80, color: this.chartColors.purple });
        }
        // Linear Interpolation: Show interpolated vs actual points
        else if (cat.includes('linear interpol')) {
            titleEl.textContent = 'Linear Interpolation: y = 2x + 1';
            const points = [];
            // Known points (larger)
            [0, 1, 2, 3].forEach(x => {
                points.push({ x: x * 20 - 30, y: (2*x + 1) * 10 - 35, z: 0, color: this.chartColors.green, size: 8 });
            });
            // Interpolated points (smaller)
            [0.5, 1.5, 2.5].forEach(x => {
                points.push({ x: x * 20 - 30, y: (2*x + 1) * 10 - 35, z: 0, color: this.chartColors.cyan, size: 5 });
            });
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Cubic Spline: Sin function interpolation
        else if (cat.includes('spline') || cat.includes('cubic')) {
            titleEl.textContent = 'Cubic Spline: sin(x) Interpolation';
            const points = [];
            // Data points
            for (let i = 0; i <= 10; i++) {
                const x = i * Math.PI / 10;
                points.push({
                    x: (x - Math.PI/2) * 30,
                    y: Math.sin(x) * 40,
                    z: 0,
                    color: this.chartColors.green,
                    size: 6
                });
            }
            // Interpolated curve
            for (let x = 0; x <= Math.PI; x += 0.05) {
                points.push({
                    x: (x - Math.PI/2) * 30,
                    y: Math.sin(x) * 40,
                    z: 5,
                    color: this.chartColors.cyan,
                    size: 2
                });
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Reference Frame Transformations
        else if (cat.includes('frame') || cat.includes('reference')) {
            titleEl.textContent = 'Frame Transform: 30° Z-Rotation';
            const points = [];
            const angle = Math.PI / 6; // 30 degrees as in test
            const r = 35;

            // Draw X-Y plane grid for context
            for (let x = -40; x <= 40; x += 20) {
                for (let y = -40; y <= 40; y += 20) {
                    points.push({ x, y, z: 0, color: this.chartColors.dim, size: 1 });
                }
            }

            // Original vector (red) - along X axis in rotating frame
            // Draw as line from origin
            for (let t = 0; t <= 1; t += 0.1) {
                points.push({ x: t * r, y: 0, z: 0, color: this.chartColors.red, size: t === 1 ? 10 : 4 });
            }

            // Transformed vector (green) - rotated 30° in inertial frame
            for (let t = 0; t <= 1; t += 0.1) {
                points.push({
                    x: t * r * Math.cos(angle),
                    y: t * r * Math.sin(angle),
                    z: 0,
                    color: this.chartColors.green,
                    size: t === 1 ? 10 : 4
                });
            }

            // Rotation arc showing the 30° angle
            for (let a = 0; a <= 30; a += 2) {
                const rad = a * Math.PI / 180;
                points.push({ x: 20 * Math.cos(rad), y: 20 * Math.sin(rad), z: 0, color: this.chartColors.cyan, size: 3 });
            }

            // Origin point
            points.push({ x: 0, y: 0, z: 0, color: this.chartColors.cyan, size: 8 });

            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Modified Equinoctial Elements
        else if (cat.includes('equinoctial')) {
            titleEl.textContent = 'Modified Equinoctial ↔ Keplerian';
            // Orbit from test case (a=10000km, e=0.1, i=30°)
            const points = [];
            const a = 40, e = 0.1, inc = 30 * Math.PI / 180;
            for (let nu = 0; nu <= 360; nu += 5) {
                const nuRad = nu * Math.PI / 180;
                const r = a * (1 - e * e) / (1 + e * Math.cos(nuRad));
                const x = r * Math.cos(nuRad);
                const y = r * Math.sin(nuRad) * Math.cos(inc);
                const z = r * Math.sin(nuRad) * Math.sin(inc);
                points.push({ x, y, z });
            }
            this.render3DChart({ type: 'trajectory', points, scale: 80, color: this.chartColors.orange });
        }
        // Statistics: Sample mean and variance
        else if (cat.includes('statistic')) {
            titleEl.textContent = 'Statistics: Mean & Variance';
            // Show data distribution {1,2,3,4,5} as in test
            const data = [1, 2, 3, 4, 5];
            const mean = 3;
            const points = data.map((v, i) => ({
                x: (v - mean) * 20,
                y: (i - 2) * 15,
                z: 0,
                color: v === mean ? this.chartColors.green : this.chartColors.cyan,
                size: 8
            }));
            // Add mean line
            points.push({ x: 0, y: -40, z: 0, color: this.chartColors.green, size: 4 });
            points.push({ x: 0, y: 40, z: 0, color: this.chartColors.green, size: 4 });
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Spherical Harmonics
        else if (cat.includes('spherical harmon')) {
            titleEl.textContent = 'Spherical Harmonics: Geodesy Pₙᵐ';
            this.render3DChart({ type: 'spherical', scale: 80, l: 2, m: 0 });
        }
        // Resource Paths (WASM)
        else if (cat.includes('resource') || cat.includes('path')) {
            titleEl.textContent = 'WASM Virtual Filesystem';
            // Show directory tree structure
            const points = [
                { x: 0, y: 30, z: 0, color: this.chartColors.cyan, size: 12 },     // /tudat_data
                { x: -30, y: 0, z: 0, color: this.chartColors.green, size: 8 },    // ephemeris
                { x: -10, y: 0, z: 0, color: this.chartColors.green, size: 8 },    // earth_orientation
                { x: 10, y: 0, z: 0, color: this.chartColors.green, size: 8 },     // spice_kernels
                { x: 30, y: 0, z: 0, color: this.chartColors.green, size: 8 },     // gravity_models
            ];
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Emscripten Environment
        else if (cat.includes('emscripten') || cat.includes('wasm')) {
            titleEl.textContent = 'WASM Environment';
            this.render3DChart({
                type: 'bars3d',
                scale: 100,
                bars: [
                    { value: 100, color: this.chartColors.green }, // Runtime initialized
                    { value: 100, color: this.chartColors.green }  // FS available
                ]
            });
        }
        // Linear Algebra: Cross product
        else if (cat.includes('linear algebra') || cat.includes('algebra')) {
            titleEl.textContent = 'Cross Product: x × y = z';
            const points = [
                { x: 40, y: 0, z: 0, color: this.chartColors.red, size: 10 },   // x-axis
                { x: 0, y: 40, z: 0, color: this.chartColors.green, size: 10 }, // y-axis
                { x: 0, y: 0, z: 40, color: this.chartColors.cyan, size: 10 },  // z-axis (result)
            ];
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // CR3BP Propagation
        else if (cat.includes('cr3bp') || cat.includes('three-body') || cat.includes('3bp')) {
            titleEl.textContent = 'CR3BP: Normalized Trajectory';
            // Generate trajectory similar to test initial state
            const points = [];
            let x = 0.994, y = 0.853, z = 0.312;
            let vx = 0.195, vy = -0.211, vz = 0.15;
            const mu = 2.528e-5;
            const dt = 0.1;
            for (let t = 0; t < 200; t++) {
                points.push({ x: x * 40, y: y * 40, z: z * 80 });
                // Simplified CR3BP equations
                const r1 = Math.sqrt((x + mu) ** 2 + y ** 2 + z ** 2);
                const r2 = Math.sqrt((x - 1 + mu) ** 2 + y ** 2 + z ** 2);
                const ax = 2 * vy + x - (1 - mu) * (x + mu) / (r1 ** 3) - mu * (x - 1 + mu) / (r2 ** 3);
                const ay = -2 * vx + y - (1 - mu) * y / (r1 ** 3) - mu * y / (r2 ** 3);
                const az = -(1 - mu) * z / (r1 ** 3) - mu * z / (r2 ** 3);
                vx += ax * dt; vy += ay * dt; vz += az * dt;
                x += vx * dt; y += vy * dt; z += vz * dt;
            }
            this.render3DChart({ type: 'trajectory', points, scale: 80, color: this.chartColors.purple });
        }
        // Custom State Propagation
        else if (cat.includes('custom state')) {
            titleEl.textContent = 'Custom Propagation: dS/dt = -0.02';
            // Show linear decay from 500 to 480
            const points = [];
            for (let t = 0; t <= 1000; t += 50) {
                const S = 500 - 0.02 * t;
                points.push({
                    x: t / 1000 * 60 - 30,
                    y: (S - 490) * 4,
                    z: 0,
                    color: this.chartColors.cyan,
                    size: 4
                });
            }
            this.render3DChart({ type: 'trajectory', points, scale: 80, color: this.chartColors.cyan });
        }
        // Mass Propagation
        else if (cat.includes('mass')) {
            titleEl.textContent = 'Mass Propagation: dm/dt = -0.01';
            // Fuel consumption visualization
            this.render3DChart({
                type: 'bars3d',
                scale: 100,
                bars: [
                    { value: 100, color: this.chartColors.green },  // t=0: 500kg
                    { value: 80, color: this.chartColors.green },   // t=200
                    { value: 60, color: this.chartColors.cyan },    // t=400
                    { value: 40, color: this.chartColors.orange },  // t=600
                    { value: 20, color: this.chartColors.red }      // t=800
                ]
            });
        }
        // Two-Body / Kepler Propagation
        else if (cat.includes('two-body') || cat.includes('two body') || cat.includes('kepler propag')) {
            titleEl.textContent = 'Two-Body: Kepler Propagation';
            const points = [];
            const a = 40, e = 0.1;
            for (let nu = 0; nu <= 360; nu += 5) {
                const nuRad = nu * Math.PI / 180;
                const r = a * (1 - e * e) / (1 + e * Math.cos(nuRad));
                points.push({
                    x: r * Math.cos(nuRad),
                    y: r * Math.sin(nuRad),
                    z: 0
                });
            }
            // Add central body
            const orbit = { type: 'trajectory', points, scale: 80, color: this.chartColors.green };
            this.render3DChart(orbit);
        }
        // TLE/SGP4
        else if (cat.includes('tle') || cat.includes('sgp4') || cat.includes('vallado')) {
            titleEl.textContent = 'TLE/SGP4: Satellite Orbit';
            const points = [];
            const r = 35, inc = 51.6 * Math.PI / 180; // ISS-like
            for (let nu = 0; nu <= 360; nu += 5) {
                const nuRad = nu * Math.PI / 180;
                points.push({
                    x: r * Math.cos(nuRad),
                    y: r * Math.sin(nuRad) * Math.cos(inc),
                    z: r * Math.sin(nuRad) * Math.sin(inc)
                });
            }
            this.render3DChart({ type: 'trajectory', points, scale: 80, color: this.chartColors.cyan });
        }
        // SPICE Interface
        else if (cat.includes('spice')) {
            titleEl.textContent = 'SPICE: Time & Frame Conversions';
            // Show J2000 -> ECLIP rotation
            const points = [];
            const obliquity = 23.4 * Math.PI / 180;
            // Ecliptic plane
            for (let a = 0; a <= 360; a += 10) {
                const rad = a * Math.PI / 180;
                points.push({
                    x: 40 * Math.cos(rad),
                    y: 40 * Math.sin(rad) * Math.cos(obliquity),
                    z: 40 * Math.sin(rad) * Math.sin(obliquity),
                    color: this.chartColors.orange,
                    size: 3
                });
            }
            // Equatorial plane
            for (let a = 0; a <= 360; a += 10) {
                const rad = a * Math.PI / 180;
                points.push({
                    x: 40 * Math.cos(rad),
                    y: 40 * Math.sin(rad),
                    z: 0,
                    color: this.chartColors.cyan,
                    size: 3
                });
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Third-Body Perturbation (Gravitation)
        else if (cat.includes('third') && cat.includes('body') || testName.toLowerCase().includes('third body')) {
            titleEl.textContent = 'Third-Body Perturbation: Sun-Earth-Moon';
            // Show perturbed orbit with 3 bodies
            const points = [];
            // Central body (Earth)
            points.push({ x: 0, y: 0, z: 0, color: this.chartColors.cyan, size: 12 });
            // Perturbing body (Sun - far away, shown at scaled distance)
            points.push({ x: -60, y: 0, z: 0, color: this.chartColors.yellow, size: 15 });
            // Moon orbit (perturbed)
            const moonDist = 25;
            for (let nu = 0; nu <= 360; nu += 5) {
                const nuRad = nu * Math.PI / 180;
                // Add small perturbation oscillation from Sun's gravity
                const perturbation = 2 * Math.sin(2 * nuRad);
                const r = moonDist + perturbation;
                points.push({
                    x: r * Math.cos(nuRad),
                    y: r * Math.sin(nuRad),
                    z: 0,
                    color: this.chartColors.dim,
                    size: 2
                });
            }
            // Moon position
            points.push({ x: moonDist, y: 0, z: 0, color: this.chartColors.purple, size: 8 });
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Libration Points (Lagrange Points)
        else if (cat.includes('libration') || testName.toLowerCase().includes('libration') || testName.toLowerCase().includes('lagrange')) {
            titleEl.textContent = 'Lagrange Points: L1, L2, L3, L4, L5';
            const points = [];
            // Primary body (larger mass - Sun)
            const mu = 0.01; // Mass ratio
            points.push({ x: -mu * 60, y: 0, z: 0, color: this.chartColors.yellow, size: 14 });
            // Secondary body (smaller mass - Earth)
            points.push({ x: (1 - mu) * 60, y: 0, z: 0, color: this.chartColors.cyan, size: 10 });
            // L1 (between bodies)
            const L1x = (1 - mu) * 60 - 10;
            points.push({ x: L1x, y: 0, z: 0, color: this.chartColors.red, size: 6 });
            // L2 (beyond secondary)
            const L2x = (1 - mu) * 60 + 10;
            points.push({ x: L2x, y: 0, z: 0, color: this.chartColors.red, size: 6 });
            // L3 (opposite side)
            points.push({ x: -65, y: 0, z: 0, color: this.chartColors.red, size: 6 });
            // L4 (leading triangle)
            points.push({ x: 30 * Math.cos(Math.PI/3), y: 30 * Math.sin(Math.PI/3), z: 0, color: this.chartColors.green, size: 6 });
            // L5 (trailing triangle)
            points.push({ x: 30 * Math.cos(Math.PI/3), y: -30 * Math.sin(Math.PI/3), z: 0, color: this.chartColors.green, size: 6 });
            // Draw orbital path
            for (let a = 0; a <= 360; a += 15) {
                const rad = a * Math.PI / 180;
                points.push({ x: 60 * Math.cos(rad), y: 60 * Math.sin(rad), z: 0, color: this.chartColors.dim, size: 1 });
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Jacobi Energy (CR3BP integral of motion)
        else if (cat.includes('jacobi') || testName.toLowerCase().includes('jacobi')) {
            titleEl.textContent = 'Jacobi Integral: Energy Conservation';
            // Show zero-velocity curves
            const points = [];
            const mu = 0.01;
            // Generate points on zero-velocity surface
            for (let x = -80; x <= 80; x += 8) {
                for (let y = -80; y <= 80; y += 8) {
                    const r1 = Math.sqrt((x/60 + mu) ** 2 + (y/60) ** 2);
                    const r2 = Math.sqrt((x/60 - 1 + mu) ** 2 + (y/60) ** 2);
                    if (r1 > 0.1 && r2 > 0.1) {
                        // Simplified effective potential
                        const U = -0.5 * ((x/60) ** 2 + (y/60) ** 2) - (1 - mu) / r1 - mu / r2;
                        const z = Math.max(-40, Math.min(40, U * 20));
                        points.push({
                            x, y, z,
                            color: z > 0 ? this.chartColors.red : this.chartColors.cyan,
                            size: 3
                        });
                    }
                }
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 100 });
        }
        // Spherical Harmonics Gravity
        else if ((cat.includes('gravity') || cat.includes('gravitation')) && (testName.toLowerCase().includes('spherical') || testName.toLowerCase().includes('harmon'))) {
            titleEl.textContent = 'Gravity Field: J2 Oblateness';
            // Show gravity potential with J2 perturbation
            const points = [];
            const J2 = 0.3; // Exaggerated for visualization
            for (let theta = 10; theta <= 170; theta += 10) {
                for (let phi = 0; phi <= 360; phi += 15) {
                    const t = theta * Math.PI / 180;
                    const p = phi * Math.PI / 180;
                    // r = r0 * (1 - J2 * P2(cos(theta))) where P2 = (3cos²θ - 1)/2
                    const P2 = (3 * Math.cos(t) ** 2 - 1) / 2;
                    const r = 35 * (1 - J2 * P2);
                    points.push({
                        x: r * Math.sin(t) * Math.cos(p),
                        y: r * Math.sin(t) * Math.sin(p),
                        z: r * Math.cos(t),
                        color: P2 > 0 ? this.chartColors.orange : this.chartColors.cyan,
                        size: 3
                    });
                }
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Central Gravity Model
        else if ((cat.includes('gravity') || cat.includes('gravitation')) && (testName.toLowerCase().includes('central') || testName.toLowerCase().includes('point'))) {
            titleEl.textContent = 'Central Gravity: 1/r² Force Field';
            // Show radial force vectors
            const points = [];
            // Central body
            points.push({ x: 0, y: 0, z: 0, color: this.chartColors.yellow, size: 15 });
            // Force vectors at different distances
            for (let r = 15; r <= 45; r += 10) {
                for (let a = 0; a < 360; a += 45) {
                    const rad = a * Math.PI / 180;
                    const x = r * Math.cos(rad);
                    const y = r * Math.sin(rad);
                    // Force magnitude ∝ 1/r²
                    const forceMag = 400 / (r * r);
                    points.push({ x, y, z: 0, color: this.chartColors.cyan, size: 4 + forceMag });
                    // Arrow toward center
                    const dx = -x * 0.2;
                    const dy = -y * 0.2;
                    points.push({ x: x + dx, y: y + dy, z: 0, color: this.chartColors.red, size: 2 });
                }
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Gravitation category (general)
        else if (cat.includes('gravitation') || cat.includes('gravity')) {
            titleEl.textContent = 'Gravitation: Field Models';
            // Show gravity well visualization
            const points = [];
            for (let r = 5; r <= 50; r += 5) {
                for (let a = 0; a < 360; a += 30) {
                    const rad = a * Math.PI / 180;
                    const x = r * Math.cos(rad);
                    const y = r * Math.sin(rad);
                    const z = -100 / r; // Potential well
                    points.push({
                        x, y, z,
                        color: r < 20 ? this.chartColors.red : this.chartColors.cyan,
                        size: 3
                    });
                }
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Exponential Atmosphere
        else if (cat.includes('exponential') || (cat.includes('aerodynamic') && testName.toLowerCase().includes('exponential'))) {
            titleEl.textContent = 'Exponential Atmosphere: ρ = ρ₀e^(-h/H)';
            // Show density profile vs altitude
            const points = [];
            const H = 7.2; // Scale height in km (normalized)
            for (let h = 0; h <= 100; h += 2) {
                const rho = Math.exp(-h / H);
                points.push({
                    x: rho * 50 - 25,
                    y: h - 50,
                    z: 0,
                    color: h < 50 ? this.chartColors.cyan : this.chartColors.purple,
                    size: 4
                });
            }
            // Add altitude markers
            [0, 50, 100].forEach(h => {
                points.push({ x: -40, y: h - 50, z: 0, color: this.chartColors.dim, size: 2 });
                points.push({ x: 30, y: h - 50, z: 0, color: this.chartColors.dim, size: 2 });
            });
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // NRLMSISE-00 Atmosphere
        else if (cat.includes('nrlmsise') || (cat.includes('aerodynamic') && testName.toLowerCase().includes('nrlmsise'))) {
            titleEl.textContent = 'NRLMSISE-00: Species Densities';
            // Show atmospheric composition by altitude
            this.render3DChart({
                type: 'bars3d',
                scale: 100,
                bars: [
                    { value: 100, color: this.chartColors.cyan },   // He
                    { value: 85, color: this.chartColors.green },    // O
                    { value: 70, color: this.chartColors.purple },   // N2
                    { value: 40, color: this.chartColors.orange },   // O2
                    { value: 25, color: this.chartColors.red }       // Ar
                ]
            });
        }
        // Aerodynamics category (general)
        else if (cat.includes('aerodynamic')) {
            titleEl.textContent = 'Aerodynamics: Atmospheric Models';
            // Show density profile
            const points = [];
            for (let h = 0; h <= 100; h += 3) {
                const rho = Math.exp(-h / 8.5);
                points.push({
                    x: rho * 40 - 20,
                    y: h - 50,
                    z: 0,
                    color: this.chartColors.cyan,
                    size: 3
                });
            }
            this.render3DChart({ type: 'trajectory', points, scale: 80, color: this.chartColors.cyan });
        }
        // Lambert Targeting (Mission Segments)
        else if (cat.includes('lambert') || cat.includes('mission segment') || testName.toLowerCase().includes('lambert')) {
            titleEl.textContent = 'Lambert Problem: Transfer Arc';
            const points = [];
            // Departure position (from test: 2R_E, 0, 0)
            const r1 = { x: 40, y: 0, z: 0 };
            // Arrival position (from test: 2R_E, 2√3 R_E, 0)
            const r2 = { x: 40, y: 40 * Math.sqrt(3), z: 0 };
            // Draw departure point
            points.push({ ...r1, color: this.chartColors.green, size: 10 });
            // Draw arrival point
            points.push({ x: r2.x * 0.5, y: r2.y * 0.5, z: 0, color: this.chartColors.red, size: 10 });
            // Central body
            points.push({ x: 0, y: 0, z: 0, color: this.chartColors.cyan, size: 12 });
            // Transfer arc (elliptical)
            for (let t = 0; t <= 1; t += 0.05) {
                const angle = t * Math.PI / 3; // 60 degree transfer
                const r = 40 * (1 + 0.2 * Math.sin(angle * 2));
                points.push({
                    x: r * Math.cos(angle),
                    y: r * Math.sin(angle),
                    z: 0,
                    color: this.chartColors.orange,
                    size: 3
                });
            }
            // Velocity vectors (departure and arrival)
            points.push({ x: r1.x + 5, y: r1.y + 12, z: 0, color: this.chartColors.green, size: 4 }); // V_dep direction
            points.push({ x: r2.x * 0.5 - 3, y: r2.y * 0.5 + 8, z: 0, color: this.chartColors.red, size: 4 }); // V_arr direction
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // RK78 Integrator
        else if (cat.includes('rk78') || cat.includes('runge-kutta 78') || (cat.includes('integrat') && testName.toLowerCase().includes('rk78'))) {
            titleEl.textContent = 'RK78: Adaptive Step Propagation';
            // Show actual RK78 integration of a harmonic oscillator
            const points = [];
            // Integrate d²x/dt² = -x (harmonic oscillator)
            let x = 40, v = 0, t = 0;
            const dt = 0.15;
            for (let step = 0; step < 150; step++) {
                points.push({
                    x: t * 3 - 40,
                    y: x,
                    z: v * 2,
                    color: this.chartColors.cyan,
                    size: 3
                });
                // RK4 step for visualization
                const k1v = -x;
                const k1x = v;
                const k2v = -(x + k1x * dt/2);
                const k2x = v + k1v * dt/2;
                const k3v = -(x + k2x * dt/2);
                const k3x = v + k2v * dt/2;
                const k4v = -(x + k3x * dt);
                const k4x = v + k3v * dt;
                x += (k1x + 2*k2x + 2*k3x + k4x) * dt / 6;
                v += (k1v + 2*k2v + 2*k3v + k4v) * dt / 6;
                t += dt;
            }
            this.render3DChart({ type: 'trajectory', points, scale: 80, color: this.chartColors.cyan });
        }
        // Bulirsch-Stoer Integrator
        else if (cat.includes('bulirsch') || (cat.includes('integrat') && testName.toLowerCase().includes('bulirsch'))) {
            titleEl.textContent = 'Bulirsch-Stoer: Extrapolation Steps';
            // Show integration with extrapolation sequence visualization
            const points = [];
            // Integrate exponential growth dy/dt = y
            let y = 10;
            for (let t = 0; t <= 1; t += 0.02) {
                const exact = 10 * Math.exp(t);
                points.push({
                    x: t * 60 - 30,
                    y: (exact - 15) * 1.5,
                    z: 0,
                    color: this.chartColors.green,
                    size: 4
                });
            }
            // Show extrapolation convergence as vertical bars at key points
            [0.25, 0.5, 0.75].forEach(t => {
                const exact = 10 * Math.exp(t);
                for (let n = 1; n <= 4; n++) {
                    const approx = exact * (1 + 0.1 / n); // Simulate convergence
                    points.push({
                        x: t * 60 - 30,
                        y: (approx - 15) * 1.5,
                        z: n * 8,
                        color: this.chartColors.purple,
                        size: 3
                    });
                }
            });
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Integrators category (general) - actual propagation visualization
        else if (cat.includes('integrat')) {
            titleEl.textContent = 'Numerical Integrators: ODE Solutions';
            // Show multiple integration methods on same problem
            const points = [];
            // Harmonic oscillator: x'' = -x, with x(0)=1, v(0)=0
            // Analytical: x(t) = cos(t)
            for (let t = 0; t <= 2 * Math.PI; t += 0.1) {
                // Exact solution
                points.push({
                    x: t * 10 - 30,
                    y: Math.cos(t) * 30,
                    z: -Math.sin(t) * 30, // velocity
                    color: this.chartColors.cyan,
                    size: 4
                });
            }
            this.render3DChart({ type: 'trajectory', points, scale: 80, color: this.chartColors.cyan });
        }
        // Propagation Termination
        else if (cat.includes('termination') || testName.toLowerCase().includes('termination')) {
            titleEl.textContent = 'Termination Conditions: Event Detection';
            // Show orbit with termination event
            const points = [];
            const a = 40, e = 0.3;
            let terminated = false;
            for (let nu = 0; nu <= 360 && !terminated; nu += 5) {
                const nuRad = nu * Math.PI / 180;
                const r = a * (1 - e * e) / (1 + e * Math.cos(nuRad));
                const x = r * Math.cos(nuRad);
                const y = r * Math.sin(nuRad);
                // Terminate when crossing a threshold (simulating altitude check)
                if (nu > 180 && y < -20) {
                    terminated = true;
                    points.push({ x, y, z: 0, color: this.chartColors.red, size: 12 }); // Termination point
                } else {
                    points.push({
                        x, y, z: 0,
                        color: nu < 180 ? this.chartColors.cyan : this.chartColors.orange,
                        size: 3
                    });
                }
            }
            // Add threshold line
            for (let x = -50; x <= 50; x += 10) {
                points.push({ x, y: -20, z: 0, color: this.chartColors.dim, size: 2 });
            }
            // Central body
            points.push({ x: 0, y: 0, z: 0, color: this.chartColors.green, size: 10 });
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
        // Default fallback
        else {
            titleEl.textContent = category;
            const points = [];
            for (let i = 0; i < 30; i++) {
                const theta = i / 30 * 2 * Math.PI;
                points.push({
                    x: 30 * Math.cos(theta),
                    y: 30 * Math.sin(theta),
                    z: i - 15,
                    color: this.chartColors.cyan,
                    size: 4
                });
            }
            this.render3DChart({ type: 'scatter3d', points, scale: 80 });
        }
    }

    // ==================== Chart-Only Visualizations (Python Example Ports) ====================

    showChartOnlyVisualization(category, testName) {
        // Get containers
        const globeContainer = document.querySelector('.globe-container');
        const trajectoryContainer = document.querySelector('.trajectory-container');
        const chartPanel = document.getElementById('chart-panel');
        const chartContainer = document.getElementById('d3-chart');
        const chartTitle = document.getElementById('chart-title');
        const orbitSelectorPanel = document.getElementById('orbit-selector-panel');
        const odModelPanel = document.getElementById('od-model-panel');

        if (!chartPanel || !chartContainer) {
            this.log('Chart container not available', 'error');
            return;
        }

        // Hide the Cesium globe completely
        if (globeContainer) globeContainer.style.display = 'none';

        // Make trajectory container take FULL space (remove max-height constraint)
        if (trajectoryContainer) {
            trajectoryContainer.style.flex = '1';
            trajectoryContainer.style.maxHeight = 'none';
            trajectoryContainer.style.height = '100%';
        }

        // Hide other panels, show chart panel that fills entire space
        if (orbitSelectorPanel) orbitSelectorPanel.style.display = 'none';
        if (odModelPanel) odModelPanel.style.display = 'none';

        // Chart panel fills the trajectory container
        chartPanel.style.display = 'flex';
        chartPanel.style.flexDirection = 'column';
        chartPanel.style.height = '100%';
        chartPanel.style.flex = '1';

        // Chart container fills the chart panel
        chartContainer.style.flex = '1';
        chartContainer.style.height = '100%';
        chartContainer.style.minHeight = '0';
        chartContainer.style.overflow = 'auto';

        // Clear previous chart content
        chartContainer.innerHTML = '';

        // Clear the Cesium globe visualization (not needed for chart-only)
        this.clearOrbitEntities();

        // Logging helper
        const log = (msg, level) => this.log(msg, level);

        // Route to the appropriate example visualization
        // Map category names to show functions
        const exampleMap = {
            'Keplerian Orbit': { title: 'Keplerian Orbit (Two-Body Problem)', fn: showKeplerianOrbitExample },
            'Perturbed Orbit': { title: 'Perturbed Orbit (J2 vs Full Force)', fn: showPerturbedOrbitExample },
            'Re-entry Trajectory': { title: 'Atmospheric Re-entry Trajectory', fn: showReentryTrajectoryExample },
            'Solar System': { title: 'Solar System Propagation', fn: showSolarSystemExample },
            'Thrust Satellite': { title: 'Low-Thrust Orbit Transfer', fn: showThrustSatelliteExample },
            'Two-Stage Rocket': { title: 'Two-Stage Rocket Ascent', fn: showTwoStageRocketExample },
            'Linear Sensitivity': { title: 'Linear Sensitivity Analysis', fn: showLinearSensitivityExample },
            'Coupled Dynamics': { title: 'Coupled Orbit-Attitude Dynamics', fn: showCoupledDynamicsExample },
            'CR3BP Manifolds': { title: 'CR3BP Halo Orbit Manifolds', fn: showCR3BPManifoldsExample },
            'Differential Drag': { title: 'Differential Drag Maneuver', fn: showDifferentialDragExample },
            'JUICE Flybys': { title: 'JUICE Jovian Moon Flybys', fn: showJuiceFlybysExample },
            'Earth-Mars Transfer': { title: 'Earth-Mars Transfer Window', fn: showEarthMarsTransferExample },
            'MGA Trajectory': { title: 'Multi-Gravity Assist Trajectory', fn: showMGATrajectoryExample },
            'Hohmann Transfer': { title: 'Hohmann Transfer (LEO to GEO)', fn: showHohmannTransferExample },
            'Gravity Assist': { title: 'Planetary Gravity Assist', fn: showGravityAssistExample },
            'Covariance Propagation': { title: 'Covariance Propagation', fn: showCovariancePropagationExample },
            'Full Estimation': { title: 'Full Parameter Estimation', fn: showFullEstimationExample },
            'Galilean Moons Estimation': { title: 'Galilean Moons State Estimation', fn: showGalileanMoonsEstimationExample },
            'Himmelblau Optimization': { title: 'Himmelblau Function Optimization', fn: showHimmelblauOptimizationExample },
            'Asteroid Orbit Optimization': { title: 'Asteroid Mission Optimization', fn: showAsteroidOrbitOptimizationExample },
            'Cassini MGA': { title: 'Cassini MGA Trajectory Optimization', fn: showCassiniMGAExample },
            'Low-Thrust Porkchop': { title: 'Low-Thrust Transfer Windows', fn: showLowThrustPorkchopExample },
            'Earth-Moon Thrust': { title: 'Earth-Moon Low-Thrust Transfer', fn: showEarthMoonThrustExample },
            'Estimation Dynamical Models': { title: 'Mars Express Estimation (Model Mismatch)', fn: showEstimationDynamicalModelsExample },
            'MPC Asteroid Estimation': { title: 'Asteroid Orbit from MPC Data', fn: showMPCAsteroidEstimationExample },
            'Hodographic Shaping MGA': { title: 'Low-Thrust MGA Optimization', fn: showHodographicShapingMGAExample }
        };

        const example = exampleMap[category];
        if (example) {
            if (chartTitle) chartTitle.textContent = example.title;
            example.fn(chartContainer, log);
        } else {
            // Unknown chart-only visualization
            if (chartTitle) chartTitle.textContent = category;
            chartContainer.innerHTML = `<div style="padding: 20px; color: var(--text-secondary);">Chart-only visualization for "${category}" not implemented yet.</div>`;
        }
    }

    // Restore normal visualization layout (globe + chart)
    restoreGlobeLayout() {
        const globeContainer = document.querySelector('.globe-container');
        const trajectoryContainer = document.querySelector('.trajectory-container');
        const chartPanel = document.getElementById('chart-panel');
        const chartContainer = document.getElementById('d3-chart');

        if (globeContainer) globeContainer.style.display = '';
        if (trajectoryContainer) {
            trajectoryContainer.style.flex = '';
            trajectoryContainer.style.height = '';
            trajectoryContainer.style.maxHeight = '';
        }
        if (chartPanel) {
            chartPanel.style.flex = '';
        }
        if (chartContainer) {
            chartContainer.style.minHeight = '';
            chartContainer.style.overflow = '';
        }
    }

    clearOrbitEntities() {
        // Stop camera tracking if active
        this.stopCameraTracking();

        // Clear tracked entity
        if (this.viewer) {
            this.viewer.trackedEntity = undefined;
        }

        // Clear OD-specific data
        this.odSatellite = null;
        this.odObservations = null;

        this.orbitEntities.forEach(entity => {
            this.viewer.entities.remove(entity);
        });
        this.orbitEntities = [];
    }

    generateOrbitalData() {
        // Keep for compatibility but not used anymore
        this.orbitalData = {};
    }

    // ==================== D3.js 3D Charts ====================

    setupCharts() {
        this.chartColors = {
            cyan: '#00f0ff',
            purple: '#8b5cf6',
            green: '#00ff9d',
            red: '#ff3366',
            orange: '#ff9f1c',
            yellow: '#ffd93d',
            dim: '#4a6066',
            bg: '#0d1620'
        };

        // D3 chart container (may not exist if using orbit selector instead)
        const container = document.getElementById('d3-chart');
        if (!container) {
            // Orbit selector is being used instead of D3 charts
            return;
        }

        this.d3Container = d3.select('#d3-chart');
        this.d3Rotation = { x: -20, y: 30, z: 0 };
        this.isDragging = false;

        // Setup mouse drag rotation
        container.addEventListener('mousedown', (e) => {
            this.isDragging = true;
            this.lastMouse = { x: e.clientX, y: e.clientY };
        });
        document.addEventListener('mousemove', (e) => {
            if (this.isDragging && this.lastMouse) {
                const dx = e.clientX - this.lastMouse.x;
                const dy = e.clientY - this.lastMouse.y;
                this.d3Rotation.y += dx * 0.5;
                this.d3Rotation.x += dy * 0.5;
                this.lastMouse = { x: e.clientX, y: e.clientY };
                if (this.currentD3Data) {
                    this.render3DChart(this.currentD3Data);
                }
            }
        });
        document.addEventListener('mouseup', () => {
            this.isDragging = false;
        });
    }

    // Project 3D point to 2D using rotation
    project3D(x, y, z, width, height, scale) {
        const rx = this.d3Rotation.x * Math.PI / 180;
        const ry = this.d3Rotation.y * Math.PI / 180;

        // Rotate around X axis
        let y1 = y * Math.cos(rx) - z * Math.sin(rx);
        let z1 = y * Math.sin(rx) + z * Math.cos(rx);

        // Rotate around Y axis
        let x2 = x * Math.cos(ry) + z1 * Math.sin(ry);
        let z2 = -x * Math.sin(ry) + z1 * Math.cos(ry);

        // Perspective projection
        const perspective = 500;
        const factor = perspective / (perspective + z2);

        return {
            x: width / 2 + x2 * scale * factor,
            y: height / 2 - y1 * scale * factor,
            z: z2
        };
    }

    render3DChart(data) {
        this.currentD3Data = data;
        const container = document.getElementById('d3-chart');
        if (!container) return;
        const width = container.clientWidth;
        const height = container.clientHeight;

        // Clear previous (both D3 and raw DOM like canvas elements)
        if (this.d3Container) {
            this.d3Container.selectAll('*').remove();
        }
        container.innerHTML = '';

        const svg = this.d3Container.append('svg')
            .attr('width', width)
            .attr('height', height)
            .style('background', this.chartColors.bg);

        // Draw based on chart type
        if (data.type === 'scatter3d') {
            this.draw3DScatter(svg, data, width, height);
        } else if (data.type === 'surface') {
            this.draw3DSurface(svg, data, width, height);
        } else if (data.type === 'trajectory') {
            this.draw3DTrajectory(svg, data, width, height);
        } else if (data.type === 'bars3d') {
            this.draw3DBars(svg, data, width, height);
        } else if (data.type === 'spherical') {
            this.draw3DSpherical(svg, data, width, height);
        }

        // Draw axes
        this.draw3DAxes(svg, data.scale || 100, width, height);
    }

    draw3DAxes(svg, scale, width, height) {
        const axes = [
            { start: [0, 0, 0], end: [scale, 0, 0], color: this.chartColors.red, label: 'X' },
            { start: [0, 0, 0], end: [0, scale, 0], color: this.chartColors.green, label: 'Y' },
            { start: [0, 0, 0], end: [0, 0, scale], color: this.chartColors.cyan, label: 'Z' }
        ];

        axes.forEach(axis => {
            const p1 = this.project3D(...axis.start, width, height, 0.8);
            const p2 = this.project3D(...axis.end, width, height, 0.8);

            svg.append('line')
                .attr('x1', p1.x).attr('y1', p1.y)
                .attr('x2', p2.x).attr('y2', p2.y)
                .attr('stroke', axis.color)
                .attr('stroke-width', 1)
                .attr('opacity', 0.5);

            svg.append('text')
                .attr('x', p2.x + 5).attr('y', p2.y)
                .attr('fill', axis.color)
                .attr('font-size', '10px')
                .attr('font-family', 'Share Tech Mono')
                .text(axis.label);
        });
    }

    draw3DScatter(svg, data, width, height) {
        const points = data.points;
        const scale = data.scale || 100;

        // Sort by z for proper depth rendering
        const projected = points.map(p => ({
            ...this.project3D(p.x, p.y, p.z, width, height, scale / Math.max(...points.map(pt => Math.abs(pt.x)), ...points.map(pt => Math.abs(pt.y)), ...points.map(pt => Math.abs(pt.z))) * 0.8),
            color: p.color || this.chartColors.cyan,
            size: p.size || 4
        })).sort((a, b) => a.z - b.z);

        projected.forEach(p => {
            svg.append('circle')
                .attr('cx', p.x).attr('cy', p.y)
                .attr('r', p.size)
                .attr('fill', p.color)
                .attr('opacity', 0.8);
        });
    }

    draw3DTrajectory(svg, data, width, height) {
        const points = data.points;
        const maxVal = Math.max(...points.flatMap(p => [Math.abs(p.x), Math.abs(p.y), Math.abs(p.z)]));
        const scaleFactor = (data.scale || 100) / maxVal * 0.8;

        const projected = points.map(p => this.project3D(p.x, p.y, p.z, width, height, scaleFactor));

        // Draw trajectory line
        const line = d3.line()
            .x(d => d.x)
            .y(d => d.y)
            .curve(d3.curveCardinal);

        svg.append('path')
            .datum(projected)
            .attr('fill', 'none')
            .attr('stroke', data.color || this.chartColors.cyan)
            .attr('stroke-width', 2)
            .attr('d', line);

        // Draw start/end markers
        if (projected.length > 0) {
            svg.append('circle')
                .attr('cx', projected[0].x).attr('cy', projected[0].y)
                .attr('r', 6)
                .attr('fill', this.chartColors.green);

            svg.append('circle')
                .attr('cx', projected[projected.length - 1].x)
                .attr('cy', projected[projected.length - 1].y)
                .attr('r', 6)
                .attr('fill', this.chartColors.red);
        }
    }

    draw3DSurface(svg, data, width, height) {
        const grid = data.grid;
        const scale = data.scale || 100;
        const rows = grid.length;
        const cols = grid[0].length;

        // Create quads
        const quads = [];
        for (let i = 0; i < rows - 1; i++) {
            for (let j = 0; j < cols - 1; j++) {
                const x1 = (j / cols - 0.5) * scale * 2;
                const x2 = ((j + 1) / cols - 0.5) * scale * 2;
                const y1 = (i / rows - 0.5) * scale * 2;
                const y2 = ((i + 1) / rows - 0.5) * scale * 2;
                const z1 = grid[i][j] * scale * 0.5;
                const z2 = grid[i][j + 1] * scale * 0.5;
                const z3 = grid[i + 1][j + 1] * scale * 0.5;
                const z4 = grid[i + 1][j] * scale * 0.5;

                const avgZ = (z1 + z2 + z3 + z4) / 4;
                const p1 = this.project3D(x1, y1, z1, width, height, 0.8);
                const p2 = this.project3D(x2, y1, z2, width, height, 0.8);
                const p3 = this.project3D(x2, y2, z3, width, height, 0.8);
                const p4 = this.project3D(x1, y2, z4, width, height, 0.8);

                quads.push({
                    points: [p1, p2, p3, p4],
                    avgZ: (p1.z + p2.z + p3.z + p4.z) / 4,
                    value: avgZ
                });
            }
        }

        // Sort by depth and render
        quads.sort((a, b) => a.avgZ - b.avgZ);

        const colorScale = d3.scaleSequential(d3.interpolateViridis)
            .domain([d3.min(quads, d => d.value), d3.max(quads, d => d.value)]);

        quads.forEach(quad => {
            const pathData = `M${quad.points[0].x},${quad.points[0].y} L${quad.points[1].x},${quad.points[1].y} L${quad.points[2].x},${quad.points[2].y} L${quad.points[3].x},${quad.points[3].y} Z`;
            svg.append('path')
                .attr('d', pathData)
                .attr('fill', colorScale(quad.value))
                .attr('stroke', this.chartColors.dim)
                .attr('stroke-width', 0.5)
                .attr('opacity', 0.9);
        });
    }

    draw3DBars(svg, data, width, height) {
        const bars = data.bars;
        const scale = data.scale || 100;
        const barWidth = scale / bars.length * 1.5;

        // Create 3D bars with faces
        const allFaces = [];
        bars.forEach((bar, i) => {
            const x = (i / bars.length - 0.5) * scale * 2;
            const h = bar.value * scale * 0.01;
            const w = barWidth * 0.4;
            const color = bar.color || this.chartColors.cyan;

            // Top face
            const top = [
                this.project3D(x - w, h, -w, width, height, 0.8),
                this.project3D(x + w, h, -w, width, height, 0.8),
                this.project3D(x + w, h, w, width, height, 0.8),
                this.project3D(x - w, h, w, width, height, 0.8)
            ];
            allFaces.push({ points: top, z: d3.mean(top, p => p.z), color: d3.color(color).brighter(0.5) });

            // Front face
            const front = [
                this.project3D(x - w, 0, w, width, height, 0.8),
                this.project3D(x + w, 0, w, width, height, 0.8),
                this.project3D(x + w, h, w, width, height, 0.8),
                this.project3D(x - w, h, w, width, height, 0.8)
            ];
            allFaces.push({ points: front, z: d3.mean(front, p => p.z), color: color });

            // Right face
            const right = [
                this.project3D(x + w, 0, -w, width, height, 0.8),
                this.project3D(x + w, 0, w, width, height, 0.8),
                this.project3D(x + w, h, w, width, height, 0.8),
                this.project3D(x + w, h, -w, width, height, 0.8)
            ];
            allFaces.push({ points: right, z: d3.mean(right, p => p.z), color: d3.color(color).darker(0.3) });
        });

        // Sort by depth and render
        allFaces.sort((a, b) => a.z - b.z);
        allFaces.forEach(face => {
            const pathData = `M${face.points[0].x},${face.points[0].y} L${face.points[1].x},${face.points[1].y} L${face.points[2].x},${face.points[2].y} L${face.points[3].x},${face.points[3].y} Z`;
            svg.append('path')
                .attr('d', pathData)
                .attr('fill', face.color)
                .attr('stroke', this.chartColors.bg)
                .attr('stroke-width', 1);
        });
    }

    draw3DSpherical(svg, data, width, height) {
        // Draw spherical harmonic visualization as colored points on a sphere
        const points = [];
        const scale = data.scale || 80;
        const l = data.l || 2;
        const m = data.m || 0;

        for (let theta = 0; theta <= 180; theta += 10) {
            for (let phi = 0; phi <= 360; phi += 10) {
                const t = theta * Math.PI / 180;
                const p = phi * Math.PI / 180;

                // Simplified spherical harmonic (real part of Y_l^m)
                let Y = Math.pow(Math.cos(t), l) * Math.cos(m * p);
                const r = scale * (0.5 + 0.5 * Math.abs(Y));

                const x = r * Math.sin(t) * Math.cos(p);
                const y = r * Math.sin(t) * Math.sin(p);
                const z = r * Math.cos(t);

                points.push({ x, y, z, value: Y });
            }
        }

        const projected = points.map(pt => ({
            ...this.project3D(pt.x, pt.y, pt.z, width, height, 0.8),
            value: pt.value
        })).sort((a, b) => a.z - b.z);

        const colorScale = d3.scaleDiverging(d3.interpolateRdYlBu)
            .domain([-1, 0, 1]);

        projected.forEach(p => {
            svg.append('circle')
                .attr('cx', p.x).attr('cy', p.y)
                .attr('r', 3)
                .attr('fill', colorScale(p.value))
                .attr('opacity', 0.8);
        });
    }

    resetCharts() {
        this.d3Container.selectAll('*').remove();
        this.currentD3Data = null;
    }

    sleep(ms) {
        return new Promise(resolve => setTimeout(resolve, ms));
    }
}

// Initialize
window.addEventListener('DOMContentLoaded', () => {
    window.tudatRunner = new TudatTestRunner();
});
