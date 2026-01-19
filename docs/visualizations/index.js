// Tudat WASM Visualizations
// Central export for all visualization modules

// Shared utilities
export * from './shared/utils.js';
export * from './shared/spice-utils.js';

// Individual visualizations
export * from './cr3bp/index.js';
export * from './libration-points/index.js';
export * from './atmospheric-drag/index.js';
export * from './reference-frames/index.js';
export * from './geostationary/index.js';
export * from './j2-vs-full-force/index.js';
export * from './omm-vs-j2/index.js';
export * from './orbit-determination/index.js';

// Python example ports (chart-only, no 3D globe)
export * from './examples/index.js';

// Visualization registry for dynamic loading
// Order: J2 vs Full Force, OMM vs Two-Body, CR3BP, then ported examples, then other demos
export const visualizationRegistry = {
    // ============ Primary Demos ============
    'J2 vs Full Force': {
        module: './j2-vs-full-force/index.js',
        showFunction: 'showJ2vsFullForceVisualization',
        description: 'J2-only vs full force model'
    },
    'OMM vs Two-Body': {
        module: './omm-vs-j2/index.js',
        showFunction: 'showOMMvsJ2Visualization',
        description: 'OMM secular drift vs Keplerian'
    },
    'CR3BP Dynamics': {
        module: './cr3bp/index.js',
        showFunction: 'showCR3BPVisualization',
        description: 'Circular Restricted 3-Body Problem'
    },

    // ============ Python Example Ports (chart-only) ============
    // === Propagation ===
    'Keplerian Orbit': {
        module: './examples/keplerian-orbit.js',
        showFunction: 'showKeplerianOrbitExample',
        description: 'Two-body satellite orbit propagation',
        chartOnly: true
    },
    'Perturbed Orbit': {
        module: './examples/perturbed-orbit.js',
        showFunction: 'showPerturbedOrbitExample',
        description: 'Multi-body perturbed orbit with drag & SRP',
        chartOnly: true
    },
    'Re-entry Trajectory': {
        module: './examples/reentry-trajectory.js',
        showFunction: 'showReentryTrajectoryExample',
        description: 'Atmospheric re-entry with heating analysis',
        chartOnly: true
    },
    'Solar System': {
        module: './examples/solar-system-propagation.js',
        showFunction: 'showSolarSystemExample',
        description: 'Multi-body planetary propagation',
        chartOnly: true
    },
    'Thrust Satellite': {
        module: './examples/thrust-satellite.js',
        showFunction: 'showThrustSatelliteExample',
        description: 'Low-thrust orbit transfer',
        chartOnly: true
    },
    'Two-Stage Rocket': {
        module: './examples/two-stage-rocket.js',
        showFunction: 'showTwoStageRocketExample',
        description: 'Rocket ascent with staging',
        chartOnly: true
    },
    'Linear Sensitivity': {
        module: './examples/linear-sensitivity.js',
        showFunction: 'showLinearSensitivityExample',
        description: 'State transition matrix analysis',
        chartOnly: true
    },
    'Coupled Dynamics': {
        module: './examples/coupled-dynamics.js',
        showFunction: 'showCoupledDynamicsExample',
        description: 'Coupled orbit and attitude propagation',
        chartOnly: true
    },
    'CR3BP Manifolds': {
        module: './examples/cr3bp-manifolds.js',
        showFunction: 'showCR3BPManifoldsExample',
        description: 'Halo orbit manifolds in Earth-Moon system',
        chartOnly: true
    },
    'Differential Drag': {
        module: './examples/differential-drag.js',
        showFunction: 'showDifferentialDragExample',
        description: 'Satellite separation using differential drag',
        chartOnly: true
    },
    'JUICE Flybys': {
        module: './examples/juice-flybys.js',
        showFunction: 'showJuiceFlybysExample',
        description: 'JUICE mission Jovian moon flybys',
        chartOnly: true
    },
    'Earth-Moon Thrust': {
        module: './examples/earth-moon-thrust.js',
        showFunction: 'showEarthMoonThrustExample',
        description: 'Low-thrust Earth-Moon transfer',
        chartOnly: true
    },
    // === Mission Design ===
    'Earth-Mars Transfer': {
        module: './examples/earth-mars-transfer.js',
        showFunction: 'showEarthMarsTransferExample',
        description: 'Porkchop plot for interplanetary transfer',
        chartOnly: true
    },
    'MGA Trajectory': {
        module: './examples/mga-trajectory.js',
        showFunction: 'showMGATrajectoryExample',
        description: 'Multi-gravity assist (Cassini-like)',
        chartOnly: true
    },
    'Hohmann Transfer': {
        module: './examples/hohmann-transfer.js',
        showFunction: 'showHohmannTransferExample',
        description: 'Classical orbit transfer (LEO to GEO)',
        chartOnly: true
    },
    'Gravity Assist': {
        module: './examples/gravity-assist.js',
        showFunction: 'showGravityAssistExample',
        description: 'Planetary flyby mechanics',
        chartOnly: true
    },
    'Cassini MGA': {
        module: './examples/cassini-mga.js',
        showFunction: 'showCassiniMGAExample',
        description: 'Cassini trajectory optimization',
        chartOnly: true
    },
    'Low-Thrust Porkchop': {
        module: './examples/low-thrust-porkchop.js',
        showFunction: 'showLowThrustPorkchopExample',
        description: 'Low-thrust transfer windows',
        chartOnly: true
    },
    // === Estimation ===
    'Covariance Propagation': {
        module: './examples/covariance-propagation.js',
        showFunction: 'showCovariancePropagationExample',
        description: 'Uncertainty propagation over time',
        chartOnly: true
    },
    'Full Estimation': {
        module: './examples/full-estimation.js',
        showFunction: 'showFullEstimationExample',
        description: 'Parameter estimation with Doppler observations',
        chartOnly: true
    },
    'Galilean Moons Estimation': {
        module: './examples/galilean-moons-estimation.js',
        showFunction: 'showGalileanMoonsEstimationExample',
        description: 'Jupiter moon state estimation',
        chartOnly: true
    },
    'Estimation Dynamical Models': {
        module: './examples/estimation-dynamical-models.js',
        showFunction: 'showEstimationDynamicalModelsExample',
        description: 'Mars Express estimation with model mismatch',
        chartOnly: true
    },
    'MPC Asteroid Estimation': {
        module: './examples/mpc-asteroid-estimation.js',
        showFunction: 'showMPCAsteroidEstimationExample',
        description: 'Asteroid orbit estimation from MPC data',
        chartOnly: true
    },
    // === Optimization (PyGMO) ===
    'Himmelblau Optimization': {
        module: './examples/himmelblau-optimization.js',
        showFunction: 'showHimmelblauOptimizationExample',
        description: 'Multi-modal function optimization',
        chartOnly: true
    },
    'Asteroid Orbit Optimization': {
        module: './examples/asteroid-orbit-optimization.js',
        showFunction: 'showAsteroidOrbitOptimizationExample',
        description: 'Multi-objective mission design',
        chartOnly: true
    },
    'Hodographic Shaping MGA': {
        module: './examples/hodographic-shaping-mga.js',
        showFunction: 'showHodographicShapingMGAExample',
        description: 'Low-thrust MGA trajectory optimization',
        chartOnly: true
    },

    // ============ Other Demos ============
    'Libration Points': {
        module: './libration-points/index.js',
        showFunction: 'showLibrationPointsVisualization',
        description: 'Earth-Moon L1-L5 Lagrange points'
    },
    'Atmospheric Drag': {
        module: './atmospheric-drag/index.js',
        showFunction: 'showAtmosphericDragVisualization',
        description: 'NRLMSISE-00 density model'
    },
    'Reference Frames': {
        module: './reference-frames/index.js',
        showFunction: 'showReferenceFramesVisualization',
        description: 'J2000 / ECLIPJ2000 rotations'
    },
    'Geostationary Orbit': {
        module: './geostationary/index.js',
        showFunction: 'showGeostationaryVisualization',
        description: 'GEO period and mean motion'
    },
    'Orbit Determination': {
        module: './orbit-determination/index.js',
        showFunction: 'showOrbitDeterminationVisualization',
        description: 'Batch least squares with OMM/Full Force toggle'
    }
};
