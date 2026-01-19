// Tudat Python Examples - WASM Ports
// 2D chart-based visualizations that mirror the official tudatpy examples
// No 3D globe - pure computational examples with D3.js charts

// Propagation examples
export * from './keplerian-orbit.js';
export * from './perturbed-orbit.js';
export * from './reentry-trajectory.js';
export * from './solar-system-propagation.js';
export * from './thrust-satellite.js';
export * from './two-stage-rocket.js';
export * from './linear-sensitivity.js';
export * from './coupled-dynamics.js';
export * from './cr3bp-manifolds.js';
export * from './differential-drag.js';
export * from './juice-flybys.js';
export * from './earth-moon-thrust.js';

// Mission design examples
export * from './earth-mars-transfer.js';
export * from './mga-trajectory.js';
export * from './hohmann-transfer.js';
export * from './gravity-assist.js';
export * from './cassini-mga.js';
export * from './low-thrust-porkchop.js';

// Estimation examples
export * from './covariance-propagation.js';
export * from './full-estimation.js';
export * from './galilean-moons-estimation.js';
export * from './estimation-dynamical-models.js';
export * from './mpc-asteroid-estimation.js';

// Optimization examples (PyGMO ports)
export * from './himmelblau-optimization.js';
export * from './asteroid-orbit-optimization.js';
export * from './hodographic-shaping-mga.js';

// Example registry for the visualization list
export const exampleRegistry = {
    // === Propagation ===
    'Keplerian Orbit': {
        module: './examples/keplerian-orbit.js',
        showFunction: 'showKeplerianOrbitExample',
        description: 'Two-body satellite orbit propagation',
        category: 'Propagation'
    },
    'Perturbed Orbit': {
        module: './examples/perturbed-orbit.js',
        showFunction: 'showPerturbedOrbitExample',
        description: 'Multi-body perturbed orbit with drag & SRP',
        category: 'Propagation'
    },
    'Re-entry Trajectory': {
        module: './examples/reentry-trajectory.js',
        showFunction: 'showReentryTrajectoryExample',
        description: 'Atmospheric re-entry with heating analysis',
        category: 'Propagation'
    },
    'Solar System': {
        module: './examples/solar-system-propagation.js',
        showFunction: 'showSolarSystemExample',
        description: 'Multi-body planetary propagation',
        category: 'Propagation'
    },
    'Thrust Satellite': {
        module: './examples/thrust-satellite.js',
        showFunction: 'showThrustSatelliteExample',
        description: 'Low-thrust orbit transfer',
        category: 'Propagation'
    },
    'Two-Stage Rocket': {
        module: './examples/two-stage-rocket.js',
        showFunction: 'showTwoStageRocketExample',
        description: 'Rocket ascent with staging',
        category: 'Propagation'
    },
    'Linear Sensitivity': {
        module: './examples/linear-sensitivity.js',
        showFunction: 'showLinearSensitivityExample',
        description: 'State transition matrix analysis',
        category: 'Propagation'
    },
    'Coupled Dynamics': {
        module: './examples/coupled-dynamics.js',
        showFunction: 'showCoupledDynamicsExample',
        description: 'Coupled orbit and attitude propagation',
        category: 'Propagation'
    },
    'CR3BP Manifolds': {
        module: './examples/cr3bp-manifolds.js',
        showFunction: 'showCR3BPManifoldsExample',
        description: 'Halo orbit manifolds in Earth-Moon system',
        category: 'Propagation'
    },
    'Differential Drag': {
        module: './examples/differential-drag.js',
        showFunction: 'showDifferentialDragExample',
        description: 'Satellite separation using differential drag',
        category: 'Propagation'
    },
    'JUICE Flybys': {
        module: './examples/juice-flybys.js',
        showFunction: 'showJuiceFlybysExample',
        description: 'JUICE mission Jovian moon flybys',
        category: 'Propagation'
    },
    'Earth-Moon Thrust': {
        module: './examples/earth-moon-thrust.js',
        showFunction: 'showEarthMoonThrustExample',
        description: 'Low-thrust Earth-Moon transfer',
        category: 'Propagation'
    },

    // === Mission Design ===
    'Earth-Mars Transfer': {
        module: './examples/earth-mars-transfer.js',
        showFunction: 'showEarthMarsTransferExample',
        description: 'Porkchop plot for interplanetary transfer',
        category: 'Mission Design'
    },
    'MGA Trajectory': {
        module: './examples/mga-trajectory.js',
        showFunction: 'showMGATrajectoryExample',
        description: 'Multi-gravity assist (Cassini-like)',
        category: 'Mission Design'
    },
    'Hohmann Transfer': {
        module: './examples/hohmann-transfer.js',
        showFunction: 'showHohmannTransferExample',
        description: 'Classical orbit transfer (LEO to GEO)',
        category: 'Mission Design'
    },
    'Gravity Assist': {
        module: './examples/gravity-assist.js',
        showFunction: 'showGravityAssistExample',
        description: 'Planetary flyby mechanics',
        category: 'Mission Design'
    },
    'Cassini MGA': {
        module: './examples/cassini-mga.js',
        showFunction: 'showCassiniMGAExample',
        description: 'Cassini trajectory optimization',
        category: 'Mission Design'
    },
    'Low-Thrust Porkchop': {
        module: './examples/low-thrust-porkchop.js',
        showFunction: 'showLowThrustPorkchopExample',
        description: 'Low-thrust transfer windows',
        category: 'Mission Design'
    },

    // === Estimation ===
    'Covariance Propagation': {
        module: './examples/covariance-propagation.js',
        showFunction: 'showCovariancePropagationExample',
        description: 'Uncertainty propagation over time',
        category: 'Estimation'
    },
    'Full Estimation': {
        module: './examples/full-estimation.js',
        showFunction: 'showFullEstimationExample',
        description: 'Parameter estimation with Doppler observations',
        category: 'Estimation'
    },
    'Galilean Moons Estimation': {
        module: './examples/galilean-moons-estimation.js',
        showFunction: 'showGalileanMoonsEstimationExample',
        description: 'Jupiter moon state estimation',
        category: 'Estimation'
    },
    'Estimation Dynamical Models': {
        module: './examples/estimation-dynamical-models.js',
        showFunction: 'showEstimationDynamicalModelsExample',
        description: 'Mars Express estimation with model mismatch',
        category: 'Estimation'
    },
    'MPC Asteroid Estimation': {
        module: './examples/mpc-asteroid-estimation.js',
        showFunction: 'showMPCAsteroidEstimationExample',
        description: 'Asteroid orbit estimation from MPC data',
        category: 'Estimation'
    },

    // === Optimization (PyGMO) ===
    'Himmelblau Optimization': {
        module: './examples/himmelblau-optimization.js',
        showFunction: 'showHimmelblauOptimizationExample',
        description: 'Multi-modal function optimization',
        category: 'Optimization'
    },
    'Asteroid Orbit Optimization': {
        module: './examples/asteroid-orbit-optimization.js',
        showFunction: 'showAsteroidOrbitOptimizationExample',
        description: 'Multi-objective mission design',
        category: 'Optimization'
    },
    'Hodographic Shaping MGA': {
        module: './examples/hodographic-shaping-mga.js',
        showFunction: 'showHodographicShapingMGAExample',
        description: 'Low-thrust MGA trajectory optimization',
        category: 'Optimization'
    }
};
