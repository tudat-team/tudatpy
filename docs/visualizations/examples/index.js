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

// Mission design examples
export * from './earth-mars-transfer.js';
export * from './mga-trajectory.js';
export * from './hohmann-transfer.js';
export * from './gravity-assist.js';

// Estimation examples
export * from './covariance-propagation.js';

// Optimization examples (PyGMO ports)
export * from './himmelblau-optimization.js';
export * from './asteroid-orbit-optimization.js';

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

    // === Estimation ===
    'Covariance Propagation': {
        module: './examples/covariance-propagation.js',
        showFunction: 'showCovariancePropagationExample',
        description: 'Uncertainty propagation over time',
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
    }
};
