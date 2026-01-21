/**
 * SPICE Utilities for Visualizations
 *
 * Provides easy-to-use functions for accessing SPICE ephemeris data
 * in browser-based visualizations.
 */

// Global reference to the Tudat module and SPICE state
let _tudatModule = null;
let _spiceReady = false;
let _ephemerisAvailable = false;  // Whether SPK ephemeris data is loaded (requires binary kernels)
let _precomputedEphemerisLoaded = new Set();  // Track which precomputed ephemerides are loaded

/**
 * Initialize SPICE support with a loaded Tudat module
 * @param {Object} tudatModule - The loaded tudatpy WASM module
 * @param {boolean} ready - Whether SPICE text kernels are loaded (time conversions, constants)
 * @param {boolean} ephemerisAvailable - Whether SPK ephemeris data is available (planetary positions)
 *                                       NOTE: In WASM/browser, this is always false because binary
 *                                       SPK kernels cannot be loaded due to FORTRAN I/O limitations.
 */
export function initSpice(tudatModule, ready = false, ephemerisAvailable = false) {
    _tudatModule = tudatModule;
    _spiceReady = ready;
    _ephemerisAvailable = ephemerisAvailable;
}

/**
 * Check if SPICE is ready for ephemeris queries (planetary positions/states)
 * NOTE: In WASM/browser, this always returns false because binary SPK kernels
 * cannot be loaded. Use isHighAccuracyEphemerisAvailable() or just call getBodyState()
 * which automatically uses the best available data source.
 * @returns {boolean}
 */
export function isSpiceReady() {
    return _ephemerisAvailable && _tudatModule !== null;
}

/**
 * Check if high-accuracy ephemeris is available for a body
 * Returns true if precomputed, CALCEPH, or SPICE ephemeris is available.
 * This is the recommended check before using getBodyState() when you need accurate data.
 * @param {string} target - Target body name
 * @param {string} observer - Observer body name (default: 'Sun')
 * @param {string} frame - Reference frame (default: 'J2000')
 * @returns {boolean}
 */
export function isHighAccuracyEphemerisAvailable(target, observer = 'Sun', frame = 'J2000') {
    // Check precomputed ephemeris first (from pre-converted SPK data via JSON)
    if (isPrecomputedEphemerisAvailable(target, observer, frame)) {
        return true;
    }
    // Check CALCEPH (direct binary SPK reading)
    if (isCalcephEphemerisAvailable(target, observer, frame)) {
        return true;
    }
    // Check if SPICE binary kernels are loaded (rare in WASM)
    return isSpiceReady();
}

/**
 * Check if SPICE text kernels are loaded (for time conversions, constants)
 * @returns {boolean}
 */
export function isSpiceTextKernelsLoaded() {
    return _spiceReady && _tudatModule !== null;
}

/**
 * Get planetary state at a given epoch
 * Priority order:
 * 1. Precomputed ephemeris (from pre-converted SPK data via JSON) - highest accuracy
 * 2. CALCEPH (direct binary SPK reading) - high accuracy, binary files
 * 3. SPICE (if binary kernels somehow loaded) - not available in WASM due to f2c issues
 * 4. Analytical JPL ephemeris - fallback for major planets
 *
 * @param {string} target - Target body (e.g., 'Earth', 'Mars', 'Jupiter')
 * @param {string} observer - Observer body (e.g., 'Sun', 'Earth')
 * @param {number} epoch - Ephemeris time (seconds since J2000)
 * @param {string} frame - Reference frame (default: 'J2000')
 * @returns {Object|null} {x, y, z, vx, vy, vz} in meters and m/s
 */
export function getBodyState(target, observer, epoch, frame = 'J2000') {
    // Priority 1: Try precomputed ephemeris (from pre-converted SPK files via JSON)
    if (isPrecomputedEphemerisAvailable(target, observer, frame)) {
        const state = getPrecomputedState(target, observer, epoch, frame);
        if (state) {
            return state;
        }
    }

    // Priority 2: Try CALCEPH (direct binary SPK reading in WASM)
    if (isCalcephEphemerisAvailable(target, observer, frame)) {
        const state = getCalcephState(target, observer, epoch, frame);
        if (state) {
            return state;
        }
    }

    // Priority 3: Try SPICE if ephemeris data is available (not available in WASM)
    if (isSpiceReady()) {
        try {
            const state = _tudatModule.interface_spice_get_body_cartesian_state_at_epoch(
                target,
                observer,
                frame,
                'NONE',
                epoch
            );

            const result = {
                x: state.get(0),
                y: state.get(1),
                z: state.get(2),
                vx: state.get(3),
                vy: state.get(4),
                vz: state.get(5)
            };

            // Clean up
            state.delete();

            return result;
        } catch (error) {
            console.error(`SPICE query failed for ${target}:`, error);
            // Fall through to analytical ephemeris
        }
    }

    // Priority 4: Use Tudat's analytical JPL ephemeris (no SPICE kernels required)
    return getBodyStateAnalytical(target, observer, epoch);
}

/**
 * Get planetary state using Tudat's analytical JPL ephemeris
 * This function does NOT require SPICE kernels - it uses built-in orbital elements.
 *
 * NOTE: Returns state in ECLIPJ2000 frame relative to Sun.
 * For other observer bodies, computes relative state.
 *
 * @param {string} target - Target body (Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto)
 * @param {string} observer - Observer body (typically 'Sun')
 * @param {number} epoch - Ephemeris time (seconds since J2000)
 * @returns {Object|null} {x, y, z, vx, vy, vz} in meters and m/s
 */
export function getBodyStateAnalytical(target, observer, epoch) {
    if (!_tudatModule) {
        console.warn('Tudat module not loaded');
        return null;
    }

    // Check if analytical ephemeris function exists
    if (!_tudatModule.ephemeris_get_planet_state) {
        console.warn('Analytical ephemeris not available in this build');
        return null;
    }

    // Validate planet name
    if (!_tudatModule.ephemeris_is_valid_planet(target)) {
        console.warn(`Invalid planet name for analytical ephemeris: ${target}`);
        return null;
    }

    try {
        // Get target state relative to Sun (in ECLIPJ2000)
        const targetState = _tudatModule.ephemeris_get_planet_state(target, epoch);

        // Check for null return (indicates error)
        if (!targetState) {
            console.warn(`Analytical ephemeris returned null for ${target} at epoch ${epoch}`);
            return null;
        }

        // If observer is Sun, we're done
        if (observer === 'Sun' || observer === 'SSB' || observer === 'Solar System Barycenter') {
            return {
                x: targetState[0],
                y: targetState[1],
                z: targetState[2],
                vx: targetState[3],
                vy: targetState[4],
                vz: targetState[5]
            };
        }

        // For other observers, get observer state and compute relative
        if (_tudatModule.ephemeris_is_valid_planet(observer)) {
            const observerState = _tudatModule.ephemeris_get_planet_state(observer, epoch);
            if (!observerState) {
                console.warn(`Analytical ephemeris returned null for observer ${observer}`);
                return null;
            }
            return {
                x: targetState[0] - observerState[0],
                y: targetState[1] - observerState[1],
                z: targetState[2] - observerState[2],
                vx: targetState[3] - observerState[3],
                vy: targetState[4] - observerState[4],
                vz: targetState[5] - observerState[5]
            };
        }

        // Unknown observer - return Sun-centered state
        console.warn(`Unknown observer '${observer}', returning Sun-centered state`);
        return {
            x: targetState[0],
            y: targetState[1],
            z: targetState[2],
            vx: targetState[3],
            vy: targetState[4],
            vz: targetState[5]
        };
    } catch (error) {
        console.error(`Analytical ephemeris failed for ${target}:`, error);
        return null;
    }
}

/**
 * Check if analytical ephemeris is available for a body
 * @param {string} bodyName - Body name to check
 * @returns {boolean} True if analytical ephemeris can compute state for this body
 */
export function isAnalyticalEphemerisAvailable(bodyName) {
    if (!_tudatModule || !_tudatModule.ephemeris_is_valid_planet) {
        return false;
    }
    return _tudatModule.ephemeris_is_valid_planet(bodyName);
}

// ============================================================================
// Precomputed Ephemeris Support (from pre-converted SPK data)
// ============================================================================

/**
 * Load precomputed ephemeris from a JSON file
 * JSON format:
 * {
 *   "metadata": { "target": "Earth", "observer": "Sun", "frame": "J2000", ... },
 *   "states": [[epoch, x, y, z, vx, vy, vz], ...]
 * }
 *
 * @param {string} url - URL to the JSON ephemeris file
 * @returns {Promise<boolean>} True if loaded successfully
 */
export async function loadPrecomputedEphemerisFromUrl(url) {
    if (!_tudatModule) {
        console.warn('[EPHEMERIS] Tudat module not loaded');
        return false;
    }

    try {
        console.log(`[EPHEMERIS] Loading precomputed ephemeris from ${url}`);
        const response = await fetch(url);
        if (!response.ok) {
            console.warn(`[EPHEMERIS] Failed to fetch ${url}: ${response.status}`);
            return false;
        }

        const data = await response.json();
        return loadPrecomputedEphemerisFromJson(data);
    } catch (error) {
        console.error(`[EPHEMERIS] Error loading from ${url}:`, error);
        return false;
    }
}

/**
 * Load precomputed ephemeris from parsed JSON data
 * @param {Object} data - Parsed JSON with metadata and states arrays
 * @returns {boolean} True if loaded successfully
 */
export function loadPrecomputedEphemerisFromJson(data) {
    if (!_tudatModule || !_tudatModule.ephemeris_load_precomputed) {
        console.warn('[EPHEMERIS] Precomputed ephemeris functions not available');
        return false;
    }

    const { metadata, states } = data;
    if (!metadata || !states || !Array.isArray(states)) {
        console.error('[EPHEMERIS] Invalid JSON format: missing metadata or states');
        return false;
    }

    const { target, observer, frame } = metadata;
    if (!target || !observer || !frame) {
        console.error('[EPHEMERIS] Invalid metadata: missing target, observer, or frame');
        return false;
    }

    // Extract epochs and flatten states
    const epochs = new Float64Array(states.length);
    const flatStates = new Float64Array(states.length * 6);

    for (let i = 0; i < states.length; i++) {
        const state = states[i];
        epochs[i] = state[0];  // First element is epoch
        flatStates[i * 6 + 0] = state[1];  // x
        flatStates[i * 6 + 1] = state[2];  // y
        flatStates[i * 6 + 2] = state[3];  // z
        flatStates[i * 6 + 3] = state[4];  // vx
        flatStates[i * 6 + 4] = state[5];  // vy
        flatStates[i * 6 + 5] = state[6];  // vz
    }

    // Pass to C++
    console.log(`[EPHEMERIS] Calling ephemeris_load_precomputed with ${epochs.length} epochs, ${flatStates.length} state values`);
    console.log(`[EPHEMERIS] First epoch: ${epochs[0]}, last epoch: ${epochs[epochs.length-1]}`);

    const success = _tudatModule.ephemeris_load_precomputed(
        target, observer, frame,
        Array.from(epochs),
        Array.from(flatStates)
    );

    if (success) {
        const key = `${target.toLowerCase()}_${observer.toLowerCase()}_${frame.toLowerCase()}`;
        _precomputedEphemerisLoaded.add(key);
        console.log(`[EPHEMERIS] Loaded precomputed ephemeris: ${target} relative to ${observer} (${states.length} states)`);

        // Verify it's now available
        const available = _tudatModule.ephemeris_is_precomputed_available(target, observer, frame);
        console.log(`[EPHEMERIS] Verification - isAvailable(${target}, ${observer}, ${frame}): ${available}`);

        // Check time bounds
        if (_tudatModule.ephemeris_get_precomputed_time_bounds) {
            const bounds = _tudatModule.ephemeris_get_precomputed_time_bounds(target, observer, frame);
            console.log(`[EPHEMERIS] Time bounds: ${bounds[0]} to ${bounds[1]} seconds`);
        }
    } else {
        console.error(`[EPHEMERIS] Failed to load precomputed ephemeris for ${target} relative to ${observer}`);
    }

    return success;
}

/**
 * Check if precomputed ephemeris is available for a body
 * @param {string} target - Target body name
 * @param {string} observer - Observer body name
 * @param {string} frame - Reference frame (default: 'J2000')
 * @returns {boolean}
 */
export function isPrecomputedEphemerisAvailable(target, observer, frame = 'J2000') {
    if (!_tudatModule || !_tudatModule.ephemeris_is_precomputed_available) {
        return false;
    }
    return _tudatModule.ephemeris_is_precomputed_available(target, observer, frame);
}

/**
 * Get state from precomputed ephemeris
 * @param {string} target - Target body
 * @param {string} observer - Observer body
 * @param {number} epoch - Ephemeris time (seconds since J2000)
 * @param {string} frame - Reference frame (default: 'J2000')
 * @returns {Object|null} {x, y, z, vx, vy, vz} in meters and m/s
 */
export function getPrecomputedState(target, observer, epoch, frame = 'J2000') {
    if (!isPrecomputedEphemerisAvailable(target, observer, frame)) {
        return null;
    }

    // Check if epoch is within bounds (if available)
    const bounds = getPrecomputedTimeBounds(target, observer, frame);
    if (bounds && (epoch < bounds.start || epoch > bounds.end)) {
        console.warn(`[EPHEMERIS] Epoch ${epoch} outside bounds [${bounds.start}, ${bounds.end}] for ${target}`);
        return null;
    }

    try {
        const state = _tudatModule.ephemeris_get_precomputed_state(target, observer, frame, epoch);
        if (!state || state.length < 6) {
            console.error(`[EPHEMERIS] Invalid state returned for ${target} at epoch ${epoch}`);
            return null;
        }
        return {
            x: state[0],
            y: state[1],
            z: state[2],
            vx: state[3],
            vy: state[4],
            vz: state[5]
        };
    } catch (error) {
        console.error(`[EPHEMERIS] Error getting precomputed state for ${target} at epoch ${epoch}:`, error);
        return null;
    }
}

/**
 * Get time bounds for precomputed ephemeris
 * @param {string} target - Target body
 * @param {string} observer - Observer body
 * @param {string} frame - Reference frame (default: 'J2000')
 * @returns {Object|null} {start, end} epochs in seconds since J2000
 */
export function getPrecomputedTimeBounds(target, observer, frame = 'J2000') {
    if (!_tudatModule || !_tudatModule.ephemeris_get_precomputed_time_bounds) {
        return null;
    }

    try {
        const bounds = _tudatModule.ephemeris_get_precomputed_time_bounds(target, observer, frame);
        return { start: bounds[0], end: bounds[1] };
    } catch (error) {
        return null;
    }
}

/**
 * List all loaded precomputed ephemerides
 * @returns {string[]} Array of loaded ephemeris keys (e.g., "earth_sun_j2000")
 */
export function listPrecomputedEphemerides() {
    if (!_tudatModule || !_tudatModule.ephemeris_list_precomputed) {
        return [];
    }
    return _tudatModule.ephemeris_list_precomputed();
}

/**
 * Clear all precomputed ephemerides
 */
export function clearPrecomputedEphemerides() {
    if (_tudatModule && _tudatModule.ephemeris_clear_precomputed) {
        _tudatModule.ephemeris_clear_precomputed();
    }
    _precomputedEphemerisLoaded.clear();
}

/**
 * Get body position only (no velocity)
 * @param {string} target - Target body
 * @param {string} observer - Observer body
 * @param {number} epoch - Ephemeris time
 * @param {string} frame - Reference frame
 * @returns {Object|null} {x, y, z} in meters
 */
export function getBodyPosition(target, observer, epoch, frame = 'J2000') {
    const state = getBodyState(target, observer, epoch, frame);
    if (!state) return null;
    return { x: state.x, y: state.y, z: state.z };
}

/**
 * Get gravitational parameter from SPICE
 * @param {string} body - Body name
 * @returns {number|null} GM in m^3/s^2
 */
export function getGM(body) {
    if (!isSpiceTextKernelsLoaded()) {
        console.warn('SPICE text kernels not loaded');
        return null;
    }

    return _tudatModule.interface_spice_get_body_gravitational_parameter(body);
}

/**
 * Get average radius from SPICE
 * @param {string} body - Body name
 * @returns {number|null} Radius in meters
 */
export function getRadius(body) {
    if (!isSpiceTextKernelsLoaded()) {
        console.warn('SPICE text kernels not loaded');
        return null;
    }

    return _tudatModule.interface_spice_get_average_radius(body);
}

/**
 * Convert Julian Date to Ephemeris Time
 * @param {number} jd - Julian Date
 * @returns {number|null} Ephemeris time (seconds since J2000)
 */
export function jdToEt(jd) {
    if (!isSpiceTextKernelsLoaded()) {
        console.warn('SPICE text kernels not loaded');
        return null;
    }

    return _tudatModule.interface_spice_convert_julian_date_to_ephemeris_time(jd);
}

/**
 * Convert Ephemeris Time to Julian Date
 * @param {number} et - Ephemeris time (seconds since J2000)
 * @returns {number|null} Julian Date
 */
export function etToJd(et) {
    if (!isSpiceTextKernelsLoaded()) {
        console.warn('SPICE text kernels not loaded');
        return null;
    }

    return _tudatModule.interface_spice_convert_ephemeris_time_to_julian_date(et);
}

/**
 * Convert Date object or string to Ephemeris Time
 * @param {Date|string} date - JavaScript Date or ISO string
 * @returns {number|null} Ephemeris time
 */
export function dateToEt(date) {
    if (!isSpiceTextKernelsLoaded()) {
        console.warn('SPICE text kernels not loaded');
        return null;
    }

    if (date instanceof Date) {
        // Convert to Julian Date first
        const jd = dateToJd(date);
        return jdToEt(jd);
    } else if (typeof date === 'string') {
        return _tudatModule.interface_spice_convert_date_string_to_ephemeris_time(date);
    }
    throw new Error('Invalid date format');
}

/**
 * Convert Date to Julian Date
 * @param {Date} date - JavaScript Date
 * @returns {number} Julian Date
 */
export function dateToJd(date) {
    const year = date.getUTCFullYear();
    const month = date.getUTCMonth() + 1;
    const day = date.getUTCDate();
    const hour = date.getUTCHours();
    const minute = date.getUTCMinutes();
    const second = date.getUTCSeconds() + date.getUTCMilliseconds() / 1000;

    const dayFraction = (hour + minute / 60 + second / 3600) / 24;

    // Julian Date formula
    let a = Math.floor((14 - month) / 12);
    let y = year + 4800 - a;
    let m = month + 12 * a - 3;

    let jd = day + Math.floor((153 * m + 2) / 5) + 365 * y + Math.floor(y / 4) -
        Math.floor(y / 100) + Math.floor(y / 400) - 32045;

    return jd + dayFraction - 0.5;
}

/**
 * Get planetary positions for solar system visualization
 * @param {number} epoch - Ephemeris time
 * @param {string[]} planets - List of planet names
 * @returns {Object} Map of planet name to position
 */
export function getSolarSystemPositions(epoch, planets = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']) {
    const positions = {};

    for (const planet of planets) {
        const state = getBodyState(planet, 'Sun', epoch);
        if (state) {
            positions[planet] = {
                x: state.x,
                y: state.y,
                z: state.z
            };
        }
    }

    return positions;
}

/**
 * Get Moon position relative to Earth
 * @param {number} epoch - Ephemeris time
 * @returns {Object|null} Position {x, y, z} in meters relative to Earth
 */
export function getMoonPosition(epoch) {
    return getBodyPosition('Moon', 'Earth', epoch);
}

/**
 * Planetary constants (fallback when SPICE not available)
 */
export const PLANETARY_GM = {
    Sun: 1.32712440018e20,
    Mercury: 2.2032e13,
    Venus: 3.24859e14,
    Earth: 3.986004418e14,
    Mars: 4.282837e13,
    Jupiter: 1.26687e17,
    Saturn: 3.7931e16,
    Uranus: 5.7940e15,
    Neptune: 6.8365e15,
    Moon: 4.9028e12
};

export const PLANETARY_RADIUS = {
    Sun: 6.9634e8,
    Mercury: 2.4397e6,
    Venus: 6.0518e6,
    Earth: 6.371e6,
    Mars: 3.3895e6,
    Jupiter: 6.9911e7,
    Saturn: 5.8232e7,
    Uranus: 2.5362e7,
    Neptune: 2.4622e7,
    Moon: 1.7374e6
};

// Semi-major axes in meters (fallback for analytical ephemeris)
export const PLANETARY_SMA = {
    Mercury: 5.79e10,
    Venus: 1.082e11,
    Earth: 1.496e11,
    Mars: 2.279e11,
    Jupiter: 7.786e11,
    Saturn: 1.4335e12,
    Uranus: 2.8725e12,
    Neptune: 4.4951e12,
    Moon: 3.844e8  // from Earth
};

// Orbital periods in seconds
export const PLANETARY_PERIOD = {
    Mercury: 87.97 * 86400,
    Venus: 224.7 * 86400,
    Earth: 365.25 * 86400,
    Mars: 687 * 86400,
    Jupiter: 4333 * 86400,
    Saturn: 10759 * 86400,
    Uranus: 30687 * 86400,
    Neptune: 60190 * 86400,
    Moon: 27.32 * 86400
};

// ============================================================================
// CALCEPH Binary SPK Support (direct binary file reading in WASM)
// ============================================================================

// Track which binary SPK files have been loaded via CALCEPH
let _calcephLoaded = new Set();

/**
 * Check if CALCEPH binary SPK support is available in this build
 * @returns {boolean} True if CALCEPH functions are available
 */
export function isCalcephAvailable() {
    return _tudatModule && typeof _tudatModule.calceph_load_spk === 'function';
}

/**
 * Load a binary SPK file directly using CALCEPH.
 * This allows reading binary SPK/BSP files in WASM without CSPICE.
 *
 * The file must be in the Emscripten virtual filesystem. You can write
 * files to the VFS using Module.FS.writeFile() after fetching them.
 *
 * @param {string} spkPath - Path to the SPK file in the virtual filesystem
 * @param {string} target - Target body name (e.g., "Earth")
 * @param {string} observer - Observer body name (e.g., "Sun")
 * @param {string} frame - Reference frame (default: "J2000")
 * @returns {boolean} True if loaded successfully
 */
export function loadBinarySpk(spkPath, target, observer, frame = 'J2000') {
    if (!isCalcephAvailable()) {
        console.warn('[CALCEPH] CALCEPH not available in this build');
        return false;
    }

    const success = _tudatModule.calceph_load_spk(spkPath, target, observer, frame);
    if (success) {
        const key = `${target.toLowerCase()}_${observer.toLowerCase()}_${frame.toLowerCase()}`;
        _calcephLoaded.add(key);
        console.log(`[CALCEPH] Loaded binary SPK: ${target} relative to ${observer} from ${spkPath}`);
    } else {
        console.error(`[CALCEPH] Failed to load ${spkPath} for ${target} relative to ${observer}`);
    }
    return success;
}

/**
 * Load a binary SPK file using NAIF IDs directly
 * @param {string} spkPath - Path to the SPK file
 * @param {number} targetId - NAIF ID of target body
 * @param {number} observerId - NAIF ID of observer body
 * @param {string} frame - Reference frame
 * @returns {boolean}
 */
export function loadBinarySpkByNaifId(spkPath, targetId, observerId, frame = 'J2000') {
    if (!isCalcephAvailable()) {
        console.warn('[CALCEPH] CALCEPH not available in this build');
        return false;
    }

    return _tudatModule.calceph_load_spk_by_naif_id(spkPath, targetId, observerId, frame);
}

/**
 * Check if CALCEPH ephemeris is available for a target/observer pair
 * @param {string} target - Target body name
 * @param {string} observer - Observer body name
 * @param {string} frame - Reference frame
 * @returns {boolean}
 */
export function isCalcephEphemerisAvailable(target, observer, frame = 'J2000') {
    if (!isCalcephAvailable()) return false;
    return _tudatModule.calceph_is_available(target, observer, frame);
}

/**
 * Get state from a loaded binary SPK file via CALCEPH
 * @param {string} target - Target body
 * @param {string} observer - Observer body
 * @param {number} epoch - Ephemeris time (seconds since J2000)
 * @param {string} frame - Reference frame
 * @returns {Object|null} {x, y, z, vx, vy, vz} in meters and m/s
 */
export function getCalcephState(target, observer, epoch, frame = 'J2000') {
    if (!isCalcephEphemerisAvailable(target, observer, frame)) {
        return null;
    }

    try {
        const state = _tudatModule.calceph_get_state(target, observer, frame, epoch);
        if (!state || state.length < 6) {
            return null;
        }
        return {
            x: state[0],
            y: state[1],
            z: state[2],
            vx: state[3],
            vy: state[4],
            vz: state[5]
        };
    } catch (error) {
        console.error(`[CALCEPH] Error getting state for ${target}:`, error);
        return null;
    }
}

/**
 * Get time bounds for a loaded binary SPK file
 * @param {string} target - Target body
 * @param {string} observer - Observer body
 * @param {string} frame - Reference frame
 * @returns {Object|null} {start, end} epochs in seconds since J2000
 */
export function getCalcephTimeBounds(target, observer, frame = 'J2000') {
    if (!isCalcephAvailable()) return null;

    try {
        const bounds = _tudatModule.calceph_get_time_bounds(target, observer, frame);
        return { start: bounds[0], end: bounds[1] };
    } catch (error) {
        return null;
    }
}

/**
 * List all SPK files loaded via CALCEPH
 * @returns {string[]}
 */
export function listCalcephLoaded() {
    if (!isCalcephAvailable()) return [];
    return _tudatModule.calceph_list_loaded();
}

/**
 * Clear all CALCEPH-loaded SPK files
 */
export function clearCalceph() {
    if (isCalcephAvailable()) {
        _tudatModule.calceph_clear_all();
    }
    _calcephLoaded.clear();
}

/**
 * Fetch a binary SPK file from a URL and load it via CALCEPH.
 * This is a convenience function that handles fetching the file,
 * writing it to the Emscripten virtual filesystem, and loading it.
 *
 * @param {string} url - URL to fetch the SPK file from
 * @param {string} target - Target body name
 * @param {string} observer - Observer body name
 * @param {string} frame - Reference frame
 * @returns {Promise<boolean>} True if loaded successfully
 */
export async function fetchAndLoadBinarySpk(url, target, observer, frame = 'J2000') {
    if (!isCalcephAvailable()) {
        console.warn('[CALCEPH] CALCEPH not available in this build');
        return false;
    }

    if (!_tudatModule.FS) {
        console.error('[CALCEPH] Emscripten FS not available');
        return false;
    }

    try {
        console.log(`[CALCEPH] Fetching binary SPK from ${url}`);
        const response = await fetch(url);
        if (!response.ok) {
            console.error(`[CALCEPH] Failed to fetch ${url}: ${response.status}`);
            return false;
        }

        const arrayBuffer = await response.arrayBuffer();
        const data = new Uint8Array(arrayBuffer);

        // Write to Emscripten virtual filesystem
        const filename = url.split('/').pop() || 'ephemeris.bsp';
        const vfsPath = `/tmp/${filename}`;

        // Create /tmp if it doesn't exist
        try {
            _tudatModule.FS.mkdir('/tmp');
        } catch (e) {
            // Directory may already exist
        }

        _tudatModule.FS.writeFile(vfsPath, data);
        console.log(`[CALCEPH] Wrote ${data.length} bytes to ${vfsPath}`);

        // Load the file via CALCEPH
        return loadBinarySpk(vfsPath, target, observer, frame);
    } catch (error) {
        console.error(`[CALCEPH] Error fetching/loading ${url}:`, error);
        return false;
    }
}

/**
 * Fetch a binary SPK file and load ephemeris for multiple bodies.
 * This is useful when a single SPK file (like de430.bsp) contains data for multiple planets.
 *
 * @param {string} url - URL to fetch the SPK file from
 * @param {Array<{target: string, observer: string}>} bodies - Array of target/observer pairs to load
 * @param {string} frame - Reference frame (default: 'J2000')
 * @returns {Promise<{success: boolean, loaded: string[], failed: string[]}>}
 */
export async function fetchAndLoadBinarySpkForBodies(url, bodies, frame = 'J2000') {
    if (!isCalcephAvailable()) {
        console.warn('[CALCEPH] CALCEPH not available in this build');
        return { success: false, loaded: [], failed: bodies.map(b => `${b.target}/${b.observer}`) };
    }

    if (!_tudatModule.FS) {
        console.error('[CALCEPH] Emscripten FS not available');
        return { success: false, loaded: [], failed: bodies.map(b => `${b.target}/${b.observer}`) };
    }

    try {
        console.log(`[CALCEPH] Fetching binary SPK from ${url}`);
        const response = await fetch(url);
        if (!response.ok) {
            console.error(`[CALCEPH] Failed to fetch ${url}: ${response.status}`);
            return { success: false, loaded: [], failed: bodies.map(b => `${b.target}/${b.observer}`) };
        }

        const arrayBuffer = await response.arrayBuffer();
        const data = new Uint8Array(arrayBuffer);

        // Write to Emscripten virtual filesystem
        const filename = url.split('/').pop() || 'ephemeris.bsp';
        const vfsPath = `/tmp/${filename}`;

        // Create /tmp if it doesn't exist
        try {
            _tudatModule.FS.mkdir('/tmp');
        } catch (e) {
            // Directory may already exist
        }

        _tudatModule.FS.writeFile(vfsPath, data);
        console.log(`[CALCEPH] Wrote ${data.length} bytes to ${vfsPath}`);

        // Load ephemeris for each body pair
        const loaded = [];
        const failed = [];

        for (const { target, observer } of bodies) {
            const success = loadBinarySpk(vfsPath, target, observer, frame);
            if (success) {
                loaded.push(`${target}/${observer}`);
            } else {
                failed.push(`${target}/${observer}`);
            }
        }

        return {
            success: loaded.length > 0,
            loaded,
            failed
        };
    } catch (error) {
        console.error(`[CALCEPH] Error fetching/loading ${url}:`, error);
        return { success: false, loaded: [], failed: bodies.map(b => `${b.target}/${b.observer}`) };
    }
}

/**
 * Load planetary ephemeris from a JPL DE kernel.
 * Convenience function that loads common planet ephemerides from a DE kernel like de430.bsp.
 *
 * @param {string} url - URL to the DE kernel (e.g., './data/de430.bsp')
 * @param {string[]} planets - Array of planet names to load (default: inner + outer planets)
 * @param {string} observer - Observer body (default: 'Sun')
 * @param {string} frame - Reference frame (default: 'J2000')
 * @returns {Promise<{success: boolean, loaded: string[], failed: string[]}>}
 */
export async function loadPlanetaryEphemeris(url, planets = null, observer = 'Sun', frame = 'J2000') {
    // Default to loading all major planets
    const defaultPlanets = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'];
    const planetsToLoad = planets || defaultPlanets;

    const bodies = planetsToLoad.map(target => ({ target, observer }));
    return fetchAndLoadBinarySpkForBodies(url, bodies, frame);
}

/**
 * Convert body name to NAIF ID
 * @param {string} name - Body name
 * @returns {number} NAIF ID or -1 if unknown
 */
export function bodyNameToNaifId(name) {
    if (!isCalcephAvailable()) {
        // Fallback mapping for common bodies
        const mapping = {
            'ssb': 0, 'solar system barycenter': 0,
            'sun': 10,
            'mercury barycenter': 1, 'mercury': 199,
            'venus barycenter': 2, 'venus': 299,
            'earth-moon barycenter': 3, 'emb': 3, 'earth': 399,
            'mars barycenter': 4, 'mars': 499,
            'jupiter barycenter': 5, 'jupiter': 599,
            'saturn barycenter': 6, 'saturn': 699,
            'uranus barycenter': 7, 'uranus': 799,
            'neptune barycenter': 8, 'neptune': 899,
            'pluto barycenter': 9, 'pluto': 999,
            'moon': 301, 'luna': 301
        };
        return mapping[name.toLowerCase()] || -1;
    }
    return _tudatModule.calceph_body_name_to_naif_id(name);
}

/**
 * Convert NAIF ID to body name
 * @param {number} naifId - NAIF ID
 * @returns {string} Body name
 */
export function naifIdToBodyName(naifId) {
    if (!isCalcephAvailable()) {
        const mapping = {
            0: 'SSB', 10: 'Sun',
            1: 'Mercury Barycenter', 199: 'Mercury',
            2: 'Venus Barycenter', 299: 'Venus',
            3: 'Earth-Moon Barycenter', 399: 'Earth',
            4: 'Mars Barycenter', 499: 'Mars',
            5: 'Jupiter Barycenter', 599: 'Jupiter',
            6: 'Saturn Barycenter', 699: 'Saturn',
            7: 'Uranus Barycenter', 799: 'Uranus',
            8: 'Neptune Barycenter', 899: 'Neptune',
            9: 'Pluto Barycenter', 999: 'Pluto',
            301: 'Moon'
        };
        return mapping[naifId] || `Body${naifId}`;
    }
    return _tudatModule.calceph_naif_id_to_body_name(naifId);
}
