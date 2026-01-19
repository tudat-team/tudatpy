/**
 * SPICE Utilities for Visualizations
 *
 * Provides easy-to-use functions for accessing SPICE ephemeris data
 * in browser-based visualizations.
 */

// Global reference to the Tudat module and SPICE state
let _tudatModule = null;
let _spiceReady = false;

/**
 * Initialize SPICE support with a loaded Tudat module
 * @param {Object} tudatModule - The loaded tudatpy WASM module
 * @param {boolean} ready - Whether SPICE kernels are loaded
 */
export function initSpice(tudatModule, ready = false) {
    _tudatModule = tudatModule;
    _spiceReady = ready;
}

/**
 * Check if SPICE is ready for queries
 * @returns {boolean}
 */
export function isSpiceReady() {
    return _spiceReady && _tudatModule !== null;
}

/**
 * Get planetary state at a given epoch
 * @param {string} target - Target body (e.g., 'Earth', 'Mars', 'Jupiter')
 * @param {string} observer - Observer body (e.g., 'Sun', 'Earth')
 * @param {number} epoch - Ephemeris time (seconds since J2000)
 * @param {string} frame - Reference frame (default: 'J2000')
 * @returns {Object|null} {x, y, z, vx, vy, vz} in meters and m/s
 */
export function getBodyState(target, observer, epoch, frame = 'J2000') {
    if (!isSpiceReady()) {
        console.warn('SPICE not initialized');
        return null;
    }

    try {
        const state = _tudatModule.interface_.spice.get_body_cartesian_state_at_epoch(
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
        return null;
    }
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
    if (!isSpiceReady()) {
        console.warn('SPICE not initialized');
        return null;
    }

    try {
        return _tudatModule.interface_.spice.get_body_gravitational_parameter(body);
    } catch (error) {
        console.error(`Failed to get GM for ${body}:`, error);
        return null;
    }
}

/**
 * Get average radius from SPICE
 * @param {string} body - Body name
 * @returns {number|null} Radius in meters
 */
export function getRadius(body) {
    if (!isSpiceReady()) {
        console.warn('SPICE not initialized');
        return null;
    }

    try {
        return _tudatModule.interface_.spice.get_average_radius(body);
    } catch (error) {
        console.error(`Failed to get radius for ${body}:`, error);
        return null;
    }
}

/**
 * Convert Julian Date to Ephemeris Time
 * @param {number} jd - Julian Date
 * @returns {number|null} Ephemeris time (seconds since J2000)
 */
export function jdToEt(jd) {
    if (!isSpiceReady()) {
        // Fallback to analytical conversion
        const J2000_JD = 2451545.0;
        return (jd - J2000_JD) * 86400.0;
    }

    try {
        return _tudatModule.interface_.spice.convert_julian_date_to_ephemeris_time(jd);
    } catch (error) {
        // Fallback
        const J2000_JD = 2451545.0;
        return (jd - J2000_JD) * 86400.0;
    }
}

/**
 * Convert Ephemeris Time to Julian Date
 * @param {number} et - Ephemeris time (seconds since J2000)
 * @returns {number} Julian Date
 */
export function etToJd(et) {
    if (!isSpiceReady()) {
        // Fallback to analytical conversion
        const J2000_JD = 2451545.0;
        return J2000_JD + et / 86400.0;
    }

    try {
        return _tudatModule.interface_.spice.convert_ephemeris_time_to_julian_date(et);
    } catch (error) {
        // Fallback
        const J2000_JD = 2451545.0;
        return J2000_JD + et / 86400.0;
    }
}

/**
 * Convert Date object or string to Ephemeris Time
 * @param {Date|string} date - JavaScript Date or ISO string
 * @returns {number} Ephemeris time
 */
export function dateToEt(date) {
    if (date instanceof Date) {
        // Convert to Julian Date first
        const jd = dateToJd(date);
        return jdToEt(jd);
    } else if (typeof date === 'string') {
        if (isSpiceReady()) {
            try {
                return _tudatModule.interface_.spice.convert_date_string_to_ephemeris_time(date);
            } catch (error) {
                // Parse manually
                const d = new Date(date);
                return dateToEt(d);
            }
        }
        // Parse manually
        const d = new Date(date);
        return dateToEt(d);
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
