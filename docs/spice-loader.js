/**
 * SPICE Kernel Loader for Browser
 *
 * This utility fetches real SPICE kernel files and mounts them to
 * Emscripten's virtual filesystem, enabling accurate ephemeris queries.
 *
 * Usage:
 *   import { SpiceKernelLoader } from './spice-loader.js';
 *
 *   const loader = new SpiceKernelLoader(tudatModule);
 *   await loader.loadStandardKernels();
 *
 *   // Now you can use SPICE functions
 *   const state = tudat.interface_.spice.get_body_cartesian_state_at_epoch(...);
 */

export class SpiceKernelLoader {
    constructor(tudatModule) {
        this.module = tudatModule;
        this.loadedKernels = new Set();
        this.kernelBasePath = '/spice_kernels';
        this._directoryCreated = false;

        // Standard kernel URLs - these can be hosted locally or fetched from NAIF
        // For local development, copy kernels to tests/wasm/web/data/spice/
        // NOTE: Order matters! Leap seconds kernel must be loaded first.
        this.standardKernels = {
            // Leap seconds kernel (required for time conversions) - MUST BE FIRST
            leapSeconds: {
                name: 'naif0012.tls',
                localPath: './data/spice/naif0012.tls',
                naifUrl: 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls',
                size: 5409  // bytes
            },
            // Planetary constants
            planetaryConstants: {
                name: 'pck00010.tpc',
                localPath: './data/spice/pck00010.tpc',
                naifUrl: 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc',
                size: 71311
            },
            // DE430 ephemeris (small version for common date range)
            de430Small: {
                name: 'de430_mar097_small.bsp',
                localPath: './data/spice/de430_mar097_small.bsp',
                // This is a reduced-size version of DE430 (~12MB)
                size: 12276736
            },
            // Frame kernel for rotating frames
            frameKernel: {
                name: 'moon_assoc_pa.tf',
                localPath: './data/spice/moon_assoc_pa.tf',
                size: 8749
            }
        };

        // Don't create directory in constructor - wait until we actually need it
        // This avoids issues if FS isn't fully ready yet
    }

    /**
     * Create directory in Emscripten's virtual filesystem
     */
    _ensureDirectory(path) {
        if (this._directoryCreated) {
            return true;
        }

        // Check if FS is available
        if (!this.module || !this.module.FS || typeof this.module.FS.mkdir !== 'function') {
            console.error('[SPICE] FS not available - cannot create directory');
            return false;
        }

        const parts = path.split('/').filter(p => p);
        let currentPath = '';

        for (const part of parts) {
            currentPath += '/' + part;
            try {
                this.module.FS.mkdir(currentPath);
            } catch (e) {
                // Directory may already exist (EEXIST) - that's OK
                if (e.code !== 'EEXIST' && e.errno !== 20) {
                    console.warn(`[SPICE] Could not create directory ${currentPath}:`, e.message || e);
                }
            }
        }

        this._directoryCreated = true;
        return true;
    }

    /**
     * Fetch a kernel file and mount it to the virtual filesystem
     */
    async loadKernel(kernelConfig, options = {}) {
        const { name, localPath, naifUrl } = kernelConfig;
        const { forceReload = false, preferLocal = true } = options;

        // Check if already loaded
        if (this.loadedKernels.has(name) && !forceReload) {
            console.log(`[SPICE] Kernel ${name} already loaded`);
            return true;
        }

        // Ensure FS is available before proceeding
        if (!this.module || !this.module.FS || typeof this.module.FS.writeFile !== 'function') {
            console.error(`[SPICE] Cannot load kernel ${name}: Emscripten FS not available`);
            return false;
        }

        // Ensure directory exists (lazy creation)
        if (!this._ensureDirectory(this.kernelBasePath)) {
            console.error(`[SPICE] Cannot load kernel ${name}: Failed to create directory`);
            return false;
        }

        const virtualPath = `${this.kernelBasePath}/${name}`;

        try {
            // Try to fetch the kernel data
            let data;
            const urls = preferLocal ? [localPath, naifUrl] : [naifUrl, localPath];

            for (const url of urls.filter(Boolean)) {
                try {
                    console.log(`[SPICE] Fetching ${name} from ${url}...`);
                    const response = await fetch(url);

                    if (!response.ok) {
                        console.warn(`[SPICE] Failed to fetch from ${url}: ${response.status}`);
                        continue;
                    }

                    data = new Uint8Array(await response.arrayBuffer());
                    console.log(`[SPICE] Downloaded ${name}: ${data.length} bytes`);
                    break;
                } catch (e) {
                    console.warn(`[SPICE] Could not fetch from ${url}:`, e.message);
                }
            }

            if (!data) {
                throw new Error(`Could not fetch kernel ${name} from any source`);
            }

            // Write to virtual filesystem
            this.module.FS.writeFile(virtualPath, data);
            console.log(`[SPICE] Written ${name} to ${virtualPath}`);

            // Load the kernel via SPICE (use flat function name from Embind)
            this.module.interface_spice_load_kernel(virtualPath);
            this.loadedKernels.add(name);

            const totalKernels = this.module.interface_spice_get_total_count_of_kernels_loaded();
            console.log(`[SPICE] Kernel ${name} loaded. Total kernels: ${totalKernels}`);

            return true;
        } catch (error) {
            console.error(`[SPICE] Failed to load kernel ${name}:`, error);
            return false;
        }
    }

    /**
     * Load all standard kernels needed for basic ephemeris queries
     * Kernels are loaded in order - leap seconds kernel must be first!
     */
    async loadStandardKernels(options = {}) {
        const { onProgress } = options;

        // Get kernels in explicit order - leap seconds MUST be first
        const kernelOrder = ['leapSeconds', 'planetaryConstants', 'de430Small', 'frameKernel'];
        const kernels = kernelOrder
            .filter(key => this.standardKernels[key])
            .map(key => this.standardKernels[key]);

        let loaded = 0;

        console.log(`[SPICE] Loading ${kernels.length} standard kernels...`);

        for (let i = 0; i < kernels.length; i++) {
            const kernel = kernels[i];
            const success = await this.loadKernel(kernel, options);

            if (success) {
                loaded++;
            } else if (i === 0) {
                // If leap seconds kernel fails, we cannot continue with other kernels
                // as SPICE requires it for time conversions used by other kernels
                console.error('[SPICE] Leap seconds kernel failed to load - cannot load remaining kernels');
                break;
            }

            if (onProgress) {
                onProgress({
                    loaded,
                    total: kernels.length,
                    currentKernel: kernel.name,
                    success
                });
            }
        }

        console.log(`[SPICE] Loaded ${loaded}/${kernels.length} standard kernels`);
        return loaded === kernels.length;
    }

    /**
     * Load a custom kernel from a URL
     */
    async loadCustomKernel(name, url) {
        return this.loadKernel({ name, localPath: url, naifUrl: url });
    }

    /**
     * Clear all loaded kernels
     */
    clearKernels() {
        this.module.interface_spice_clear_kernels();
        this.loadedKernels.clear();
        console.log('[SPICE] All kernels cleared');
    }

    /**
     * Get count of loaded kernels
     */
    getLoadedCount() {
        return this.module.interface_spice_get_total_count_of_kernels_loaded();
    }

    /**
     * Check if kernels are loaded and ready
     */
    isReady() {
        return this.loadedKernels.size > 0;
    }

    /**
     * Test SPICE functionality with a simple query
     */
    testSpice() {
        try {
            // Try a simple time conversion
            const jd = 2451545.0;  // J2000 epoch
            const et = this.module.interface_spice_convert_julian_date_to_ephemeris_time(jd);
            console.log(`[SPICE] Test: JD ${jd} -> ET ${et}`);

            // Try to get Earth state (requires ephemeris kernel)
            if (this.loadedKernels.has('de430_mar097_small.bsp')) {
                const state = this.module.interface_spice_get_body_cartesian_state_at_epoch(
                    'Earth',
                    'Sun',
                    'J2000',
                    'NONE',
                    et
                );
                console.log('[SPICE] Test: Earth position at J2000:', state);
            }

            return true;
        } catch (error) {
            console.error('[SPICE] Test failed:', error);
            return false;
        }
    }
}

/**
 * Utility function to get planetary state using SPICE
 */
export function getBodyState(module, body, observer, epoch, frame = 'J2000') {
    const state = module.interface_.spice.get_body_cartesian_state_at_epoch(
        body,
        observer,
        frame,
        'NONE',
        epoch
    );

    // Convert from Vector6d to plain object
    return {
        x: state.get(0),
        y: state.get(1),
        z: state.get(2),
        vx: state.get(3),
        vy: state.get(4),
        vz: state.get(5)
    };
}

/**
 * Convert date to ephemeris time
 */
export function dateToEt(module, date) {
    if (typeof date === 'string') {
        return module.interface_.spice.convert_date_string_to_ephemeris_time(date);
    } else if (date instanceof Date) {
        // Convert JS Date to ISO string
        const isoString = date.toISOString().replace('Z', '');
        return module.interface_.spice.convert_date_string_to_ephemeris_time(isoString);
    } else if (typeof date === 'number') {
        // Assume Julian date
        return module.interface_.spice.convert_julian_date_to_ephemeris_time(date);
    }
    throw new Error('Invalid date format');
}

/**
 * Get gravitational parameter from SPICE
 */
export function getGM(module, body) {
    return module.interface_.spice.get_body_gravitational_parameter(body);
}

/**
 * Get average radius from SPICE
 */
export function getRadius(module, body) {
    return module.interface_.spice.get_average_radius(body);
}

export default SpiceKernelLoader;
