/**
 * Tudat WASM TypeScript Definitions
 *
 * This file provides TypeScript type definitions for the Tudat WASM module.
 * The API mirrors the Python tudatpy library.
 *
 * @packageDocumentation
 */

// =============================================================================
// Vector Types
// =============================================================================

/**
 * 3-dimensional vector wrapper for Eigen::Vector3d
 */
export interface Vector3d {
  /** Get the value at index i (0-2) */
  get(i: number): number;
  /** Set the value at index i (0-2) */
  set(i: number, value: number): void;
  /** Get the size (always 3) */
  size(): number;
  /** Convert to JavaScript array */
  toArray(): [number, number, number];
  /** Free the underlying memory */
  delete(): void;
}

/**
 * 6-dimensional vector wrapper for Eigen::Vector6d (state vector)
 */
export interface Vector6d {
  /** Get the value at index i (0-5) */
  get(i: number): number;
  /** Set the value at index i (0-5) */
  set(i: number, value: number): void;
  /** Get the size (always 6) */
  size(): number;
  /** Convert to JavaScript array */
  toArray(): [number, number, number, number, number, number];
  /** Free the underlying memory */
  delete(): void;
}

/**
 * Dynamic-size vector wrapper for Eigen::VectorXd
 */
export interface VectorXd {
  /** Get the value at index i */
  get(i: number): number;
  /** Set the value at index i */
  set(i: number, value: number): void;
  /** Get the size */
  size(): number;
  /** Resize the vector */
  resize(newSize: number): void;
  /** Convert to JavaScript array */
  toArray(): number[];
  /** Free the underlying memory */
  delete(): void;
}

/**
 * 3x3 matrix wrapper for Eigen::Matrix3d
 */
export interface Matrix3d {
  /** Get the value at row i, column j */
  get(i: number, j: number): number;
  /** Set the value at row i, column j */
  set(i: number, j: number, value: number): void;
  /** Get number of rows (always 3) */
  rows(): number;
  /** Get number of columns (always 3) */
  cols(): number;
  /** Free the underlying memory */
  delete(): void;
}

// =============================================================================
// Constants Module
// =============================================================================

/**
 * Physical and mathematical constants
 * @namespace
 */
export namespace constants {
  /** Gravitational constant [m^3 kg^-1 s^-2] */
  export const GRAVITATIONAL_CONSTANT: number;
  /** Speed of light in vacuum [m/s] */
  export const SPEED_OF_LIGHT: number;
  /** Astronomical Unit [m] */
  export const ASTRONOMICAL_UNIT: number;
  /** Julian day in seconds [s] */
  export const JULIAN_DAY: number;
  /** Julian year in seconds [s] */
  export const JULIAN_YEAR: number;
  /** Sidereal day in seconds [s] */
  export const SIDEREAL_DAY: number;
  /** Sidereal year in seconds [s] */
  export const SIDEREAL_YEAR: number;
  /** Standard gravitational acceleration [m/s^2] */
  export const SEA_LEVEL_GRAVITATIONAL_ACCELERATION: number;
  /** Earth equatorial radius [m] */
  export const EARTH_EQUATORIAL_RADIUS: number;
  /** Earth flattening factor */
  export const EARTH_FLATTENING_FACTOR: number;
  /** Earth gravitational parameter [m^3/s^2] */
  export const EARTH_GRAVITATIONAL_PARAMETER: number;
  /** Sun gravitational parameter [m^3/s^2] */
  export const SUN_GRAVITATIONAL_PARAMETER: number;
  /** Moon gravitational parameter [m^3/s^2] */
  export const MOON_GRAVITATIONAL_PARAMETER: number;
  /** Pi */
  export const PI: number;
}

// =============================================================================
// Math Module
// =============================================================================

/**
 * Mathematical utilities
 * @namespace
 */
export namespace math {
  /**
   * Interpolation utilities
   */
  export namespace interpolators {
    /** Available lookup schemes for interpolation */
    export enum AvailableLookupScheme {
      huntingAlgorithm,
      binarySearch
    }

    /** Boundary handling for interpolation */
    export enum BoundaryInterpolationType {
      throw_exception_at_boundary,
      use_boundary_value,
      use_boundary_value_with_warning,
      extrapolate_at_boundary,
      extrapolate_at_boundary_with_warning
    }

    /** Lagrange interpolator types */
    export enum AvailableLagrangeInterpolatorTypes {
      cubic,
      sixth_order,
      eighth_order
    }

    /** Base class for interpolator settings */
    export interface InterpolatorSettings {
      delete(): void;
    }

    /** Create linear interpolation settings */
    export function linear_interpolation(
      lookupScheme?: AvailableLookupScheme,
      boundaryHandling?: BoundaryInterpolationType
    ): InterpolatorSettings;

    /** Create cubic spline interpolation settings */
    export function cubic_spline_interpolation(
      boundaryHandling?: BoundaryInterpolationType
    ): InterpolatorSettings;

    /** Create Lagrange interpolation settings */
    export function lagrange_interpolation(
      numberOfPoints: number,
      boundaryHandling?: BoundaryInterpolationType
    ): InterpolatorSettings;

    /** Create piecewise constant interpolation settings */
    export function piecewise_constant_interpolation(
      lookupScheme?: AvailableLookupScheme,
      boundaryHandling?: BoundaryInterpolationType
    ): InterpolatorSettings;

    /** Create Hermite spline interpolation settings */
    export function hermite_spline_interpolation(
      boundaryHandling?: BoundaryInterpolationType
    ): InterpolatorSettings;
  }

  /**
   * Numerical integration utilities
   */
  export namespace numerical_integrators {
    /** Available integrator types */
    export enum AvailableIntegrators {
      rungeKutta4,
      euler,
      rungeKuttaVariableStepSize,
      adamsBashforthMoulton,
      bulirschStoer
    }

    /** Runge-Kutta coefficient sets */
    export enum CoefficientSets {
      rungeKutta4Classic,
      rungeKuttaFehlberg45,
      rungeKuttaFehlberg56,
      rungeKuttaFehlberg78,
      rungeKutta87DormandPrince
    }

    /** Base class for integrator settings */
    export interface IntegratorSettings {
      delete(): void;
    }
  }

  /**
   * Root finding utilities
   */
  export namespace root_finders {
    /** Available root finder types */
    export enum AvailableRootFinders {
      bisection,
      newtonRaphson,
      secant,
      halley
    }

    /** Base class for root finder settings */
    export interface RootFinderSettings {
      delete(): void;
    }
  }
}

// =============================================================================
// Astro Module
// =============================================================================

/**
 * Astrodynamics utilities
 * @namespace
 */
export namespace astro {
  /**
   * Orbital element conversion functions
   */
  export namespace element_conversion {
    /** Keplerian element indices */
    export enum KeplerianElementIndices {
      semi_major_axis_index,
      eccentricity_index,
      inclination_index,
      argument_of_periapsis_index,
      longitude_of_ascending_node_index,
      true_anomaly_index
    }

    /** Cartesian element indices */
    export enum CartesianElementIndices {
      x_cartesian_position_index,
      y_cartesian_position_index,
      z_cartesian_position_index,
      x_cartesian_velocity_index,
      y_cartesian_velocity_index,
      z_cartesian_velocity_index
    }

    /** Spherical state indices */
    export enum SphericalStateIndices {
      radius_index,
      latitude_index,
      longitude_index,
      speed_index,
      flight_path_angle_index,
      heading_angle_index
    }

    /**
     * Convert Keplerian elements to Cartesian state
     * @param keplerianElements - [a, e, i, ω, Ω, θ]
     * @param gravitationalParameter - Central body GM [m³/s²]
     * @returns Cartesian state [x, y, z, vx, vy, vz] [m, m/s]
     */
    export function keplerian_to_cartesian(
      keplerianElements: Vector6d,
      gravitationalParameter: number
    ): Vector6d;

    /**
     * Convert Cartesian state to Keplerian elements
     * @param cartesianState - [x, y, z, vx, vy, vz] [m, m/s]
     * @param gravitationalParameter - Central body GM [m³/s²]
     * @returns Keplerian elements [a, e, i, ω, Ω, θ]
     */
    export function cartesian_to_keplerian(
      cartesianState: Vector6d,
      gravitationalParameter: number
    ): Vector6d;

    /**
     * Convert Cartesian state to spherical state
     * @param cartesianState - [x, y, z, vx, vy, vz]
     * @returns Spherical state [r, lat, lon, v, γ, χ]
     */
    export function cartesian_to_spherical(cartesianState: Vector6d): Vector6d;

    /**
     * Convert spherical state to Cartesian state
     * @param sphericalState - [r, lat, lon, v, γ, χ]
     * @returns Cartesian state [x, y, z, vx, vy, vz]
     */
    export function spherical_to_cartesian(sphericalState: Vector6d): Vector6d;

    /**
     * Convert mean anomaly to eccentric anomaly
     * @param eccentricity - Orbital eccentricity
     * @param meanAnomaly - Mean anomaly [rad]
     * @returns Eccentric anomaly [rad]
     */
    export function mean_to_eccentric_anomaly(
      eccentricity: number,
      meanAnomaly: number
    ): number;

    /**
     * Convert eccentric anomaly to true anomaly
     * @param eccentricity - Orbital eccentricity
     * @param eccentricAnomaly - Eccentric anomaly [rad]
     * @returns True anomaly [rad]
     */
    export function eccentric_to_true_anomaly(
      eccentricity: number,
      eccentricAnomaly: number
    ): number;

    /**
     * Convert true anomaly to eccentric anomaly
     * @param eccentricity - Orbital eccentricity
     * @param trueAnomaly - True anomaly [rad]
     * @returns Eccentric anomaly [rad]
     */
    export function true_to_eccentric_anomaly(
      eccentricity: number,
      trueAnomaly: number
    ): number;

    /**
     * Convert eccentric anomaly to mean anomaly
     * @param eccentricity - Orbital eccentricity
     * @param eccentricAnomaly - Eccentric anomaly [rad]
     * @returns Mean anomaly [rad]
     */
    export function eccentric_to_mean_anomaly(
      eccentricity: number,
      eccentricAnomaly: number
    ): number;
  }

  /**
   * Reference frame conversion utilities
   */
  export namespace frame_conversion {
    /**
     * Get rotation matrix from inertial to RSW frame
     * @param state - Cartesian state in inertial frame
     * @returns 3x3 rotation matrix
     */
    export function inertial_to_rsw_rotation_matrix(state: Vector6d): Matrix3d;

    /**
     * Get rotation matrix from inertial to TNW frame
     * @param state - Cartesian state in inertial frame
     * @returns 3x3 rotation matrix
     */
    export function inertial_to_tnw_rotation_matrix(state: Vector6d): Matrix3d;
  }

  /**
   * Time representation utilities
   */
  export namespace time_representation {
    /**
     * DateTime class for calendar date/time operations
     */
    export class DateTime {
      /**
       * Create a DateTime from calendar components
       * @param year - Year (e.g., 2024)
       * @param month - Month (1-12)
       * @param day - Day (1-31)
       * @param hour - Hour (0-23)
       * @param minute - Minute (0-59)
       * @param second - Second (0-59.999...)
       */
      constructor(
        year: number,
        month: number,
        day: number,
        hour?: number,
        minute?: number,
        second?: number
      );

      /**
       * Get seconds since J2000 epoch (2000-01-01 12:00:00 TDB)
       * @returns Seconds since J2000
       */
      epoch(): number;

      /** Get the year */
      year(): number;
      /** Get the month (1-12) */
      month(): number;
      /** Get the day of month (1-31) */
      day(): number;
      /** Get the hour (0-23) */
      hour(): number;
      /** Get the minute (0-59) */
      minute(): number;
      /** Get the second (0-59.999...) */
      second(): number;

      /** Free the underlying memory */
      delete(): void;
    }

    /**
     * Get Julian date from seconds since J2000
     * @param secondsSinceJ2000 - Seconds since J2000 epoch
     * @returns Julian date
     */
    export function julian_day_from_epoch(secondsSinceJ2000: number): number;

    /**
     * Get seconds since J2000 from Julian date
     * @param julianDay - Julian date
     * @returns Seconds since J2000
     */
    export function epoch_from_julian_day(julianDay: number): number;
  }

  /**
   * Two-body dynamics utilities
   */
  export namespace two_body_dynamics {
    /**
     * Compute Keplerian orbital period
     * @param semiMajorAxis - Semi-major axis [m]
     * @param gravitationalParameter - Central body GM [m³/s²]
     * @returns Orbital period [s]
     */
    export function compute_kepler_orbit_period(
      semiMajorAxis: number,
      gravitationalParameter: number
    ): number;

    /**
     * Propagate Keplerian orbit
     * @param initialState - Initial Keplerian elements
     * @param propagationTime - Time to propagate [s]
     * @param gravitationalParameter - Central body GM [m³/s²]
     * @returns Final Keplerian elements
     */
    export function propagate_kepler_orbit(
      initialState: Vector6d,
      propagationTime: number,
      gravitationalParameter: number
    ): Vector6d;
  }
}

// =============================================================================
// Dynamics Module
// =============================================================================

/**
 * Dynamics simulation and propagation
 * @namespace
 */
export namespace dynamics {
  /**
   * Environment models (Body, SystemOfBodies)
   */
  export namespace environment {
    /**
     * A celestial body with physical properties
     */
    export interface Body {
      /** Get the body's state at a given time */
      state_in_base_frame_from_ephemeris(time: number): Vector6d;
      /** Get the body's position at a given time */
      position(time: number): Vector3d;
      /** Get the body's velocity at a given time */
      velocity(time: number): Vector3d;
      /** Get the body's gravitational parameter */
      gravitational_parameter: number;
      /** Free the underlying memory */
      delete(): void;
    }

    /**
     * Collection of bodies forming the simulation environment
     */
    export interface SystemOfBodies {
      /** Get a body by name */
      get(bodyName: string): Body;
      /** Get all body names */
      bodies(): string[];
      /** Add a new body */
      create_empty_body(bodyName: string): void;
      /** Free the underlying memory */
      delete(): void;
    }
  }

  /**
   * Environment setup functions
   */
  export namespace environment_setup {
    /**
     * Settings for creating body models
     */
    export interface BodySettings {
      delete(): void;
    }

    /**
     * Settings for all bodies in the simulation
     */
    export interface SystemOfBodySettings {
      /** Get settings for a specific body */
      get(bodyName: string): BodySettings;
      /** Add empty settings for a new body */
      add_empty_settings(bodyName: string): void;
      delete(): void;
    }

    /**
     * Get default settings for standard solar system bodies
     * @param bodies - Names of bodies to create settings for
     * @param frameOrigin - Origin of the reference frame
     * @param frameOrientation - Orientation of the reference frame
     * @returns Settings object for all requested bodies
     */
    export function get_default_body_settings(
      bodies: string[],
      frameOrigin: string,
      frameOrientation?: string
    ): SystemOfBodySettings;

    /**
     * Create bodies from settings
     * @param bodySettings - Settings for all bodies
     * @returns System of bodies ready for simulation
     */
    export function create_system_of_bodies(
      bodySettings: SystemOfBodySettings
    ): environment.SystemOfBodies;

    /**
     * Atmosphere model types
     */
    export namespace atmosphere {
      export enum AtmosphereModelType {
        exponential_atmosphere,
        tabulated_atmosphere,
        nrlmsise00
      }

      export interface AtmosphereSettings {
        delete(): void;
      }

      export function exponential_atmosphere(
        surfaceDensity: number,
        scaleHeight: number
      ): AtmosphereSettings;

      export function nrlmsise00(): AtmosphereSettings;
    }

    /**
     * Gravity field model types
     */
    export namespace gravity_field {
      export enum GravityFieldType {
        central,
        central_spice,
        spherical_harmonic
      }

      export interface GravityFieldSettings {
        delete(): void;
      }

      export function central(
        gravitationalParameter: number
      ): GravityFieldSettings;

      export function central_spice(): GravityFieldSettings;

      export function spherical_harmonic(
        gravitationalParameter: number,
        referenceRadius: number,
        cosineCoefficients: number[][],
        sineCoefficients: number[][],
        associatedReferenceFrame: string
      ): GravityFieldSettings;
    }

    /**
     * Ephemeris model types
     */
    export namespace ephemeris {
      export enum EphemerisType {
        approximate_jpl,
        direct_spice,
        tabulated,
        constant
      }

      export interface EphemerisSettings {
        delete(): void;
      }

      export function direct_spice(
        frameOrigin?: string,
        frameOrientation?: string
      ): EphemerisSettings;

      export function tabulated(
        bodyStateHistory: Map<number, Vector6d>,
        frameOrigin: string,
        frameOrientation: string
      ): EphemerisSettings;

      export function constant(
        constantState: Vector6d,
        frameOrigin: string,
        frameOrientation: string
      ): EphemerisSettings;
    }
  }

  /**
   * Propagation setup
   */
  export namespace propagation_setup {
    /**
     * Numerical integrator settings
     */
    export namespace integrator {
      export enum AvailableIntegrators {
        rungeKutta4,
        euler,
        rungeKuttaVariableStepSize,
        adamsBashforthMoulton,
        bulirschStoer
      }

      export enum CoefficientSets {
        rungeKutta4Classic,
        rungeKuttaFehlberg45,
        rungeKuttaFehlberg56,
        rungeKuttaFehlberg78,
        rungeKutta87DormandPrince
      }

      export interface IntegratorSettings {
        delete(): void;
      }

      /**
       * Create fixed step size Runge-Kutta integrator
       * @param timeStep - Integration step size [s]
       * @param coefficientSet - RK coefficient set to use
       * @returns Integrator settings
       */
      export function runge_kutta_fixed_step_size(
        timeStep: number,
        coefficientSet?: CoefficientSets
      ): IntegratorSettings;

      /**
       * Create variable step size Runge-Kutta integrator
       * @param initialTimeStep - Initial step size [s]
       * @param coefficientSet - RK coefficient set to use
       * @param minimumStepSize - Minimum allowed step [s]
       * @param maximumStepSize - Maximum allowed step [s]
       * @param relativeErrorTolerance - Relative error tolerance
       * @param absoluteErrorTolerance - Absolute error tolerance
       * @returns Integrator settings
       */
      export function runge_kutta_variable_step_size(
        initialTimeStep: number,
        coefficientSet: CoefficientSets,
        minimumStepSize: number,
        maximumStepSize: number,
        relativeErrorTolerance?: number,
        absoluteErrorTolerance?: number
      ): IntegratorSettings;

      /**
       * Create Bulirsch-Stoer integrator
       * @param initialTimeStep - Initial step size [s]
       * @param extrapolationSequence - Extrapolation sequence type
       * @param maximumNumberOfSteps - Max steps per integration
       * @param minimumStepSize - Minimum allowed step [s]
       * @param maximumStepSize - Maximum allowed step [s]
       * @param relativeErrorTolerance - Relative error tolerance
       * @param absoluteErrorTolerance - Absolute error tolerance
       * @returns Integrator settings
       */
      export function bulirsch_stoer(
        initialTimeStep: number,
        extrapolationSequence: number,
        maximumNumberOfSteps: number,
        minimumStepSize: number,
        maximumStepSize: number,
        relativeErrorTolerance?: number,
        absoluteErrorTolerance?: number
      ): IntegratorSettings;
    }

    /**
     * Acceleration model types
     */
    export namespace acceleration {
      export enum AvailableAcceleration {
        undefined_acceleration,
        point_mass_gravity,
        spherical_harmonic_gravity,
        aerodynamic,
        radiation_pressure,
        thrust_acceleration,
        relativistic_correction,
        empirical_acceleration
      }

      export interface AccelerationSettings {
        delete(): void;
      }

      /**
       * Point mass gravity acceleration
       * @returns Acceleration settings
       */
      export function point_mass_gravity(): AccelerationSettings;

      /**
       * Spherical harmonic gravity acceleration
       * @param maximumDegree - Maximum degree of expansion
       * @param maximumOrder - Maximum order of expansion
       * @returns Acceleration settings
       */
      export function spherical_harmonic_gravity(
        maximumDegree: number,
        maximumOrder: number
      ): AccelerationSettings;

      /**
       * Aerodynamic acceleration
       * @returns Acceleration settings
       */
      export function aerodynamic(): AccelerationSettings;

      /**
       * Cannonball radiation pressure acceleration
       * @returns Acceleration settings
       */
      export function radiation_pressure(): AccelerationSettings;
    }

    /**
     * Propagator settings
     */
    export namespace propagator {
      export enum TranslationalPropagatorType {
        undefined_translational_propagator,
        cowell,
        encke,
        gauss_keplerian,
        gauss_modified_equinoctial,
        unified_state_model_quaternions,
        unified_state_model_modified_rodrigues_parameters,
        unified_state_model_exponential_map
      }

      export enum RotationalPropagatorType {
        undefined_rotational_propagator,
        quaternions,
        modified_rodrigues_parameters,
        exponential_map
      }

      export interface PropagationTerminationSettings {
        delete(): void;
      }

      export interface PropagatorSettings {
        delete(): void;
      }

      /**
       * Create time-based termination condition
       * @param terminationTime - Simulation end time [s since J2000]
       * @returns Termination settings
       */
      export function time_termination(
        terminationTime: number
      ): PropagationTerminationSettings;

      /**
       * Create CPU time termination condition
       * @param cpuTimeLimit - Maximum CPU time [s]
       * @returns Termination settings
       */
      export function cpu_time_termination(
        cpuTimeLimit: number
      ): PropagationTerminationSettings;

      /**
       * Create translational state propagator settings
       * @param centralBodies - Central bodies for each propagated body
       * @param accelerationModelMap - Acceleration models
       * @param bodiesToPropagate - Bodies to propagate
       * @param initialStates - Initial states
       * @param terminationSettings - Termination condition
       * @param propagatorType - Type of propagator to use
       * @param dependentVariables - Variables to save
       * @param printInterval - Console output interval [s]
       * @returns Propagator settings
       */
      export function translational(
        centralBodies: string[],
        accelerationModelMap: any,
        bodiesToPropagate: string[],
        initialStates: VectorXd,
        terminationSettings: PropagationTerminationSettings,
        propagatorType?: TranslationalPropagatorType,
        dependentVariables?: any[],
        printInterval?: number
      ): PropagatorSettings;
    }

    /**
     * Create acceleration models from settings
     * @param bodies - System of bodies
     * @param accelerationSettings - Acceleration settings per body
     * @param bodiesToPropagate - Bodies being propagated
     * @param centralBodies - Central bodies for propagation
     * @returns Acceleration model map
     */
    export function create_acceleration_models(
      bodies: environment.SystemOfBodies,
      accelerationSettings: any,
      bodiesToPropagate: string[],
      centralBodies: string[]
    ): any;
  }

  /**
   * Simulator classes
   */
  export namespace simulator {
    export interface SimulationResults {
      delete(): void;
    }

    export interface SingleArcSimulator {
      /** Run the simulation */
      integrate_equations_of_motion(): void;
      /** Get simulation results */
      propagation_results(): SimulationResults;
      /** Check if integration succeeded */
      integration_completed_successfully(): boolean;
      delete(): void;
    }

    /**
     * Create a dynamics simulator
     * @param bodies - System of bodies
     * @param integratorSettings - Integrator settings
     * @param propagatorSettings - Propagator settings
     * @returns Dynamics simulator
     */
    export function create_dynamics_simulator(
      bodies: environment.SystemOfBodies,
      integratorSettings: propagation_setup.integrator.IntegratorSettings,
      propagatorSettings: propagation_setup.propagator.PropagatorSettings
    ): SingleArcSimulator;
  }
}

// =============================================================================
// Interface Module
// =============================================================================

/**
 * External interfaces (SPICE, etc.)
 * @namespace
 */
export namespace interface_ {
  /**
   * SPICE ephemeris interface
   */
  export namespace spice {
    /**
     * Load a SPICE kernel file
     * @param kernelPath - Path to kernel file
     */
    export function load_kernel(kernelPath: string): void;

    /**
     * Clear all loaded SPICE kernels
     */
    export function clear_kernels(): void;

    /**
     * Get number of loaded kernels
     * @returns Number of kernels
     */
    export function get_total_count_of_kernels_loaded(): number;

    /**
     * Convert Julian date to ephemeris time (TDB)
     * @param julianDate - Julian date
     * @returns Ephemeris time [s since J2000]
     */
    export function convert_julian_date_to_ephemeris_time(
      julianDate: number
    ): number;

    /**
     * Convert ephemeris time (TDB) to Julian date
     * @param ephemerisTime - Ephemeris time [s since J2000]
     * @returns Julian date
     */
    export function convert_ephemeris_time_to_julian_date(
      ephemerisTime: number
    ): number;

    /**
     * Convert date string to ephemeris time
     * @param dateString - Date string (e.g., "2024-01-15 12:00:00")
     * @returns Ephemeris time [s since J2000]
     */
    export function convert_date_string_to_ephemeris_time(
      dateString: string
    ): number;

    /**
     * Get body Cartesian state at epoch
     * @param targetBody - Name of target body
     * @param observerBody - Name of observer body
     * @param referenceFrame - Reference frame name
     * @param aberrationCorrections - Aberration correction type
     * @param ephemerisTime - Epoch [s since J2000]
     * @returns Cartesian state [x, y, z, vx, vy, vz]
     */
    export function get_body_cartesian_state_at_epoch(
      targetBody: string,
      observerBody: string,
      referenceFrame: string,
      aberrationCorrections: string,
      ephemerisTime: number
    ): Vector6d;

    /**
     * Get body gravitational parameter from SPICE
     * @param body - Body name
     * @returns Gravitational parameter [m³/s²]
     */
    export function get_body_gravitational_parameter(body: string): number;

    /**
     * Get body average radius from SPICE
     * @param body - Body name
     * @returns Average radius [m]
     */
    export function get_average_radius(body: string): number;
  }
}

// =============================================================================
// Data Module
// =============================================================================

/**
 * Data paths and file utilities
 * @namespace
 */
export namespace data {
  /** Get path to tudat resources directory */
  export function get_resource_path(): string;
  /** Get path to ephemeris data files */
  export function get_ephemeris_path(): string;
  /** Get path to Earth orientation data files */
  export function get_earth_orientation_path(): string;
  /** Get path to SPICE kernel files */
  export function get_spice_kernel_path(): string;
  /** Get path to atmosphere tables */
  export function get_atmosphere_tables_path(): string;
  /** Get path to gravity field models */
  export function get_gravity_models_path(): string;
  /** Get path to space weather data */
  export function get_space_weather_path(): string;
}

// =============================================================================
// Trajectory Design Module
// =============================================================================

/**
 * Trajectory design utilities
 * @namespace
 */
export namespace trajectory_design {
  /**
   * Transfer trajectory utilities
   */
  export namespace transfer_trajectory {
    export enum TransferLegTypes {
      unpowered_unperturbed_leg_type,
      dsm_position_based_leg_type,
      dsm_velocity_based_leg_type,
      spherical_shaping_low_thrust_leg,
      hodographic_low_thrust_leg
    }

    export interface TransferLeg {
      state_along_trajectory(timeSinceLegBeginning: number): Vector6d;
      delete(): void;
    }

    export interface TransferTrajectory {
      evaluate(nodeParameters: number[], legParameters: number[]): void;
      delta_v(): number;
      time_of_flight(): number;
      delete(): void;
    }

    export interface TransferNodeSettings {
      delete(): void;
    }

    export interface TransferLegSettings {
      delete(): void;
    }

    export function swingby_node(minimumPeriapsisRadius: number): TransferNodeSettings;
    export function departure_node(
      departureSemiMajorAxis: number,
      departureEccentricity: number
    ): TransferNodeSettings;
    export function capture_node(
      captureSemiMajorAxis: number,
      captureEccentricity: number
    ): TransferNodeSettings;

    export function unpowered_leg(): TransferLegSettings;
    export function dsm_position_based_leg(): TransferLegSettings;
    export function dsm_velocity_based_leg(): TransferLegSettings;
  }
}

// =============================================================================
// Module Initialization
// =============================================================================

/**
 * Options for initializing the Tudat WASM module
 */
export interface TudatModuleOptions {
  /** Callback for standard output */
  print?: (text: string) => void;
  /** Callback for error output */
  printErr?: (text: string) => void;
  /** Pre-run callbacks for filesystem setup */
  preRun?: Array<(module: TudatModule) => void>;
  /** Post-run callbacks */
  postRun?: Array<(module: TudatModule) => void>;
  /** Locator for .wasm file */
  locateFile?: (path: string, prefix: string) => string;
}

/**
 * The main Tudat WASM module interface
 */
export interface TudatModule {
  /** Filesystem API for mounting data files */
  FS: any;

  /** Memory management */
  _malloc(size: number): number;
  _free(ptr: number): void;

  /** Call a C function by name */
  ccall(
    ident: string,
    returnType: string,
    argTypes: string[],
    args: any[]
  ): any;

  /** Get a wrapped C function */
  cwrap(
    ident: string,
    returnType: string,
    argTypes: string[]
  ): (...args: any[]) => any;

  // Vector constructors
  Vector3d: new () => Vector3d;
  Vector6d: new () => Vector6d;
  VectorXd: new () => VectorXd;
  Matrix3d: new () => Matrix3d;

  // All namespaces available on Module
  constants: typeof constants;
  math: typeof math;
  astro: typeof astro;
  dynamics: typeof dynamics;
  interface: typeof interface_;
  data: typeof data;
  trajectory_design: typeof trajectory_design;
}

/**
 * Create and initialize the Tudat WASM module
 * @param options - Module initialization options
 * @returns Promise resolving to the initialized module
 *
 * @example
 * ```typescript
 * import createTudatModule from '@tudat/tudatpy-wasm';
 *
 * async function main() {
 *   const tudat = await createTudatModule();
 *
 *   // Convert Keplerian elements to Cartesian
 *   const kepler = new tudat.Vector6d();
 *   kepler.set(0, 7000000);  // semi-major axis [m]
 *   kepler.set(1, 0.01);     // eccentricity
 *   kepler.set(2, 0.5);      // inclination [rad]
 *   kepler.set(3, 0);        // argument of periapsis [rad]
 *   kepler.set(4, 0);        // RAAN [rad]
 *   kepler.set(5, 0);        // true anomaly [rad]
 *
 *   const GM = 3.986004418e14;  // Earth GM
 *   const cartesian = tudat.astro.element_conversion.keplerian_to_cartesian(kepler, GM);
 *
 *   console.log('Position:', cartesian.get(0), cartesian.get(1), cartesian.get(2));
 *   console.log('Velocity:', cartesian.get(3), cartesian.get(4), cartesian.get(5));
 *
 *   // Clean up
 *   kepler.delete();
 *   cartesian.delete();
 * }
 *
 * main();
 * ```
 */
declare function createTudatModule(
  options?: TudatModuleOptions
): Promise<TudatModule>;

export default createTudatModule;
