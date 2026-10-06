/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_statistics.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/basics/basicTypedefs.h>

#include "tudat/astro/system_models/timingSystem.h"
#include "tudat/math/statistics/allanVariance.h"

namespace py = pybind11;

namespace ts = tudat::statistics;
namespace tsm = tudat::system_models;

namespace tudatpy
{

void expose_statistics( py::module& m )
{
    m.def( "calculate_allan_variance_of_dataset",
           &ts::calculateAllanVarianceOfTimeDataSet,
           py::arg( "timing_errors" ),
           py::arg( "time_step_size" ),
           R"doc(

         Calculate Allan variance from uniformly spaced clock ``timing_errors`` in seconds.

         ``time_step_size`` is the sample spacing in seconds.

         Parameters
         ----------
         timing_errors : list[float]
             Uniformly sampled clock timing errors, in seconds.
         time_step_size : float
             Interval between successive timing-error samples, in seconds.

         Returns
         -------
         dict[float, float]
             Allan variance values keyed by averaging interval in seconds.

      )doc" );

    m.def( "convert_allan_variance_amplitudes_to_phase_noise_amplitudes",
           &tsm::convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes,
           py::arg( "allan_variance_amplitudes" ),
           py::arg( "frequency_domain_cutoff_frequency" ),
           py::arg( "is_inverse_square_term_flicker_phase_noise" ) = 0,
           R"doc(

         Convert an Allan variance power-law model to phase-noise power-law amplitudes.

         ``allan_variance_amplitudes`` maps integer powers of averaging time to their amplitudes.
         ``frequency_domain_cutoff_frequency`` is the high-frequency cutoff in hertz.
         ``is_inverse_square_term_flicker_phase_noise`` selects flicker phase noise for the inverse-square term;
         otherwise white phase noise is used.

         Parameters
         ----------
         allan_variance_amplitudes : dict[int, float]
             Amplitudes of the Allan variance power-law terms, keyed by integer powers of averaging time.
         frequency_domain_cutoff_frequency : float
             High-frequency cutoff of the phase-noise model, in hertz.
         is_inverse_square_term_flicker_phase_noise : bool, optional
             Whether the inverse-square Allan variance term represents flicker phase noise instead of white phase noise.

         Returns
         -------
         dict[int, float]
             Phase-noise amplitudes keyed by integer powers of frequency.

      )doc" );

#if ( TUDAT_BUILD_WITH_FFTW3 )
    m.def( "generate_noise_from_allan_deviation",
           &tsm::generateClockNoise,
           py::arg( "allan_variance_amplitudes" ),
           py::arg( "start_time" ),
           py::arg( "end_time" ),
           py::arg( "number_of_time_steps" ),
           py::arg( "is_inverse_square_term_flicker_phase_noise" ) = 0,
           py::arg( "seed" ) = ts::defaultRandomSeedGenerator->getRandomVariableValue( ),
           R"doc(

         Generate clock timing-noise samples from the supplied Allan variance power-law amplitudes.

         ``allan_variance_amplitudes`` maps averaging-time powers to amplitudes. ``start_time`` and ``end_time``
         delimit the interval in seconds and ``number_of_time_steps`` selects the number of samples.
         ``is_inverse_square_term_flicker_phase_noise`` selects the inverse-square noise interpretation; ``seed``
         controls random generation. Available
         when FFTW support is enabled.

         Parameters
         ----------
         allan_variance_amplitudes : dict[int, float]
             Amplitudes of the Allan variance power-law terms, keyed by integer powers of averaging time.
         start_time : float
             Start of the noise-generation interval, in seconds.
         end_time : float
             End of the noise-generation interval, in seconds.
         number_of_time_steps : int
             Number of samples in the generated noise series.
         is_inverse_square_term_flicker_phase_noise : bool, optional
             Whether the inverse-square Allan variance term represents flicker phase noise instead of white phase noise.
         seed : float, optional
             Seed used by the random number generator.

         Returns
         -------
         tuple[list[float], float]
             Pair containing clock timing-error samples in seconds and their sampling interval in seconds.

      )doc" );

    m.def( "generate_colored_clock_noise",
           &tsm::generateColoredClockNoise,
           py::arg( "allan_variance_amplitudes" ),
           py::arg( "variance_type" ),
           py::arg( "start_time" ),
           py::arg( "end_time" ),
           py::arg( "number_of_time_steps" ),
           py::arg( "seed" ) = ts::defaultRandomSeedGenerator->getRandomVariableValue( ),
           R"doc(

         Generate colored clock timing noise from a set of Allan or Hadamard deviation nodes.

         Despite its name, ``allan_variance_amplitudes`` maps averaging times to deviation values. ``variance_type``
         selects Allan or Hadamard variance. ``start_time``, ``end_time`` and ``number_of_time_steps`` define the
         sample grid; ``seed`` controls random generation. Available when FFTW support is enabled.

         Parameters
         ----------
         allan_variance_amplitudes : dict[int, float]
             Clock deviation values keyed by averaging time in seconds.
         variance_type : str
             Clock-stability variance to use: "Allan" or "Hadamard".
         start_time : float
             Start of the noise-generation interval, in seconds.
         end_time : float
             End of the noise-generation interval, in seconds.
         number_of_time_steps : int
             Number of samples in the generated noise series.
         seed : float, optional
             Seed used by the random number generator.

         Returns
         -------
         tuple[list[float], float]
             Pair containing clock timing-error samples in seconds and their sampling interval in seconds.

      )doc" );

    m.def( "get_clock_noise_interpolator",
           &tsm::getClockNoiseInterpolator,
           py::arg( "allan_variance_amplitudes" ),
           py::arg( "start_time" ),
           py::arg( "end_time" ),
           py::arg( "time_step" ),
           py::arg( "is_inverse_square_term_flicker_phase_noise" ) = 0,
           py::arg( "seed" ) = ts::defaultRandomSeedGenerator->getRandomVariableValue( ),
           R"doc(

         Generate and interpolate clock timing noise with the specified Allan variance power-law amplitudes.

         ``start_time``, ``end_time`` and ``time_step`` define the sampling interval in seconds.
         ``is_inverse_square_term_flicker_phase_noise`` selects the inverse-square noise interpretation and ``seed``
         controls random generation. Available when FFTW support is enabled.

         Parameters
         ----------
         allan_variance_amplitudes : dict[int, float]
             Amplitudes of the Allan variance power-law terms, keyed by integer powers of averaging time.
         start_time : float
             Start of the noise-generation interval, in seconds.
         end_time : float
             End of the noise-generation interval, in seconds.
         time_step : float
             Sampling interval, in seconds.
         is_inverse_square_term_flicker_phase_noise : bool, optional
             Whether the inverse-square Allan variance term represents flicker phase noise instead of white phase noise.
         seed : float, optional
             Seed used by the random number generator.

         Returns
         -------
         Callable[[float], float]
             Function returning the clock timing error in seconds at the supplied epoch in seconds.

      )doc" );

    m.def( "get_colored_clock_noise_interpolator",
           &tsm::getColoredClockNoiseInterpolator,
           py::arg( "allan_variance_nodes" ),
           py::arg( "variance_type" ),
           py::arg( "start_time" ),
           py::arg( "end_time" ),
           py::arg( "time_step" ),
           py::arg( "seed" ) = ts::defaultRandomSeedGenerator->getRandomVariableValue( ),
           R"doc(

         Generate and interpolate colored clock timing noise from ``allan_variance_nodes``, mapping averaging times to
         deviation values.

         ``variance_type`` selects Allan or Hadamard variance. ``start_time``, ``end_time`` and ``time_step`` define
         the sample grid in seconds and ``seed`` controls random generation.
         Available when FFTW support is enabled.

         Parameters
         ----------
         allan_variance_nodes : dict[int, float]
             Clock deviation values keyed by averaging time in seconds.
         variance_type : str
             Clock-stability variance to use: "Allan" or "Hadamard".
         start_time : float
             Start of the noise-generation interval, in seconds.
         end_time : float
             End of the noise-generation interval, in seconds.
         time_step : float
             Sampling interval, in seconds.
         seed : float, optional
             Seed used by the random number generator.

         Returns
         -------
         Callable[[float], float]
             Function returning the clock timing error in seconds at the supplied epoch in seconds.

      )doc" );
#endif
};

}  // namespace tudatpy
