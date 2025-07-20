import typing
__all__ = ['calculate_allan_variance_of_dataset', 'convert_allan_variance_amplitudes_to_phase_noise_amplitudes']

def calculate_allan_variance_of_dataset(timing_errors: list[float], time_step_size: float) -> dict[float, float]:
    """No documentation found."""

def convert_allan_variance_amplitudes_to_phase_noise_amplitudes(allan_variance_amplitudes: dict[int, float], frequency_domain_cutoff_frequency: float, is_inverse_square_term_flicker_phase_noise: bool=0) -> dict[int, float]:
    """No documentation found."""