import numpy as np
import pytest

from tudatpy.math import statistics


@pytest.mark.skipif(
    not hasattr(statistics, "generate_colored_clock_noise"),
    reason="Clock noise generation requires FFTW support",
)
@pytest.mark.parametrize("variance_type", ["Allan", "Hadamard"])
def test_colored_clock_noise_deviation_nodes_keyword(variance_type):
    nodes = {1: 1.0e-11, 10: 3.0e-12, 100: 1.0e-12}
    samples, time_step = statistics.generate_colored_clock_noise(
        deviation_nodes=nodes,
        variance_type=variance_type,
        start_time=1000.0,
        end_time=1128.0,
        number_of_time_steps=128,
        seed=42.0,
    )
    assert samples
    assert time_step == pytest.approx(1.0)
    assert np.all(np.isfinite(samples))

    interpolator = statistics.get_colored_clock_noise_interpolator(
        deviation_nodes=nodes,
        variance_type=variance_type,
        start_time=1000.0,
        end_time=1128.0,
        time_step=time_step,
        seed=42.0,
    )
    # With the same seed and grid, interpolation reproduces the generated series.
    for index in (0, len(samples) // 2, len(samples) - 1):
        assert interpolator(1000.0 + index * time_step) == pytest.approx(
            samples[index], rel=1.0e-12, abs=1.0e-25
        )
