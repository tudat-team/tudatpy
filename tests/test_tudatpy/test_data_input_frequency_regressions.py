from tudatpy.astro.time_representation import Time
from tudatpy.data_input.tracking_data import FrequencyRamp, RampedFrequencySupplementaryData


def test_ramp_epochs_preserve_time_precision_and_accept_floats():
    start = Time(194444, 1600.0) + 1.0e-8
    ramps = RampedFrequencySupplementaryData()
    ramps.add_frequency_ramp(start, start + 60.0, 8.4e9, 0.1)
    ramp = ramps.frequency_ramps[0]
    assert float(ramp.start_time - start) == 0.0
    assert float(ramp.end_time - (start + 60.0)) == 0.0
    assert float(FrequencyRamp(start, start + 60.0, 8.4e9, 0.1).start_time - start) == 0.0
    ramps.add_frequency_ramp(1.0, 2.0, 8.4e9, 0.0)
    assert float(ramps.frequency_ramps[1].start_time) == 1.0
    ramp.start_time = 3.0
    assert float(ramp.start_time) == 3.0
