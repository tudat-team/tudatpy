from tudatpy.astro.time_representation import Time
from tudatpy.data_input.tracking_data import FrequencyRamp, RampedFrequencySupplementaryData


def test_ramp_epochs_preserve_time_precision_and_accept_floats():
    """Preserve extended ramp epochs while accepting float inputs and independent record edits."""
    start = Time(194444, 1600.0) + 1.0e-8
    ramps = RampedFrequencySupplementaryData()
    ramps.add_frequency_ramp(start, start + 60.0, 8.4e9, 0.1)
    ramp = ramps.frequency_ramps[0]
    # Container and direct-constructor round trips retain a sub-double epoch increment.
    assert float(ramp.start_time - start) == 0.0
    assert float(ramp.end_time - (start + 60.0)) == 0.0
    assert float(FrequencyRamp(start, start + 60.0, 8.4e9, 0.1).start_time - start) == 0.0
    # Ordinary float arguments remain accepted by the container and record setter.
    ramps.add_frequency_ramp(1.0, 2.0, 8.4e9, 0.0)
    assert float(ramps.frequency_ramps[1].start_time) == 1.0
    ramp.start_time = 3.0
    assert float(ramp.start_time) == 3.0
    # Editing the returned record must not update the container's stored epoch.
    assert float(ramps.frequency_ramps[0].start_time - start) == 0.0


def test_frequency_ramp_records_survive_container_growth_and_deletion():
    """Keep retrieved ramp records valid and writable after container growth and deletion."""
    ramps = RampedFrequencySupplementaryData()
    ramps.add_frequency_ramp(1.0, 2.0, 8.4e9, 0.1)
    record = ramps.frequency_ramps[0]
    # Verify copy independence before deliberately growing the underlying vector.
    record.start_frequency = 9.0e9
    assert ramps.frequency_ramps[0].start_frequency == 8.4e9

    # Exercise reallocation and destruction of the original owner.
    for index in range(64):
        ramps.add_frequency_ramp(3.0 + index, 4.0 + index, 8.4e9, 0.0)
    del ramps

    # All fields survive, and the independent record remains writable.
    assert float(record.start_time) == 1.0
    assert float(record.end_time) == 2.0
    assert record.start_frequency == 9.0e9
    assert record.frequency_rate == 0.1
    record.end_time = 2.5
    assert float(record.end_time) == 2.5
