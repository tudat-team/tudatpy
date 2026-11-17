"""Check that transmitted-frequency schedules keep exact times and remain safe to edit."""

from tudatpy.astro.time_representation import Time
from tudatpy.data_input.tracking_data import FrequencyRamp, RampedFrequencySupplementaryData


def test_ramp_epochs_preserve_time_precision_and_accept_floats():
    """Keep an extra 0.00000001 second, while still accepting ordinary second values."""
    start = Time(194444, 1600.0) + 1.0e-8
    ramps = RampedFrequencySupplementaryData()
    ramps.add_frequency_ramp(start, start + 60.0, 8.4e9, 0.1)
    ramp = ramps.frequency_ramps[0]
    # The start and end retain the small extra fraction that was supplied.
    assert float(ramp.start_time - start) == 0.0
    assert float(ramp.end_time - (start + 60.0)) == 0.0
    assert float(FrequencyRamp(start, start + 60.0, 8.4e9, 0.1).start_time - start) == 0.0
    # Times supplied as ordinary numbers are still accepted and can be changed.
    ramps.add_frequency_ramp(1.0, 2.0, 8.4e9, 0.0)
    assert float(ramps.frequency_ramps[1].start_time) == 1.0
    ramp.start_time = 3.0
    assert float(ramp.start_time) == 3.0


def test_frequency_ramp_records_survive_container_growth_and_deletion():
    """Save one ramp, add 64 more, delete the schedule, and check the saved ramp still works."""
    ramps = RampedFrequencySupplementaryData()
    ramps.add_frequency_ramp(1.0, 2.0, 8.4e9, 0.1)
    record = ramps.frequency_ramps[0]
    # Changing the saved ramp must not change the ramp kept in the full schedule.
    record.start_frequency = 9.0e9
    assert ramps.frequency_ramps[0].start_frequency == 8.4e9

    # Add enough ramps to force the schedule to grow, and then delete that schedule.
    for index in range(64):
        ramps.add_frequency_ramp(3.0 + index, 4.0 + index, 8.4e9, 0.0)
    del ramps

    # The saved ramp keeps every original value and remains editable on its own.
    assert float(record.start_time) == 1.0
    assert float(record.end_time) == 2.0
    assert record.start_frequency == 9.0e9
    assert record.frequency_rate == 0.1
    record.end_time = 2.5
    assert float(record.end_time) == 2.5
