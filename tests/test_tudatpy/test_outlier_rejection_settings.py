import pytest

from tudatpy.estimation import estimation_analysis
from tudatpy.estimation.observable_models_setup import links, model_settings
from tudatpy.estimation.observations import ObservationSelectionCondition, observation_query

CONDITION_TYPE = type(ObservationSelectionCondition.all().condition_type)

ALGORITHMS = [
    estimation_analysis.carpino_outlier_rejection_settings,
    estimation_analysis.simple_outlier_rejection_settings,
]


def _range_and_angles():
    return (
        observation_query.observable_type == model_settings.one_way_range_type,
        observation_query.observable_type == model_settings.angular_position_type,
    )


@pytest.mark.parametrize(
    "create, settings_type",
    [
        (
            estimation_analysis.carpino_outlier_rejection_settings,
            estimation_analysis.CarpinoOutlierRejectionSettings,
        ),
        (
            estimation_analysis.simple_outlier_rejection_settings,
            estimation_analysis.SimpleOutlierRejectionSettings,
        ),
    ],
)
def test_single_pair_applies_to_all_observations(create, settings_type):
    settings = create((8.0, 9.0), 2)
    assert isinstance(settings, settings_type)
    assert isinstance(settings, estimation_analysis.OutlierRejectionSettings)
    assert settings.first_iteration_with_rejection == 2

    thresholds = settings.thresholds_per_condition
    assert len(thresholds) == 1
    condition, pair = thresholds[0]
    assert condition.condition_type == CONDITION_TYPE.all
    assert pair == (8.0, 9.0)


@pytest.mark.parametrize("create", ALGORITHMS)
def test_thresholds_per_condition_round_trip(create):
    range_condition, angle_condition = _range_and_angles()
    settings = create([(range_condition, (24.0, 25.0)), (angle_condition, (28.0, 29.0))])
    assert settings.first_iteration_with_rejection == 1

    thresholds = settings.thresholds_per_condition
    assert [pair for _, pair in thresholds] == [(24.0, 25.0), (28.0, 29.0)]
    assert [condition.condition_type for condition, _ in thresholds] == [
        CONDITION_TYPE.observable_type,
        CONDITION_TYPE.observable_type,
    ]
    assert thresholds[0][0].observable_type_value == model_settings.one_way_range_type
    assert thresholds[1][0].observable_type_value == model_settings.angular_position_type


def test_simple_thresholds_per_link_end():
    from_gaia = observation_query.receiver == links.body_origin_link_end_id("Gaia")
    settings = estimation_analysis.simple_outlier_rejection_settings(
        [(from_gaia, (5.0e-9, 1.0e-8)), (~from_gaia, (5.0e-6, 1.0e-5))],
        first_iteration_with_rejection=0,
    )
    thresholds = settings.thresholds_per_condition
    assert thresholds[1][0].condition_type == CONDITION_TYPE.not_condition
    assert thresholds[1][1] == (5.0e-6, 1.0e-5)
    assert settings.first_iteration_with_rejection == 0


@pytest.mark.parametrize("create", ALGORITHMS)
@pytest.mark.parametrize("valid_pair", [(0.0, 0.0), (0.0, 9.0), (9.0, 9.0)])
def test_valid_thresholds_are_accepted(create, valid_pair):
    """Both algorithms accept any pair with 0 <= recovery <= rejection."""
    range_condition, _ = _range_and_angles()
    assert create(valid_pair).thresholds_per_condition[0][1] == valid_pair
    assert create([(range_condition, valid_pair)]).thresholds_per_condition[0][1] == valid_pair


@pytest.mark.parametrize("create", ALGORITHMS)
@pytest.mark.parametrize("invalid_pair", [(9.0, 8.0), (-1.0, 9.0), (0.0, -1.0)])
def test_invalid_thresholds_are_rejected(create, invalid_pair):
    range_condition, angle_condition = _range_and_angles()
    with pytest.raises(RuntimeError):
        create(invalid_pair)
    with pytest.raises(RuntimeError):
        create([(range_condition, (1.0, 2.0)), (angle_condition, invalid_pair)])


@pytest.mark.parametrize("create", ALGORITHMS)
def test_empty_list_of_conditions_is_rejected(create):
    with pytest.raises(RuntimeError):
        create([])
