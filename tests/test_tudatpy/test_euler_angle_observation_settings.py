from tudatpy.estimation.observable_models_setup import links, model_settings


def test_euler_angles_313_settings():
    link_definition = links.LinkDefinition(
        {links.observed_body: links.body_origin_link_end_id("Earth")}
    )

    settings = model_settings.euler_angles_313(link_definition)

    assert isinstance(settings, model_settings.ObservationModelSettings)
    assert model_settings.get_observable_size(model_settings.euler_angle_313_observable_type) == 3
