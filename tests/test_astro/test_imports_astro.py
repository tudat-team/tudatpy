import datetime


def test_import_from():
    from tudatpy import astro

    cartesian = astro.element_conversion.keplerian_to_cartesian_elementwise(
        gravitational_parameter=3.986004418e14,
        semi_major_axis=6.99276221e06,
        eccentricity=4.03294322e-03,
        inclination=1.71065169e00,
        argument_of_periapsis=1.31226971e00,
        longitude_of_ascending_node=3.82958313e-01,
        true_anomaly=3.07018490e00,
    )
    _ = astro.frame_conversion.inertial_to_rsw_rotation_matrix(cartesian)
    _ = astro.fundamentals.compute_shadow_function(
        [1, 2, 3], 10, [0, 1, 2], 2, [1, 1, 1]
    )
    _ = astro.gravitation.legendre_normalization_factor(2, 2)
    # Replace with a call if possible
    astro.polyhedron_utilities.surface_area
    #######################################
    _ = astro.time_conversion.datetime_to_tudat(
        datetime.datetime.fromisoformat("2023-06-20T00:05:23.281765")
    )
    _ = astro.two_body_dynamics.compute_escape_or_capture_delta_v(
        3.986004418e14,
        6.99276221e06,
        4.03294322e-03,
        1.71065169e00,
    )


def test_import_as():
    import tudatpy.astro as astro

    cartesian = astro.element_conversion.keplerian_to_cartesian_elementwise(
        gravitational_parameter=3.986004418e14,
        semi_major_axis=6.99276221e06,
        eccentricity=4.03294322e-03,
        inclination=1.71065169e00,
        argument_of_periapsis=1.31226971e00,
        longitude_of_ascending_node=3.82958313e-01,
        true_anomaly=3.07018490e00,
    )
    _ = astro.frame_conversion.inertial_to_rsw_rotation_matrix(cartesian)
    _ = astro.fundamentals.compute_shadow_function(
        [1, 2, 3], 10, [0, 1, 2], 2, [1, 1, 1]
    )
    _ = astro.gravitation.legendre_normalization_factor(2, 2)
    # Replace with a call if possible
    astro.polyhedron_utilities.surface_area
    #######################################
    _ = astro.time_conversion.datetime_to_tudat(
        datetime.datetime.fromisoformat("2023-06-20T00:05:23.281765")
    )
    _ = astro.two_body_dynamics.compute_escape_or_capture_delta_v(
        3.986004418e14,
        6.99276221e06,
        4.03294322e-03,
        1.71065169e00,
    )


def test_import_submodules_from():

    from tudatpy.astro import (
        element_conversion,
        frame_conversion,
        fundamentals,
        gravitation,
        polyhedron_utilities,
        time_conversion,
        two_body_dynamics,
    )

    cartesian = element_conversion.keplerian_to_cartesian_elementwise(
        gravitational_parameter=3.986004418e14,
        semi_major_axis=6.99276221e06,
        eccentricity=4.03294322e-03,
        inclination=1.71065169e00,
        argument_of_periapsis=1.31226971e00,
        longitude_of_ascending_node=3.82958313e-01,
        true_anomaly=3.07018490e00,
    )
    _ = frame_conversion.inertial_to_rsw_rotation_matrix(cartesian)
    _ = fundamentals.compute_shadow_function(
        [1, 2, 3], 10, [0, 1, 2], 2, [1, 1, 1]
    )
    _ = gravitation.legendre_normalization_factor(2, 2)
    # Replace with a call if possible
    polyhedron_utilities.surface_area
    #######################################
    _ = time_conversion.datetime_to_tudat(
        datetime.datetime.fromisoformat("2023-06-20T00:05:23.281765")
    )
    _ = two_body_dynamics.compute_escape_or_capture_delta_v(
        3.986004418e14,
        6.99276221e06,
        4.03294322e-03,
        1.71065169e00,
    )


def test_import_submodules_as():
    """Import each submodule with alias and call some function from it"""

    import tudatpy.astro.element_conversion as element_conversion
    import tudatpy.astro.frame_conversion as frame_conversion
    import tudatpy.astro.fundamentals as fundamentals
    import tudatpy.astro.gravitation as gravitation
    import tudatpy.astro.polyhedron_utilities as polyhedron_utilities
    import tudatpy.astro.time_conversion as time_conversion
    import tudatpy.astro.two_body_dynamics as two_body_dynamics

    cartesian = element_conversion.keplerian_to_cartesian_elementwise(
        gravitational_parameter=3.986004418e14,
        semi_major_axis=6.99276221e06,
        eccentricity=4.03294322e-03,
        inclination=1.71065169e00,
        argument_of_periapsis=1.31226971e00,
        longitude_of_ascending_node=3.82958313e-01,
        true_anomaly=3.07018490e00,
    )
    _ = frame_conversion.inertial_to_rsw_rotation_matrix(cartesian)
    _ = fundamentals.compute_shadow_function(
        [1, 2, 3], 10, [0, 1, 2], 2, [1, 1, 1]
    )
    _ = gravitation.legendre_normalization_factor(2, 2)
    # Replace with a call if possible
    polyhedron_utilities.surface_area
    #######################################
    _ = time_conversion.datetime_to_tudat(
        datetime.datetime.fromisoformat("2023-06-20T00:05:23.281765")
    )
    _ = two_body_dynamics.compute_escape_or_capture_delta_v(
        3.986004418e14,
        6.99276221e06,
        4.03294322e-03,
        1.71065169e00,
    )
