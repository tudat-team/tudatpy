def test_imports_astro() -> None:

    # Importing from astro doesn't raise errors [PYBIND11]
    from tudatpy import astro
    from tudatpy.astro import (
        element_conversion,
        ephemerides,
        frame_conversion,
        fundamentals,
        gravitation,
        polyhedron_utilities,
        time_conversion,
        two_body_dynamics,
    )

    return None


def test_imports_constants() -> None:

    from tudatpy import constants

    return None


def test_imports_data() -> None:

    from tudatpy import data
    from tudatpy.data import (
        horizons,
        mpc,
        sbdb,
    )

    return None


def test_imports_interface() -> None:

    from tudatpy import interface
    from tudatpy.interface import spice

    return None


class TestImportsMath:

    def test_import_from(self):
        from tudatpy import math

    def test_import_as(self):
        import tudatpy.math as math

    def test_import_submodules_from(self):
        from tudatpy.math import (
            geometry,
            interpolators,
            numerical_integrators,
            root_finders,
            statistics,
        )

    def test_import_submodules_as(self):
        """Import each submodule with alias and call some function from it"""

        # Geometry
        import tudatpy.math.geometry as geometry
        import tudatpy.math.interpolators as interpolators
        import tudatpy.math.numerical_integrators as numerical_integrators
        import tudatpy.math.root_finders as root_finders
        import tudatpy.math.statistics as statistics

        _ = geometry.Capsule(1, 2, 3, 0.5, 0.2)
        _ = interpolators.hermite_spline_interpolation()


def test_imports_math() -> None:

    from tudatpy import math

    from tudatpy import math
    from tudatpy.math import (
        geometry,
        interpolators,
        numerical_integrators,
        root_finders,
        statistics,
    )

    return None


def test_imports_numerical_simulation() -> None:

    from tudatpy import numerical_simulation
    from tudatpy.numerical_simulation import (
        environment,
        environment_setup,
        propagation_setup,
        propagation,
        estimation,
        estimation_setup,
    )

    from tudatpy.numerical_simulation.estimation_setup import (
        observation,
        parameter,
    )

    from tudatpy.numerical_simulation.environment_setup import (
        aerodynamic_coefficients,
        atmosphere,
        ephemeris,
        gravity_field,
        gravity_field_variation,
        ground_station,
        radiation_pressure,
        rigid_body,
        rotation_model,
        shape,
        shape_deformation,
        vehicle_systems,
    )

    from tudatpy.numerical_simulation.propagation_setup import (
        acceleration,
        dependent_variable,
        integrator,
        mass_rate,
        propagator,
        thrust,
        torque,
    )

    return None


def test_imports_plotting() -> None:

    from tudatpy import plotting

    return None


def test_import_trajectory_design() -> None:

    from tudatpy import trajectory_design

    from tudatpy.trajectory_design import (
        porkchop,
        shape_based_thrust,
        transfer_trajectory,
    )

    return None


def test_imports_util() -> None:

    from tudatpy import util

    return None


def test_import_all() -> None:

    from tudatpy import (
        astro,
        constants,
        data,
        interface,
        math,
        numerical_simulation,
        plotting,
        trajectory_design,
        util,
    )
    from tudatpy.astro import (
        element_conversion,
        ephemerides,
        frame_conversion,
        fundamentals,
        gravitation,
        polyhedron_utilities,
        time_conversion,
        two_body_dynamics,
    )
    from tudatpy.interface import spice
    from tudatpy.math import (
        geometry,
        interpolators,
        numerical_integrators,
        root_finders,
        statistics,
    )
    from tudatpy.numerical_simulation import (
        environment,
        environment_setup,
        propagation_setup,
        propagation,
        estimation,
        estimation_setup,
    )

    from tudatpy.numerical_simulation.estimation_setup import (
        observation,
        parameter,
    )

    from tudatpy.numerical_simulation.environment_setup import (
        aerodynamic_coefficients,
        atmosphere,
        ephemeris,
        gravity_field,
        gravity_field_variation,
        ground_station,
        radiation_pressure,
        rigid_body,
        rotation_model,
        shape,
        shape_deformation,
        vehicle_systems,
    )

    from tudatpy.numerical_simulation.propagation_setup import (
        acceleration,
        dependent_variable,
        integrator,
        mass_rate,
        propagator,
        thrust,
        torque,
    )
    from tudatpy.trajectory_design import (
        porkchop,
        shape_based_thrust,
        transfer_trajectory,
    )

    return None
