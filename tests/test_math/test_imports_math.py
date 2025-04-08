def test_import_from():
    from tudatpy import math

    # Call a function from each submodule to check if the import works
    _ = math.geometry.Capsule(1, 2, 3, 0.5, 0.2)
    _ = math.interpolators.hermite_spline_interpolation()
    _ = math.root_finders.NewtonRaphsonCore(1e-2, 10)
    _ = math.statistics.calculate_allan_variance_of_dataset([0.4, 0.2], 0.1)


def test_import_as():
    import tudatpy.math as math

    # Call a function from each submodule to check if the import works
    _ = math.geometry.Capsule(1, 2, 3, 0.5, 0.2)
    _ = math.interpolators.hermite_spline_interpolation()
    _ = math.root_finders.NewtonRaphsonCore(1e-2, 10)
    _ = math.statistics.calculate_allan_variance_of_dataset([0.4, 0.2], 0.1)


def test_import_submodules_from():
    from tudatpy.math import (
        geometry,
        interpolators,
        numerical_integrators,
        root_finders,
        statistics,
    )

    # Call a function from each submodule to check if the import works
    _ = geometry.Capsule(1, 2, 3, 0.5, 0.2)
    _ = interpolators.hermite_spline_interpolation()
    _ = root_finders.NewtonRaphsonCore(1e-2, 10)
    _ = statistics.calculate_allan_variance_of_dataset([0.4, 0.2], 0.1)


def test_import_submodules_as():
    """Import each submodule with alias and call some function from it"""

    # Geometry
    import tudatpy.math.geometry as geometry
    import tudatpy.math.interpolators as interpolators
    import tudatpy.math.root_finders as root_finders
    import tudatpy.math.statistics as statistics

    # Call a function from each submodule to check if the import works
    _ = geometry.Capsule(1, 2, 3, 0.5, 0.2)
    _ = interpolators.hermite_spline_interpolation()
    _ = root_finders.NewtonRaphsonCore(1e-2, 10)
    _ = statistics.calculate_allan_variance_of_dataset([0.4, 0.2], 0.1)
