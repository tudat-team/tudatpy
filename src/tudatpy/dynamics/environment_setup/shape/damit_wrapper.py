from tudatpy.data_input.environment_data.damit import DamitModel
from .. import shape


def damit_polyhedron(damit_model: DamitModel, *args, **kwargs) -> shape.PolyhedronBodyShapeSettings:
    """Factory function to create polyhedron body shape model settings from a DAMIT model.

    Creates polyhedron body shape model settings (see
    :func:`~tudatpy.dynamics.environment_setup.shape.polyhedron`) with the vertices and facets of a
    :class:`~tudatpy.data_input.environment_data.damit.DamitModel`. Use together with
    :func:`~tudatpy.dynamics.environment_setup.rotation_model.damit` for the rotation model.

    This function wraps :class:`~tudatpy.data_input.environment_data.damit.DamitModel`.

    Parameters
    ----------
    damit_model : DamitModel
        DAMIT model from which the shape model settings are created.
    *args, **kwargs
        Remaining inputs of :func:`~tudatpy.dynamics.environment_setup.shape.polyhedron`.

    Returns
    -------
    PolyhedronBodyShapeSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings` derived
        :class:`~tudatpy.dynamics.environment_setup.shape.PolyhedronBodyShapeSettings` class.

    Examples
    --------
    >>> eros = DamitModel.query(mpc_number=433)
    >>> body_settings.get("Eros").shape_settings = environment_setup.shape.damit_polyhedron(eros)
    """
    return shape.polyhedron(damit_model.vertices, damit_model.facets, *args, **kwargs)


def damit_best_fitting_ellipsoid(
    damit_model: DamitModel,
) -> shape.TriAxialEllipsoidBodyShapeSettings:
    """Factory function to create tri-axial ellipsoid body shape model settings from a DAMIT model.

    Creates tri-axial ellipsoid body shape model settings (see
    :func:`~tudatpy.dynamics.environment_setup.shape.triaxial_ellipsoid`) for the ellipsoid with the same moments of
    inertia as the polyhedron shape of a :class:`~tudatpy.data_input.environment_data.damit.DamitModel` (see
    :meth:`~tudatpy.data_input.environment_data.damit.DamitModel.best_fitting_ellipsoid`). Use together with
    :func:`~tudatpy.dynamics.environment_setup.rotation_model.damit` for the rotation model.

    This function wraps :class:`~tudatpy.data_input.environment_data.damit.DamitModel`.

    Parameters
    ----------
    damit_model : DamitModel
        DAMIT model from which the shape model settings are created.

    Returns
    -------
    TriAxialEllipsoidBodyShapeSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings` derived
        :class:`~tudatpy.dynamics.environment_setup.shape.TriAxialEllipsoidBodyShapeSettings` class.
    """
    return shape.triaxial_ellipsoid(*damit_model.best_fitting_ellipsoid())
