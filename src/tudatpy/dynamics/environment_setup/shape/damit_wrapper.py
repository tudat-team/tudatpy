from tudatpy.data_input.environment_data.damit import DamitModel
from .. import shape


def damit_polyhedron(damit_model: DamitModel, *args, **kwargs) -> shape.PolyhedronBodyShapeSettings:
    """
    Factory function to create polyhedron shape model settings from a DAMIT shape model.

    This function wraps :class:`~tudatpy.data_input.environment_data.damit.DamitModel`.

    Parameters
    ----------
    damit_model : DamitModel
        DAMIT model to create the shape settings for
    args, kwargs
        Remaining inputs are passed to :func:`~tudatpy.dynamics.environment_setup.shape.polyhedron`

    Returns
    -------
    PolyhedronBodyShapeSettings
        Polyhedron shape settings for the DAMIT model

    Examples
    --------
    >>> eros = DamitModel.query(mpc_number=433)
    >>> body_settings.get("Eros").shape_settings = environment_setup.shape.damit_polyhedron(eros)
    """
    return shape.polyhedron(damit_model.vertices, damit_model.facets, *args, **kwargs)


def damit_best_fitting_ellipsoid(
    damit_model: DamitModel,
) -> shape.TriAxialEllipsoidBodyShapeSettings:
    """
    Factory function to create the best fitting triaxial ellipsoid shape model settings from a DAMIT shape model.

    Creates a triaxial ellipsoid shape model for an ellipsoid with equivalent moment of inertia to the DAMIT
    polyhedron shape model. Use together with
    :func:`~tudatpy.dynamics.environment_setup.rotation_model.damit_wrapper.damit_best_fitting_ellipsoid` for the
    rotation model, to correctly align the ellipsoid.

    This function wraps :class:`~tudatpy.data_input.environment_data.damit.DamitModel`.

    Parameters
    ----------
    damit_model : DamitModel
        DAMIT model to create the shape settings for

    Returns
    -------
    TriAxialEllipsoidBodyShapeSettings
        Shape settings for the best fitting triaxial ellipsoid
    """
    semi_axes, _ = damit_model.best_fitting_ellipsoid()
    return shape.triaxial_ellipsoid(*semi_axes)
