import numpy as np
from tudatpy.data_input.environment_data.damit import DamitModel
from .. import rotation_model


def _damit_simple_rotation(
    damit_model: DamitModel, base_frame: str, target_frame: str, extra_rotation: float = 0.0
) -> rotation_model.RotationModelSettings:
    """Simple rotation model settings of a DAMIT model, with an optional extra offset around the z-axis"""
    return rotation_model.simple(
        base_frame=base_frame,
        target_frame=target_frame,
        initial_orientation=damit_model.initial_orientation(base_frame, extra_rotation),
        initial_time=damit_model.reference_epoch,
        rotation_rate=2 * np.pi / damit_model.period,
    )


def damit_polyhedron(
    damit_model: DamitModel, base_frame: str, target_frame: str
) -> rotation_model.RotationModelSettings:
    """
    Factory function to create rotation model settings of a DAMIT polyhedron shape model.

    Rotation is assumed to be around a fixed axis in inertial space with a fixed rate of rotation.

    This function wraps :class:`~tudatpy.data_input.environment_data.damit.DamitModel`.

    Parameters
    ----------
    damit_model : DamitModel
        DAMIT model to create the rotation model settings for
    base_frame : str
        Base frame used in the environment. Options: "J2000", "ECLIPJ2000"
    target_frame : str
        Name of the body-fixed frame

    Returns
    -------
    RotationModelSettings
        Simple rotation model settings for the DAMIT polyhedron shape model

    Examples
    --------
    >>> eros = DamitModel.query(mpc_number=433)
    >>> body_settings.get("Eros").rotation_model_settings = environment_setup.rotation_model.damit_polyhedron(
    ...     eros, "ECLIPJ2000", "Eros_fixed"
    ... )
    """
    return _damit_simple_rotation(damit_model, base_frame, target_frame)


def damit_best_fitting_ellipsoid(
    damit_model: DamitModel, base_frame: str, target_frame: str
) -> rotation_model.RotationModelSettings:
    """
    Factory function to create rotation model settings of the best fitting ellipsoid of a DAMIT shape model.

    Creates the rotation model settings of the best fitting ellipsoid with equivalent moment of inertia to the DAMIT
    polyhedron shape model (see
    :func:`~tudatpy.dynamics.environment_setup.shape.damit_wrapper.damit_best_fitting_ellipsoid`). Rotation is
    assumed to be around a fixed axis in inertial space with a fixed rate of rotation.

    NOTE: this rotation model is **not** the same as computed by
    :func:`~tudatpy.dynamics.environment_setup.rotation_model.damit_wrapper.damit_polyhedron`,
    because an extra offset to the initial phase of rotation is necessary for the ellipsoid model.

    This function wraps :class:`~tudatpy.data_input.environment_data.damit.DamitModel`.

    Parameters
    ----------
    damit_model : DamitModel
        DAMIT model to create the rotation model settings for
    base_frame : str
        Base frame used in the environment. Options: "J2000", "ECLIPJ2000"
    target_frame : str
        Name of the body-fixed frame

    Returns
    -------
    RotationModelSettings
        Simple rotation model settings for the best fitting ellipsoid
    """
    _, principal_axis_offset = damit_model.best_fitting_ellipsoid()
    return _damit_simple_rotation(damit_model, base_frame, target_frame, principal_axis_offset)
