import numpy as np
from tudatpy.data_input.environment_data.damit import DamitModel
from .. import rotation_model


def damit(
    damit_model: DamitModel, base_frame: str, target_frame: str
) -> rotation_model.RotationModelSettings:
    """Factory function to create rotation model settings from a DAMIT model.

    Creates simple rotation model settings (see
    :func:`~tudatpy.dynamics.environment_setup.rotation_model.simple`), with the rotation axis, rotation rate and
    initial orientation of a :class:`~tudatpy.data_input.environment_data.damit.DamitModel`. The body-fixed frame is
    that of the model, in which both its polyhedron shape and its best fitting ellipsoid are defined (see
    :func:`~tudatpy.dynamics.environment_setup.shape.damit_polyhedron` and
    :func:`~tudatpy.dynamics.environment_setup.shape.damit_best_fitting_ellipsoid`).

    This function wraps :class:`~tudatpy.data_input.environment_data.damit.DamitModel`.

    Parameters
    ----------
    damit_model : DamitModel
        DAMIT model from which the rotation model settings are created.
    base_frame : str
        Base frame orientation of the environment: ``"J2000"`` or ``"ECLIPJ2000"``.
    target_frame : str
        Name of the body-fixed frame.

    Returns
    -------
    RotationModelSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings` derived
        :class:`~tudatpy.dynamics.environment_setup.rotation_model.SimpleRotationModelSettings` class.

    Examples
    --------
    >>> eros = DamitModel.query(mpc_number=433)
    >>> body_settings.get("Eros").rotation_model_settings = environment_setup.rotation_model.damit(
    ...     eros, "ECLIPJ2000", "Eros_fixed"
    ... )
    """
    return rotation_model.simple(
        base_frame=base_frame,
        target_frame=target_frame,
        initial_orientation=damit_model.initial_orientation(base_frame),
        initial_time=damit_model.reference_epoch,
        rotation_rate=2 * np.pi / damit_model.period,
    )
