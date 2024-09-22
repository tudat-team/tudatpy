import numpy
import typing
__all__ = ['BodyPanelGeometrySettings', 'BodyPanelSettings', 'FullPanelledBodySettings', 'body_panel_settings', 'body_tracking_panel_geometry', 'box_wing_panelled_body_settings', 'frame_fixed_panel_geometry', 'full_panelled_body_settings', 'time_varying_panel_geometry']

class BodyPanelGeometrySettings:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class BodyPanelSettings:
    """
		"""
    reflection_law_settings: ...

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class FullPanelledBodySettings:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

def body_panel_settings(panel_geometry: BodyPanelGeometrySettings, panel_reflection_law: ...=None, panel_type_id: str='') -> BodyPanelSettings:
    ...

def body_tracking_panel_geometry(body_to_track: str, towards_tracked_body: bool, area: float, frame_orientation: str='') -> BodyPanelGeometrySettings:
    ...

def box_wing_panelled_body_settings(length: float, width: float, height: float, solar_array_area: float, box_specular_reflectivity: float, box_diffuse_reflectivity: float, solar_array_specular_reflectivity: float, solar_array_diffuse_reflectivity: float, box_instantaneous_reradiation: bool=True, solar_array_instantaneous_reradiation: bool=True) -> FullPanelledBodySettings:
    ...

def frame_fixed_panel_geometry(surface_normal: numpy.ndarray, area: float, frame_orientation: str='') -> BodyPanelGeometrySettings:
    ...

def full_panelled_body_settings(panel_settings: list[BodyPanelSettings], part_rotation_model_settings: dict[str, ...]={}) -> FullPanelledBodySettings:
    ...

def time_varying_panel_geometry(surface_normal_function: typing.Callable[[], numpy.ndarray], area: float, frame_orientation: str) -> BodyPanelGeometrySettings:
    ...