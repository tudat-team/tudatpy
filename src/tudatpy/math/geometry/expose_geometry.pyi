import typing
__all__ = ['Capsule', 'CompositeSurfaceGeometry', 'SurfaceGeometry']

class Capsule(CompositeSurfaceGeometry):

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, nose_radius: float, middle_radius: float, rear_length: float, rear_angle: float, side_radius: float) -> None:
        ...

    @property
    def length(self) -> float:
        ...

    @property
    def middle_radius(self) -> float:
        ...

    @property
    def volume(self) -> float:
        ...

class CompositeSurfaceGeometry(SurfaceGeometry):

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class SurfaceGeometry:

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...