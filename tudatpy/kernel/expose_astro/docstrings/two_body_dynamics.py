from math import sqrt


class BaseState:

    def __init__(self, gravitational_parameter):
        self._gravitational_parameter = gravitational_parameter

    @property
    def state_vector(self):
        raise NotImplementedError("Not implemented in base class.")

    @property
    def mean_motion(self):
        return sqrt(self.gravitational_parameter / abs(self.to_keplerian().semi_major_axis ** 3))

    @property
    def gravitational_parameter(self):
        return self._gravitational_parameter

    @classmethod
    def to_cartesian(cls):
        raise NotImplementedError("Not implemented in base class.")

    @classmethod
    def to_keplerian(cls):
        raise NotImplementedError("Not implemented in base class.")

    @classmethod
    def to_equinoctial(cls):
        raise NotImplementedError("Not implemented in base class.")


class CartesianState(BaseState):

    def __init__(self, gravitation_parameter, position_vector, velocity_vector):
        super().__init__(gravitation_parameter)
        self._position_vector = position_vector
        self._velocity_vector = velocity_vector

    @property
    def position_vector(self):
        return self._position_vector

    @property
    def velocity_vector(self):
        return self._velocity_vector
