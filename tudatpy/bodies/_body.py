from ..core.simulation_setup import Body as _Body
import numpy as np


class Body(_Body):

    def __init__(self, name=None, state=np.zeros(6)):
        """

        Parameters
        ----------
        state
        """
        super().__init__(state)
        self._name = name

    @property
    def state(self):
        """
        Current state of the body.

        Returns
        -------
        np.ndarray[dim=6]
        """
        return self.get_state()

    @state.setter
    def state(self, x):
        self.set_state(x)

    @property
    def inertia_tensor(self):
        return self.get_body_inertia_tensor()

    @inertia_tensor.setter
    def inertia_tensor(self, x):
        self.set_body_inertia_tensor(x)

    @property
    def atmosphere_model(self):
        return self.get_atmosphere_model()

    @atmosphere_model.setter
    def atmosphere_model(self, x):
        self.set_atmosphere_model(x)
