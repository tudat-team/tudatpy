__all__ = ["satellite_propagator"]


class _BaseApplication:

    def __init__(self):
        pass


class SSP(_BaseApplication):

    def __init__(self):
        super(SSP, self).__init__()
