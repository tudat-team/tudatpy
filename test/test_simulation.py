from simulation.environment.body import Body
import unittest
import numpy as np


class TestBody(unittest.TestCase):

    def test1(self):
        test = Body(np.array([0,0,0,0,0,0]))
        print(test.get_ephemeris_frame_to_base_frame())