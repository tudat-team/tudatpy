from simulation.environment import body
import unittest
import numpy as np


class TestBody(unittest.TestCase):

    def test1(self):
        test = body(np.array([0,0,0,0,0,0]))
        print(test.get_ephemeris_frame_to_base_frame())