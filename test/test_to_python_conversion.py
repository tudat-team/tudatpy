from .test_module import python_wrapper_for_string_vector_solar_system
from .test_module import python_wrapper_for_ints_to_ten
import unittest

class TestConversions(unittest.TestCase):

    def test_string_std_vector_to_python_list(self):
        test = python_wrapper_for_string_vector_solar_system()
        assert(test == ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"])

    def test_string_std_vector_to_python_list(self):
        test = python_wrapper_for_ints_to_ten()
        assert(test == [1,2,3,4,5,6,7,8,9,10])
