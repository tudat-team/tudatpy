# THIS FILE IS OUT OF DATE
# import unittest
# import numpy as np
# from tudatpy.apps import satellite_propagator
# from tudatpy.kernel import constants
#
#
# class TestTutorials(unittest.TestCase):
#
#     def test_tutorial_1(self):
#         test_output = satellite_propagator.single(
#             start_epoch=0.0,
#             fixed_step_size=10.0,
#             end_epoch=constants.JULIAN_DAY,
#             parent_body="Earth",
#             frame_orientation="ECLIPJ2000",
#             frame_origin="SSB",
#             satellite_name="Delfi-C3",
#             sat_sma=7500.0E3,
#             sat_ecc=0.1,
#             sat_inc=np.deg2rad(85.3),
#             sat_raan=np.deg2rad(23.4),
#             sat_argp=np.deg2rad(235.7),
#             sat_nu=np.deg2rad(139.87),
#             return_output=True,
#             print_output=True,
#             output_type='array'
#         )
#         np.testing.assert_almost_equal(
#             actual=test_output,
#             desired=np.array(
#                 [[+7.037484001334376e+03, -4.560454309804188e+03],
#                  [+3.238059017921996e+03, -1.438318384929149e+03],
#                  [+2.150724187502574e+03, +5.973990914838605e+03],
#                  [-1.465657627242447e+00, -4.550213315831223e+00],
#                  [-4.095839489293598e-02, -2.412544138773537e+00],
#                  [+6.622797609422496e+00, -4.950630568956334e+00]]
#             ).T
#
#         )
#
#     def test_tutorial_2(self):
#         test_output = satellite_propagator.single(
#             start_epoch=0.0,
#             fixed_step_size=10.0,
#             end_epoch=constants.JULIAN_DAY,
#             parent_body="Earth",
#             frame_orientation="ECLIPJ2000",
#             frame_origin="SSB",
#             satellite_name="Delfi-C3",
#             sat_sma=7500.0E3,
#             sat_ecc=0.1,
#             sat_inc=np.deg2rad(85.3),
#             sat_raan=np.deg2rad(23.4),
#             sat_argp=np.deg2rad(235.7),
#             sat_nu=np.deg2rad(139.87),
#             return_output=True,
#             print_output=True,
#             output_type='array'
#         )
#         np.testing.assert_almost_equal(
#             actual=test_output,
#             desired=np.array(
#                 [[+7.037484001334376e+03, -4.560454309804188e+03],
#                  [+3.238059017921996e+03, -1.438318384929149e+03],
#                  [+2.150724187502574e+03, +5.973990914838605e+03],
#                  [-1.465657627242447e+00, -4.550213315831223e+00],
#                  [-4.095839489293598e-02, -2.412544138773537e+00],
#                  [+6.622797609422496e+00, -4.950630568956334e+00]]
#             ).T
#
#         )
