import os

if bool(os.getenv("READTHEDOCS")):
    os.environ['PROJ_LIB'] = "/home/docs/checkouts/readthedocs.org/user_builds/tudatpy/share/proj"
else:
    os.environ['PROJ_LIB'] = os.environ['CONDA_PREFIX'] + '/share/proj'

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap

#
# def plot_blue_marble_ground_track(lon: np.ndarray, lat: np.ndarray, lat_0: float = 0, lon_0: float = 0):
#     """Plots a blue marbel ground track.
#
#     Parameters
#     ----------
#     lon : np.ndarray
#         Array of longitude values for plotting the ground track.
#     lat : np.ndarray
#         Array of latitude values for plotting the ground track.
#     lat_0 : float, optional, default=0.0
#         Parallel line (line of latitude) about which the plot is centered.
#     lon_0 : float, optional, default=0.0
#         Meridian line (line of longitude) about which the plot is centered.
#
#     Examples
#     --------
#     .. code-block:: python
#
#         from tudatpy.plotting import *
#         plot_blue_marble_ground_track(lon, lat)
#
#     .. image:: _static/blue-marble-groundtrack.png
#        :width: 400
#        :align: center
#
#     Returns
#     -------
#     None
#
#     """
#     plt.figure(figsize=(8, 8))
#     m = Basemap(projection='ortho', resolution=None, lat_0=lat_0, lon_0=lon_0)
#     m.bluemarble(scale=0.5)
#     x, y = m(lat, lon)
#     m.plot(x, y, color="red", latlon=False, marker='.', linestyle='None')
#
#
# def plot_miller_ground_track(lon, lat, lon_0=0):
#     """Plots a Miller ground track.
#
#     Parameters
#     ----------
#     lon : np.ndarray
#         Array of longitude values for plotting the ground track.
#     lat : np.ndarray
#         Array of latitude values for plotting the ground track.
#     lon_0 : float, optional, default=0.0
#         Meridian line (line of longitude) about which the plot is centered.
#
#     Examples
#     --------
#     .. code-block:: python
#
#         from tudatpy.plotting import *
#         plot_miller_ground_track(lon, lat)
#
#     .. image:: _static/miller-groundtrack.png
#        :width: 400
#        :align: center
#
#     Returns
#     -------
#     None
#
#     """
#     f = plt.figure(figsize=(10, 7.5))
#     m = Basemap(projection="mill", lon_0=lon_0)
#     m.drawcoastlines()
#     m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
#     m.drawmeridians(np.arange(-180, 181, 60), labels=[0, 0, 0, 1])
#     x, y = m(lat, lon)
#     m.plot(x, y, color="red", latlon=False, marker='.', linestyle='None')
