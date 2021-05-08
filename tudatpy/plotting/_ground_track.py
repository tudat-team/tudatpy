import os

os.environ['PROJ_LIB'] = os.environ['CONDA_PREFIX'] + '/share/proj'  # Required fix for Jupyter
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def plot_blue_marble_ground_track(lon, lat, lat_0=0, lon_0=0):
    plt.figure(figsize=(8, 8))
    m = Basemap(projection='ortho', resolution=None, lat_0=lat_0, lon_0=lon_0)
    m.bluemarble(scale=0.5)
    x, y = m(lat, lon)
    m.plot(x, y, color="red", latlon=False, marker='.', linestyle='None')


def plot_miller_ground_track(lon, lat, lon_0=0):
    f = plt.figure(figsize=(10, 7.5))
    m = Basemap(projection="mill", lon_0=lon_0)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 60), labels=[0, 0, 0, 1])
    x, y = m(lat, lon)
    m.plot(x, y, color="red", latlon=False, marker='.', linestyle='None')
