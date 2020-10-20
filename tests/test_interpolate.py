from tudatpy.kernel.math import interpolators
import numpy as np

x = np.array([0, 1, 2, 3, 4])
y = np.sin(x * 2 * np.pi / np.max(x))

f_linear = interpolators.interp1d(x, y)
f_lagrange = interpolators.interp1d(x, y, kind='lagrange')
# f_hermite_spline = interpolators.interp1d(x, y, kind="hermite_spline")
f_cubic_spline = interpolators.interp1d(x, y, kind="cublic_spline")

plt.plot(x_interp, f_lagrange(x_interp).T, label="tudat_lagrange")
# plt.plot(x_interp, f_hermite_spline(x_interp).T, label="tudat_hermite_spline")
plt.plot(x_interp, f_cubic_spline(x_interp).T, label="tudat_cubic_spline")
plt.plot(x_interp, f_linear(x_interp).T, label="tudat_linear")

from scipy import interpolate

f_linear_scipy = interpolate.interp1d(x, y)
plt.plot(x_interp, f_linear_scipy(x_interp).T, label="scipy_linear")
plt.legend()
plt.show()
