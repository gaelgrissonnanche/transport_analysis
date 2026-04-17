import numpy as np
from scipy.interpolate import PchipInterpolator
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Load coefficients of thermocouples for E = f(T_Celsius) where E is in Volts
coeff_E_below_273K = np.array([0, 5.8665508708E1, 4.5410977124E-2, -7.7998048686E-4,
                            -2.5800160843E-5, -5.9452583057E-7, -9.3214058667E-9,
                            -1.0287605534E-10, -8.0370123621E-13, -4.3979497391E-15,
                            -1.6414776355E-17, -3.9673619516E-20, -5.5827328721E-23,
                            -3.4657842013E-26])[::-1]*1e-6
coeff_E_above_273K = np.array([0, 5.8665508710E1, 4.5032275582E-2, 2.8908407212E-5,
                            -3.3056896652E-7, 6.5024403270E-10, -1.9197495504E-13,
                            -1.2536600497E-15, 2.1489217569E-18, -1.4388041782E-21,
                            3.5960899481E-25])[::-1]*1e-6

def Sther(T_Kelvin):
    """
    This function returns the Seebeck coefficient of the thermocouple
    concerned (by default type "E") at a certain temperature. The input of the
    function is a temperature in Kelvin, but the coefficients below are for a
    polynomial function with T in Celsius. The output is S in [V / K]
    """
    # Convert T_Kelvin to Celsius
    x = T_Kelvin - 273.15
    # T < 273 K
    E_below = np.poly1d(coeff_E_below_273K)  # is a poly1d object in Volt
    S_below = np.polyder(E_below)  # is a poly1d object in Volt / Celsius
    # T > 273 K
    E_above = np.poly1d(coeff_E_above_273K)  # is a poly1d object in Volt
    S_above = np.polyder(E_above)  # is a poly1d object in Volt / Celsius
    S_values = np.piecewise(x, [x <= 0, x > 0], [S_below, S_above]) # is in Volt / K
    return S_values


def Ether(T_Kelvin):
    """
    This function returns the integrated Seebeck coefficient in temperature
    of the thermocouple concerned (by default type "E") at a certain
    temperature. The input of the function is a temperature in Kelvin,
    but the coefficients below are for a polynomial function with T
    in Celsius. The output is E in Volts
    """
    # Convert T_Kelvin to Celsius
    x = T_Kelvin - 273.15
    # Create a piecewise function of E(T_Celsius) in Volts
    E_values = np.piecewise(x, [x <= 0, x > 0],
                    [lambda x: np.polyval(coeff_E_below_273K,x),
                    lambda x: np.polyval(coeff_E_above_273K,x)])
    return E_values


def Ether_inv(y, niter=3):
    """
    Compute the inverse of the monotonic function Ether(y) for scalar or array inputs.

    The inverse is obtained in two steps: first, a monotonic spline interpolator
    provides an initial guess on a predefined x-range; then a few Newton iterations
    refine the result for improved accuracy.

    Parameters
    ----------
    y : float or array-like
        Value(s) of Ether(x) for which to compute the inverse.
    niter : int, optional
        Number of Newton iterations used to refine the interpolated initial guess.
        Default is 3.

    Returns
    -------
    x : float or ndarray
        Approximate value(s) such that Ether(x) = y.

    Notes
    -----
    This function assumes that Ether(x) is monotonic on the interval [0.0, 10.0].
    Values of y outside the corresponding range of Ether(x) will return invalid
    results because extrapolation is disabled.
    """
    T_min, T_max = 0.0, 663.0 # Kelvin
    T_grid = np.linspace(T_min, T_max, 10000)
    E_grid = Ether(T_grid)

    index = np.argsort(E_grid)
    inv_guess = PchipInterpolator(E_grid[index], T_grid[index], extrapolate=False)

    x = inv_guess(y)
    for _ in range(niter):
        x -= (Ether(x) - y) / Sther(x)
    return x