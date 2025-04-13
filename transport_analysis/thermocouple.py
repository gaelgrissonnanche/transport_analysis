import numpy as np
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

def S_ther(T_Kelvin):
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

def E_ther(T_Kelvin):
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