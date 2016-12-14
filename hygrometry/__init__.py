"""
hygrometry module.

To install:

`pip install -e git+https://github.com/bluthen/hygrometry.git@1.0.0`

Sources and Citations:

.. [BriceHalls] Tim Brice, Todd Halls. Javascript wet bulb calculator. http://www.srh.noaa.gov/epz/?n=wxcalc_rh
.. [SensHumGl] Sensirion. Humidity at a Glance. https://www.sensirion.com/fileadmin/user_upload/customers/sensirion/Dokumente/Humidity_Sensors/Sensirion_Humidity_Sensors_Introduction_to_Relative_Humidity_V2.pdf
.. [SensIntroHum] Sensirion. Introduction to Humidity. https://www.sensirion.com/fileadmin/user_upload/customers/sensirion/Dokumente/Humidity_Sensors/Sensirion_Humidity_Sensors_at_a_Glance_V1.pdf
.. [NOAAHeatIndex] The Heat Index Equation. http://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml
.. [PlanetAbsHum] PlanetCalc. Relative humidity to absolute humidity and vise versa calculators. http://planetcalc.com/2167/

:Example:
    >>> import hygrometry
    >>> hygrometry.conv_f2c(70.0)
    21.111128

"""

import math


def conv_f2c(f):
    """
    Convert fahrenheit to Celsius

    :param f: Temperature in Fahrenheit
    :type f: float
    :return: Temperature in Celsius
    :rtype: float
    :Example:

        >>> import hygrometry
        >>> hygrometry.conv_f2c(70.0)
        21.111128

    """
    return (f - 32.0) * 0.555556


def conv_c2f(c):
    """
    Convert Celsius to Fahrenheit

    :param c: Temperature in Celsius
    :type c: float
    :return: Temperature in Fahrenheit
    :rtype: float
    :Example:

        >>> import hygrometry
        >>> hygrometry.conv_c2f(21.111128)
        70.0000304

    """
    return c*1.8+32.0


def calc_es_v_dew(t_c, rh):
    """
    Calculates vapor pressure, saturation vapor pressure, and dew point in celcius. See [BriceHalls]_

    :param t_c: Temperature in Celcius.
    :type t_c: float
    :param rh: Relative humdity 0-100
    :type rh: float
    :return: [vapor pressure, saturation vapor pressure, dew point]
    :rtype: [float, float, float]
    :Example:

        >>> import hygrometry
        >>> hygrometry.calc_es_v_dew(20.1, 50.3)
        [23.514683799663736, 11.827885951230858, 9.451033779734948]

    """
    t_c = float(t_c)
    rh = float(rh)
    es = 6.112 * math.exp(17.67 * t_c / (t_c + 243.5))
    v = es * rh / 100.0
    if v != 0:
        l = math.log(v / 6.112)
        dew_c = 243.5 * l / (17.67 - l)
    else:
        dew_c = float('nan')
    return [es, v, dew_c]


def calc_wb(e_difference, t_w_guess, c_temp, mb_pressure, e2, previoussign, incr):
    """
    Incremental wetbulb calculation. See [BriceHalls]_

    Recommend not to use directly, use wetbulb() instead
    """
    # global debug_counts
    count = 0

    # Make sure everything is a float
    e_difference = float(e_difference)
    t_w_guess = float(t_w_guess)
    c_temp = float(c_temp)
    mb_pressure = float(mb_pressure)
    e2 = float(e2)
    previoussign = previoussign
    incr = float(incr)

    while math.fabs(e_difference) > 0.0005:
        e_w_guess = 6.112 * math.exp((17.67 * t_w_guess) / (t_w_guess + 243.5))
        e_guess = e_w_guess - mb_pressure * (c_temp - t_w_guess) * 0.00066 * (1.0 + (0.00115 * t_w_guess))
        e_difference = e2 - e_guess

        # Had to change this from e_difference == 0
        if math.fabs(e_difference) < 0.0005:
            break
        else:
            if e_difference < 0.0:
                cursign = -1
                if cursign != previoussign:
                    previoussign = float(cursign)
                    incr /= 10.0
                else:
                    incr = incr
            else:
                cursign = 1
                if cursign != previoussign:
                    previoussign = cursign
                    incr /= 10.0
                else:
                    incr = incr

        t_w_guess = float(t_w_guess) + float(incr) * float(previoussign)
        count += 1
    return t_w_guess


def wetbulb(t_c, rh, p):
    """
    Calculate Wet bulb in Celcius given temperature, relative humidity, and pressure.
    See: [BriceHalls]_

    :param t_c: Temperature in Celcius.
    :type t_c: float
    :param rh: Relative humdity 0-100
    :type rh: float
    :param p: Pressure in hPa
    :type p: float
    :Example:

        >>> import hygrometry
        >>> hygrometry.wetbulb(40, 50, 3)
        27.62960000000001
 
    """
    rh = float(rh)
    if rh <= 0:
        rh = 0.
    [es, v, dew_c] = calc_es_v_dew(float(t_c), rh)
    return calc_wb(1., 0., float(t_c), float(p), float(v), 1, 10.)


def humidity_adjust_temp(rh1, t_c1, t_c2):
    """
    Gives you would the relative humidity would be if just the temperature changed. See: [SensIntroHum]_

    :param rh1: Initial relative humidity 0-100
    :type rh1: float
    :param t_c1: Initial temperature in Celsius.
    :type t_c1: float
    :param t_c2: The temperature to find the new RH at.
    :type t_c2: float
    :return: The adjusted RH (0-100) at Temperature t_c2
    :rtype: float
    :Example:

        >>> import hygrometry
        >>> hygrometry.humidity_adjust_temp(60, 25, 30)
        44.784059201238314

    """
    rh2 = rh1*math.exp(4283.78*(t_c1-t_c2)/(243.12+t_c1)/(243.12+t_c2))
    return rh2


def dew(t_c, rh):
    """
    Gives you the Dew point given a RH at a Temperature. See [SensIntroHum]_

    :param t_c: Temperature in Celsius
    :type t_c: float
    :param rh: Relative Humidity 0-100
    :type rh: float
    :return: Dew point in Celsius
    :rtype: float
    :Example:
        >>> import hygrometry
        >>> hygrometry.dew(25, 60)
        16.693149006198954

    """
    t_n = 243.12  # C
    m = 17.62
    h = math.log(rh/100.0) + (m*t_c)/(t_n+t_c)
    dew_c = t_n * h / (m-h)
    return dew_c


def absolute_humidity(t, rh):
    """
    Gives you the mass of water vapor in volume of dry air. Units in g/m^3 See [SensIntroHum]_

    Different pressure seem to affect absolute humidity slightly. For a more accurate calculation that uses pressure,
    see [PlanetAbsHum]_.

    :param t: Temperature in Celsius.
    :type t: float
    :param rh: Relative Humidity 0-100
    :return: Absolute humidity g/m^3
    :rtype: float
    :Example:

        >>> import hygrometry
        >>> hygrometry.absolute_humidity(25, 60)
        13.780667458722558

    """
    t_n = 243.12  # C
    m = 17.62
    a = 6.112  # hPa
    
    dv = 216.7*(rh/100.0*a*math.exp(m*t/(t_n+t))/(273.15+t))
    return dv


def mixing_ratio(t, rh, p):
    """
    Gives you the mixing ratio in g/kg. See [SensIntroHum]_

    :param t: Temperature in Celsius
    :type t: float
    :param rh: Relative humidity 0-100
    :type rh: float
    :param p: Barometric Air pressure in hPa
    :type p: float
    :return: Mixing ratio g/kg
    :rtype: float

    :Example:

        >>> import hygrometry
        >>> hygrometry.mixing_ratio(30, 80, 980)
        22.266502023175242

    """
    t_n = 243.12  # C
    m = 17.62
    a = 6.112  # hPa

    e = rh/100.0*a*math.exp(m*t/(t_n+t))
    r = 622.0*e/(p-e)
    return r


def heat_index(t, rh):
    """
    Gives you the heat index in Celsius. See [NOAAHeatIndex]_

    :param t: Temperature in Celsius
    :type t: float
    :param rh: Relative humidity 0-100
    :type rh: float
    :return: Heat index in Celsius.
    :rtype: float
 
    :Example:

        >>> import hygrometry
        >>> hygrometry.heat_index(25, 80)
        25.644464960000008

    """
    t_f = conv_c2f(t)
    h_i = -42.379 + 2.04901523 * t_f + 10.14333127 * rh - .22475541 * t_f * rh - .00683783 * t_f * t_f - .05481717 * \
          rh * rh + .00122874 * t_f * t_f * rh + .00085282 * t_f * rh * rh - .00000199 * t_f * t_f * rh * rh
    if rh < 13 and t_f > 80 and t_f < 112:
        adj = ((13.0-rh)/4.0)*math.sqrt((17.0-math.fabs(t_f-95.))/17.0)
        h_i -= adj
    elif rh > 85 and t_f > 80 and t_f < 87:
        adj = ((rh-85.0)/10.0) * ((87.0-t_f)/5.0)
        h_i += adj
    elif t_f < 80:
        h_i = 0.5 * (t_f + 61.0 + ((t_f-68.0)*1.2) + (rh*0.094))
    return conv_f2c(h_i)

