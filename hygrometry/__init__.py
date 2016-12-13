"""
hygrometry module.

Sources and Citations:

.. [BriceHalls] Tim Brice, Todd Halls. Javascript wet bulb calculator. http://www.srh.noaa.gov/epz/?n=wxcalc_rh
.. [SensHumGl] Sensirion. Humidity at a Glance. https://www.sensirion.com/fileadmin/user_upload/customers/sensirion/Dokumente/Humidity_Sensors/Sensirion_Humidity_Sensors_Introduction_to_Relative_Humidity_V2.pdf
.. [SensIntroHum] Sensirion. Introduction to Humidity. https://www.sensirion.com/fileadmin/user_upload/customers/sensirion/Dokumente/Humidity_Sensors/Sensirion_Humidity_Sensors_at_a_Glance_V1.pdf
.. [NOAAHeatIndex] The Heat Index Equation. http://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml
.. [PlanetAbsHum] PlanetCalc. Relative humidity to absolute humidity and vise versa calculators. http://planetcalc.com/2167/
:Example:
>>> import hygrometry

"""

import math


def conv_F2C(f):
    """
    Convert fahrenheit to Celsius

    :param f: Temperature in Fahrenheit
    :type f: float
    :return: Temperature in Celsius
    :rtype: float
    :Example:

    >>> hygrometry.conv_F2C(70.0)
    21.111128

    """
    return (f - 32.0) * 0.555556


def conv_C2F(c):
    """
    Convert Celsius to Fahrenheit

    :param c: Temperature in Celsius
    :type c: float
    :return: Temperature in Fahrenheit
    :rtype: float
    :Example:

    >>> hygrometry.conv_C2F(21.111128)
    70.0000304

    """
    return c*1.8+32.0


# Python implementation of wetbulb calculation by Tim Brice and Todd Hall
# http://www.srh.noaa.gov/epz/?n=wxcalc_rh

def Td(Tc, RH):
    """
    Calculates vapor pressure, saturation vapor pressure, and dew point in celcius. See [BriceHalls]_

    :param Tc: Temperature in Celcius.
    :type Tc: float
    :param RH: Relative humdity 0-100
    :type RH: float
    :return: [vapor pressure, saturation vapor pressure, dew point]
    :rtype: [float, float, float]
    :Example:

    >>> hygrometry.Td(20.1, 50.3)
    [23.514683799663736, 11.827885951230858, 9.451033779734948]

    """
    Tc = float(Tc)
    RH = float(RH)
    es = 6.112 * math.exp(17.67 * Tc / (Tc + 243.5))
    V = es * RH / 100.0
    if V != 0:
        l = math.log(V / 6.112)
        Td = 243.5 * l / ( 17.67 - l)
    else:
        Td = float('nan')
    return [es, V, Td]


#debug_counts=[]
def calc_wb(Edifference, Twguess, Ctemp, MBpressure, E2, previoussign, incr):
    """
    Incremental wetbulb calculation. See [BriceHalls]_

    Recommend not use directly, use calc_sWB() instead
    """
    #global debug_counts
    count = 0

    #Make sure everything is a float
    Edifference = float(Edifference)
    Twguess = float(Twguess)
    Ctemp = float(Ctemp)
    MBpressure = float(MBpressure)
    E2 = float(E2)
    previoussign = previoussign
    incr = float(incr)

    while math.fabs(Edifference) > 0.0005:
        Ewguess = 6.112 * math.exp((17.67 * Twguess) / (Twguess + 243.5))
        Eguess = Ewguess - MBpressure * (Ctemp - Twguess) * 0.00066 * (1.0 + (0.00115 * Twguess))
        Edifference = E2 - Eguess
        #print Edifference

        #print Twguess, Ewguess, E2, Eguess, Edifference, count

        #Had to change this from Edifference == 0
        if math.fabs(Edifference) < 0.0005:
            break
        else:
            if Edifference < 0.0:
                cursign = -1
                if cursign != previoussign:
                    previoussign = float(cursign)
                    incr = incr/10.0
                else:
                    incr = incr
            else:
                cursign = 1
                if cursign != previoussign:
                    previoussign = cursign
                    incr = incr/10.0
                else:
                    incr = incr

        Twguess = float(Twguess) + float(incr) * float(previoussign)
        count += 1
        #if count > 15:
        #    break
    wetbulb = Twguess
    #print wetbulb
    #debug_counts.append(count)
    #print "Count %d" % (count,)
    return wetbulb


def calc_sWB(Tc, RH, P):
    """
    Calculate Wet bulb in Celcius given temperature, relative humidity, and pressure.
    See: [BriceHalls]_

    :param Tc: Temperature in Celcius.
    :type Tc: float
    :param RH: Relative humdity 0-100
    :type RH: float
    :param P: Pressure in hPa
    :type P: float
    :Example:

    >>> hygrometry.calc_sWB(40, 50, 3)
    27.62960000000001
 
    """
    RH = float(RH)
    if RH <= 0:
        RH = 0.
    [es, V, Tdv] = Td(float(Tc), RH)
    return calc_wb(1., 0., float(Tc), float(P), float(V), 1, 10.)


def humidity_adjust_temp(RH_1, Tc_1, Tc_2):
    """
    Gives you would the relative humidity would be if just the temperature changed. See: [SensIntroHum]_

    :param RH_1: Initial relative humidity 0-100
    :type RH_1: float
    :param Tc_1: Initial temperature in Celsius.
    :type Tc_1: float
    :param Tc_2: The temperature to find the new RH at.
    :type Tc_2: float
    :return: The adjusted RH (0-100) at Temperature Tc_2
    :rtype: float
    :Example:

    >>> hygrometry.humidity_adjust_temp(60, 25, 30)
    44.784059201238314

    """
    RH_2 = RH_1*math.exp(4283.78*(Tc_1-Tc_2)/(243.12+Tc_1)/(243.12+Tc_2));
    return RH_2

def dew(Tc, RH):
    """
    Gives you the Dew point given a RH at a Temperature. See [SensIntroHum]_

    :param Tc: Temperature in Celsius
    :type Tc: float
    :param RH: Relative Humidity 0-100
    :type RH: float
    :return: Dew point in Celsius
    :rtype: float
    :Example:
    >>> hygrometry.dew(25, 60)
    16.693149006198954

    """
    Tn = 243.12 # C
    m = 17.62
    H = math.log(RH/100.0) + (m*Tc)/(Tn+Tc)
    Td = Tn * H / (m-H)
    return Td

def absolute_humidity(t, RH):
    """
    Gives you the mass of water vapor in volume of dry air. Units in g/m^3 See [SensIntroHum]_

    Different pressure seem to affect absolute humidity slightly. For a more accurate calculation that uses pressure, see [PlanetAbsHum]_.

    :param t: Temperature in Celsius.
    :type t: float
    :param RH: Relative Humidity 0-100
    :return: Absolute humidity g/m^3
    :rtype: float
    :Example:

    >>> hygrometry.absolute_humidity(25, 60)
    13.780667458722558

    """
    Tn = 243.12 # C
    m = 17.62
    A = 6.112 # hPa
    
    dv = 216.7*(RH/100.0*A*math.exp(m*t/(Tn+t))/(273.15+t));
    return dv

def mixing_ratio(t, RH, p):
    """
    Gives you the mixing ratio in g/kg. See [SensIntroHum]_

    :param t: Temperature in Celsius
    :type t: float
    :param RH: Relative humidity 0-100
    :type RH: float
    :param p: Barometric Air pressure in hPa
    :type p: float
    :return: Mixing ratio g/kg
    :rtype: float

    :Example:

    >>> hygrometry.mixing_ratio(30, 80, 980)
    22.266502023175242

    """
    Tn = 243.12 # C
    m = 17.62
    A = 6.112 # hPa

    e = RH/100.0*A*math.exp(m*t/(Tn+t))
    r = 622.0*e/(p-e)
    return r

def heat_index(t,RH):
    """
    Gives you the heat index in Celsius. See [NOAAHeatIndex]_

    :param t: Temperature in Celsius
    :type t: float
    :param RH: Relative humidity 0-100
    :type RH: float
    :return: Heat index in Celsius.
    :rtype: float
 
    :Example:

     >>> hygrometry.heat_index(25, 80)
     25.644464960000008

    """
    tF = conv_C2F(t)
    HI = -42.379 + 2.04901523*tF + 10.14333127*RH - .22475541*tF*RH - .00683783*tF*tF - .05481717*RH*RH + .00122874*tF*tF*RH + .00085282*tF*RH*RH - .00000199*tF*tF*RH*RH
    if RH < 13 and tF > 80 and tF < 112:
        adj = ((13.0-RH)/4.0)*math.sqrt((17.0-math.fabs(tF-95.))/17.0)
        HI = HI - adj
    elif RH > 85 and tF > 80 and tF < 87:
        adj = ((RH-85.0)/10.0) * ((87.0-tF)/5.0)
        HI = HI + adj
    elif tF < 80:
        HI = 0.5 * (tF + 61.0 + ((tF-68.0)*1.2) + (RH*0.094))
    return conv_F2C(HI)

