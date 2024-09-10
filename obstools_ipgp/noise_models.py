#!/usr/bin/env python3
import scipy.signal as ssig
"""
Various noise models:
- Sismo (low and high).  My version (S_low and S_high) is Peterson et al.
  down to 0.1 s (Hz), Wolin and McNamara above XX Hz 
- Pressure (low and high) hand picked from Brown et al., 2014 out to 10 s,
  Webb 1998 for longer periods
- SHOULD ADD 
- 2-pole high-pass filters corresponding to geophone and hydrophone cutoffs
"""
import numpy as np

#        Period    Level       Slope
LPAB = [[0.10,    -162.36,     5.64],
        [0.17,    -166.7,      0.00],
        [0.40,    -170.0,     -8.30],
        [0.80,    -166.4,     28.90],
        [1.24,    -168.60,    52.48],
        [2.40,    -159.98,    29.81],
        [4.30,    -141.10,     0.00],
        [5.00,     -71.36,   -99.77],
        [6.00,     -97.26,   -66.49],
        [10.00,   -132.18,   -31.57],
        [12.00,   -205.27,    36.16],
        [15.60,    -37.65,  -104.33],
        [21.90,   -114.37,   -47.10],
        [31.60,   -160.58,   -16.28],
        [45.00,   -187.50,     0.00],
        [70.00,   -216.47,    15.70],
        [101.00,  -185.00,     0.00],
        [154.00,  -168.34,    -7.61],
        [328.00,  -217.43,    11.90],
        [600.00,  -258.28,    26.60],
        [10000.0, -346.88,    48.75],
        [100000,  -346.88,    48.75]]

HPAB = [[0.10,    -108.73,   -17.23],
        [0.22,    -150.34,   -80.50],
        [0.32,    -122.31,   -23.87],
        [0.80,    -116.85,    32.51],
        [3.80,    -108.48,    18.08],
        [4.60,     -74.66,   -32.95],
        [6.30,       0.66,  -127.18],
        [7.90,     -93.37,   -22.42],
        [15.40,     73.54,  -162.98],
        [20.00,   -151.52,    10.01],
        [354.80,  -206.66,    31.63],
        [100000,  -206.66,    31.63]]

# Seismo signal levels
#        Period    Level
S_low = [[0.01,    -140],
         [0.05,    -150],
         [0.10,    -168],
         [0.17,    -166.7],
         [0.40,    -166.7],
         [0.80,    -169.2],
         [1.24,    -163.7],
         [2.40,    -148.6],
         [4.30,    -141.1],
         [5.00,    -141.1],
         [6.00,    -149.0],
         [10.00,   -163.75],
         [12.00,   -166.25],
         [15.60,   -162.1],
         [21.90,   -177.5],
         [31.60,   -185.0],
         [45.00,   -187.5],
         [70.00,   -187.5],
         [101.00,  -185.00],
         [154.00,  -185.0],
         [328.00,  -187.5],
         [600.00,  -184.4],
         [10000.0, -151.9],
         [100000,  -103.1]]

S_high = [[0.01,     -90],
          [0.03,     -90],
          [0.06,     -88],
          [0.10,     -91.5],
          [0.22,     -97.4],
          [0.32,    -110.5],
          [0.80,    -120.0],
          [3.80,     -98.0],
          [4.60,     -96.5],
          [6.30,    -101.0],
          [7.90,    -113.5],
          [15.40,   -120.0],
          [20.00,   -138.5],
          [354.80,  -126.0],
          [100000,   -48.5]]


# Pressure signal levels, from Brown et al. 2014 (<=10 s) and Webb 1998 (>10 s)
#        Period    PSD (re uPa^2/Hz)
P_high = [[0.01,      82],
          [0.02,      89],
          [0.05,      95],
          [0.10,      95],
          [0.20,     100],
          [0.50,     105],
          [1.0,      118],
          [2.0,      127],
          [5.0,      154],
          [10.0,     140],
          [20.0,     150],
          [50.0,     165],
          [100.0,    170],
          [1000,     180]]

P_low = [[0.01,     64],
         [0.02,     72.5],
         [0.05,     71],
         [0.10,     71],
         [0.20,     72],
         [0.50,     90],
         [1.0,      95],
         [2.0,     113],
         [5.0,     130],
         [7.0,     111],
         [10.0,    114],
         [20.0,     90],
         [50.0,    100],
         [100.0,   130],
         [1000,    150]]


def PetersonNoiseModel(periods, as_freqs=False):
    """
    Return Peterson low and high seismological noise models, extended by
    Wolin and McNamara from >10 to 100 Hz

    Args:
        periods (list): periods to use (should be increasing)
        as_freqs (bool): interpret "periods" as frequencies instead
    Returns:
        models (list): Noise models in  dB ref 1 (m/s^2)^2/Hz
            low_noise (list)
            high_noise (list)
    """
    return _from_NoiseModel(S_low, S_high, periods, as_freqs)


def PressureNoiseModel(periods, as_freqs=False):
    """
    Return Brown low and high presure noise models in dB ref to 1 Pa^2/Hz,
    extended by Webb 1998 for periods > 10s

    Args:
        periods (list): periods to use (should be increasing)
        as_freqs (bool): interpret "periods" as frequencies instead
    """
    return _from_NoiseModel([[x[0], x[1]-120] for x in P_low],
                            [[x[0], x[1]-120] for x in P_high],
                            periods, as_freqs)


def two_pole_HP(periods, corner_freq, noise_level):
    """
    Return high-pass filter effect on data

    (multipling noise by the correction factor)

    Modified from obspy.paz_to_freq_resp
    """
    f = [1./x for x in periods]
    zeros = [0., 0.]
    p = corner_freq*2*np.pi
    poles = [p, p]
    b, a = ssig.zpk2tf(zeros, poles, noise_level)
    if not isinstance(a, np.ndarray) and a == 1.0:
        a = [1.0]
    _w, h = ssig.freqs(b, a, f * 2 * np.pi)
    return [20*np.log10(1./x) for x in h]


def _from_NoiseModel(low_base, high_base, periods, as_freqs):
    """
    Args:
        periods (list): periods to use (should be increasing)
        as_freqs (bool): interpret "periods" as frequencies instead
    """
    if not as_freqs:
        lownoise = _fit_points(periods, low_base)
        highnoise = _fit_points(periods, high_base)
    else:
        periods = np.power(periods[::-1], -1)
        lownoise = _fit_points(periods, low_base)
        highnoise = _fit_points(periods, high_base)
        lownoise = lownoise[::-1]
        highnoise = highnoise[::-1]
    return lownoise, highnoise


def _fit_points(periods, model):
    """
    Fit points to a noise model

    Args:
        periods (list): periods in increasing order
        model (list of lists): list of [period, value, slope]
    """
    x = np.log10(periods)
    xp = np.log10([x[0] for x in model])
    yp = [x[1] for x in model]
    assert np.all(np.diff(x) > 0), 'x is not increasing'
    assert np.all(np.diff(xp) > 0), 'xp is not increasing'
    return np.interp(x, xp, yp, left=np.nan, right=np.nan)


if __name__ == '__main__':
    print('Printing Peterson as non-line model')
    print('Peterson low')
    print([x[1] + x[2]*np.log10(x[0]) for x in LPAB])
    print('Peterson high')
    print([x[1] + x[2]*np.log10(x[0]) for x in HPAB])
