#! /usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
from obspy.signal.util import next_pow_2
from scipy.fftpack import fft, ifft



USE_FFTW3 = False

def xcorrf(data1, data2, shift=None, shift_zero=0, oneside=False,
           demean=True, window=0, ndat1d=0, ndat2d=0, N1=None, N2=None,
           normalize=True,
           freq_domain=False, transform_back=True,
           stdev1=None, stdev2=None):
    """
    Cross-correlation of numpy arrays data1 and data2 in frequency domain.
    
    We define cross-corelation as:
    xcorr[i] = sum_j (tr1[i+j-shift_zero] * tr2[j])
    The data is demeaned before cross-correlation and the result normalized
    after cross-correlation.
    data1, data2: data
    shift:    maximum samples to shift
              (window for i in the above formula)
    shift_zero: shift tr1 before cross-correlation by this amount of samples to
              the right (this means correlation function is shifted to the
              right or better: the window of what you get of the function
              is shifted to the left)
    oneside:  if True only the right/positive side of the correlation function
              is returned. Overrides parameter shift_zero.
    demean:   if True demean data beforehand
    normalize: if True normalize correlation function
              (1 means perfect correlation)
    window:   Use only data in this window for demeaning and normalizing
              0: window = min(ndat1, ndat2)
              >0: window = this parameter
    ndat1d, ndat2d: If >0 use different values for the length of the arrays when
              calculating the mean (defaults to window parameter)
    return:   numpy array with correlation function of length 2*shift+1 for
              oneside=False and of length shift+1 for oneside=True
    """
    if freq_domain and not transform_back:
        return data1 * np.conjugate(data2)
    elif freq_domain:
        min_size = max(2 * shift + 1 + abs(shift_zero),
                       (N1 + N2) // 2 + shift + abs(shift_zero))
        if len(data1) < min_size:
            raise ValueError('NFFT was not large enough to cover the desired '
                          'xcorr!\nnfft: %d, required minimum: %d' %
                          (len(data1), min_size))
        ret = (ifft(data1 * np.conjugate(data2))).real
    else:
        complex_result = (data1.dtype == np.complex or
                          data2.dtype == np.complex)
        N1 = len(data1)
        N2 = len(data2)

        #if isinstance(data1[0], np.integer) or isinstance(data2[0], np.integer):
        data1 = data1.astype('float64')
        data2 = data2.astype('float64')
        #if (N1-N2)%2==1:
        #    raise ValueError('(N1-N2)%2 has to be 0')
        if window == 0:
            window = min(N1, N2)
        if ndat1d == 0:
            ndat1d = window
        if ndat2d == 0:
            ndat2d = window

        # determine indices for demeaning and normalization
        ind1 = max(0, (N1 - window) // 2)
        ind2 = min(N1, (N1 + window) // 2)
        ind3 = max(0, (N2 - window) // 2)
        ind4 = min(N2, (N2 + window) // 2)


        # demean and normalize data
        if demean:
            data1 -= np.sum(data1[ind1:ind2]) / ndat1d
            data2 -= np.sum(data2[ind3:ind4]) / ndat2d
        if normalize:
            data1 /= np.max(data1[ind1:ind2])
            data2 /= np.max(data2[ind3:ind4])

        # Always use 2**n-sized FFT, perform xcorr
        size = max(2 * shift + 1 + abs(shift_zero),
                   (N1 + N2) // 2 + shift + abs(shift_zero))

        nfft = next_pow_2(size)
        IN1 = fft(data1, nfft)

        if USE_FFTW3:
            IN1 = IN1.copy()
        IN1 *= np.conjugate(fft(data2, nfft))
        ret = ifft(IN1)
        if not USE_FFTW3:
            del IN1
        if not complex_result:
            ret = ret.real
    # shift data for time lag 0 to index 'shift'

    ret = np.roll(ret, -(N1 - N2) // 2 + shift + shift_zero)[:2 * shift + 1]
    # normalize xcorr
    if normalize:
        if not freq_domain:
            stdev1 = (np.sum(data1[ind1:ind2] ** 2)) ** 0.5
            stdev2 = (np.sum(data2[ind3:ind4] ** 2)) ** 0.5
#            stdev1 = (np.sum(data1 ** 2)) ** 0.5
#            stdev2 = (np.sum(data2 ** 2)) ** 0.5
        if stdev1 == 0 or stdev2 == 0:
            log.warning('Data is zero!!')
            ret[:] = 0.
        else:
            ret /= stdev1 * stdev2
    if oneside:
        ret = ret[shift:]
    return np.copy(ret)