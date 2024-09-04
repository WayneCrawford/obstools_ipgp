#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
import glob
import os
import time
from obspy import Stream
from obspy import read_inventory
import multiprocessing
import numpy as np
from matplotlib import pyplot as plt
from obspy.core import read, UTCDateTime
from obspy.signal.filter import highpass, lowpass
from scipy import signal
from obspy.signal.util import next_pow_2
from xcorr_trichter import xcorrf # xcorr function of Tom Richter (available on github)
from msnoise.whiten import whiten # whiten function from MSNoise
from datelist_fct import datelist

# Preprocessing and cross-correlating



# ###################### zerophase_chebychev_lowpass_filter #############


def zerophase_chebychev_lowpass_filter(trace, freqmax):
    """
    Custom Chebychev type two zerophase lowpass filter useful for
    decimation filtering.
    This filter is stable up to a reduction in frequency with a factor of
    10. If more reduction is desired, simply decimate in steps.
    Partly based on a filter in ObsPy.
    :param trace: The trace to be filtered.
    :param freqmax: The desired lowpass frequency.
    Will be replaced once ObsPy has a proper decimation filter.

    This code is from LASIF repository (Lion Krischer).
    """
    # rp - maximum ripple of passband, rs - attenuation of stopband
    rp, rs, order = 1, 96, 1e99
    ws = freqmax / (trace.stats.sampling_rate * 0.5)  # stop band frequency
    wp = ws  # pass band frequency

    while True:
        if order <= 12:
            break
        wp *= 0.99
        order, wn = signal.cheb2ord(wp, ws, rp, rs, analog=0)

    b, a = signal.cheby2(order, rs, wn, btype="low", analog=0, output="ba")

    # Apply twice to get rid of the phase distortion.
    trace.data = signal.filtfilt(b, a, trace.data)

# ###################### decimate_trace #############


def decimate_trace(tr, dt):
    """
    Decimate ObsPy Trace (tr) with dt as delta (1/sampling_rate).

    This code is from LASIF repository (Lion Krischer) with some
    minor modifications.
    """
    while True:
        decimation_factor = int(dt / tr.stats.delta)
        # Decimate in steps for large sample rate reductions.
        if decimation_factor > 5:
            decimation_factor = 5
        if decimation_factor > 1:
            new_nyquist = tr.stats.sampling_rate / 2.0 / decimation_factor
            zerophase_chebychev_lowpass_filter(tr, new_nyquist)
            tr.decimate(factor=decimation_factor, no_filter=True)
        else:
            return tr

# ##################### check_par_jobs ######################################


def check_par_jobs(jobs, sleep_time=1):
    """
    check whether all the parallel jobs are finished or not
    :param jobs:
    :param sleep_time:
    :return:
    """
    pp_flag = True
    while pp_flag:
        for proc in jobs:
            if proc.is_alive():
                time.sleep(sleep_time)
                pp_flag = True
                break
            else:
                pp_flag = False
    if not pp_flag:
        print('\n\nAll %s processes are finished...\n' % len(jobs))


def tr_preproc_iterator(date_net_sta_cha_list, overlap, rs_freq, w_freq1,
                        w_freq2, BP_filter, ic, savedir_tr, starti, endi):

    for i in range(starti, endi):
        date = UTCDateTime(date_net_sta_cha_list[i].split('_')[0])
        net = date_net_sta_cha_list[i].split('_')[1]
        sta = date_net_sta_cha_list[i].split('_')[2]
        cha = date_net_sta_cha_list[i].split('_')[3]
        filename = '/path/to/continuous/waveforms/%s/' \
                   'continuous%s/%s.%s.00.%s' % (date.strftime('%Y'),
                                                 date.strftime('%Y-%m-%d'),
                                                 net, sta, cha)
        tr_preproc(filename=filename, overlap=overlap, date=date,
                   rs_freq=rs_freq, w_freq1=w_freq1, w_freq2=w_freq2,
                   BP_filter=BP_filter, ic=ic, savedir_tr=savedir_tr)


def tr_preproc(filename, overlap, date, rs_freq, w_freq1, w_freq2, BP_filter,
               ic, savedir_tr):

    tr_preproc_savedir = savedir_tr + '%s/continuous%s/' \
                         % (date.strftime('%Y'), date.strftime('%Y-%m-%d'))
    scnl = os.path.split(filename)[1].split('.')
    net = scnl[0]
    sta = scnl[1]
    cha = scnl[3]
    tr_preproc_savename = os.path.join(tr_preproc_savedir, '%s_%s_%s.npy'
                                       % (net, sta, cha))
    tr_preproc_savename_windows = os.path.join(tr_preproc_savedir,
                                               '%s_%s_%s_windows.npy'
                                               % (net, sta, cha))

    if os.path.exists(filename) and not
            (os.path.exists(tr_preproc_savename) and
             os.path.exists(tr_preproc_savename_windows))):

        print(filename.split('/')[-2][10:], filename.split('/')[-1])
        st = read(filename, 'MSEED')

        # hourly windows for preprocessing -> 47 windows per day
        r_overlap = 1. - overlap
        times = datelist(date, date+23*3600, r_overlap*3600)
        n_win = len(times)

        tr_all = np.zeros((n_win, int(rs_freq)*3600))
        window_number = []
        window_counter = 0

        # if stream contains gaps, interpolate those
        # with length <= 500 samples (~5/10 sec, depending on sampling_rate)
        if len(st) > 1:
            len_st = len(st)
            gaps = st.get_gaps()
            g_list = []
            for g in xrange(0, len(gaps)):
                if gaps[g][-1] > 500.0:
                    g_list.append(g)
            if len(g_list) == 0:
                st.merge(method=1, fill_value='interpolate')
            else:
                st_merge = Stream()
                g_list = [x+1 for x in g_list]
                g_list.insert(0, 0)
                g_list.append(len_st)
                for g2 in xrange(0, len(g_list)-1):
                    tr_merge = st[g_list[g2]:g_list[g2+1]].merge(
                            method=1, fill_value='interpolate')
                    st_merge.extend(tr_merge)
                st = st_merge
        st_copy = st.copy()

        # preprocessing of each hourly window
        for timeW in times:
            st = st_copy.copy()
            st.trim(timeW, timeW+3600)

            if len(st) == 1:
                if len(st[0]) == st[0].stats.sampling_rate*3600. + 1:
                    st[0].data = st[0].data[:-1]
                # necessary, if some data of the trace hour are missing:
                if len(st[0]) == st[0].stats.sampling_rate*3600.:

                    tr = st[0]
                    sta = tr.stats.station
                    cha = tr.stats.channel
                    net = tr.stats.network

                    # demean & detrend
                    tr.detrend('demean')
                    tr.detrend('linear')

                    # remove response
                    if ic == 'stationXML':
                        stxmldir = '/path/to/STXML/' \
                                   'STXML.%s.%s.00.xml' % (net, sta)

                        stxml = read_inventory(stxmldir, format='stationXML')
                        try:
                            tr.attach_response(stxml)
                        except Exception:
                            break

                        tr.remove_response(output='VEL', water_level=60,
                                           pre_filt=BP_filter, zero_mean=True,
                                           taper=True, taper_fraction=0.05)

                        # demean & detrend
                        tr.detrend('demean')
                        tr.detrend('linear')

                    # bandpass filtering --> use low- and highpass,
                    # since this is more stable than bandpass
                    tr.taper(type='blackman', max_percentage=0.05)
                    tr.data = lowpass(tr.data, freq=10.,
                                      df=tr.stats.sampling_rate, corners=4,
                                      zerophase=True)
                    tr.data = highpass(tr.data, freq=0.01,
                                       df=tr.stats.sampling_rate, corners=4,
                                       zerophase=True)

                    # resampling
                    tr.resample(rs_freq)

                    # clipping to 2*std of amplitudes
                    std = np.std(tr.data)
                    fact = 2.
                    tr.data = np.clip(tr.data, -fact*std, fact*std)

                    # spectral whitening (using whiten function from MSNoise)
                    tr_data_fft = whiten(tr.data, Nfft=next_pow_2(len(tr.data)),
                                         delta=tr.stats.delta, freqmin=w_freq1,
                                         freqmax=w_freq2, plot=False)
                    tr.data = np.real(np.fft.ifft(tr_data_fft)[0:len(tr.data)])

                    # filtering again
                    tr.taper(type='blackman', max_percentage=0.05)
                    tr.data = lowpass(tr.data, freq=w_freq2,
                                      df=tr.stats.sampling_rate, corners=4,
                                      zerophase=True)
                    tr.data = highpass(tr.data, freq=w_freq1,
                                       df=tr.stats.sampling_rate, corners=4,
                                       zerophase=True)

                    # 1-bit normalization
                    tr.data = np.sign(tr.data)

                    tr_all[window_counter, :] = tr.data
                    window_number.append(window_counter)

                window_counter += 1

            else:
                window_counter += 1

        # store the preprocessed hourly windows + information if a window is available
        np.save(tr_preproc_savename, tr_all)
        np.save(tr_preproc_savename_windows, np.asarray(window_number))


def CCF_creator_iterator(date_station_channel_pairs, savedir_tr, savedir_CCF,
                         starti, endi):

    for i in range(starti, endi):
        date = UTCDateTime(date_station_channel_pairs[i].split('+')[0])
        stationpair = date_station_channel_pairs[i].split('+')[1]
        channelpair = date_station_channel_pairs[i].split('+')[2]


        shift_len = 40000 # -> length of CCF: 80001 samples
        saveDir_CCF = savedir_CCF + 'F%s/%s/%s/%s/' \
            % (whiten_freq, date.strftime('%Y'),
               date.strftime('%Y-%m-%d'), channelpair)
        CCF_creator(stationpair=stationpair, channelpair=channelpair,
                    date=date, shift_len=shift_len, savedir_tr=savedir_tr,
                    saveDir_CCF=saveDir_CCF)


def CCF_creator(stationpair, channelpair, date, shift_len, savedir_tr,
                saveDir_CCF):

    savename = os.path.join(saveDir_CCF, '%s_%s.npy' % (stationpair,
                                                        channelpair))
    if not os.path.exists(savename):
        print('%s %s %s' % (stationpair, channelpair, date))

        # check if the waveforms of both stations are available for 
        # cross-correlating
        filedirs = []
        for sc in xrange(0, 2):
            sta_name = stationpair.split('_')[sc]
            cha_name = channelpair[sc]

            filedir = glob.glob(savedir_tr + '%s/continuous%s/*_%s_??%s.npy'
                                % (date.strftime('%Y'),
                                   date.strftime('%Y-%m-%d'), sta_name,
                                   cha_name))

            if len(filedir) > 1:
                raise NameError('Path to file is not unique!')
            if len(filedir) == 1:
                filedirs.append(filedir[0])

        if len(filedirs) == 2:
            # determine the windows that are available in both stations
            for filename in filedirs:
                net = filename.split('/')[-1].split('_')[0]
                sta = filename.split('/')[-1].split('_')[1]
                cha = filename.split('/')[-1].split('_')[2][:3]

                tr_preproc_dir = savedir_tr + '%s/continuous%s/' \
                    % (date.strftime('%Y'),
                       date.strftime('%Y-%m-%d'))
                tr_preproc_name = os.path.join(tr_preproc_dir, '%s_%s_%s.npy'
                                               % (net, sta, cha))
                tr_preproc_name_windows = os.path.join(tr_preproc_dir,
                                                       '%s_%s_%s_windows.npy'
                                                       % (net, sta, cha))

                if filedirs.index(filename) == 0:
                    tr_all1 = np.load(tr_preproc_name)
                    window_number1 = np.load(tr_preproc_name_windows).tolist()

                if filedirs.index(filename) == 1:
                    tr_all2 = np.load(tr_preproc_name)
                    window_number2 = np.load(tr_preproc_name_windows).tolist()

            windows = set(window_number1).intersection(window_number2)
            CCF_all = np.zeros((tr_all1.shape[0], shift_len*2+1))

            # cross-correlation of each window
            for i in windows:
                ccf = xcorrf(tr_all2[i], tr_all1[i], shift=shift_len)
                CCF_all[i, :] = ccf

            # daily CCF as averaged of the hourly windows of 1 day
            CCF = np.sum(CCF_all, axis=0)/len(windows)
            np.save(savename, CCF)

        else:
            # save CCFs with zeros when the waveforms of one or both
            # stations are missing
            CCF = np.zeros(shift_len*2+1)
            np.save(savename, CCF)


if __name__ == '__main__':

    t1 = time.time()

    # list of stationpairs
    stationpairs = ['STA1_STA2', 'STA1_STA3']

    # list of stations
    stations = ['STA1', 'STA2', 'STA3']

    # list of 16 component pairs
    channels = ['1', '2', 'Z', 'H']
    channelpairs = []
    for cha1 in channels:
        for cha2 in channels:
            channelpairs.append('%s%s' % (cha1, cha2))
    channelpairs.sort()

    savedir_tr = '/path/to/preprocessed/traces/PREPROC/'
    savedir_CCF = '/path/to/CCFs/'

    # PARAMETERS
    whiten_freq = '005_05' # freq of 0.05-0.5 Hz used for whitening
    rs_freq = 50. # resampling frequency
    BP_filter = (0.005, 0.01, 9., 10.) # bandpass filter
    date_start = UTCDateTime(2012, 9, 29)
    date_end = UTCDateTime(2013, 11, 29)
    overlap = 0.5 # windows overlap


    # number of parallel processes for trace preprocessing
    tr_np = 30
    # number of parallel processes for CCF calculation
    ccf_np = 30

    ic = 'stationXML' # instrument correction based on stationXML file

    # calculate frequency values from whiten_freq
    f = whiten_freq.split('_')
    if f[0][0] == '0':  
        w_freq1 = '0.' + f[0][1:]
    else:
        w_freq1 = f[0][:] + '.0'

    if f[1][0] == '0':
        w_freq2 = '0.' + f[1][1:]
    else:
        w_freq2 = f[1] + '.0'

    w_freq1 = float(w_freq1)
    w_freq2 = float(w_freq2)


    # make list of all dates
    dates = datelist(date_start, date_end, 24.*3600.)

    # determine all traces to preprocess
    date_net_sta_cha_list = []              
    for c in channels:
        for s in stations:
            net = 'YV'
            for date in dates:
                if c != 'H':
                    date_net_sta_cha = '%s_%s_%s_BH%s' \
                                        % (date.strftime('%Y-%m-%d'), net, s, c,)
                if c == 'H':
                    date_net_sta_cha = '%s_%s_%s_BD%s' \
                                        % (date.strftime('%Y-%m-%d'), net, s, c,)
                date_net_sta_cha_list.append(date_net_sta_cha)
    date_net_sta_cha_list.sort()
           

    # PREPROCESSING
    if tr_np:
        print('Preprocessing all traces...')
        start = 0
        end = len(date_net_sta_cha_list)
        step = (end - start) / tr_np + 1
        step = int(step)

        jobs = []
        for index in xrange(0, tr_np):
            starti = start + index * step
            endi = min(start + (index + 1) * step, end)
            print(index, starti, endi)
            starti = int(starti)
            endi = int(endi)
            p = multiprocessing.Process(target=tr_preproc_iterator,
                                        args=(date_net_sta_cha_list,
                                              overlap, rs_freq, w_freq1,
                                              w_freq2, BP_filter, ic,
                                              savedir_tr, starti, endi))

            jobs.append(p)
        for i in range(len(jobs)):
            jobs[i].start()
        check_par_jobs(jobs)






    # determine all possible station/channel pairs
    date_station_channel_pairs = []
    for s in xrange(0, len(stationpairs)):
        for c in xrange(0, len(channelpairs)):
            for date in dates:
                d_s_c = date.strftime('%Y-%m-%d') + '+' + stationpairs[s] + \
                        '+' + channelpairs[c]
                date_station_channel_pairs.append(d_s_c)
    date_station_channel_pairs.sort()

    for d in dates:
        for c in channelpairs:
            saveDir = savedir_CCF + 'F%s/%s/%s/%s/' % (whiten_freq,
                                                       d.strftime('%Y'),
                                                       d.strftime('%Y-%m-%d'), c)
    # CROSS-CORRELATION
    if ccf_np:
        print('Cross-correlating...')
        start = 0
        end = len(date_station_channel_pairs)
        step = (end - start) / ccf_np + 1
        step = int(step)

        jobs = []
        for index in xrange(0, ccf_np):
            starti = start + index * step
            endi = min(start + (index + 1) * step, end)
            print(index, starti, endi)
            starti = int(starti)
            endi = int(endi)
            p = multiprocessing.Process(target=CCF_creator_iterator,
                                        args=(date_station_channel_pairs,
                                              savedir_tr, savedir_CCF,
                                              starti, endi))
            jobs.append(p)
        for i in range(len(jobs)):
            jobs[i].start()
        check_par_jobs(jobs)

    t2 = time.time()
    print('Elapsed time:', t2-t1, 'sec')
