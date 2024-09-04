#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import glob
from datelist_fct import datelist
from obspy.core import UTCDateTime
from obspy.signal.cross_correlation import xcorr, xcorr_max
from xcorr_trichter import xcorrf
import os
import time
import multiprocessing
from scipy.interpolate import InterpolatedUnivariateSpline
import cPickle as pickle


# For each station pair and component pair, a pickle (pck) file is 
# created. This file contains information about the clock errors (+ CC)
# of this station/component pair together with RCF, SNR information etc.


# ##################### getSNR ##############################################


def getSNR(signal, sig_start, sig_end, noise_start, noise_end, timevec):
    """
    Determine SNR between signal and noise windows of CCF:
    max(abs(sig))/var(noise)
    signal : 1D numpy array containing correlation or stack trace
    sig_start, sig_end : start and end of signal window, in seconds
    noise_start, noise_end : start and end of noise window, in seconds
    timevec : 1D numpy array with timesteps corresponding to signal trace
    """
    if signal.shape != timevec.shape:
        raise ValueError("signal (%s) and time vector (%s) are not the same \
                         size" % (signal.shape, timevec.shape))

    # find indices timevec corresponding to start and end of signal window
    win_sig = (range(find_nearest(timevec, -sig_end),
                     find_nearest(timevec, -sig_start)) +
               range(find_nearest(timevec, sig_start),
                     find_nearest(timevec, sig_end)))
    # find indices of timevec corresponding to start and end of noise window
    win_noise = (range(find_nearest(timevec, -noise_end),
                       find_nearest(timevec, -noise_start)) +
                 range(find_nearest(timevec, noise_start),
                       find_nearest(timevec, noise_end)))

    noisemax = np.std(signal[win_noise])
    sigmax = max(abs(signal[win_sig]))

    if noisemax == 0:
        SNR = 0
    else:
        SNR = sigmax/noisemax
    return SNR


# ##################### find_nearest ########################################


def find_nearest(array, value):
    """
    find index for nearest value in array
    from:
    http://stackoverflow.com/a/2566508
    """
    idx = (np.abs(array-value)).argmin()
    return idx


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


# ##################### pickle_creator_iterator #############################


def pickle_creator_iterator(station_channel_pairs, filedir_ccf, filedir_pck,
                            starti, endi):

    for i in range(starti, endi):
        stationpair = station_channel_pairs[i].split('+')[0]
        channelpair = station_channel_pairs[i].split('+')[1]
        pickle_creator(stationpair=stationpair, channelpair=channelpair,
                       filedir_ccf=filedir_ccf, filedir_pck=filedir_pck)


# ##################### pickle_creator ######################################


def pickle_creator(stationpair, channelpair, filedir_ccf, filedir_pck):
        
        # directory structure:
        # '/path_to_ccfs/freq/year/date/channelpair/sta1_sta2_channelpair.npy'
        filedirs = glob.glob(filedir_ccf + '/*/*/*/%s_%s.npy' %
                                     (stationpair, channelpair))
        filedirs.sort()
        # get startdate and enddate
        start_filedate = filedirs[0].split('/')[-3]
        end_filedate = filedirs[-1].split('/')[-3]
        # list of dates
        datevec = datelist(UTCDateTime(start_filedate),
                           UTCDateTime(end_filedate), 24.*3600.)
        days = len(filedirs)
        

        # PARAMETERS
        # start and end index for RCF calculation
        # (can be changed for special cases)
        start_idx = 0 
        end_idx = datevec.index(UTCDateTime(2013, 11, 29))

        # largest delay that could be expected
        expec_delay = 1000

        # length of the CCFs
        corrlen = 80001

        # time vector calculation for CCFs
        fs = 50. # sampling rate
        corrlag = np.floor(corrlen/2) / fs
        timevec = np.linspace(-corrlag, corrlag, corrlen)

        s_len = 10 # length of stack, if 1: daily CCFs are used, if 10: 10-day stacks
        stackshift = 1  # shift each stack by stackshift days
        stacklength = s_len  # redundant definition


        sig_start = 0   # start of signal window, in seconds
        sig_end = 400   # end of signal window, in seconds

        noise_start = 650   # start of noise window, in seconds
        noise_end = 800   # end of noise window, in seconds

        SNR_threshold = 4. # SNR threshold


        SNR_list = []
        SNR_list_test = []
        norm_fact_list = []
        avail_list = []
        stack = np.zeros((days, corrlen))
        stack_norm = np.zeros((days, corrlen))
        stack_it = np.zeros((days, corrlen))

        

        for filedir in filedirs:
            CCF = np.load(filedir)
            f_idx = filedirs.index(filedir)

            norm_fact = np.max(np.abs(CCF))
            norm_fact_list.append(norm_fact)
            if norm_fact != 0.:
                CCF_norm = CCF/norm_fact
            else:
                CCF_norm = CCF.copy()

            # SNR calculation
            SNR = getSNR(CCF_norm, sig_start, sig_end, noise_start,
                         noise_end, timevec)
            SNR_list.append(SNR)
            SNR_list_test.append(SNRtest)

            
            if SNR >= SNR_threshold:
                stack_norm[f_idx] = CCF_norm
                stack[f_idx] = CCF
                avail_list.append(1)
            else:
                avail_list.append(0)


        # RCF (ref) calculation 
        ref = (np.sum(stack[start_idx:end_idx+1], axis=0) /
               np.sum(avail_list[start_idx:end_idx+1]))

        ref_norm_fact = np.max(np.abs(ref))
        ref_norm = ref/ref_norm_fact
        ref_startdate = datevec[start_idx]
        ref_enddate = datevec[end_idx]


        # if stacks of several days (with a sliding window of one day)
        # are used instead of daily CCFs
        if s_len != 1:
            Nstack = int(np.floor((stack.shape[0]-stacklength)/stackshift))+1
            slidestack_norm = np.zeros((Nstack, stack.shape[1]))

            avail_list_slide = []
            norm_fact_list_slide = []
            SNR_list_slide = []
            stackmiddate = []
            for iStack in range(Nstack):
                stackstart = iStack * stackshift
                stackend = iStack * stackshift + stacklength
                stackmiddate.append(datevec[int(stackstart)] +
                                    (stacklength/2.)*3600.*24.)
                n_ccfs = np.sum(avail_list[stackstart:stackend])

                if n_ccfs == 0:
                    avail_list_slide.append(0)
                    norm_fact_list_slide.append(0)
                    SNR_list_slide.append(0)
                else:
                    CCF_slide = (np.sum(stack[stackstart:stackend, :],
                                        axis=0) / n_ccfs)
                    norm_fact_slide = np.max(np.abs(CCF_slide))
                    norm_fact_list_slide.append(norm_fact_slide)
                    CCF_norm_slide = CCF_slide/norm_fact_slide

                    # SNR calculation
                    SNR = getSNR(CCF_norm_slide, sig_start, sig_end,
                                 noise_start, noise_end, timevec)

                    SNR_list_slide.append(SNR)

                    if SNR >= SNR_threshold:
                        slidestack_norm[iStack] = CCF_norm_slide
                        avail_list_slide.append(n_ccfs)
                    else:
                        avail_list_slide.append(0)

            datevec = stackmiddate
            norm_fact_list = norm_fact_list_slide
            avail_list = avail_list_slide
            SNR_list = SNR_list_slide
            stack_norm = slidestack_norm

        
        # calculation of the delays (= clock errors)
        delays = []
        CCs = []
        
        for i in xrange(0, stack_norm.shape[0]):
            # delay and CC are both 0, if ref_norm, stack_norm[i] or both are zero
            delay, CC = xcorr(ref_norm, stack_norm[i], shift_len=expec_delay)
            # exclude clock errors with negative CC
            if CC < 0.0:    
                delay = 0.0
                CC = 0.0

            delays.append(delay)
            CCs.append(CC)


        # clock errors resulting in delayed waveforms (-> station clock
        # faster than GPS clock) are defined as positive here
        delays = [-x for x in delays]

        savename = os.path.join(filedir_pck, '%s_%s' % (stationpair, channelpair) +
                                '.pck')
        # store as pickle file
        pickle.dump([stationpair, channelpair, datevec, timevec, fs,
                     delays, CCs, ref_norm, ref_norm_fact, ref_startdate,
                     ref_enddate, norm_fact_list, SNR_list,
                     SNR_threshold, avail_list, stacklength, stackshift],
                    open(savename, "wb"))

 
if __name__ == '__main__':

    t1 = time.time()

    stationpairs = ['STA1_STA2', 'STA1_STA3']

    # list of 16 component pairs
    channels = ['1', '2', 'Z', 'H']
    channelpairs = []
    for cha1 in channels:
        for cha2 in channels:
            channelpairs.append('%s%s' % (cha1, cha2))
    channelpairs.sort()

    filedir_ccf = '/path/to/stored/ccfs'
    filedir_pck = '/path/to/pickle/files'
    
    if not os.path.exists(filedir_pck):
        os.makedirs(filedir_pck)

    # number of parallel processes for pickle calculation
    pickle_np = 32

    # all station pairs and component pairs
    station_channel_pairs = []
    for channelpair in channelpairs:
        for stationpair in stationpairs:
            sta_cha_pair = '%s+%s' % (stationpair, channelpair)
            station_channel_pairs.append(sta_cha_pair)


    # create one pck file for each station pair and component pair
    print('Creating pickles...')
    start = 0
    end = len(station_channel_pairs)
    step = (end - start) / pickle_np + 1
    step = int(step)

    jobs = []
    for index in xrange(0, pickle_np):
        starti = start + index * step
        endi = min(start + (index + 1) * step, end)
        print(index, starti, endi)
        starti = int(starti)
        endi = int(endi)
        p = multiprocessing.Process(target=pickle_creator_iterator,
                                    args=(station_channel_pairs,
                                          filedir_ccf, filedir_pck,
                                          starti, endi))
        jobs.append(p)
    for i in range(len(jobs)):
        jobs[i].start()
    check_par_jobs(jobs)

    t2 = time.time()
    print 'Elapsed time:', t2-t1, 'sec'
