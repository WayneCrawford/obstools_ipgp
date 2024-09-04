#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import cPickle as pickle

# For each stationpair, the clock errors (delays) and correlation
# coefficients (CCs)) are averaged over the 16 component pairs
# according to the formulas given by Hobiger et al., 2012. The 
# averaged clock errors are stored as npz

stationpairs = ['STA1_STA2', 'STA1_STA3']

# list of 16 component pairs
channels = ['1', '2', 'Z', 'H']
channelpairs = []
for cha1 in channels:
    for cha2 in channels:
        channelpairs.append('%s%s' % (cha1, cha2))
channelpairs.sort()

# number of days
n_days = 418

filedir_pck = '/path/to/stored/pickle/files'
filedir_npz = '/path/to/npz/files'

for stationpair in stationpairs:
    delays_all = np.zeros((16, n_days))
    CCs_all = np.zeros((16, n_days))

    sta1 = stationpair.split('_')[0]
    sta2 = stationpair.split('_')[1]


    for channelpair in channelpairs:

        # path to pickle files
        filename_pck = os.path.join(filedir_pck, '%s_%s.pck' %(stationpair, channelpair))

        # load pickle
        [stationpair, channelpair, datevec, timevec, fs, delays, CCs,
         ref_norm, ref_norm_fact, ref_startdate, ref_enddate, norm_fact_list,
         SNR_list, SNR_threshold, avail_list, stacklength, stackshift] = \
            pickle.load(open(filename_pck, "rb"))

        # delays and CCs as numpy array
        delays = np.asarray(delays)
        CCs = np.asarray(CCs)

        i = channelpairs_all.index(channelpair)
        delays_all[i] = delays
        CCs_all[i] = CCs


    delays_it = []  
    CCs_it = [] 

    # reject CCFs with CC<0.85*CCav
    # CCav is the average CC of entire deployment period for a given
    # station pair and component pair
    for j1 in xrange(0, delays_all.shape[0]):
        # as delays_comp isn't a copy of delays_all[j1],
        # changes will also be applied to delays_all
        delays_comp = delays_all[j1] 
        CCs_comp = CCs_all[j1]
        counter = 0
        for j2 in xrange(0, len(CCs_comp)):
            if CCs_comp[j2] != 0.0:
                counter = counter + 1
        CC_avv = np.sum(CCs_comp)/counter
        for j3 in xrange(0, len(CCs_comp)):
            if CCs_comp[j3] < 0.85*CC_avv:
                delays_comp[j3] = 0.0
                CCs_comp[j3] = 0.0
    
    # average over 16 component pairs according to the formulas given
    # by Hobiger et al., 2012, the same formula can be used to average
    # over different station pairs (replace 16 by number of station pairs)
    n_comps = []
    for i1 in xrange(0, n_days, 1):
        add1 = 0.0
        add2 = 0.0
        add3 = 0.0
        n_comps.append(len(CCs_all[:,i1][CCs_all[:,i1]!=0.0]))
        for i2 in xrange(0, 16, 1):
            add1 = add1 + (np.power(CCs_all[i2][i1], 2) * delays_all[i2][i1])
            add2 = add2 + np.power(CCs_all[i2][i1], 2)
            add3 = add3 + np.power(CCs_all[i2][i1], 3)
        if add2 != 0:
            delays_it.append(add1/add2)
            CCs_it.append(add3/add2)
        else:
            delays_it.append(0)
            CCs_it.append(0)

    # saving the average delays and CCs
    savename = os.path.join(filedir_npz, '%s.npz' % (stationpair))
    np.savez(savename, delays=delays_it, CCs=CCs_it, n_comps=n_comps)