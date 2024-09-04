Calculation of clock errors using noise cross-correlations functions (CCFs)

1. create_CCFs.py
performs the preprocessing of the individual traces and the calculation
of the CCFs. For the spectral whitening a function from MSNoise is used.
For the CCF calculation the function 'xcorrf' from Tom Richter (available
on github) is used, as the obspy 'xcorr' function needs too long for a long
shift length.

2. create_pickle.py
determines the clock errors (named delays) of one station and component pair
and stores them together with RCF, SNR and other parameters as pck.


3. create_npz.py
determines the clock errors of one station pair (-> performs an average
over all component pairs) and stores the obtained clock errors together
with the corresponding correlation coefficients as npz. A similar averaging
procedure can be used to average over multiple station pairs.