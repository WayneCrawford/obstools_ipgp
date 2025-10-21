import math
from pathlib import Path
import re
import argparse
import sys

from obspy.clients.filesystem.sds import Client
from obspy.core.inventory import read_inventory
from obspy.signal.spectral_estimation import PPSD
from matplotlib import pyplot as plt

from .noise_models import PetersonNoiseModel, PressureNoiseModel, two_pole_HP

MED_PSD_DIR = 'plot_medPSDs'
PPSD_DIR = 'plot_PPSDs'


def main(save_ppsds=False, plot_spectrograms=False, verbose=True):

    args = _parse_input()

    client = Client(args.sds_dir)
    inv = read_inventory(args.inv_file)

    # Get lists of "NET.STA.LOC.CHAN"
    client_nslcs = _client_get_nslc(client, args)
    inv_nslcs = _inv_get_nslc(inv, args)  # list of 4-tuples, too
    shared_nslcs = _get_shared_nslcs(client_nslcs, inv_nslcs)

    if len(shared_nslcs)== 0:
        print('there are no shared seed_ids between the data and the inventory. Quitting...')
        return

    # Calculate PPSDs for each available NET.STA.LOC.CHAN
    median_psds, periods = {}, {}
    Path(PPSD_DIR).mkdir(exist_ok=True)
    print('CALCULATING PPSDS')
    for nslc in shared_nslcs:
        ppsd, datestr = _makePPSD(client, inv, nslc, args)
        if ppsd is None:
            print(f'{ppsd=}, skipping...')
            continue
        ppsd.plot(filename=f'{PPSD_DIR}/{nslc}_{datestr}_PPSD.png')
        if save_ppsds is True:
            ppsd.save_npz(f'{PPSD_DIR}/{nslc}_{datestr}_PPSD.npz')
        if plot_spectrograms is True:
            ppsd.plot_spectrogram(clim=[ppsd.db_bin_edges[0], ppsd.db_bin_edges[-1]],
                                  filename=f'{PPSD_DIR}/{nslc}_{datestr}_spectrogram.png')
        periods[nslc], median_psds[nslc] = ppsd.get_percentile(50)

    # For each unique channel name, plot the median PSDs together
    unique_channels = list(set([x.split('.')[3] for x in shared_nslcs]))
    unique_components = list(set([x.split('.')[3][-1] for x in shared_nslcs]))
    Path(MED_PSD_DIR).mkdir(exist_ok=True)
    print('PLOTTING median PSDS')
    for x in unique_channels:
        _plot_compare_channel(periods, median_psds, x, MED_PSD_DIR, verbose)
    for x in unique_components:
        _plot_compare_component(periods, median_psds, x, MED_PSD_DIR, verbose)


def _parse_input():
    parser = argparse.ArgumentParser(
        prog='plotPSDs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Calculates and compares Power Spectral Densities for an '
                    'SDS database and StationXML inventory')
    parser.add_argument('sds_dir', help='Path of SDS directory')
    parser.add_argument('inv_file', help='Path of StationXML file')
    parser.add_argument('-m', '--maxdays', default=1., type=float,
                        help='maximum # of days over which to calculate PPSDs')
    parser.add_argument('-k', '--skipdays', default=1., type=float,
                        help='# of days to skip before calculating PPSDs')
    parser.add_argument('-s', '--stations', default='*',
                        help='only stations matching the string '
                             '("*" and "?" wildcards accepted)')
    parser.add_argument('-c', '--channels', default='*',
                        help='only channels matching the string '
                             '("*" and "?" wildcards accepted)')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='verbose')
    return parser.parse_args()


def _get_shared_nslcs(client_nslcs, inv_nslcs, verbose=True):
    """
    Separate nslcs into client-only, inventory-only, and shared

    Return the "both"
    """
    shared = []
    client_only = []
    inv_only = []

    nslcs = list(set(inv_nslcs + client_nslcs))  # all unique nslcs
    for nslc in nslcs:
        if nslc in client_nslcs:
            if nslc in inv_nslcs:
                shared.append(nslc)
            else:
                client_only.append(nslc)
        elif nslc in inv_nslcs:
            inv_only.append(nslc)
        else:
            raise ValueError('{nslc=} not found in Client nor Inventory: '
                             'SHOULD BE IMPOSSIBLE!')
    assert (len(shared)+len(client_only)+len(inv_only)) == len(nslcs)
    if verbose:
        print(f'Inventory had {len(inv_nslcs)} station-channels')
        print(f'Client had {len(client_nslcs)} station-channels')
        print(f'{len(shared)} shared station-channels: {shared}')
        print(f'{len(client_only)} client-only station-channels: {client_only}')
        print(f'{len(inv_only)} inventory-only station-channels: {inv_only}')
    return list(set(shared))


def _inv_get_nslc(inv, args):
    """
    Return nslcs as a list of 'NET.STA.LOC.CHAN'

    Only return items matching the station and channel search strings
    """
    nslcs = []
    for net in inv:
        for sta in net:
            if re.match(_fdsn_to_regex(args.stations), sta.code) is None:
                continue
            for chan in sta:
                if re.match(_fdsn_to_regex(args.channels), chan.code) is None:
                    continue
                nslcs.append('.'.join((net.code, sta.code, chan.location_code,
                                       chan.code)))
    return nslcs


def _client_get_nslc(client, args):
    """
    Return client nslcs as a list of 'NET.STA.LOC.CHAN'

    Only return items matching the station and channel search strings
    """
    nslcs = []
    for x in client.get_all_nslc():
        if re.match(_fdsn_to_regex(args.stations), x[1]) is None:
            continue
        if re.match(_fdsn_to_regex(args.channels), x[3]) is None:
            continue
        nslcs.append('.'.join(x))
    return nslcs


def _fdsn_to_regex(x):
    """Convert FDSN search strings to regex"""
    x = x.replace('*', '.*')
    x = x.replace('?', '.')
    return '^' + x + '$'


def _inv_nslc_date_range(inv, nslc, verbose=False):
    n, s, l, c = nslc.split('.')
    nslc_inv = inv.select(network=n, station=s, channel=c, location=l)
    if len(nslc_inv) == 0:
        raise ValueError(f'{nslc=} not found in inventory!')
    for net in nslc_inv:
        for sta in net:
            for cha in sta:
                if cha.start_date is not None:
                    start_date = cha.start_date
                    if verbose:
                        print(f'{nslc} Channel start_date = {start_date}')
                elif sta.start_date is not None:
                    print(f'{nslc} Channel start_date not found, returning '
                          'Station start_date')
                    start_date = sta.start_date
                elif sta.start_date is not None:
                    print(f'{nslc} Channel and Station start_date not found, '
                          'returning Network start_date')
                    start_date = net.start_date
                if cha.end_date is not None:
                    end_date = cha.end_date
                    if verbose:
                        print(f'{nslc} Channel end_date = {end_date}')
                elif sta.end_date is not None:
                    print(f'{nslc} Channel end_date not found, returning '
                          'Station end_date')
                    end_date = sta.end_date
                elif sta.start_date is not None:
                    print(f'{nslc} Channel and Station end_date not found, '
                          'returning Network end_date')
                    end_date = net.end_date
    return start_date, end_date


def _makePPSD(client, inv, nslc, args):
    max_secs = args.maxdays*86400
    skip_secs = args.skipdays*86400

    ppsd = None
    n, s, l, c = nslc.split('.')
    start_date, end_date = _inv_nslc_date_range(inv, nslc)
    print(f'{nslc=}, dates={start_date.isoformat()} - {end_date.isoformat()}')
    inv_secs = end_date - start_date
    # Adjust start_date
    if inv_secs > skip_secs + 86400:
        # If there is at least one day of data after args.skipdays
        start_date += skip_secs
    elif inv_secs < skip_secs:
        # If there is less data than the requested skip time
        print(f"    Skip days ({arg.skipdays}) > inv days ({inv_secs/86400.:.1f})",
               end='')
        start_try = end_date - max_secs
        if start_try < start_date:
            print(" and max days (args.maxdays) > data days: skipping nothing")
        else:
            print(": only skipping {(start_try-start_date)/8600:.1f} days")
            start_date = start_try
    
    delta = 86400 if max_secs >= 86400 else int(max_secs)

    process_secs = end_date - start_date
    if max_secs is None:
        max_secs = math.ceil(process_secs)
    elif max_secs > inv_secs:
        print(f'    Reducing max_days ({args.maxdays}) to match inv_days ({process_secs/86400:.0f})')
        max_secs = math.ceil(process_secs)
    max_days = int(math.ceil(max_secs/86400))
    print(f'    Calculating {max_days}-day PPSD starting on {start_date.isoformat()} for {nslc=}')
    for day in range(max_days):
        stime = start_date + day*86400.
        etime = stime + delta
        st = client.get_waveforms(n, s, l, c, stime, etime, merge=1)
        range_text = f'nslc={n}.{s}.{l}.{c}, {stime.isoformat()}-{etime.isoformat()}'
        if len(st) == 0:
            print(f'    No data found for {range_text}, skipping...')
            continue
        elif len(st) > 1:
            print(f'    More than one stream ({len(st)}) for {range_text}, using first one')
        tr = st[0]
        if ppsd is None:
            if c[-1] in ('H', 'G', 'O'):
                # Pressure channel needs special handling
                ppsd = PPSD(tr.stats, metadata=inv,
                            special_handling='hydrophone',
                            db_bins=(-60, 60, 1.0))
            else:
                ppsd = PPSD(tr.stats, metadata=inv)
        ppsd.add(tr)
    return ppsd, f"{start_date.strftime('%Y%j')}-{end_date.strftime('%Y%j')}"


def _plot_compare_channel(periods, psds, channel, path, verbose=False):
    """
    Args:
        periods (dict of np.array): key=nslc, value = periods
        psds (np.array):  key=nslc, value = median psd value
        channel (str): channel name to use
        path (str or Path): where to write to
    """
    fig, ax = plt.subplots()
    key_list = list(psds.keys())
    used_keys = []
    for key in key_list:
        if key.split('.')[3] == channel:
            used_keys.append(key)
            ax.semilogx(periods[key], psds[key], label=key)
    if verbose:
        print(f'{channel=} keys={", ".join(used_keys)}')
    ax.legend(fontsize='xx-small', loc='best')
    ax.set_xlabel('Period(s)')
    T = periods[key_list[0]]  # Use first PPSD's periods for Peterson models
    if channel[-1] in ('H', 'G', 'O'):
        ln, hn = PressureNoiseModel(T)
        ax.semilogx(T, ln, '--')
        ax.semilogx(T, hn, '--')
        ax.set_ylabel('Power Spectral Density (dB ref 1 Pa/sqrt(Hz))')
    else:
        ln, hn = PetersonNoiseModel(T)
        ax.semilogx(T, ln, '--')
        ax.semilogx(T, hn, '--')
        ax.set_ylabel('Power Spectral Density (dB ref 1 m/s^2/sqrt(Hz))')
    plt.savefig(f'{path}/{channel}.PSDs.png')

def _plot_compare_component(periods, psds, component, path, verbose=False):
    """
    Args:
        periods (dict of np.array): key=nslc, value = periods
        psds (np.array):  key=nslc, value = median psd value
        component (str): component name to use
        path (str or Path): where to write to
    """
    fig, ax = plt.subplots()
    key_list = list(psds.keys())
    used_keys = []
    for key in key_list:
        if key.split('.')[3][-1] == component:
            used_keys.append(key)
            ax.semilogx(periods[key], psds[key], label=key)
    if verbose:
        print(f'{component=} keys={", ".join(used_keys)}')
    ax.legend(fontsize='xx-small', loc='best')
    ax.set_xlabel('Period(s)')
    T = periods[key_list[0]]  # Use first PPSD's periods for Peterson models
    if component in ('H', 'G', 'O'):
        ln, hn = PressureNoiseModel(T)
        ax.semilogx(T, ln, '--')
        ax.semilogx(T, hn, '--')
        ax.set_ylabel('Power Spectral Density (dB ref 1 Pa/sqrt(Hz))')
    else:
        ln, hn = PetersonNoiseModel(T)
        ax.semilogx(T, ln, '--')
        ax.semilogx(T, hn, '--')
        ax.set_ylabel('Power Spectral Density (dB ref 1 m/s^2/sqrt(Hz))')
    plt.savefig(f'{path}/{component}.PSDs.png')


if __name__ == '__main__':
    main()
