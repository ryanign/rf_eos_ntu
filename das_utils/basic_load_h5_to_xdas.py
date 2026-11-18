"""
Ryan Pranantyo
EOS, 18 November 2026

Load h5 DAS data, then convert to xdas format
"""
import os
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse
import xdas as xd
from xdas.coordinates import SampledCoordinate
from glob import glob

#------------------
# FUNCTION
#------------------
def read_generic_h5_data(args):
    """
    main goal is to convert generic h5 data into an XDAS data_array format for future processing
    """
    input_folder = args.input_folder
    event = args.name
    fs = args.sampling_frequency
    var2read = args.variable2read
    start_time = np.datetime64(args.start_time)
    outdir = args.outdir
    dx = args.dx
    os.makedirs(outdir, exist_ok = True)

    print(f'-'*60)
    print(f' {event}')
    print(f' Loading data from  = {input_folder}')
    print(f' Sampling frequency = {fs} Hz')
    print(f' Start time         = {start_time}')
    print(f'-'*60)

    files = glob(os.path.join(input_folder, '*.h5'))
    files.sort()
    data = []
    for fname in files:
        with h5py.File(fname, 'r') as f:
            d = f[var2read][:]
        data.append(d)
    data = np.concatenate(data, axis=0)

    nt, nx = data.shape

    #end_time = start_time + np.timedelta64(int(round(1e9 / fs)), 'ns') * (nt -1)

    dt = np.timedelta64(int(round(1e9 / fs)), 'ns')
    time = SampledCoordinate({
        'tie_values' : [start_time],
        'tie_lengths': [nt],
        'sampling_interval': dt,
        })

    da = xd.DataArray(
            data=data,
            coords={
                'time': time,
                #'time': {'tie_indices': [0, nt-1], 'tie_values': [start_time, end_time]},
                'distance': dx * np.arange(nx),
                },
            )

    return da

#------------------
# argument
#------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description = 'Read decimated generic h5 DAS data, convert to XDAS format, and do a quick plot',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=__doc__,
            )
    parser.add_argument(
            '--input_folder',
            default = '/home/ignatius.pranantyo/DAS/00_exercise/distributed_acoustic_sensing_data/2023p152354',
            help = 'path to decimated DAS data from data source'
            )
    parser.add_argument(
            '--name',
            default = '2023p152354',
            help = 'name of the event or data',
            )
    parser.add_argument(
            '--sampling_frequency',
            default = '50',
            type = float,
            help = 'sampling frequency of the data',
            )
    parser.add_argument(
            '--dx',
            default = 1.,
            help = 'channel distance',
            )
    parser.add_argument(
            '--variable2read',
            default = 'DAS',
            type = str,
            help = 'variable inside data to read',
            )
    parser.add_argument(
            '--start_time',
            default = '2023-02-26T08:53:00.122',
            help = 'data starting time in UTC',
            )
    parser.add_argument(
            '--outdir',
            default = '/home/ignatius.pranantyo/DAS/00_exercise/figures',
            help = 'output folder to save any processed file',
            )

    args = parser.parse_args()

    da = read_generic_h5_data(args)
    print(da)


    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # quick plotting
    print('-'*40)
    print(f'Quick plotting ...')
    da = da - da.median()
    vmax = float(np.percentile(np.abs(da.values), 99))

    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1,1,1)
    da.plot(ax=ax, yincrease=False, vmin=-vmax, vmax=vmax, cmap='seismic')
    fout = os.path.join(outdir, f'Data_xdas__{args.name}.png')
    fig.savefig(fout, dpi=300)
    plt.close()
    print(f'  check : {fout}')

    # simple signal processing
    print(f'-'*40)
    print(f'simple xfft and fk')
    import xdas.fft as xfft
    import xdas.signal as xs

    fk = xs.taper(da, dim='distance')
    fk = xs.taper(fk, dim='time')
    fk = xfft.rfft(fk, dim={'time': 'frequency'})
    fk = xfft.fft(fk, dim={'distance': 'wavenumber'})
    amp = np.abs(fk)
    fk_db = 20 * np.log10(amp / float(amp.values.max()) + 1e-12)
    
    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1,1,1)
    fk_db.plot(ax=ax, xlim=(-0.5, 0.5), vmin=-60, vmax=0, cmap='viridis', interpolation='none')
    fout = os.path.join(outdir, f'FK_xdas__{args.name}.png')
    fig.savefig(fout, dpi=300)
    plt.close()
    print(f'  check : {fout}')

    #stft ?
    print(f'-'*40)
    print(f'simple STFT')
    sub = da.isel(distance=slice(100, None, 10))
    s = xs.stft(sub, nperseg=256, noverlap=192, dim={'time': 'frequency'}, scaling='psd')

    p = np.abs(s.values) ** 2
    t0 = da['time'].values[0]
    tt = (s['time'].values - t0) / np.timedelta64(1, 's')
    ff = s['frequency'].values
    xx = s['distance'].values

    fig, axs = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    # spectrogram at one channel
    ich = int(np.argmin(np.abs(xx - 2000)))           # nearest to channel 2000
    db = 10 * np.log10(p[:, ich, :].T / p[:, ich, :].max() + 1e-12)
    im = axs[0].pcolormesh(tt, ff, db, shading='auto', cmap='viridis', vmin=-60, vmax=0)
    axs[0].set(xlabel='Time (s)', ylabel='Frequency (Hz)', title=f'Spectrogram, ch {xx[ich]:.0f}')
    fig.colorbar(im, ax=axs[0], label='dB')

    # mean spectrum per channel (whole record)
    m = p.mean(axis=0)                                # (distance, frequency)
    db2 = 10 * np.log10(m.T / m.max() + 1e-12)
    im2 = axs[1].pcolormesh(xx, ff, db2, shading='auto', cmap='viridis', vmin=-60, vmax=0)
    axs[1].set(xlabel='Channel', ylabel='Frequency (Hz)', title='Mean spectrum per channel')
    fig.colorbar(im2, ax=axs[1], label='dB')
    
    fig.savefig(os.path.join(outdir, f'STFT_xdas__{args.name}.png'), dpi=200)
    plt.close(fig)
    print(f'  check : {fout}')

    #bandpass
    print(f'-'*40)
    print(f'Bandpass and checking earthquake signal')
    #1. drop bad channels
    nch = da.shape[1]
    bad = set(range(100))
    for c in range(2930, 5620):
        bad.update(range(c - 3, c + 4))
    keep = np.array([i for i in range(nch) if i not in bad])
    eq = da.isel(distance=keep)

    #2. detrend, taper, bandpass
    fmin, fmax = 1.0, 5.0
    eq = xs.detrend(eq, dim='time')
    eq = xs.taper(eq, dim='time')
    eq = xs.filter(eq, [fmin, fmax], 'bandpass', corners=4, zerophase=True, dim='time')

    #3. normalise each channel
    eq = eq / eq.std('time')

    #4. zoom on the event window
    t0 = da['time'].values[0]
    win = eq.sel(time=slice(t0 + np.timedelta64(100, 's'), t0 + np.timedelta64(200, 's')))

    vmax = float(np.percentile(np.abs(win.values), 99))
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    win.plot(ax=ax, yincrease=False, vmin=-vmax, vmax=vmax, cmap='seismic')
    ax.set_title(f'Band-pass {fmin}–{fmax} Hz, normalised per channel')
    fig.savefig(os.path.join(outdir, f'EQ_{fmin}-{fmax}Hz__{args.name}.png'), dpi=200)
    plt.close(fig)
    print(f'  check : {fout}')




