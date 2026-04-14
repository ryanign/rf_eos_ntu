"""
Ryan Pranantyo
EOS, 13 April 2026

FFT and Wavelet analysis

Useage:

"""
import os
import sys
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cmcrameri.cm as cm
import argparse
import xarray as xr
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from pathlib import Path

#----------------------------------------------------
# BASIC FFT
#----------------------------------------------------
def compute_fft(df, data_freq, param):
    """
    compute FFT
    df        = data
    data_freq = sampling interval in seconds
    param     = variable to be analysed
    """
    print(f'\n  Compute FFT')
    
    print(f'   {param}')
    signal = df[param].values
    n      = len(signal)
    fs     = 1.0 / data_freq
        
    # FFT
    fft_vals = np.fft.rfft(signal)
    freqs    = np.fft.rfftfreq(n, d=data_freq)
    
    # Power spectrum
    power = (np.abs(fft_vals) ** 2) / n
        
    # convert frequency to period in minutes
    with np.errstate(divide='ignore'):
        periods_min = 1.0 / (freqs * 60.)

    # remove DC component (freq=0)
    mask = freqs > 0
    result = pd.DataFrame({
            'frequency_hz' : freqs[mask],
            'period_min'   : periods_min[mask],
            'power'        : power[mask],
            })

    return result

#----------------------------------------------------
# BASIC WAVELET ANALYSIS
#----------------------------------------------------
def morlet_wavelet(t, s, w0=6.0):
    """ Morlet wavelet function"""
    x = t / s
    return (np.pi**-0.25) * np.exp(1j * w0 * x) * np.exp(-0.5 * x**2)

def _morlet_cwt(signal, scales, dt, w0=6.0):
    """
    Manual CWT using Morlet wavelet via convulation in frequency domain.
    """
    n          = len(signal)
    freqs      = np.fft.fftfreq(n, d=dt)
    signal_fft = np.fft.fft(signal)
    coef = np.zeros((len(scales), n), dtype=complex)

    for i, s in enumerate(scales):
        # Morlet in frequency domain
        psi_fft = (np.sqrt(2 * np.pi * s /dt) *
                   (np.pi**-0.25) *
                   np.exp(-0.5 * (2 * np.pi * freqs * s - w0)**2))
        coef[i] = np.fft.ifft(signal_fft * np.conj(psi_fft))

    return coef

def compute_wavelet(df, data_freq, param, period_min=2., period_max=180., n_periods=100, morlet_w=6.0):
    print(f'\n  Compute wavelet')
    signal = df[param].values
    dt     = float(data_freq)

    periods_min = np.logspace(
            np.log10(period_min),
            np.log10(period_max),
            n_periods
            )
    periods_sec = periods_min * 60.0
    scales      = (morlet_w * periods_sec) / (2 * np.pi)

    coef  = _morlet_cwt(signal, scales, dt, w0=morlet_w)
    power = np.abs(coef) ** 2 / scales[:, np.newaxis]
    power_log2 = np.log2(power + 1e-20)

    return power_log2, periods_min, df['time'].values

#----------------------------------------------------
# SAVING
#----------------------------------------------------
def saving_fft(df, name, output_dir):
    """
    saving FFT analysis
    name = station name
    """
    print(f'\n  Saving FFT')
    fout = os.path.join(output_dir, f'{name}__FFT_analysis.csv')
    df.to_csv(fout, index=False)
    print(f'   output file : {fout}')

def saving_wavelet(power, periods, times, output_dir, name):
    """
    saving wavelet analysis result
    name = station name
    """
    print(f'\n  Saving wavelet')

    ds = xr.Dataset(
            {
                'power' : (['period_min', 'time'], power,
                           {'long_name' : 'Wavelet power spectrum',
                            'units'     : 'log2(m^2)',
                            'wavelet'   : 'Morlet'}),
            },
            coords = {
                'time' : times,
                'period_min' : periods,
                },
            attrs = {
                'title':f'Wavelet analysis - {name}',
                'description':'CWT using Morlet wavelet',
                'misc':'script by ignatius.pranantyo@ntu.edu.sg',
                }
            )

    fout = os.path.join(output_dir, f'{name}__WAVELET_analysis.nc')
    ds.to_netcdf(fout)
    print(f'   output file : {fout}')
    return ds

#----------------------------------------------------
# NICE PLOTTING
#----------------------------------------------------
def nice_plot(df, fft, ds, name, output_dir, eq_time, start, end):
    print(f'\n  Nice plotting')
    # colourmap for wavelet
    cmap = cm.batlowK
    clrs = np.linspace(ds['power'].min().values, ds['power'].max().values, int(cmap.N/4))
    norm = colors.BoundaryNorm(clrs, cmap.N)

    what_to_plot = ['tide', 'tsunami', 'WAVELET', 'FFT']
    fig = plt.figure(figsize=(7,8))
    fig.suptitle(name, fontsize=10)
    gs = GridSpec(nrows=5, ncols=1, 
                  height_ratios=[0.7, 0.7, 1.3, 0.3, 1.3],
                  width_ratios=[1.],
                  left=0.1, right=0.9,
                  bottom=0.08, top=0.95,
                  hspace=0.09,
                  wspace=0.05,)
    
    t0 = eq_time+pd.Timedelta(f'{start}h')
    t1 = eq_time+pd.Timedelta(f'{end}h')
    df = df[(df['time'] >= t0) & (df['time'] <= t1)]

    time_ticks = pd.date_range(
            start = t0.ceil('1h'),
            end = t1,
            freq = '3h')

    ### plot tide and tsunami
    for ii, var in enumerate(['z_m', 'z_filtered_m']):
        ax = fig.add_subplot(gs[ii,0])
        ax.plot(df['time'], df[var], c='black', linewidth=0.8,)
        ax.vlines(eq_time, ymin=df[var].min()*2, ymax=df[var].max()*2, colors='red', linestyles='--')
        ax.set_xlim(t0, t1)
        ax.set_ylim(df[var].min()*1.1,  df[var].max()*1.1)
        ax.tick_params(labelsize=8, labelbottom=False)
        ax.set_xticks(time_ticks)
        if ii == 1:
            ax.set_ylabel('Tsunami (m)', fontsize=10)
            #ax.tick_params(labelbottom=True)
            #ax.set_xlabel('Time (UTC)', fontsize=10)
        else:
            ax.set_ylabel('Data (m)', fontsize=10)

    ### plot wavelet
    ax = fig.add_subplot(gs[2,0])
    WV = ax.imshow(ds['power'].values,
                   origin='lower', aspect='auto', interpolation=None,
                   cmap=cmap, norm=norm,
                   extent=[0, len(ds['time'].values)-1, 0, len(ds['period_min'].values)-1],
                   )

    eq_idx = np.argmin(np.abs(pd.to_datetime(ds['time'].values) - eq_time))
    ax.axvline(x=eq_idx, color='red', linestyle='--')
    
    xtick_indices = [np.argmin(np.abs(pd.to_datetime(ds['time'].values) - t)) for t in time_ticks]
    ax.set_xticks(xtick_indices)
    ax.set_xticklabels(
            [pd.Timestamp(time).strftime('%H:%M\n%Y-%m-%d') for time in time_ticks], 
            fontsize=8)

    n_yticks  = 8
    y_indices = np.linspace(0, len(ds['period_min'].values)-1, n_yticks, dtype=int)
    ax.set_yticks(y_indices)
    ax.set_yticklabels(
            [f'{ds["period_min"][i]:.2f}' for i in y_indices],
            fontsize=8)
    ax.set_xlabel('Time (UTC)', fontsize=10)
    ax.set_ylabel('Period (min)', fontsize=10)
    
    cax = ax.inset_axes([0.97, 0., 0.03, 1.])
    cb  = fig.colorbar(WV, cax=cax,
                       orientation='vertical',
                       ticks=np.linspace(ds['power'].values.min(), ds['power'].values.max(), 5, dtype=int),
                       )
    cb.set_label('Power ($log_2(m^2)$)')

    ### plot FFT
    ax = fig.add_subplot(gs[4,0])
    ax.plot(fft['frequency_hz'], fft['power'], c='black', linewidth=1.)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1./(3*60*60), 1/(5*60))
    ax.tick_params(labelsize=8)
    ax.set_xlabel('Frequency (Hz)', fontsize=10)
    ax.set_ylabel('Power', fontsize=10)
    # saving
    fout = os.path.join(output_dir, f'{name}__analysis.png')
    fig.savefig(fout, dpi=300)
    plt.close()



#----------------------------------------------------
# MAIN
#----------------------------------------------------
def main(args):
    input_data = Path(args.input_data)
    eq_time = pd.to_datetime(args.eq_time)
    start   = args.time_window[0]
    end     = args.time_window[1]

    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    station_name = args.station_name

    print(f'='*60)
    print(f' SPECTRAL ANALYSIS')
    print(f'='*60)
    print(f' station name    = {station_name}')
    print(f' input data      = {input_data}')
    print(f' earthquake time = {eq_time}')
    print(f' time window     = {start} h -- {end} h')
    print(f'='*60)

    df = pd.read_csv(input_data)
    df['time'] = pd.to_datetime(df['time'])
    
    # bring it to MSL
    df['z_m']  = df['z_m'] - df['z_m'].mean()

    t0 = eq_time + pd.Timedelta(f'{start}h')
    t1 = eq_time + pd.Timedelta(f'{end}h')

    print(f'\n  Cutting data: {t0} - {t1}')

    df = df[(df['time'] >= t0) & (df['time'] <= t1)]

    # start to do FFT analysis
    df_fft = compute_fft(df, args.data_sampling, args.column_name)
    
    # saving fft
    saving_fft(df_fft, station_name, output_dir)

    # wavelet analysis
    power, periods, times = compute_wavelet(
            df, 
            param = args.column_name,
            data_freq  = args.data_sampling,
            period_min = args.period_min,
            period_max = args.period_max,
            n_periods  = args.n_periods,
            morlet_w   = args.morlet_w,
            )
    
    # saving wavelet
    ds = saving_wavelet(power, periods, times, output_dir, station_name)
    
    # plotting
    nice_plot(df, df_fft, ds, station_name, output_dir, eq_time, start, end)

    return df, df_fft, ds

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description = 'Spectral analysis',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=__doc__,
            )
    parser.add_argument(
            '--input_data',
            default='/home/ryan/OneDrive_NTU_Projects/202604XX__HalmaheraDoubleSubduction_EQ-Tsunami/DATA_tsunami/COMPILED/BIG__BITG_Bitung__20260324__20260407__filtered.csv',
            help='Resulted from detide_extract_tsunami_signal.py. \
                    Other input file is accepted as long as follow the standard',
            )
    parser.add_argument(
            '--station_name',
            default='BIG__BITG_Bitung',
            help='Station observation name',
            )
    parser.add_argument(
            '--data_sampling',
            default=60,
            help='Sampling interval in seconds',
            )
    parser.add_argument(
            '--column_name',
            default='z_filtered_m',
            help='Column name inside input data to be analysed. "z_filtered_m" or "z_m"'
            )
    parser.add_argument(
            '--output_dir',
            default='/home/ryan/OneDrive_NTU_Projects/202604XX__HalmaheraDoubleSubduction_EQ-Tsunami/Tsunami_SpectralAnalysis/',
            help='Directory to save output file',
            )
    parser.add_argument(
            '--eq_time',
            default='2026-04-01T22:48:00',
            help='Earthquake time (UTC)',
            )
    parser.add_argument(
            '--time_window',
            nargs='+',
            default=[-6, 12],
            help='Time window for analysis, hours before and after the earthquake time',
            )
    parser.add_argument(
            '--period_min',
            type=float,
            default=2.0,
            help='Minimum period for wavelet in minutes (default: 2)',
            )
    parser.add_argument(
            '--period_max',
            type=float,
            default=180.,
            help='Maximum period for wavelet in minutes (default: 180)',
            )
    parser.add_argument(
            '--n_periods',
            type=int,
            default=50,
            help='Numer of period scales log-spaced (default: 50)',
            )
    parser.add_argument(
            '--morlet_w',
            type=float,
            default=6.0,
            help='Morlet omega0 parameter (default: 6.0)',
            )
    args = parser.parse_args()

    df, df_fft, ds = main(args)
