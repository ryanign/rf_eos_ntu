"""
Ryan Pranantyo
EOS, 8 April 2026

Description:
    Basic script to extract tsunami signal from tide gauge data using 
    Butterworth Band pass filter

Outputs:

Usage:
    python detide_extract_tsunami_signal.py \
            --input_data \
            --data_frequency \
            --data_source \
            --output_dir \
            --freqs_to_cut \
            --eq_time \
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import butter, filtfilt

#----------------------------------------------
# Step 1: loading data
#----------------------------------------------
def loading_data(input_data, data_source, freq):
    """
    freq = sampling interval in seconds (e.g., 60)
    """
    print(f'\n  Loading data...')
    df = pd.read_csv(input_data)
    if data_source == 'IOC-UNESCO':
        print(f'  data from IOC-UNESCO')
        df.rename(columns={'Time_UTC':'time', 'batV':'bat', 'rad_m':'z_m'}, inplace=True)
    if data_source == 'BIG':
        print(f'  data from BIG')
        df.rename(columns={'data':'z_m',
                           'residu': 'residu_m'}, 
                  inplace=True)
        df['z_m'] = df['z_m'] / 100.
        df['residu_m'] = df['residu_m'] / 100.

    df['time'] = pd.to_datetime(df['time'])
    df = df.set_index('time')

    print(f'  Resampling or interpolating NaN data...')
    df = df.resample(f'{freq}s').mean()
    df['z_m'] = df['z_m'].interpolate(method='time')
    df['z_m'] = df['z_m'].ffill().bfill()

    return df

#----------------------------------------------
# Step 2: band pass filter
#----------------------------------------------
def butter_bandpass_filter(data, low_freq, high_freq, fs, order=6):
    print(f'\n  Filtering ...')
    nyq  = 0.5 * fs
    low  = low_freq / nyq
    high = high_freq / nyq

    #safety check
    if not (0 < low < high < 1):
        raise ValueError(f'\n  INVALID BANDPASS RANGE')

    b, a = butter(order, [low, high], btype='band')
    
    # ensure no NaN in the data
    if np.isnan(data).any():
        raise valueError(f'\n  INPUT DATA CONTAINS NaNs BEFORE FILTERING')
    
    # check length
    min_len = 3 * max(len(a),len(b))
    if len(data) <= min_len:
        raise ValueError(f'\n  DATA TOO SHORT FOR filtfilt')

    filt = filtfilt(b, a, data)

    return filt

#----------------------------------------------
# Step 3: saving
#----------------------------------------------
def saving_output(df, fname, output_dir, t0):
    print(f'\n  Saving output')
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    fout = os.path.join(output_dir, f'{fname[:-4]}__filtered.csv')
    df['time'] = df.index
    df.to_csv(fout, index=False)
    print(f'   fout : {fout}')

    # quick plotting
    t0 = pd.to_datetime(t0)
    print(f'   quick plotting ...')
    fig = plt.figure(figsize=(10,5), constrained_layout=True)
    ax = fig.add_subplot(2,1,1)
    ax.set_title(f'{fname[:-4]}', loc='left', pad=0.01)
    ax.plot(df.index, df['z_m'].values, linewidth=0.75, c='black')
    ax.vlines(t0, ymin=df['z_m'].min(), ymax=df['z_m'].max(), colors='red', linestyles='--')
    ax.set_xlim(t0-pd.Timedelta('1D'),df.index[-1])

    ax = fig.add_subplot(2,1,2)
    ax.set_title('Filtered', loc='left', pad=0.01)
    ax.plot(df.index, df['z_filtered_m'].values, linewidth=0.75, c='black')
    ax.set_xlim(t0-pd.Timedelta('1h'), t0+pd.Timedelta('12h'))
    ax.vlines(t0, ymin=df['z_filtered_m'].min(), ymax=df['z_filtered_m'].max(), colors='red', linestyles='--')

    fout = os.path.join(output_dir, f'{fname[:-4]}__data.png')
    fig.savefig(fout, dpi=300)
    plt.close()
    print(f'   fout : {fout}')

#----------------------------------------------
# MAIN
#----------------------------------------------
def main(args):
    input_data = Path(args.input_data)
    fname      = input_data.name

    # period to cut
    short_period = float(args.freqs_to_cut[0]) # in hour
    short_period = short_period * 3600.        # convert to seconds
    long_period  = float(args.freqs_to_cut[1]) # in hour
    long_period  = long_period  * 3600.        # convert to seconds

    # freq to cut
    low_freq  = 1./long_period          # to remove tide signal
    high_freq = 1./short_period         # to remove background ocean noise

    # data frequency in seconds
    fs = 1. / args.data_frequency

    print(f'='*60)
    print(f' Butterworth bandpass filter')
    print(f'='*60)
    print(f' Input data    : {args.input_data}')
    print(f' Low cut freq  : {low_freq:.5f} Hz ')
    print(f' High cut freq : {high_freq:.5f} Hz')
    print(f'='*60)

    df = loading_data(input_data, args.data_source, args.data_frequency)

    # start to filter
    df['z_filtered_m'] = butter_bandpass_filter(
            df['z_m'].values, 
            low_freq, 
            high_freq, 
            fs
            )

    # saving output
    saving_output(df, fname, args.output_dir, args.eq_time)

    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Butterworth band pass filter to extract tsunami signal from tide gauge data',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=__doc__,
            )
    parser.add_argument(
            '--input_data',
            default='/home/ryan/OneDrive_NTU_Projects/202604XX__HalmaheraDoubleSubduction_EQ-Tsunami/DATA_tsunami/IOC-UNESCO__PortKemaFishing_NS__20260328__20260404.csv',
            help='raw tide gauge data',
            )
    parser.add_argument(
            '--data_frequency',
            type=int,
            default=60,
            help='frequency of input data in seconds, will be used to interpolate missing data as well',
            )
    parser.add_argument(
            '--data_source',
            default='IOC-UNESCO',
            help='Currently only support IOC-UNESCO',
            )
    parser.add_argument(
            '--output_dir',
            default='/home/ryan/OneDrive_NTU_Projects/202604XX__HalmaheraDoubleSubduction_EQ-Tsunami/DATA_tsunami/filtered',
            help='Directory to save output file',
            )
    parser.add_argument(
            '--freqs_to_cut',
            nargs='+',
            default=[0.1, 3.],
            help='Low and high period frequencies to filter',
            )
    parser.add_argument(
            '--eq_time',
            default='2026-04-01T22:48:00',
            help='Earthquake time (UTC)',
            )
    args = parser.parse_args()

    df = main(args)


