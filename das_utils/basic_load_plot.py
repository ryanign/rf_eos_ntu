"""
Ryan Pranantyo
EOS, 17 November 2026

Mission:
    - load and quick plot decimated DAS data
    - data source: https://datacommons.anu.edu.au/DataCommons/rest/records/anudc:6317/data/
"""
import os
import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import argparse
from glob import glob

def plot_data(data, times, outdir, event, plot_type='original'):
    print(f'\n Plotting {plot_type} data ...')
    
    data = data - np.median(data, axis=0)
    vmax = np.percentile(np.abs(data), 99)

    times = pd.DatetimeIndex(times)
    dt = times[1] - times[0]
    t0 = mdates.date2num(times[0])
    t1 = mdates.date2num(times[-1] + dt)

    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot(1,1,1)
    ax.set_title(event)

    im = ax.imshow(
            data, aspect='auto',
            cmap='seismic',
            vmin=-vmax, vmax=vmax,
            extent=[0, data.shape[1], t1, t0],
            interpolation=None
            )

    ax.yaxis_date()
    ax.yaxis.set_major_locator(mdates.MinuteLocator())
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.yaxis.set_minor_locator(
            mdates.SecondLocator(bysecond=range(0, 60, 15)))
    ax.yaxis.set_minor_formatter(mdates.DateFormatter(':%S'))
    ax.tick_params(axis='y', which='minor', labelsize=6)

    ax.set_xlabel('Channel')
    ax.set_ylabel(f'Time (UTC), {times[0]:%Y-%m-%d}')

    fig.colorbar(im, orientation='horizontal', shrink=0.7)

    fout = os.path.join(outdir, f'Data_{plot_type}__{event}.png')
    fig.savefig(fout, dpi=300)
    plt.close()
    print(f'  check = {fout}')



def main(args):
    input_folder = args.input_folder
    event = args.name
    fs = args.sampling_frequency
    var2read = args.variable2read
    start_time = pd.to_datetime(args.start_time)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok = True)

    print(f'-'*60)
    print(f' {event}')
    print(f' Loading data from  = {input_folder}')
    print(f' Sampling frequency = {fs} Hz')
    print(f' Start time         = {start_time}')
    print(f'-'*60)

    # load data
    files = glob(os.path.join(input_folder, '*.h5'))
    files.sort()
    data, lengths = [], []
    for fname in files:
        with h5py.File(fname, 'r') as f:
            d = f[var2read][:]
        data.append(d)
        lengths.append(d.shape[0])
    data = np.concatenate(data, axis=0)

    # check num of channel in data if similar as inputed
    if data.shape[1] != args.num_of_channels:
        num_of_channels = data.shape[1]
    else:
        num_of_channels = args.num_of_channels
    print(f'\n Number of channels = {num_of_channels}')

    # check length of data
    length = data.shape[0]
    duration = length / fs
    print(f' Data duration = {duration} seconds')

    # create time_array
    times = pd.date_range(
            start = start_time,
            periods = length,
            freq = pd.Timedelta(seconds=1/fs)
            )

    # plotting original data
    plot_data(data, times, outdir, event, 'original')

    # a simple bandpass filter




    return files, data, event, outdir, times

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description = 'load decimated h5 DAS data and do a quick plot',
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
            '--num_of_channels',
            default = 7488,
            type = int,
            help = 'number of channels if known',
            )
    parser.add_argument(
            '--outdir',
            default = '/home/ignatius.pranantyo/DAS/00_exercise/figures',
            help = 'output folder to save any processed file',
            )
    args = parser.parse_args()

    files, data, event, outdir, times = main(args)

    ### sandpit, a simple bandpass filter using obspy
    print(f'  trying to do a simple bandpass filter')
    from obspy.signal.filter import bandpass, highpass, lowpass
    freqMin = 1
    freqMax = 20
    skipped_channel_num = 1
    first_channel = 0
    last_channel = data.shape[1]
    cropped_channel_list = np.arange(first_channel, last_channel, skipped_channel_num)
    data_filtered = np.empty((data.shape[0], cropped_channel_list.shape[0]), dtype='float32')
    for index, ch in enumerate(cropped_channel_list):
        data_filtered[:,index] = bandpass(data[:,ch], freqMin, freqMax, args.sampling_frequency, corners=4, zerophase=True)

    plot_data(data_filtered, times, outdir, 'filtered')

