"""
Ryan Pranantyo
EOS, 11 February 2026

A script to prepare input fault files for JAGURS for rupture propagation mode.
I use finite fault model from USGS to test. 
As long as I have a similar format to the example, the script should be OK.

Idea, with example from one subfault (A):
    We have a total slip of 5 m. A starts to rupture (trup) at t=5s with duration (trise) of 10s.
    Interval rupture duration (tau, in JAGURS) is 2 s.
    tau should be within interval of dt.
    We have to prepare an increment slip values. Hence, we should have:
        slip at t0-2 s   = 0 m
        slip at t2-4 s   = 0 m
        slip at t4-6 s   = 0.833 m (at t=5 s)
        slip at t6-8 s   = 0.833 m
        slip at t8-10 s  = 0.833 m
        slip at t10-12 s = 0.833 m
        slip at t12-14 s = 0.833 m (finish at t=15 s)
        slip at t14-XX s = 0 m
    
    Slip should be distributed over num_slip_distributed = round(trise / tau + 1)
    Slip increment = total slip / num_slip_distributed
"""
import os
import sys
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from pathlib import Path
from pyproj import Geod

def slip_increment(slip, trise, tau):
    """estimate slip increment over tau"""
    num_slip_distributed = round(trise/tau + 1).astype(int)
    slip_inc = slip/num_slip_distributed
    return slip_inc

def distribute_slip(trup, trise, tau, slip_inc, rupdur):
    t0 = trup
    slip_dist = np.zeros(len(rupdur)-1)
    for ii, tt in enumerate(rupdur[:-1]):
        if trup >= rupdur[ii] and trup < rupdur[ii+1]:
            if trup < t0+trise:
                slip_dist[ii] = slip_inc
            trup += tau
    return slip_dist

def quick_check(df, rupdur_cols, tau):
    ncols = int(np.ceil(np.sqrt(len(rupdur_cols))))
    nrows = ncols

    fig = plt.figure(figsize=(12,12), constrained_layout=True)
    for ii, tt in enumerate(rupdur_cols):
        slip2plot = df[tt]
        slip2plot = slip2plot.where(slip2plot!=0,np.nan)
        ax = fig.add_subplot(nrows, ncols, ii+1)
        ax.set_title(tt, fontsize=6)
        SLIP = ax.scatter(df['lon'], df['lat'], c=slip2plot, vmin=0, vmax=df['slip'].max(),
                          s=10, alpha=0.8, cmap='Reds_r')
        ax.set_xlim(df['lon'].min()-0.1, df['lon'].max()+0.1)
        ax.set_ylim(df['lat'].min()-0.1, df['lat'].max()+0.1)
    #fig.colorbar(SLIP)
    fig.savefig(f'quick_check__tau{tau}.png', dpi=300)
    plt.close()

def transfer_coordinates(df, length, width, strike, dip):
    """lon, lat, z_depth are centroid of subfault"""
    g = Geod(ellps='WGS84')
    for ii in df.index:
        xt, yt, _ = g.fwd(df['lon'][ii], df['lat'][ii], strike-180, 0.5*length*1000.)
        xrt,yrt,_ = g.fwd(xt, yt, strike-90, np.cos(np.radians(dip))*0.5*width*1000.)
        df.loc[ii, 'lon_tr'] = xrt
        df.loc[ii, 'lat_tr'] = yrt
        df.loc[ii, 'zt_depth'] = df['z_depth'][ii] - np.sin(np.radians(dip)) * 0.5*width
        df.loc[ii, 'length'] = length
        df.loc[ii, 'width'] = width
        df.loc[ii, 'strike'] = strike
        df.loc[ii, 'dip'] = dip
    return df

def write2jagurs_fault_rupture(df, where2save, rupdur_cols):
    df_jagurs = df[['lat_tr', 'lon_tr', 'zt_depth', 'length', 'width', 'dip', 'strike', 'rake']]
    fault_list = open(os.path.join(where2save, 'fault.list'), 'w')
    for ii, tt in enumerate(rupdur_cols):
        df_jagurs['slip'] = df[tt]
        fout = os.path.join(where2save, f'fault_rupt-{ii:04d}.txt')
        df_jagurs = df_jagurs.round(decimals=5).astype('float32')
        df_jagurs.to_csv(fout, header=None, index=None, sep=' ')
        fault_list.write(f'{fout}\n')
        print(fout)
    fault_list.close()

def write2jagurs_fault_totalslip(df, where2save):
    df_jagurs = df[['lat_tr', 'lon_tr', 'zt_depth', 'length', 'width', 'dip', 'strike', 'rake', 'slip']]
    fout = os.path.join(where2save, f'fault_rupt-totalslip.txt')
    df_jagurs.to_csv(fout, header=None, index=None, sep=' ')
    print(df_jagurs)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str, default='./examples/usgs_ffm.csv')
    parser.add_argument('--jagurs_dt', type=float, default=1.5)
    parser.add_argument('--jagurs_tau', type=float, default=1.5)
    parser.add_argument('--sf_length', type=float, default=15.0)
    parser.add_argument('--sf_width', type=float, default=10.0)
    parser.add_argument('--sf_strike', type=float, default=325)
    parser.add_argument('--sf_dip', type=float, default=11.62)
    parser.add_argument('--outdir', type=str, default='fault_list')
    args = parser.parse_args()

    print(args)

    input_file = Path(args.input_file)
    fname = input_file.name[:-4]
    # load ffm parameters
    dt = args.jagurs_dt
    tau = args.jagurs_tau
    sf_length = args.sf_length
    sf_width = args.sf_width
    sf_strike = args.sf_strike
    sf_dip = args.sf_dip

    where2save = Path(os.path.join(args.outdir + f'__{fname}__' + f'tau__{tau}' ))
    where2save.mkdir(exist_ok = True)

    df = pd.read_csv(input_file)

    #check if tau/dt is integer, if not, change tau!
    check_tau = tau%dt
    if check_tau != 0.0:
        print('tau cannot be used because not within dt interval! Change tau!')
        sys.exit()
    else:
        print('tau can be used ...')
    
    #estimate slip increment over tau
    df['slip_increment'] = slip_increment(df['slip'], df['trise'], tau)

    #estimate maximum total rupture duration
    total_duration = np.ceil(df['trup'].max()+df['trise'].max()+tau)
    
    #rupture duration arrays
    rupdur = np.arange(0, total_duration, tau)
    rupdur_cols = [f't_{dt}s' for dt in rupdur[:-1]]
    df[rupdur_cols] = 0.0

    #distribute slip
    for sf in df.index:
        slip_dist = distribute_slip(df['trup'][sf], df['trise'][sf], tau, df['slip_increment'][sf], rupdur)
        df.loc[sf,rupdur_cols] = slip_dist

    quick_check(df, rupdur_cols, tau)

    #transfer coordinates for jagurs
    df = transfer_coordinates(df, sf_length, sf_width, sf_strike, sf_dip)

    #exporting
    write2jagurs_fault_rupture(df, where2save, rupdur_cols)
    write2jagurs_fault_totalslip(df, where2save)

    """
    trup = 49.4
    trise = 9.5
    slip_inc = 0.02455
    for ii, tt in enumerate(rupdur[:-1]):
        if trup >= rupdur[ii] and trup < rupdur[ii+1]:
            if trup < 49.4 + trise:
            slips_tmp[ii] = slip_inc
        trup += tau
        print(f'{rupdur[ii]} <= t < {rupdur[ii+1]}', slips_tmp[ii])
    """



