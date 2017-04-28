# Script to compute deconvolution, cross-correlation, and coherency.
# L. Viens 04/28/2017 

import numpy as np
from obspy import read, Trace
from obspy.signal.filter import bandpass
import matplotlib.pyplot as plt
from Functions_GJI_2017 import *


def main():

    save_file =  "./foo.png" # Output file (Plot figure if: save_file = None)
    lag_dur = 100       # Duration of the signal: -100 to 100 s
    cut_per = .1        # Pass band low corner frequency of the filter in Hz
    cut_per2 = .5       # Pass band high corner frequency of the filter in Hz
    len_time_win = 48   # Divide the 1 day record into 192 segments (7.5 min segments)
    std_cont = 10       # 10 times standard deviation criteria. 
    stabilization = 10  # To stabilize the deconvolution

    directory ="./"
    
    # Methods 
    methods = ["1-bit" ,  "Cross-correlation" , "Coherency","Deconvolution"]
    # Station names
    sta_source = "station_1"
    sta_rec = "station_2"

    # Path and file names
    name_virt = directory + sta_source + "_1d.sac"
    name_rec = directory + sta_rec + "_1d.sac"                                

    # Load data
    dat_s1 = read(name_virt, debug_headers=True)
    dat_r1 = read(name_rec, debug_headers=True)
    # Get dt (in s) and delta (in Hz)
    dt = dat_s1[0].stats.delta
    delta =int(1/dt)
    # Split the 1 day data into "len_time_win" segments
    dat_s = np.split(dat_s1[0].data, len_time_win)
    dat_r = np.split(dat_r1[0].data, len_time_win)

    cnt = 1
    q=0

    for ii in range(np.size(methods)):
        cnt = 1
        q=0

        for i in range(len(dat_s)):
            # Compute the waveforms 
            if methods[ii]=='1-bit':
                corr, nodata , cnt = one_bit(dat_s[i], dat_r[i], delta, tpm=lag_dur, std_control= std_cont, count=cnt)
            elif  methods[ii]=='Deconvolution': 
                corr, nodata , cnt = deconvolution_stab(dat_s[i], dat_r[i], delta, tpm =lag_dur, std_control = std_cont, count=cnt, stab= 10)
            elif  methods[ii]=='Cross-correlation':
                corr, nodata , cnt = cross_corr(dat_s[i], dat_r[i], delta, tpm=lag_dur, std_control= std_cont, count=cnt)
            elif  methods[ii]=='Coherency':
                corr, nodata , cnt = coherency_stab(dat_s[i], dat_r[i], delta, tpm=lag_dur, std_control = std_cont, count=cnt, stab= 10)

            if q == 0:
                sumfin = corr
                q = 1
            elif q==1 and nodata ==0:
                sumfin = sumfin + corr # stack in time domain
            elif q==1 and nodata ==1:
                nodata =0
        
        # Filter data
        sumfin = bandpass(data=sumfin, freqmin = cut_per, freqmax = cut_per2, df=delta, corners=4, zerophase=True)
        # create time axis
        t = np.arange(-lag_dur, lag_dur, dt) 

        # Plot
        plt.subplot(np.size(methods), 1, ii+1)
        plt.plot(t,sumfin)
        plt.text(50, max(sumfin)/1.2,methods[ii])
        plt.ylabel('Amplitude')
        frame1 = plt.gca()
        if ii+1==1:
            plt.title('Band-pass filter: ' + str(cut_per) + '-' + str(cut_per2) + ' Hz')
        elif ii+1==np.size(methods):
            plt.xlabel('Time (s)')

        if ii+1<np.size(methods):
            plt.tick_params(axis='x', labelbottom='off')

 
    # Save data
    if  save_file:
        plt.savefig(save_file, dpi=400)
        plt.close()
        print("Figure saved as " + save_file)
    else:
        plt.show()


if __name__=="__main__":
   main()


