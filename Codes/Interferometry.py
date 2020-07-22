# Script to compute deconvolution, cross-correlation, and coherency.
# L. Viens 04/28/2017
# Code cleaned by L. Viens  22/07/2020

import numpy as np
from obspy import read
from obspy.signal.filter import bandpass
import matplotlib.pyplot as plt
from Functions_GJI_2017 import one_bit, cross_corr, deconvolution_stab, coherency_stab

save_file =  "../Figures/foo.png" # Output file 
lag_dur = 100       # Duration of the signal to plot: between -100 to 100 s
cut_per = .1        # Filter lower corner frequency (in Hz)
cut_per2 = .5       # Filter higher corner frequency (in Hz)
len_time_win = 48   # Divide the 1 day record into 48 30-min segments
std_cont = 10       # Reject windows with spikes larger than 10 times the standard deviation of the window 
stabilization = 10  # To stabilize the denominator of the deconvolution and coherency methods


# Methods 
methods = ["1-bit cross-corr" ,  "Cross-correlation" , "Coherency","Deconvolution"]

# Load virtual source and receiver data
name_virt = '../Data/station_1_1d.sac' # Virtual source data (1 day)
name_rec = '../Data/station_2_1d.sac'  # Receiver station data (1 day)                           
dat_s1 = read(name_virt, debug_headers = True)
dat_r1 = read(name_rec, debug_headers = True)

# Get sampling rate (in s) and sampling frequency (in Hz)
dt = dat_s1[0].stats.delta
delta = int(1/dt)

# create time axis
t = np.arange(-lag_dur, lag_dur, dt)

# Split the 1 day data into "len_time_win" segments
dat_s = np.split(dat_s1[0].data, len_time_win)
dat_r = np.split(dat_r1[0].data, len_time_win)

fig = plt.figure(figsize = (8, 10) )
for ii in range(np.size(methods)): # Loop over the 4 methods
    cnt = 1
    sumfin = np.zeros( len(t))
    for i in range(len(dat_s)): # Loop over the data
        if methods[ii]=='1-bit cross-corr':
            corr, nodata, cnt = one_bit(dat_s[i], dat_r[i], delta, tpm=lag_dur, std_control= std_cont, count=cnt)
        elif  methods[ii]=='Cross-correlation':
            corr, nodata, cnt = cross_corr(dat_s[i], dat_r[i], delta, tpm=lag_dur, std_control= std_cont, count = cnt)
        elif  methods[ii]=='Coherency':
            corr, nodata, cnt = coherency_stab(dat_s[i], dat_r[i], delta, tpm=lag_dur, std_control = std_cont, count=cnt, stab =10)
        elif  methods[ii]=='Deconvolution': 
            corr, nodata, cnt = deconvolution_stab(dat_s[i], dat_r[i], delta, tpm =lag_dur, std_control = std_cont, count=cnt, stab=10)
        
        sumfin += corr # linear stack in time domain
    
    # Band-pass filter the data with a 4-pole 2-pass Butterworth filter
    sumfin = bandpass(data = sumfin, freqmin = cut_per, freqmax = cut_per2, df = delta, corners = 4, zerophase = True)
    
    # Plot the data
    plt.subplot(np.size(methods), 1, ii+1)
    plt.plot(t, sumfin)
    plt.text(50, max(sumfin)/1.2, methods[ii])
    plt.ylabel('Amplitude')
    plt.grid()
    plt.xlim(t[0], t[-1])
    if ii == 0:
        plt.title('Band-pass filter: ' + str(cut_per) + '-' + str(cut_per2) + ' Hz')
    elif ii+1 == np.size(methods):
        plt.xlabel('Time (s)')

# Save plot
plt.tight_layout()
plt.plot()
fig.savefig(save_file, dpi=100)
print("Figure saved as " + save_file)


