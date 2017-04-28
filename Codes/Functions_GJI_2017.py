import numpy as np
import obspy.signal as obssig

def one_bit(dat_s, dat_r, delta, tpm, std_control, count):
    """Compute cross-correlation of 1-bit data.
        Parameters
        ----------
        Inputs
        dat_s: numpy.ndarray
            Record at the virtual source (1-h long in Viens et al. (2017, GJI))
        dat_r: numpy.ndarray
            Record at the virtual source (1-h long in Viens et al. (2017, GJI))
        delta: float
            Sampling rate in Hz (4 Hz in Viens et al. (2017, GJI))
        tpm: int
            To return the cross-correlation between -"tpm" seconds and +"tmp" seconds (tpm = 1500 s in Viens et al. (2017, GJI)).
        count: int
            When looping over many files, count the number of cross-correlations (should be equal to 1 for the first cross-correlation).
        std_control:
            If the virtual source or receiver records have spikes larger than std_control times the standard deviation of the signal, the cross-correlation is not computed (std_control = 10 in Viens et al. (2017, GJI)).
        
        -------
        Outputs
        corr: numpy.ndarray
           Cross-correlated waveform of 1-bit data in time domain between -"tpm" and +"tpm" seconds.
        nodata: int
            If the virtual source or receiver records have spikes larger than "std_control" times the standard deviation of the signal, "nodata" is equal to 1
        count: int
            Incremented "count" variable.
            """

    n = len(dat_s) 
    std_s = np.std(dat_s)
    std_r = np.std(dat_r)
    mx_s  = max(np.absolute(dat_s))
    mx_r  = max(np.absolute(dat_r))

    if mx_s< std_s*std_control and mx_r< std_r*std_control:
        dat_r[dat_r[:] >0 ]=1
        dat_r[dat_r[:] <0 ]=-1
        dat_s[dat_s[:] >0 ]=1
        dat_s[dat_s[:] <0 ]=-1

        fft_s = np.fft.fft(dat_s, n*5)
        fft_r = np.fft.fft(dat_r, n*5)
        cc_t1 = np.real(np.fft.ifft( (fft_r * np.conj(fft_s))))        
        corr2 = np.concatenate((cc_t1[int(len(cc_t1)/2):], cc_t1[:int(len(cc_t1)/2)]))
        corr  = corr2[int(len(corr2)/2)-tpm*delta:int(len(corr2)/2)+tpm*delta]
        nodata = 0
        count +=1
    else:
        corr = 0
        nodata = 1
    return corr, nodata, count



def cross_corr(dat_s, dat_r, delta, tpm, std_control, count):
    """Compute cross-correlation of raw data
        Parameters
        ----------
        Inputs
        dat_s: numpy.ndarray
            Record at the virtual source (1-h long in Viens et al. (2017, GJI))
        dat_r: numpy.ndarray
            Record at the virtual source (1-h long in Viens et al. (2017, GJI))
        delta: float
            Sampling rate in Hz (4 Hz in Viens et al. (2017, GJI))
        tpm: int
            To return the cross-correlation between -"tpm" seconds and +"tmp" seconds (tpm = 1500 s in Viens et al. (2017, GJI)).
        count: int
            When looping over many files, count the number of cross-correlations (should be equal to 1 for the first cross-correlation).
        std_control: int
            If the virtual source or receiver records have spikes larger than std_control times the standard deviation of the signal, the cross-correlation is not computed (std_control = 10 in Viens et al. (2017, GJI)).
        
        -------
        Outputs
        corr: numpy.ndarray
           Cross-correlated waveform of raw data in time domain between -tpm and +tpm seconds.
        nodata: int
            If the virtual source or receiver records have spikes larger than "std_control" times the standard deviation of the signal, "nodata" is equal to 1
        count: int
            Incremented "count" variable.
            """
    n = len(dat_s) 
    std_s = np.std(dat_s)
    std_r = np.std(dat_r)
    mx_s  = max(np.absolute(dat_s))
    mx_r  = max(np.absolute(dat_r))

    if mx_s< std_s*std_control and mx_r< std_r*std_control: 

        fft_s = np.fft.fft(dat_s, n*5)
        fft_r = np.fft.fft(dat_r, n*5)
        cc_t1 = np.real(np.fft.ifft( (fft_r * np.conj(fft_s))))        
        corr2 = np.concatenate((cc_t1[int(len(cc_t1)/2):], cc_t1[:int(len(cc_t1)/2)]))
        corr  = corr2[int(len(corr2)/2)-tpm*delta:int(len(corr2)/2)+tpm*delta]
        nodata = 0
        count +=1
    else:
        corr = 0
        nodata = 1
    return corr, nodata,  count


def deconvolution_stab(dat_s, dat_r, delta, tpm, std_control, count, stab):
    """Compute Deconvolution of raw data.
        Parameters
        ----------
        Inputs
        dat_s: numpy.ndarray
            Record at the virtual source (1-h long in Viens et al. (2017, GJI))
        dat_r: numpy.ndarray
            Record at the virtual source (1-h long in Viens et al. (2017, GJI))
        delta: float
            Sampling rate in Hz (4 Hz in Viens et al. (2017, GJI))
        tpm: int
            To return the cross-correlation between -"tpm" seconds and +"tmp" seconds (tpm = 1500 s in Viens et al. (2017, GJI)).
        std_control: int
            If the virtual source or receiver records have spikes larger than "std_control times" the standard deviation of the signal, the cross-correlation is not computed (std_control = 10 in Viens et al. (2017, GJI)).
        count: int
            When looping over many files, count the number of cross-correlations (should be equal to 1 for the first cross-correlation).
        stab: int
            To smooth the denominator term over "stab" points (stab = 10 in Viens et al. (2017, GJI)).

        -------
        Outputs
        corr: numpy.ndarray
           Cross-correlated waveform in time domain between -"tpm" and +"tpm" seconds.
        nodata: int
            If the virtual source or receiver records have spikes larger than "std_control" times the standard deviation of the signal, "nodata" is equal to 1
        count: int
            Incremented "count" variable.
            """
    n = len(dat_s)   
    std_s = np.std(dat_s)
    std_r = np.std(dat_r)
    mx_s  = max(np.absolute(dat_s))
    mx_r  = max(np.absolute(dat_r))
    if mx_s< std_s*std_control and mx_r< std_r*std_control:
        fft_s = np.fft.fft(dat_s, n*5)
        fft_r = np.fft.fft(dat_r, n*5)

        sj = obssig.util.smooth(np.absolute(fft_s), stab)

        dec_t1 =  np.real(np.fft.ifft( (fft_r * np.conj(fft_s))/ (sj**2) ))
        dec_t2 = np.concatenate((dec_t1[int(len(dec_t1)/2):], dec_t1[:int(len(dec_t1)/2)]))
        dec_t  = dec_t2[int(len(dec_t2)/2)-tpm*delta:int(len(dec_t2)/2)+tpm*delta]
        nodata = 0
        count +=1
    else:
        dec_t = 0
        nodata = 1
    return dec_t, nodata, count


def coherency_stab(dat_s, dat_r, delta, tpm, std_control, count, stab):
    """Compute coherency of raw data.
        Parameters
        ----------
        Inputs
        dat_s: numpy.ndarray
            Record at the virtual source (1-h long in Viens et al. (2017, GJI))
        dat_r: numpy.ndarray
            Record at the virtual source (1-h long in Viens et al. (2017, GJI))
        delta: float
            Sampling rate in Hz (4 Hz in Viens et al. (2017, GJI))
        tpm: int
            To return the cross-correlation between -"tpm" seconds and +"tmp" seconds (tpm = 1500 s in Viens et al. (2017, GJI)).
        std_control: int
            If the virtual source or receiver records have spikes larger than std_control times the standard deviation of the signal, the cross-correlation is not computed (std_control = 10 in Viens et al. (2017, GJI)).
        count: int
            When looping over many files, count the number of cross-correlations (should be equal to 1 for the first cross-correlation).
        stab: int
            To smooth the denominator term over "stab" points (stab = 10 in Viens et al. (2017, GJI)).
       
        -------
        Outputs
        corr: numpy.ndarray
            Coherency waveform in time domain between -"tpm" and +"tpm" seconds.
        nodata: int
            If the virtual source or receiver records have spikes larger than "std_control" times the standard deviation of the signal, "nodata" is equal to 1
        count: int
            Incremented "count" variable.
            """

    n = len(dat_s)   
    std_s = np.std(dat_s)
    std_r = np.std(dat_r)
    mx_s  = max(np.absolute(dat_s))
    mx_r  = max(np.absolute(dat_r))
    if mx_s< std_s*std_control and mx_r< std_r*std_control:
        fft_s = np.fft.fft(dat_s, n*5)
        fft_r = np.fft.fft(dat_r, n*5)

        sj = obssig.util.smooth(np.absolute(fft_s), stab) 
        si = obssig.util.smooth(np.absolute(fft_r), stab)

        coh_t1 =  np.real(np.fft.ifft( (fft_r * np.conj(fft_s))/ (si*sj) ))
        coh_t2 = np.concatenate((coh_t1[int(len(coh_t1)/2):], coh_t1[:int(len(coh_t1)/2)]))
        coh_t  = coh_t2[int(len(coh_t2)/2)-tpm*delta:int(len(coh_t2)/2)+tpm*delta]
        nodata = 0
        count +=1
    else:
        coh_t = 0
        nodata = 1
    return coh_t, nodata, count

