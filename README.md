# Viens et al. (2017, GJI)
This repository contains the functions and an example for the paper:
- Loïc Viens, Marine Denolle, Hiroe Miyake, Shin’ichi Sakai, Shigeki Nakagawa; Retrieving impulse response function amplitudes from the ambient seismic field. Geophys J Int 2017 ggx155. doi: 10.1093/gji/ggx155

Files:
- Functions_GJI_2017.py -> Functions to compute cross-correlation, deconvolution, coherency of raw data, and cross-correlation of 1-bit data.
- Interferometry.py -> runs the different techniques on 1 day of data recorded at "station_1" and "station_2". 
- station_1_1d.sac -> 1 day of record with a sampling rate of 4 Hz at station 1.
- station_2_1d.sac -> Same at station 2.
- Foo.png -> Output of Interferometry.py.

In this example, the anticausal and causal parts are strongly asymetric due to the location of the noise sources (Station_1 is closer to the ocean than Station_2).
