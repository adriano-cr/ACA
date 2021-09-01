# ACA exam project

Program to compute the FFT of a vector of 2^n elements.

For the parallel versions, currently:

**V1:**
- v1.0 only "reverse" parallelized, default schedule
- v1.1 "reverse" and "order" parallelized, default schedule
- v1.2 "reverse", "order" and "transform" parallelized, default schedule
- v1.3 "reverse", "order", "transform" and "printvec" parallelized, default schedule

**V2:**
- v2.0 same as v1.3 but writes to file instead of printing to screen
- v2.1 file read is now in its own function instead of main

The data files are randomly generated vectors of length 32 (data.txt), 1024 (data-long.txt) and 1048576 (data-huge.txt). sax.txt and violin.txt contain samples from two short audio clips, sampled at 44.1 kHz, and are both vectors of 262944 elements.

**TODO:**
- [ ] swap while loops with for loops; 
- [x] print to file instead of screen; 
- [ ] finish parallelization;
