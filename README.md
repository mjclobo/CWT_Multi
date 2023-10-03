# CWT_Multi
A MATLAB routine used to calculate nonstationary tidal amplitudes and phases, as well as water level reconstructions, using a special form of continuous wavelet transform (CWT) analysis.

## METHOD ##
CWT_Multi.m is a MATLAB routine that performs continuous wavelet
transform (CWT) analysis on an Nx1 real, scalar time series.
This is done using two different sets of filter banks.

The species analysis uses one filter per tidal species
to identify changes in tidal species energy on a short
timescale (approx. 80% of energy found within a one week window).
The filters are shorter in time and therefore have a wider
response in frequency-space.

The constituents analysis uses a suite of filters to identify changes
in major constituents within a tidal species on a still
relatively short time scale (approx. 80% of energy found within a
one week window), while also tracking minor constituents on longer
timescales (approx. 80% of energy found within a three week window).
The details of this analysis are covered in the manuscript.

CWT_Multi.m also includes a dynamic inference feature (accessed by using
the optional argument 'dynamicInference') that will use results from
six month filters of K1, P1 and S2, K2, respectively, to form
time-dependent ratios, which are then used to infer the amplitudes of
the respective constituents on the timescale of the length of the major
constituent filters from the constituents analysis.

Note that due to two-way, three-way, and higher beating interactions
between tidal constituents, the results must be interpreted partly
by understanding these interactions and their manifestations in the
results.

## CODE ##
CWT_Multi.m is loosely written in a modular sense, i.e., most parts of the script are calls
to general functions which can take arguments for either the species or
the constituent analysis.
This modular approach increases the readability and reliability of the code, as well
as benefits the code's further development.
Work still needs to be done to make the CWT_Errors.m sub-script
more modular.

The routine, and its subroutines, have comments throughout, and one
can follow the routine's work flow by first opening CWT_Multi.m,
and then following the script (which calls ``sub-scripts'').
A summary of the inputs, outputs, and optional arguments for
this routine is included in the header of CWT_Multi.m.

If the user is solely interested in using the routine, please refer
to the required (and optional) arguments/format found in the header
of the main CWT_Multi.m file.
This header also details the routine's outputs.


