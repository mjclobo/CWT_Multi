function [constits,species,ref,admittances] = cwtMulti(dtimeIn,dataIn,varargin)

% This routine performs continuous wavelet transform analysis on a scalar
% time series; in our case, water levels.
% The analysis is done using both (a) major constituents within tidal species,
% and (b) one weighted tidal frequency per tidal species.
% The latter is useful when time length scales of interest don't allow
% for resolution between 'close' constituents, e.g., storm surges.
% Weights of tidal species frequencies are hard-wired, but can be changed
% manually in cwtFreqsFilts.m file; or you may pass new species frequencies as
% an optional arg.
% 
% A couple of naming conventions: co - constits, sp - species, pf - plot flag 
%     
% *** INPUTS ***
% 
% dtimeIn: nx1 time vector in datetime format. Must align with dataIn vector,
%     temporally
% 
% dataIn: nx1 scalar series of values of interest
% 
% lon: longitude corresponding to data station
% 
% lat: latitude corresponding to data station
% 
% *** OPTIONAL ARGUMENTS ***
%
% 'performAdmittance',lon,lat: this will have the routine perform a tidal
%     admittance analysis. The default reference station is the tidal
%     potential at the lat/lon provided as inputs. The user can also
%     include the 'refStation' argument, described below.
% 
% 'spFiltLength',[length1,length2,etc.]: this changes filter lengths (hours) from default
%     for corresponding species (single, weighted) frequencies; MUST BE
%     ODD; see cwtFreqsFilts.m for default frequencies (and their order)
% 
% 'coFiltLength',[length1,length2,etc.]: filter length (hours) for
%     corresponding constituents, unless specific constituent frequencies
%     are defined below, then this is overridden; MUST BE ODD; see
%     cwtFreqsFilts.m for default frequencies (and their order)
% 
% 'spFreqs',{FREQS},[LENGTHS],[NAMES]: define your own frequencies (in hr^-1)
%     for tidal species analysis; include filter lengths in hours.
%     NAMES is vector of corresponding names, each in double quotes, i.e., "name".
% 
% 'coFreqs',{[D1 FREQS],[D2 FREQS],[D3 FREQS],[D4+ FREQS]},[LENGTHS],[NAMES]: define your own frequencies (in hr^-1) for
%     tidal species analysis; note that you must pass in corresponding filter 
%     lengths (in hours); organize FREQS, LENGTHS, and NAMES with braces around the whole list, and
%     square brackets around groups of tidal species, e.g.
%     {[Q1,O1,K1],[N2,M2,S2],[MO3,MK3],[MN4,M4,MS4]}.
%     NAMES must be strings in double quotes, not single, i.e., "name".
% 
% 'decimate',decFact: decimation factor for analysis outputs in time steps; MUST BE
%     EVEN! Default is 20 time steps.
% 
% 'refStation',refData: Use reference station data, rather than
%     tidal potential, as the denominator when calculating admittances.
%     This data set must have exact same timestamps/length as input data.
%     User must also include 'performAdmittance' argument!
% 
% 'bootstrapError',nBoot:  This performs the bootstrap error method as
%     described by Innocenti et al., 2021. nBoot is how many residual time
%     series you want to use to create confidence intervals. 1,000 is
%     robust, but computationally intensive. 
% 
% 'dynamicInference': performs dynamic inference of K2/S2 and P1/K1 using
%     6-month long filters. MUST HAVE AT LEAST 6 MONTHS OF DATA.
%
% 'dynamicInferenceCustom',[allNames],[baseFreq,auxFreqs],
%     [baseFreqLength,auxFreqsLengths]: Perform dynamic inference using a
%     custom set of filters and lengths. FREQS IN HR^-1!
% 
% 'coCutoffHi',coCutoffHi: If you want to change the cutoff frequency for
%     the highpassed constituents data. Passband is automatically two days.
% 
% 'spCutoffHi',coCutoffHi: If you want to change the cutoff frequency for
%     the highpassed species data. Passband is automatically two days.
% 
% 'pfResid': Returns plots of original/reconstructed time series
%     for constits and species; also returns plots of the residuals; ignore
%     edge effects for now
% 
% 'pfResp': Returns plot of frequency responses for all constituent
%     filters; look for overlap within tidal species
% 
% 'pfFilt': Returns plot of real and imaginary parts of example
%     filter; currently set to species semidiurnal filter
% 
% 'pfAmps': Returns two plots for amps found from constits and
%     species analyses
% 
% 'pfM4M2': Returns plot of M2 and M4 amps as (M4/M2)^2, and
%     separately
% 
% 'pfResidSpectra': Plots spectrum of residuals for both species
%     and constits results
% 
% 'pfLoD2': Plot time series of D2 amplitude(s) versus low-passed
%     data for both species and constits analysis
%
% 'alt_phase': There is the phase referenced to the center of the analysis
%     window and that referenced to the beginning of the analysis window.
%     The default is the latter, and passing this arg gives the former.
% 'pct_valid', wtCrit: This defines the percentage of data that must be
%     valid (i.e., not NaN) in order to convolve the filter with the data
%     and return a complex filter response value. Note that the window
%     length is set by the filter length. So one could have an analysis
%     where, for example, D2 amplitudes are returned but D1 amplitudes are
%     not, if the D2 filter is shorter than the D1 filter. Default is 0.9.
% 
% *** OUTPUTS ***
% 
% constits: struct type containing results from CWT analysis
%     that uses major constituents within species
% 
% species: struct type containing results from CWT analysis that
%     uses a single modified frequency per tidal species
% 
% ref: struct type containing analysis results of reference station (default is tidal
%     potential) data, for both species and constits analyses. This will be
%     empty if 'performAdmittance' argument isn't used.
% 
% admittances: struct type containing tidal admittances and absolute
%     phases. This will be empty if 'performAdmittance' argument isn't
%     used.

%% Perform CWT analysis on input data
cwtStart;               % formatting data, defining values, setting flags
cwtFreqsFilts;          % assigning filter-related variables and constructing filters
cwtConvolveFilters;     % setting up & carrying out complex filter convolution (after highpassing data)
cwtSolve;               % solve response coefficient matrices
cwtDefineOutputs;       % assigning output values

cwtErrors;              % error analysis
cwtPlots;               % plotting optional plot args

%% Dynamic inference
cwtDynamicInference;    % perform dynamic inference and plot results

%% Possibly admittance
if admitFlag==1
    
    % Perform CWT analysis on reference station data (default is tidal potential)
    cwtRefStn;              % format reference station data
    cwtConvolveFilters;
    cwtSolve;
    cwtDefineRefOutputs;

    % Define admittances
    cwtDefineAdmits;        % define admittances as dataIn/refData results
    
else
    
    ref = [];
    admittances = [];
    
end

end