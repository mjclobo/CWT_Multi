% This script creates a stationary signal with tidal and subtidal frequencies
% and then exports the data to a .mat file
% 

% define frequencies
fm2 = 1/12.4206012;         % M2
fs2 = 1/12;                 % S2
fn2 = 1/12.65834751;        % N2
fk1 = 1/23.93447213;        % K1
fm4 = 1/6.210300601;        % M4
fo1 = 1/25.81933871;        % O1
% fan = 1/(365.25*24);             % annual signal will be river flow
% fr1 = 1/(30*24);            % random monthly frequency
% fr2 = 1/(7*24);             % random weekly frequency 
% fr3 = 1/(14*24);            % random biweekly signal

% define sampling frequency/times
fs = 1/(3600);
dat.t = 0:1/fs:(365.25*24*3600);           % i.e., 1 yr of hourly data, but values are in seconds

freq1 = [fm2 fs2 fn2 fk1 fm4 fo1] / 3600;   % change frequencies to from hr^-1 to Hz

amps = [0.5 0.42 0.23 0.156 0.06 0.119];  % define amplitudes for corresponding frequencies


data = amps' .* sin((2*pi*freq1' * dat.t + 2*pi*[1, 0.5, 1.4, -0.37, -1, 2]'));   % add random phases between 0 and 2pi

dat.wl = sum(data);
dat.wlnoise = sum(data) + 0.5 * randn(size(dat.t));              % sum data column-wise, i.e. superimpose all frequencies to make one signal

dat.info = 'This file contains artificially created water level data (m) and time values (s)';

save stationary.mat dat

















