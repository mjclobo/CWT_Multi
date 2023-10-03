% This script creates a nonstationary signal with tidal and subtidal frequencies
% and then exports the data to a .mat file
% 

% define frequencies
fm2 = 1/12.4206012;         % M2
fs2 = 1/12;                 % S2
fn2 = 1/12.65834751;        % N2
fk1 = 1/23.93447213;        % K1
fm4 = 1/6.210300601;        % M4
fo1 = 1/25.81933871;        % O1
fan = 1/(365.25*24);             % annual signal will be river flow
% fr1 = 1/(30*24);            % random monthly frequency
% fr2 = 1/(7*24);             % random weekly frequency 
% fr3 = 1/(14*24);            % random biweekly signal

% define sampling frequency/times
fs = 1/3600;
dat.t = 0:1/fs:(365.25*24*3600);           % i.e., 1 yr of hourly data, values are in seconds

freq1 = [fm2 fk1 fn2 fs2 fm4 fo1 fan] / 3600;   % change frequencies to from hr^-1 to Hz

amps = [0.5 0.42 0.23 0.156 0.06 0.119 0.75];  % define amplitudes for corresponding frequencies

% add an annual modulation
ann = sin(2*pi*(fan/3600)*dat.t);
ann = ann+abs(min(ann));             % only positive modulation
annmod = ann/(max(ann));

%%
data_nonstationary = sum((amps(1:2)'*annmod).*sin(2*pi*freq1(1:2).'*dat.t +2*pi*[1,0.5]'),1) + sum(amps(3:end)*sin(2*pi*freq1(3:end).'*dat.t+2*pi*[-0.1, -0.37, -1, -0.66, 0]'),1);

% data = amps' .* sin((2*pi*freq1' * dat.t + 2*pi*[1, 0.5, -0.1, -0.37, -1, -0.66, 0]'));   % add random phases between 0 and 2pi, no phase for annual

dat.wl = data_nonstationary;
dat.wlnoise = sum(data_nonstationary + median(amps) * randn(size(dat.t)));              % sum data column-wise, i.e. superimpose all frequencies to make one signal

dat.info = 'This file contains artificially created water level data (m) and time values (s)';

save nonstationary.mat dat

















