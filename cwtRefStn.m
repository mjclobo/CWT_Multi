%% Set up reference station for CWT analysis

% NOTE that we assume there aren't any missing time values
start_time = datenum(startDate);
end_time = start_time+length(tPts)*deltaT/3600/24;      % datenum format of end date

if refStnFlag == 1
    dataIn = refData;
else
    [gmTime,poten] = phortran_2020(start_time,end_time,lon,lat,deltaT);
    dataIn = poten;
end

if size(dataIn,1) > 1 && size(dataIn,2) > 1
    error('Your data input isnt a 1D vector!')
elseif size(dataIn,2) > 1
    dataIn = dataIn.';
end
nDataIn = length(dataIn);
nDataOdd = nDataIn-mod(nDataIn+1,2);       % modified this from t_tide

% account for even-length inouts
if nDataIn ~= nDataOdd
    dataIn = dataIn(1:end-1);
end

nanInd = find(isnan(dataIn));
dataIn(nanInd) = 0.0009;
wtsIn(nanInd) = 10^-6;

if length(dataIn) ~= length(tPts)
    error("Your tidal potential (or reference station data) isn't the same length as your time array; maybe you forgot to include all timepoints, and NaN-fill missing data?")
end

if size(dataIn,1)==size(wtsIn,1)
    wtsIn = wtsIn.';
end

ref.wl = dataIn;

ref.species.names = species.names;
ref.constits.names = constits.names;
