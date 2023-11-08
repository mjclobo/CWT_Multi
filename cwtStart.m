%% #################  cwtStart  ################# %%
% 
% This file properly formats input data, sets default values for constants and
% other variables, and sets all flag parameters based on input arguments.
% 
%% default argument options

pfResid           = 0;
pfResp            = 0;
pfFilt            = 0;
pfAmps            = 0;
pfM4M2            = 0;
pfLoD2            = 0;
pfResidSpectra    = 0;
decFact           = 20;       % default decimation factor
spFiltLength      = [363,181,145,73,65,65,65,65,65,65,65,65,65,65];         % filter lengths in hours for 4d, 2d, d1, d2, d3, d4, d6, d8
coFiltLength      = [1081,1081,363,363,1081,363,1081,1081,1081,1081,1081,363,363,1081,363,1081,1081,1081,363,1081,1081,363,1081,1081,65,65,65,65,65,65,65,65];   % see OPTIONS for corresponding freqs; lengths also in hours
spFiltLength_flag = 0;
coFiltLength_flag = 0;
dynInfFlag        = 0;
dynInfCustFlag    = 0;
coOmegaFlag       = 0;
spOmegaFlag       = 0;
refStnFlag        = 0;
bootstrapFlag     = 0;
bootstrapActive   = 0;
admitFlag         = 0;
spOnlyBool        = 0;
statWindowBool    = 0;
doMonthly         = 0;
alt_phase_bool    = 0;
rel_phase_bool    = 0;

coCutoffHi = 5*24;      % Cutoff frequency of 5 days for high-passed data
spCutoffHi = 5*24;
nData = length(dataIn);

%% checking for optional args and setting appropriate bools

while ~isempty(varargin)
    switch(varargin{1})
        case 'pfResid'
            pfResid = 1;
            varargin(1) = [];
        case 'pfResp'
            pfResp = 1;
            varargin(1) = [];
        case 'pfFilt'
            pfFilt = 1;
            varargin(1) = [];
        case 'pfAmps'
            pfAmps = 1;
            varargin(1) = [];
        case 'pfM4M2'
            pfM4M2 = 1;
            varargin(1) = [];
        case 'pfResidSpectra'
            pfResidSpectra = 1;
            varargin(1) = [];
        case 'pfLoD2'
            pfLoD2 = 1;
            varargin(1) = [];
        case 'spFiltLength'
            spFiltLength = varargin{2};
            varargin(1:2) = [];
        case 'coFiltLength'
            coFiltLength = varargin{2};
            varargin(1:2) = [];
        case 'dynamicInference'
            dynInfFlag = 1;
            varargin(1) = [];
        case 'dynamicInferenceCustom'
            dynInfCustFlag = 1;
            dynInf.names   = varargin{2};
            diFreqs        = {varargin{3}};
            diFiltLength   = varargin{4};
            varargin(1:4)  = [];
        case 'wtsIn'
            wtsIn = varargin(2);
            varargin(1:2) = [];
        case 'coFreqs'
            coOmegaFlag=1;
            coFreqs = varargin{2};
            coFiltLength = varargin{3};
            constits.names = varargin{4};
            varargin(1:4) = [];
        case 'spFreqs'
            spOmegaFlag=1;
            spFreqs = varargin{2};
            spFiltLength = varargin{3};
            species.names = varargin{4};
            varargin(1:4) = [];
        case 'performMonthly'
            doMonthly = 1;
            varargin(1) = [];
        case 'decimate'
            decFact = varargin{2};
            display('Your decimation factor is: '+string(decFact))
            varargin(1:2) = [];
            if mod(decFact,2) ~= 0 && decFact ~=1
                error('Your decimation factor must be even.')
            end
        case 'refStation'
            refStnFlag = 1;
            refData = varargin{2};
            varargin(1:2) = [];
        case 'bootstrapError'
            bootstrapFlag=1;
            nBoot = varargin{2};
            varargin(1:2)=[];
        case 'coCutoffHi'
            coCutoffHi = varargin{2};
            varargin(1:2) = [];
        case 'spCutoffHi'
            spCutoffHi = varargin{2};
            varargin(1:2) = [];
        case 'statWindow'
            statWindow = varargin{2};
            statWindowBool = 1;
            varargin(1:2) = [];
        case 'performAdmittance'
            admitFlag = 1;
            lon = varargin{2};
            lat = varargin{3};
            varargin(1:3) = [];
        case 'alt_phase'
            alt_phase_bool = 1;
            varargin(1) = [];
        case 'rel_phase'
            rel_phase_bool = 1;
            rel_phase_date = varargin{2};
            varargin(1:2)  = [];
         otherwise
            error(['Unrecognized argument: ' varargin{1}]);
    end
end

%% basic data info
here = fileparts(mfilename('fullpath'));

%% defining time-related variables
dtimeIn.TimeZone = 'UTC';

startDate = datevec(dtimeIn(1));
% timeIn = convertTo(dtimeIn,'epochtime','Epoch','1970-01-01');
if rel_phase_bool==1
    timeIn = double(convertTo(dtimeIn,'epochtime','Epoch',rel_phase_date));  
    % ^^this gives seconds(dtimeIn(1)-datetime(1899,12,31,12,0,0,'TimeZone','UTC')) - timeIn(1) = 0, as expected!
    
%     midInd = floor(length(dtimeIn)/2);
%     timeIn = dtimeIn + (rel_phase_date-dtimeIn(midInd));
%     z_date = datetime(0,0,0,0,0,0); z_date.TimeZone = 'UTC';
%     timeIn = seconds(timeIn-z_date);
else
    timeIn = convertTo(dtimeIn,'epochtime','Epoch',dtimeIn(1));
end

deltaT = double(median(diff(timeIn),'omitnan'));

%% reformatting data if needed;
if size(timeIn,1) > 1 && size(timeIn,2) > 1
    error('Your time input isnt a 1D vector!')
end

if size(dataIn,1) > 1 && size(dataIn,2) > 1
    error('Your data input isnt a 1D vector!')
elseif size(dataIn,2) > 1
    dataIn = dataIn.';
end

% downweighting NaNs
wtsIn = ones(size(dataIn))';

nanInd = find(isnan(dataIn));
% dataIn(nanInd) = 0.0009;
wtsIn(nanInd) = 10^-6;

nDataIn = length(dataIn);
nDataOdd = nDataIn-mod(nDataIn+1,2);       % modified this from t_tide

tPts = deltaT*([1:nDataOdd]-ceil(nDataOdd/2)); % SHOULD THIS BE FLOOR???

% account for even-length inouts
if nDataIn ~= nDataOdd
    dataIn = dataIn(1:end-1);
    wtsIn = wtsIn(1:end-1);
    nanInd = nanInd(1:end-1);
    
    halfInd = nDataIn/2;
    centerDate = dtimeIn(halfInd);
else
    halfInd = ceil(nDataIn/2);
    centerDate = dtimeIn(halfInd);
end

% more formatting
if size(tPts,1) > 1
    tPts = tPts';
end

nData = length(dataIn);

%% setting defaults for optional args
t1 = tPts(end);

midInd = floor(length(timeIn)/2);
% t1Orig = double(timeIn(midInd));
t1Orig = timeIn(1); %seconds(rel_phase_date - datetime(0,0,0,0,0,0,'TimeZone','UTC'));

grav = 9.81;
rad2deg = 180 / pi;

ampLimit = 7;       % max valid constit amplitude
ampFloor = 10^-6;
respThresh = 10^-4;
wtCrit = 0.9;       % fraction of good data (i.e., not missing) needed to report a wavelet result

%% checking for errors
if dynInfFlag==1 && bootstrapFlag==1
    error("Please only use the boostrap error method or perform dynamic inference. Doing both at the same time doesn't work with the code's architecture.")
end

if (dynInfFlag + dynInfCustFlag)==2
    error("Either use 'dynamicInference' or 'dynamicInferenceCustom' flag, not both.")
end

% make sure data is long enough to perform dynamic inference w/ 6mo filters
if dynInfFlag==1
    if (deltaT*nDataIn)<1.6e7
        error("Your data isn't long enough to perform dynamic inference with 6 month filters.")
    end
end

if nDataIn < max(coFiltLength) && nDataIn < max(spFiltLength)
    error("Your time series is shorter than maximum species and constituents filter lengths.")
elseif nDataIn < max(coFiltLength) && nDataIn >= max(spFiltLength)
    warning("Your time series is longer than max species filters, but shorter than max constituent filters. Performing species analysis only.")
    spOnlyBool = 1;
end

