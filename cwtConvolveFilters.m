%% #################  cwtConvolve Filters  ################# %%
% 
% This file highpasses the data, convolves the data with all filters, and
% stores the outputs from the convolutions.
% The variables ampFloor and ampLimit are taken into account when
% validating the filter outputs.
% 

%% Highpass time series
if bootstrapActive==0
    [coDataInLo,coDataInHi] = idealLowPass(dataIn,wtsIn,coCutoffHi,deltaT,nanInd);
    [spDataInLo,spDataInHi] = idealLowPass(dataIn,wtsIn,spCutoffHi,deltaT,nanInd);
else
    [coDataInLo,coDataInHi] = idealLowPass(coDataIn,wtsIn,coCutoffHi,deltaT,nanInd);
    [spDataInLo,spDataInHi] = idealLowPass(spDataIn,wtsIn,spCutoffHi,deltaT,nanInd);
end

%% Convolve data with wavelet filters
if spOnlyBool==0
    [coCompConv,coMean,coDataCent,coRefTime] = waveletConv(coDataInHi,wtsIn,coMaxFiltLength,coFilt,coNorm,coN,decFact,wtCrit,nData,timeIn);
end

[spCompConv,spMean,spDataCent,spRefTime] = waveletConv(spDataInHi,wtsIn,spMaxFiltLength,spFilt,spNorm,spN,decFact,wtCrit,nData,timeIn);

%% Defines decimated time axes and unpacks CWT outputs
if spOnlyBool==0
    [coTimes,coOdd,coNDec] = defineTimes(tPts,coMaxFiltLength,decFact,nData);
end
[spTimes,spOdd,spNDec] = defineTimes(tPts,spMaxFiltLength,decFact,nData);

%% store and weight
if spOnlyBool==0
    [coCompConvTrim,coWts] = storeWeight(coCompConv,coNDec,coN,coOdd,ampFloor,ampLimit);
end

[spCompConvTrim,spWts] = storeWeight(spCompConv,spNDec,spN,spOdd,ampFloor,ampLimit);

%% do same for tidal monthly and two-weekly filters
if doMonthly==1
    % Convolve data with wavelet filters
    [moCompConv,moMean,moDataCent,moRefTime] = waveletConv(dataIn,wtsIn,moMaxFiltLength,moFilt,moNorm,moN,decFact,wtCrit,nData,timeIn);  % convolves monthly etc. filters with original data

    % Defines time axes and unpack CWT outputs
    [moTimes,moOdd,moNDec] = defineTimes(tPts,moMaxFiltLength,decFact,nData);

    % store and weight
    [moCompConvTrim,moWts] = storeWeight(moCompConv,moNDec,moN,moOdd,ampFloor,ampLimit);
end
%% Dynamic inference
if dynInfFlag==1 || dynInfCustFlag==1
    
    % Highpass time series
    [diDataInLo,diDataInHi] = idealLowPass(dataIn,wtsIn,diMaxFiltLength,deltaT,nanInd);

    % Convolve data with wavelet filters
    [diCompConv,diMean,diDataCent,diRefTime] = waveletConv(coDataInHi,wtsIn,diMaxFiltLength,diFilt,diNorm,diN,decFact,wtCrit,nData,timeIn);

    % Defines time axes and unpack CWT outputs
    [diTimes,diOdd,diNDec] = defineTimes(tPts,diMaxFiltLength,decFact,nData);

    % store and weight
    [diCompConvTrim,diWts] = storeWeight(diCompConv,diNDec,diN,diOdd,ampFloor,ampLimit);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataLo,dataHi] = idealLowPass(data,wts,maxFL,deltaT,nanInd)
    % stopband frequency is inverse of maximum filter length
    % transition band is one day wide

    % applying low-pass filter to data using optimal low-pass windowed by Kaiser fcn
    % defining filter criteria
    stop = (maxFL*3600)^-1;
    pass = ((maxFL+2*24)*3600)^-1;     % i.e., pass band of 2 cyc/day

    passPctLo = 0.1; % 5 % passband ripple
    stopDbLo = 0.005; % 47db stopband attenuation

    % defining optimal filters
    [sp_n,wn,beta,ftype] = kaiserord([pass,stop],[1,0],[passPctLo,stopDbLo],deltaT^-1);
    filtLo = fir1(sp_n,wn,ftype,kaiser(sp_n+1,beta),'noscale');
    
    % applying filters
    data(nanInd) = nanmean(data);
    
    dataLo = conv(data,filtLo,'same');
    dataHi = data - dataLo;

end

function [result,meanO,dataCent,refTime] = waveletConv(dataHi,wts,maxFL,filt,norm,N,decFact,wtCrit,nData,times)
    
    meanO = mean(dataHi.*wts.');
    dataCent = dataHi.*wts.' - meanO;
    
    nStart = 1 + cell2mat(cellit(@(it1) (maxFL-1)/2 - (length(filt{it1})-1)/2,1:N));
    nEnd = nData - nStart + 1;
    
    result = pardo2(N,@(it1) norm(it1) * NConvolved([],dataCent(nStart(it1):nEnd(it1)),filt{it1}.',wts.',decFact,wtCrit));
    
    refTime = cell2mat(pardo2(1,@(it1) NConvolved(times(nStart(it1):nEnd(it1)),dataCent(nStart(it1):nEnd(it1)),filt{it1}.',wts.',decFact,wtCrit)));

end

function [times,odd,nDec] = defineTimes(tPts,maxFL,decFact,nData)
    
    times = tPts(floor((maxFL-1)/2):decFact:nData-floor((maxFL-1)/2));      % times without edges
    
    if mod(length(times),2) == 1
        odd = 1;
    else
        times = times(1:end-1);
        odd = 0;
    end
    
    nDec = length(times);
    
end

function [compConvTrim,wts] = storeWeight(compConv,nDec,N,oddBool,ampFloor,ampLimit)
    compConvTrim = zeros(nDec,N);
    wts = zeros(nDec,N);

    for i=1:N
        if oddBool==1
            compConvTrim(:,i) = compConv{i};
            wts(abs(compConv{i})>ampFloor & abs(compConv{i})<ampLimit,i) = 1;
        else
            compConvTrim(:,i) = compConv{i}(1:end-1);
            wts(abs(compConv{i}(1:end-1))>ampFloor & abs(compConv{i}(1:end-1))<ampLimit,i) = 1;        
        end
    end
end

