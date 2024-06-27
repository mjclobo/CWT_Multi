%% #################  cwtFreqsFilts  ################# %%
% 
% This file defines frequencies for constituent and species analyses.
% It also creates the filter banks from which all filters will be convolved
% with the data.
% 

%% constituent frequencies
if coOmegaFlag==0
    
    constitsD1Freqs = [29.07266634,28.00621204, 26.868350, 25.81933871, 24.84120241, 23.93447213, 23.09848146, 22.30608083,21.578236629].^-1;
    constitsD1Names = ["Alp1","TwoQ1","Q1","O1","NO1","K1","J1","OO1","Ups1"];

    constitsD2Freqs = [13.1272668, 12.8717576, 12.65834751, 12.4206012, 12.19162085, 12,  11.75452172].^-1;
    constitsD2Names = ["Eps2","Mu2","N2","M2","L2","S2","Eta2"];

    constitsD3Freqs = [8.3863068, 8.280400802, 8.177140247, 7.9936].^-1;
    constitsD3Names = ["MO3","M3","MK3","SK3"];

    constitsD4PlusFreqs = [6.2691737, 6.21030, 6.103339, 6.0, 4.93088021306, 4.140200399, 3.52964079728,...
        3.10515029954, 2.748563985947, 2.48412023963, 2.25054027184, 2.07010019969].^-1;
    constitsD4PlusNames = ["MN4","M4","MS4","S4","MK5","M6","MK7","M8","MK9","M10","MK11","M12"];

    coFreqs = {constitsD1Freqs,constitsD2Freqs,constitsD3Freqs,constitsD4PlusFreqs};
    constits.names = horzcat(constitsD1Names,constitsD2Names,constitsD3Names,constitsD4PlusNames);

end

[coN,coND1,coND2,coND3,coND4] = deal(length(coFiltLength),length(coFreqs{1}),...
    length(coFreqs{2}),length(coFreqs{3}),length(coFreqs{4}));

if spOmegaFlag==0
    
    spFreqs = {(1/96), ...                              % 4D 
                 (1/48), ...                            % 2D 
                 (0.0417807462 * 0.96882) , ...         % D1 
                 (0.0805114007 * 1.00764), ...          % D2
                 (0.0417807462 + 0.0805114007), ...     % MK3
                 (0.0805114007 * 2) , ...               % M4
                 (0.0417807462 + 0.0805114007 * 2), ... % MK5
                 (0.0805114007 * 3) , ...               % M6 
                 (0.0805114007 * 3 + 0.0417807462), ... % MK7
                 (0.0805114007 * 4),...                 % M8
                 (0.0805114007 * 4 + 0.0417807462), ... % MK9
                 (0.0805114007 * 5),...                 % M10
                 (0.0805114007 * 5 + 0.0417807462), ... % MK11
                 (0.0805114007 * 6)};                   % M12
    
    species.names = ["low4d","low2d","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12"];

end

if spOmegaFlag == 0
    spNLow = 2;       % anything with frequency < diurnal
else
    spNLow = 0;
end

spN = length(spFiltLength);
spFreqWts = [ones(spNLow, 1); 2 .* ones(spN-spNLow,1)];   % downweight subtidal filters

coPeriods = cell2mat(coFreqs).^-1;
spPeriods = cell2mat(spFreqs).^-1;

[coOmega,coMaxFiltLength,coStartTime,coStartTimeOrig,coFilt,coNFilt,coNorm,coFiltLength] = filtBuild(coFreqs,coFiltLength,deltaT,t1,t1Orig,gabor_wavelet);
[spOmega,spMaxFiltLength,spStartTime,spStartTimeOrig,spFilt,spNFilt,spNorm,spFiltLength] = filtBuild(spFreqs,spFiltLength,deltaT,t1,t1Orig,gabor_wavelet);

[coNFiltD1,coNFiltD2] = deal(coNFilt(1:coND1),coNFilt(coND1+1:coND1+coND2));
[coNFiltD3,coNFiltD4] = deal(coNFilt(coND1+coND2+1:coND1+coND2+coND3),coNFilt(coND1+coND2+coND3+1:coN));

coOmegaFlat = cell2mat(coOmega);
spOmega = cell2mat(spOmega);

coTimeRel = @(omega,ptnum) (2 * pi * deltaT * decFact * (ptnum-1))/(3600 * coPeriods(omega));
coRadians = (2 * pi * coStartTime(1:coN)) ./ (coPeriods(1:coN) * 3600).';

% coTimeRelOrig = @(omega,ptnum) (2 * pi * deltaT * decFact * (ptnum-1))/(3600 * coPeriods(omega));
coRadiansOrig = (2 * pi * coStartTimeOrig(1:coN)) ./ (coPeriods(1:coN) * 3600).';

spTimeRel = @(omega,ptnum) (2 * pi * deltaT * decFact * (ptnum-1))/(3600 * spPeriods(omega));
spRadians = (2 * pi * spStartTime(1:spN)) ./ (spPeriods(1:spN) * 3600).';
spRadiansOrig = (2 * pi * spStartTimeOrig(1:spN)) ./ (spPeriods(1:spN) * 3600).';

% MONTHLY, TWO-WEEKLY, AND WEEKLY FILTERS
if doMonthly==1
    m2s2beat = (1/12-1/12.4206012)^-1/24;

    moFreqs = {2*m2s2beat, ...        % tidal monthly
        m2s2beat, ...      % two weekly
        m2s2beat/2};         % weekly

    moFiltLength = [4383,2017,1009];            % each are 6x period
    moN = length(moFiltLength);

    monthly.names = ["monthly","twoWeekly","weekly"];

    moPeriods = cell2mat(moFreqs).^-1;

    [moOmega,moMaxFiltLength,moStartTime,moStartTimeOrig,moFilt,moNFilt,moNorm,moFiltLength] = filtBuild(moFreqs,moFiltLength,deltaT,t1,t1Orig,gabor_wavelet);

    moOmega = cell2mat(moOmega);

    moTimeRel = @(omega,ptnum) (2 * pi * deltaT * decFact * (ptnum-1))/(3600 * moPeriods(omega));
    moRadians = (2 * pi * moStartTime(1:moN)) ./ (moPeriods((1:moN)) * 3600).';
    moRadiansOrig = (2 * pi * moStartTime(1:moN)) ./ (moPeriods((1:moN)) * 3600).';
end

if dynInfFlag==1
    
    diFreqs = {(1/12), ...        % S2
        (1/11.96723606), ...      % K2
        (1/12.01644934),...       % T2
        (1/23.93447213), ...      % K1 
        (1/24.06588766),...       % P1
        (1/24)};                  % S1
    
    diFiltLength = [4383,4383,4383,4383,4383,4383];
    diN = length(diFiltLength);
    
    dynInf.names = ["S2","K2","T2","K1","P1","S1"];
    dynInf.namesS2 = ["S2","K2","T2"];
    dynInf.namesK1 = ["K1","P1","S1"];

    diPeriods = cell2mat(diFreqs).^-1;

    [diOmega,diMaxFiltLength,diStartTime,diStartTimeOrig,diFilt,diNFilt,diNorm,diFiltLength] = filtBuild(diFreqs,diFiltLength,deltaT,t1,t1Orig,gabor_wavelet);

    diOmega = cell2mat(diOmega);

    diTimeRel = @(omega,ptnum) (2 * pi * deltaT * decFact * (ptnum-1))/(3600 * diPeriods(omega));
    diRadians = (2 * pi * diStartTime(1:diN)) ./ (diPeriods((1:diN)) * 3600).';
    diRadiansOrig = (2 * pi * diStartTimeOrig(1:diN)) ./ (diPeriods((1:diN)) * 3600).';

end
if dynInfCustFlag==1

    diN = length(diFiltLength);

    diPeriods = cell2mat(diFreqs).^-1;

    [diOmega,diMaxFiltLength,diStartTime,diStartTimeOrig,diFilt,diNFilt,diNorm,diFiltLength] = filtBuild(diFreqs,diFiltLength,deltaT,t1,t1Orig,gabor_wavelet);

    diOmega = cell2mat(diOmega);

    diTimeRel = @(omega,ptnum) (2 * pi * deltaT * decFact * (ptnum-1))/(3600 * diPeriods(omega));
    diRadians = (2 * pi * diStartTime(1:diN)) ./ (diPeriods((1:diN)) * 3600).';
    diRadiansOrig = (2 * pi * diStartTimeOrig(1:diN)) ./ (diPeriods((1:diN)) * 3600).';
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [omega,maxFiltLength,startTime,startTimeOrig,filt,nFilt,norm,filtLength] = filtBuild(fFreqs,filtLength,deltaT,t1,t1Orig,gabor_wavelet)
    
    N = length(filtLength);

    omega = cellfun(@(x) ((2*pi)/3600) * x, fFreqs,'UniformOutput',false);
    
    omegaFlat = cell2mat(omega);

    %% More parameters/optional args
    periods =  ((2 * pi)/3600) ./ omegaFlat;

    % makes sure filter lengths are odd
    filtLength = (2 * floor(0.5 * filtLength *(3600/deltaT)) + 1)';

    %% define some time-related variables
    % center times of the analysis period
    maxFiltLength = max(filtLength);

    startTime = (-t1 + deltaT * (maxFiltLength-1) / 2) .* ones(N,1);
%     startTimeOrig = (t1Orig + deltaT * (maxFiltLength-1) / 2) .* ones(N,1);
    startTimeOrig = double(double(t1Orig) + deltaT * (maxFiltLength-1) / 2) .* ones(N,1);

    
    % startTimeOrig = t1Orig .* ones(N,1); % idea here is that initial phase is just relative to start of signal, regardless of when filters start
    % startTimeOrig = t1Orig + deltaT * ((maxFiltLength-1)/2 - (filtLength-1)/2) ; % idea here is that we need to start at nStart in line 97 of cwtConvolveFilters.m
    

    %% Construct filters
    % define the halflength parameter for each tidal and subtidal frequency
    tMid = deltaT*((filtLength-1)/2);        % in sec

    % CWT definitions 
    beta=4.6;        % roll-off parameter of tidal wavelet filters
    c11 = 0.97;

    % for tidal frequencies
    tcon = 1.11504;

    %% Define wavelet functions
    tidalwavelet = @(t,om,sp,TT) om * exp(-1i * om * t) .* (besseli(0,beta * sqrt(1-((tcon*t) ./ (TT(sp))).^2)) / besseli(0,beta));

    %% Define filter properties
    filtFunc = @(om) tidalwavelet(-tMid(om) + ((1:filtLength(om))-1)*deltaT, omegaFlat(om),om,tMid); 
    
    filt = cellit(filtFunc,(1:N));
    
    %% Gabor transform test
    if gabor_wavelet==1
        gabor_filt_func = @(t,om,t0,a) exp(-(t-t0).^2./(a^2)) .* exp(-1i * om .* (t - t0));

        filtFunc = @(om) gabor_filt_func(-tMid(om) + ((1:filtLength(om))-1)*deltaT, omegaFlat(om),0,0.4*tMid(om));
        filt = cellit(filtFunc,(1:N));
    end
    %% Filter normalization
    testc = @(t,om) cos(om * t);

    cdata = @(sp,flngthlist,freqlist) testc((1:flngthlist(sp)) .* deltaT - (flngthlist(sp)+1)/2 * deltaT,freqlist(sp));    % cos wave at filter frequencies, of filter length, centered at t=0
    normwts = @(sp,flngthlist) ones(flngthlist(sp),1);

    % normalizing filter parameters
    norm = cell2mat(arrayfun(@(x) chop(x{1,1:end},24), cellit(@(it1) NConvolved([],cdata(it1,filtLength,omegaFlat),filt{it1},normwts(it1,filtLength).',1,0.9).^-1,1:N),'UniformOutput',false));
    
    % normalize filters
    nFilt = cellfun(@(x) x ,cellit(@(it1) filt{it1} .* norm(it1),1:length(norm)),'UniformOutput',false);

end
