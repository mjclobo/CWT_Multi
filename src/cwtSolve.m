%% #################  cwtSolve  ################# %%
% 
% This file solves all response coefficient matrix problems,
% leading to amplitudes and phases for all resolved constituents.
% It also creates reconstructed data from these solutions.
% 

%% Solve using responses within tidal species
if spOnlyBool==0
    coSoln = zeros(coN,coNDec);

    [coSoln(1:coND1,:),coRespD1] = invertSol(coCompConvTrim,coNDec,1,coND1,respThresh,coNFiltD1,coOmega{1},deltaT);
    [coSoln(coND1+1:coND1+coND2,:),coRespD2] = invertSol(coCompConvTrim,coNDec,coND1+1,coND1+coND2,respThresh,coNFiltD2,coOmega{2},deltaT);
    [coSoln(coND1+coND2+1:coND1+coND2+coND3,:),coRespD3] = invertSol(coCompConvTrim,coNDec,coND1+coND2+1,coND1+coND2+coND3,respThresh,coNFiltD3,coOmega{3},deltaT);
    [coSoln(coND1+coND2+coND3+1:coND1+coND2+coND3+coND4,:),coRespD4] = invertSol(coCompConvTrim,coNDec,coND1+coND2+coND3+1,coND1+coND2+coND3+coND4,respThresh,coNFiltD4,coOmega{4},deltaT);

    [coSolnWts,coSolnAmps,coSolnPhases,coTimesAll,coSolnFunc] = solnProcess(coTimes,coSoln,coN,coNDec,coRadians,coRadiansOrig,tPts,coTimeRel,rad2deg,ampFloor,ampLimit,coMaxFiltLength,decFact,nData,coRefTime,coOmegaFlat,alt_phase_bool);
    
    if bootstrapActive==0
        [coReconAll,coReconHi,coReconWts,coNRecon] = reconSoln(coSolnFunc,tPts,wtsIn,coDataInLo,coMean,cell2mat(coOmega),nData,coN,ampLimit);
    end
end


%% Solve using responses between species filters
[spSoln,spResp] = invertSol(spCompConvTrim,spNDec,1,spN,respThresh,spNFilt,spOmega,deltaT);

[spSolnWts,spSolnAmps,spSolnPhases,spTimesAll,spSolnFunc] = solnProcess(spTimes,spSoln,spN,spNDec,spRadians,spRadiansOrig,tPts,spTimeRel,rad2deg,ampFloor,ampLimit,spMaxFiltLength,decFact,nData,spRefTime,spOmega,alt_phase_bool);

if bootstrapActive==0
    [spReconAll,spReconHi,spReconWts,spNRecon] = reconSoln(spSolnFunc,tPts,wtsIn,spDataInLo,spMean,spOmega,nData,spN,ampLimit);
end

%% monthly etc. filters
if doMonthly==1
    [moSoln,moResp] = invertSol(moCompConvTrim,moNDec,1,moN,respThresh,moNFilt,moOmega,deltaT);

    [moSolnWts,moSolnAmps,moSolnPhases,moTimesAll,moSolnFunc] = solnProcess(moTimes,moSoln,moN,moNDec,moRadians,moRadiansOrig,tPts,moTimeRel,rad2deg,ampFloor,ampLimit,moMaxFiltLength,decFact,nData,moRefTime,moOmega,alt_phase_bool);
end

%% and possibly dynamic inference solution
if dynInfFlag==1 || dynInfCustFlag==1
    
    [diSoln,diResp] = invertSol(diCompConvTrim,diNDec,1,diN,respThresh,diNFilt,diOmega,deltaT);

    [diSolnWts,diSolnAmps,diSolnPhases,diTimesAll,diSolnFunc] = solnProcess(diTimes,diSoln,diN,diNDec,diRadians,diRadiansOrig,tPts,diTimeRel,rad2deg,ampFloor,ampLimit,diMaxFiltLength,decFact,nData,diRefTime,diOmega,alt_phase_bool);
    
    if bootstrapActive==0
        [diReconAll,diReconHi,diReconWts,diNRecon] = reconSoln(diSolnFunc,tPts,wtsIn,diDataInLo,diMean,diOmega,nData,diN,ampLimit);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resp = freqResp(filts,freqs,deltaT)
% defines frequency response matrix between all frequencies and filters
%

    Nfilt = length(filts);
    Nfreq = length(freqs);
    resp = zeros(Nfilt,Nfreq);

    for i =1:Nfilt
        for k=1:Nfreq
            N2 = (length(filts{i})-1)/2;
            t = (-N2:N2) * deltaT;
            c = cos(freqs(k) * t);

            % for NConvolved() data must be same dim
            resp(i,k) = real(NConvolved([],c,filts{i},ones(size(filts{i})),1,0.5));
        end
    end
    
end

function [solnAll,validResp] = invertSol(filtOutput,nVals,indStart,indStop,respThresh,cnf,cod,dt)
% finds solution for matrix inversion involving filter response matrix
% removes constituents with small responses in order to maintain stability
% 

    solnAll = zeros(indStop-indStart+1,nVals);

    for j=1:nVals
        
        dataResp = filtOutput(j,indStart:indStop);
        validVals = find(abs(dataResp)> respThresh);
        invalidVals = find(abs(dataResp)<= respThresh);
        invec = dataResp(validVals);
        validResp = freqResp(cnf(validVals),cod(validVals),dt);
        
        if isempty(validVals)
            solnAll(:,j) = 0;
        elseif length(validVals) ~= (indStop-indStart+1)
            solnAll(validVals,j) = validResp\invec.';
            solnAll(invalidVals,j) = 0.0;
        else
            solnAll(:,j) = validResp\invec.';
        end
        
    end
    
end

function [solnWts,solnAmpsAll,phaseOut,timesAll,solnFunc] = solnProcess(timesIn,soln,N,nDec,radians,radiansOrig,tPts,timeRel,rad2deg,ampFloor,ampLimit,maxFL,decFact,nData,tConv,omega,alt_phase)
% Defines most outputs of interest from the solutions to the response
% coefficient matrix problem.
% Also extrapolates out to time ends, using constant values.
% 

    solnWts = cellit(@(it2,it3) (ampFloor < abs(soln(it2,it3)) < ampLimit)*1,1:N,1:nDec);

    solnAmps = cell2mat(cellit(@(it2,it3) (solnWts{it2,it3} > 0.01) * abs(soln(it2,it3)) +...
        (solnWts{it2,it3} < 0.01) * (0),1:N,1:nDec));
    
    solnPhases = cell2mat(cellit(@(it2,it3)(solnWts{it2,it3} > 0.01) * mod(radians(it2) + timeRel(it2,it3) - angle(soln(it2,it3)),2*pi) * ...
        rad2deg + (solnWts{it2,it3} < 0.01) * (0),1:N,1:nDec));    % angle(exp(1i * timeRel(it2,it3))) = timeRel(it2,it3)
    
    solnAmpsAll = cell2mat(cellit(@(it2) [solnAmps(it2,1) * ones(length(1:decFact:floor(maxFL/2)-1),1) ; solnAmps(it2,:).';...
        solnAmps(it2,nDec) * ones(length(nData-floor((maxFL/2))+1:decFact:nData),1)] ,1:N)).';
    
    solnPhasesAll = cell2mat(cellit(@(it2) [solnPhases(it2,1) * ones(length(1:decFact:floor(maxFL/2)-1),1) ; solnPhases(it2,:).';...
        solnPhases(it2,nDec) * ones(length(nData-floor((maxFL/2))+1:decFact:nData),1)] ,1:N)).';
    
    timesAll = [tPts(1:decFact:floor(maxFL/2)-1) , timesIn, tPts(nData - floor((maxFL/2))+1:decFact:nData)];
    
    solnsBounds = cell2mat(cellit(@(it2) (ampFloor < solnAmpsAll(it2,:) < ampLimit) .* solnAmpsAll(it2,:) .* exp(-1i * solnPhasesAll(it2,:) * pi/180) +...
        (ampFloor > solnAmpsAll(it2,:)  + solnAmpsAll(it2,:) > ampLimit) * (0),1:N).');

    solnFunc = cellit(@(it2) spline(timesAll,solnsBounds(it2,:)),1:N);

%     solnPhasesUnshifted = cell2mat(cellit(@(it2,it3)(solnWts{it2,it3} > 0.01) * mod(radiansOrig(it2) + timeRel(it2,it3) - angle(soln(it2,it3)),2*pi) * rad2deg +...
%         (solnWts{it2,it3} < 0.01) * (0),1:N,1:nDec));           % phases relative to non-centered time axis (i.e., agrees with UTIDE) ; angle(exp(1i * timeRel(it2,it3)))
    
    % tConv is reftime, which is center time of analysis window (see NConvolved() definition for details).

    solnPhasesUnshifted = cell2mat(cellit(@(it2,it3)(solnWts{it2,it3} > 0.01) * mod(tConv(it3)*omega(it2) - angle(soln(it2,it3)),2*pi) * rad2deg +...
        (solnWts{it2,it3} < 0.01) * (0),1:N,1:nDec));   % THIS GIVES PHASES RELATIVE TO startDate (NEEDED FOR ARTIFICIAL DATA)
    
    solnPhasesUnshiftedAll = cell2mat(cellit(@(it2) [solnPhasesUnshifted(it2,1) * ones(length(1:decFact:floor(maxFL/2)-1),1) ; solnPhasesUnshifted(it2,:).';...
        solnPhasesUnshifted(it2,nDec) * ones(length(nData-floor((maxFL/2))+1:decFact:nData),1)] ,1:N)).';
    
    if alt_phase==1
        phaseOut = solnPhasesAll;
    else
        phaseOut = solnPhasesUnshiftedAll;
    end
    
end

function [reconAll,reconHi,reconWts,nRecon] = reconSoln(solnFunc,tPts,wtsIn,dataInLo,meanI,omega,nData,N,ampLimit)
% Builds reconstruction of data using the complex solutions to response
% coefficient matrix problems.
% This is where the low-frequency data may be added back to result.
% 
    reconMat = zeros(N,nData);

    for i=1:nData
        for j=1:N
            reconMat(j,i)= real(ppvalFast(solnFunc{j},tPts(i)) * exp(1i*(omega(j) * tPts(i))));
        end
    end

    nRecon = N;

    reconWts = cellit(@(it2,it3) (abs(reconMat(it2,it3)) < ampLimit) * 1, 1:nRecon,1:nData);

    reconHi = cell2mat(cellit(@(it2) sum(reconMat(:,it2))*wtsIn(it2),1:nData));
    
    reconAll = reconHi + dataInLo.' + meanI;

end

