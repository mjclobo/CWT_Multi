%% #################  cwtDynamicInference  ################# %%
% 
% This script uses the dynamic inference method to...

if dynInfFlag==1
    
    S2ind = find(constits.names=="S2");
    K1ind = find(constits.names=="K1");
    
    [dynInf,offsetS2,nDiffS2] = dynamicInf(dynInf.namesS2,1,diSoln,diSolnAmps,diSolnPhases,diFilt,constits,dynInf,decFact,diMaxFiltLength,coMaxFiltLength,nData,coFilt{S2ind},diOmega,deltaT);
    [dynInf,offsetK1,nDiffK1] = dynamicInf(dynInf.namesK1,4,diSoln,diSolnAmps,diSolnPhases,diFilt,constits,dynInf,decFact,diMaxFiltLength,coMaxFiltLength,nData,coFilt{K1ind},diOmega,deltaT);
    
    constits.dynamicInference = dynInf;
    
    %% plot results
    p = figure;
    if nDiffK1>0
        x=constits.decTimesAll(1+offsetK1:end-offsetK1);
    else
        x=constits.decTimesAll;
    end
    
    plot(x,dynInf.K2.shortAmps,'k-','linewidth',2.0,'DisplayName','K2 inferred from S2')
    hold on
    plot(x,dynInf.S2.shortAmps,'b-','linewidth',2.0,'DisplayName','S2 inferred from S2')
    plot(x,dynInf.T2.shortAmps,'r-','linewidth',2.0,'DisplayName','T2 inferred from S2')
    plot(x,dynInf.P1.shortAmps,'k:','linewidth',2.0,'DisplayName','P1 inferred from K1')
    plot(x,dynInf.K1.shortAmps,'b:','linewidth',2.0,'DisplayName','K1 inferred from K1')
    plot(x,dynInf.S1.shortAmps,'r:','linewidth',2.0,'DisplayName','S1 inferred from K1')
    grid on
    legend('location','northeast')
    ylabel('Tidal Amplitude')
    xlabel('Time [days]')
    xlim([min(x) max(x)])
    title('Dynamically inferred amplitudes')
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])
    
    p = figure;
    if nDiffK1>0
        x=constits.decTimesAll(1+offsetK1:end-offsetK1);
        plot(x,dynInf.K2.sixMoAmp,'k-','linewidth',2.0,'DisplayName','K2 6mo amp')
        hold on
        plot(x,dynInf.S2.sixMoAmp,'b-','linewidth',2.0,'DisplayName','S2 6mo amp')
        plot(x,dynInf.T2.sixMoAmp,'r-','linewidth',2.0,'DisplayName','T2 6mo amp')
        plot(x,dynInf.P1.sixMoAmp,'k:','linewidth',2.0,'DisplayName','P1 6mo amp')
        plot(x,dynInf.K1.sixMoAmp,'b:','linewidth',2.0,'DisplayName','K1 6mo amp')
        plot(x,dynInf.S1.sixMoAmp,'r:','linewidth',2.0,'DisplayName','S1 6mo amp')
    elseif nDiffK1==0
        x=constits.decTimesAll;
        plot(x,dynInf.K2.sixMoAmp,'k-','linewidth',2.0,'DisplayName','K2 6mo amp')
        hold on
        plot(x,dynInf.S2.sixMoAmp,'b-','linewidth',2.0,'DisplayName','S2 6mo amp')
        plot(x,dynInf.T2.sixMoAmp,'r-','linewidth',2.0,'DisplayName','T2 6mo amp')
        plot(x,dynInf.P1.sixMoAmp,'k:','linewidth',2.0,'DisplayName','P1 6mo amp')
        plot(x,dynInf.K1.sixMoAmp,'b:','linewidth',2.0,'DisplayName','K1 6mo amp')
        plot(x,dynInf.S1.sixMoAmp,'r:','linewidth',2.0,'DisplayName','S1 6mo amp')
    elseif nDiffK1<0
        x=constits.decTimesAll;
        plot(x,dynInf.K2.sixMoAmp(1+offsetK1:end-offsetK1),'k-','linewidth',2.0,'DisplayName','K2 6mo amp')
        hold on
        plot(x,dynInf.S2.sixMoAmp(1+offsetK1:end-offsetK1),'b-','linewidth',2.0,'DisplayName','S2 6mo amp')
        plot(x,dynInf.T2.sixMoAmp(1+offsetK1:end-offsetK1),'r-','linewidth',2.0,'DisplayName','T2 6mo amp')
        plot(x,dynInf.P1.sixMoAmp(1+offsetK1:end-offsetK1),'k:','linewidth',2.0,'DisplayName','P1 6mo amp')
        plot(x,dynInf.K1.sixMoAmp(1+offsetK1:end-offsetK1),'b:','linewidth',2.0,'DisplayName','K1 6mo amp')
        plot(x,dynInf.S1.sixMoAmp(1+offsetK1:end-offsetK1),'r:','linewidth',2.0,'DisplayName','S1 6mo amp')
    end
    hold off
    grid on
    legend('location','northeast')
    ylabel('Tidal Amplitude')
    xlabel('Time [days]')
    xlim([min(x) max(x)])
    title('Six-month filter amplitudes')
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])

end

%% custom dynamic inference
if dynInfCustFlag==1
    baseInd = find(constits.names==dynInf.names(1));
    [dynInf,offsetBase,nDiffBase] = dynamicInf(dynInf.names,1,diSoln,diSolnAmps,diSolnPhases,diFilt,constits,dynInf,decFact,diMaxFiltLength,coMaxFiltLength,nData,coFilt{baseInd},diOmega,deltaT);    
    constits.dynamicInferenceCust = dynInf;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dynInf,offset,nDiff] = dynamicInf(namesIn,nStart,soln,amps,phases,filt,constitsIn,dynInf,decFact,maxFLLong,maxFLShort,nData,shortFilt,freqs,deltaT)
    % define frequency response coefficients for short filter
    resp = freqResp(shortFilt,freqs(nStart:(nStart+length(namesIn)-1)),deltaT);
    
    % normalize (i.e., we need to partition energy from short filter results with these coefficients)
    % resp = resp/sum(resp);
    
    % perform dynamic inference analysis for one constituent group
    n = nStart;
    
    % define basic DI variables
    dynInf.(namesIn(1)).groupResponse6mo = zeros(size(amps(1,:))).';
    for k=namesIn
        dynInf.(k).sixMoCompSol = extendDI(soln(n,:),decFact,maxFLLong,nData);
        dynInf.(k).sixMoAmp     = amps(n,:);
        dynInf.(k).sixMoPhases  = phases(n,:);
        dynInf.(k).filter       = filt{n};
        
        % defining group response at six month level
        dynInf.(namesIn(1)).groupResponse6mo = dynInf.(namesIn(1)).groupResponse6mo + dynInf.(k).sixMoCompSol;
        
        n=n+1;
    end
    
    for i=1:length(namesIn)
        inferredCon = namesIn(i);
        dynInf.(inferredCon+'toGroup').ratio = dynInf.(inferredCon).sixMoCompSol ./ dynInf.(namesIn(1)).groupResponse6mo;
    end
    
    % extend dominant constituent to time edges
    domShortCompSol = extendDI(constitsIn.(namesIn(1)).compSol,decFact,maxFLShort,nData);
    
    % Use ratios to define time series of constituents on shorter timescales
    for k=1:length(namesIn)
        domL = length(domShortCompSol);
        inferredCon = namesIn(k);
        infL = length(dynInf.(inferredCon+'toGroup').ratio);
        
        nDiff = domL-infL;
        
        inferredCon = namesIn(k);

        if nDiff==0
            offset=0;
            ratioSoln = dynInf.(inferredCon+'toGroup').ratio;
            dynInf.(inferredCon).shortCompSol = resp(k)*ratioSoln.*domShortCompSol;
            dynInf.(inferredCon).shortAmps = abs(dynInf.(inferredCon).shortCompSol);
%             dynInf.(inferredCon).phases = ;
        elseif nDiff>0
            offset = nDiff/2; 
            ratioSoln = dynInf.(inferredCon+'toGroup').ratio;
            dynInf.(inferredCon).shortCompSol = resp(k)*ratioSoln.*domShortCompSol(1+offset:domL-offset);
            dynInf.(inferredCon).shortAmps = abs(dynInf.(inferredCon).shortCompSol);
%             dynInf.(inferredCon).phases = ;
        elseif nDiff<0
            offset = -1*nDiff/2;
            ratioSoln = dynInf.(inferredCon+'toGroup').ratio(1+offset:infL-offset);
            dynInf.(inferredCon).shortCompSol = resp(k)*ratioSoln.*domShortCompSol;
            dynInf.(inferredCon).shortAmps = abs(dynInf.(inferredCon).shortCompSol);
%             dynInf.(inferredCon).phases = ;
        end
        
    end
end


function extSoln = extendDI(diSoln,decFact,maxFL,nData)

    extSoln = [diSoln(1) * ones(length(1:decFact:floor(maxFL/2)-1),1) ; diSoln.';...
        diSoln(end) * ones(length(nData-floor((maxFL/2))+1:decFact:nData),1)];

end

function resp = freqResp(filt,freqs,deltaT)
% defines frequency response matrix between all frequencies and filters

    Nfreq = length(freqs);
    resp = zeros(Nfreq,1);

    for k=1:Nfreq
        N2 = (length(filt)-1)/2;
        t = (-N2:N2) * deltaT;
        c = cos(freqs(k) * t);

        % for NConvolved() data must be same dim
        resp(k) = real(NConvolved([],c,filt,ones(size(filt)),1,0.5));
    end

    
end