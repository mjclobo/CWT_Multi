%% #################  cwtDefineOutputs  ################# %%
% 
% This file defines the constits and species structs which will be output
% through the main cwtMulti routine.
% 

%%
if spOnlyBool==0
    constits = assignResults(constits,coReconAll,coReconHi,coFiltLength,coTimes,tPts,coTimesAll,coDataInLo,coDataInHi,...
        coSoln,coSolnAmps,coSolnPhases,coFilt,centerDate,{"D1","D2","D3","D4"},{coRespD1,coRespD2,coRespD3,coRespD4},coOmegaFlat);
end

species = assignResults(species,spReconAll,spReconHi,spFiltLength,spTimes,tPts,spTimesAll,spDataInLo,spDataInHi,spSoln,...
    spSolnAmps,spSolnPhases,spFilt,centerDate,"Sp",spResp,spOmega);

if doMonthly==1
    monthly = assignResults(monthly,[],[],moFiltLength,moTimes,tPts,moTimesAll,dataIn,[],moSoln,moSolnAmps,moSolnPhases,moFilt,centerDate,"Mo",moResp);
end

% include lower-frequency filter results
if spOnlyBool==0
    if doMonthly==1
        constits.lowerFreq = monthly;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function structIn = assignResults(structIn,reconAll,reconHi,filtLength,times,tPts,timesAll,dataInLo,dataInHi,soln,solnAmps,solnPhases,filt,centerDate,resp_names,resp,omega)

    structIn.reconstruction.reconAll = reconAll;
    structIn.reconstruction.reconHi  = reconHi;
    
    structIn.resp.names = resp_names;
    structIn.resp.mats = resp;

    structIn.input.filter_length = filtLength;

%     structIn.dectimes = datetime(times,'ConvertFrom','epochtime','Epoch',centerDate);
    structIn.alltimes = datetime(tPts,'ConvertFrom','epochtime','Epoch',centerDate);
    structIn.decTimesAll = datetime(timesAll,'ConvertFrom','epochtime','Epoch',centerDate);

    structIn.input.dataIn_lo = dataInLo;
    structIn.input.dataIn_hi = dataInHi;
    
    structIn.input.omegas = omega;

    % 
    n = 1;
    for k=structIn.names
        structIn.(k).compSol = soln(n,:);
        structIn.(k).amps    = solnAmps(n,:);
        structIn.(k).phases  = solnPhases(n,:);
        structIn.(k).filter  = filt{n};
        
%         nanInd = find(solnAmps(n,:)==0.0);   % this is where NaNs should be for decimated data
%         structIn.(k).compSol(nanInd) = NaN;
%         structIn.(k).amps(nanInd)    = NaN;
%         structIn.(k).phases(nanInd)  = NaN;
        n=n+1;
    end

end

