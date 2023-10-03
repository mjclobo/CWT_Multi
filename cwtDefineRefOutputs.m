%% #################  cwtDefineRefOutputs  ################# %%
% 
% This file defines the ref.constits and ref.species structs which will be output
% through the main cwtMulti routine.

%%
if spOnlyBool==0
    ref.constits.reconAll   = coReconAll;
    ref.constits.reconHi    = coReconHi;
    ref.constits.dectimes   = coTimes;
end

% species stuff now
ref.species.reconAll    = spReconAll;
ref.species.reconHi     = spReconHi;
ref.species.dectimes    = spTimes;

n = 1;
for k=ref.species.names
    ref.species.(k).amps    = spSolnAmps(n,:);
    ref.species.(k).phases  = spSolnPhases(n,:);
    n=n+1;
end


if spOnlyBool==0
    n = 1;
    for k=ref.constits.names
        ref.constits.(k).amps   = coSolnAmps(n,:);
        ref.constits.(k).phases = coSolnPhases(n,:);
        n=n+1;
    end
end

if refStnFlag==0
    ref.potential = poten;
end




