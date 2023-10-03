%% #################  cwtErrors  ################# %%
% 
% This file characterizes the residual error of the species and
% constitituents reconstructions through basic statistic (RMSE, bias) and
% power spectra.
% Note that this file is yet to be modularized.
% 

%%
% define which values to zero out for power spectrum
% there's a more elegant way of doing this that will be implemented in a
% future release

all_ind = 1:length(tPts);
overlap = ismember(all_ind,nanInd);
nan_change = find(abs(diff(overlap))==1)+1;

bad_window = ceil(wtCrit*max(coFiltLength));

nanEdges = [];
for k = 1:length(nan_change)
    nanEdges = [nanEdges,((nan_change(k)-bad_window):(nan_change(k)+bad_window))];
end

nanEdges(nanEdges>length(tPts)) = [];
nanEdges(nanEdges<1) = [];

% sloppy, but works for now
% bad_vals = int32(sort(unique([(1:bad_window),nanEdges,nanInd',(length(tPts)-bad_window:length(tPts))])));
bad_vals = int32(sort(unique([(1:bad_window),nanInd',(length(tPts)-bad_window:length(tPts))])));


%% designate NFFT for hi-res spectra
nfft = 2^12;

%% define power spectrum of original signal
c = coDataInHi .* wtsIn.';

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

input.pwr=pxx/3600/24;
input.freqs=f*3600*24;      % frequency in cpd

constits.input_pwr = input.pwr;
constits.input_freqs = input.freqs;

%% power spectra of reconstructed signals
if spOnlyBool==0
    % constits
    c = constits.reconHi .* wtsIn;

    d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

    [pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

    constits.recon_pwr=pxx/3600/24;
    constits.recon_freqs=f*3600*24;      % frequency in cpd
end

% species
c = species.reconHi .* wtsIn;

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

species.recon_pwr=pxx/3600/24;
species.recon_freqs=f*3600*24;      % frequency in cpd


%% define residual noise's power spectrum
hband = 0.2/(3600*24);         % frequency band of 0.2 cyc/day [cpd] in cyc/sec; half width

if spOnlyBool==0
    % for constits
    freq_bands = zeros(length(coOmegaFlat),2);

    for k=1:length(coOmegaFlat)
        freq_bands(k,1) = coOmegaFlat(k)/(2*pi) - hband;
        freq_bands(k,2) = coOmegaFlat(k)/(2*pi) + hband;
    end

    c = coDataInHi .* wtsIn.' - constits.reconHi.';
    % break pt
    % should converge to essentially zero mean
    for i=1:10
        c(bad_vals) = 0;
        c = c-mean(c);
    end

    d = floor(length(c)/2);

    [pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

    fband_var = zeros(length(coOmegaFlat),1);
    constits_pwr = zeros(length(coOmegaFlat),1);
    constits_SNR = zeros(length(coOmegaFlat),1);
    constits_modRC = zeros(length(coOmegaFlat),1);
    df_MH = zeros(length(coOmegaFlat),1);
    df_resolved = zeros(length(coOmegaFlat),1);
    factor_i = zeros(length(coOmegaFlat),1);
    fband_ci = zeros(length(coOmegaFlat),1);
    fband_ci_ph = zeros(length(coOmegaFlat),1);     % CI for phases
    df = f(2)-f(1);

    for k=1:length(coOmegaFlat)
        ind = find(f > freq_bands(k,1) & f < freq_bands(k,2));
        [~,freq_ind] = min(abs(f(ind)-coOmegaFlat(k)*((1)/(2*pi))));
        fband_pwr = pxx(ind);
        a = constits.recon_pwr(ind);
        constits_pwr(k) = a(freq_ind)/3600/24;
        fband_pwr(freq_ind) = 0;
        fband_var(k) = sum(df*fband_pwr);
        fband_ci(k) = sqrt(fband_var(k)) * 1.96;            % 95 % CI of normal dist
        A2B2 = mean(real(coSoln(k,:)).^2 + imag(coSoln(k,:)).^2);   % phase CI per David's comments
        fband_ci_ph(k) = (180/pi) * sqrt(fband_var(k)/A2B2) * 1.96;
        constits_SNR(k) = a(freq_ind)*df/(fband_var(k)/3600/24);
        constits_modRC(k) = (coFiltLength(k)*3600*sqrt(constits_SNR(k)))^-1;
    end

    % M-H criterion
    k=1;
    nT=0; % total number of constituents gone through at the start of each j-loop; needed to index flattened names from constits
    for i=1:length(constits.resp.names)
        nC = length(constits.resp.mats{i}(1,:));
        for j=1:nC
            factor_i_denom = 0;
            for l=1:nC
                factor_i_denom = factor_i_denom + (constits.resp.mats{i}(j,l)*nanmean(constits.(constits.names{nT+l}).amps))^2;
            end
            factor_i_denom = sqrt(factor_i_denom);
            factor_i(k) = (constits.resp.mats{i}(j,j)*nanmean(constits.(constits.names{k}).amps))/factor_i_denom;
            LOF = coFiltLength(k)/24;   % now in days
            df_MH(k)= 1/(LOF*sqrt(constits_SNR(k)));        
            df_resolved(k)=2*df_MH(k)*(1-factor_i(k));      % factor of 2 is bc we want half of LOF

            k=k+1;
        end
        nT = nT+nC;
    end
    
    constits.fband.bands = freq_bands*3600*24;
    constits.fband.pwrs = pxx/3600/24;
    constits.fband.freqs = f * 3600 * 24;
    constits.fband.var = fband_var;
    constits.fband.ci = fband_ci;
    constits.fband.ci_ph = fband_ci_ph;
    constits.resid = c;
    constits.fband.constits_pwr = constits_pwr;
    constits.fband.SNR = constits_SNR;
    constits.modRC = constits_modRC;
    constits.df_resolved = df_resolved;
    constits.df_MH = df_MH;
    constits.factor_i = factor_i;
    
end

% for species now
bad_window = ceil(wtCrit*max(spFiltLength));

nanEdges = [];
for k = 1:length(nan_change)
    nanEdges = [nanEdges,(nan_change(k)-bad_window:nan_change(k)+bad_window)];
end

nanEdges(nanEdges>length(tPts)) = [];
nanEdges(nanEdges<1) = [];

% bad_vals = int32(sort(unique([(1:bad_window),nanEdges,nanInd',(length(tPts)-bad_window:length(tPts))])));
bad_vals = int32(sort(unique([(1:bad_window),nanInd',(length(tPts)-bad_window:length(tPts))])));


freq_bands = zeros(length(spOmega),2);

for k=1:length(spOmega)
    freq_bands(k,1) = spOmega(k)/(2*pi) - hband;
    freq_bands(k,2) = spOmega(k)/(2*pi) + hband;
end

c = spDataInHi .* wtsIn.' - species.reconHi';

% should converge to essentially zero mean
for i=1:10
    c(bad_vals) = 0;
    c = c-mean(c);
end

d = floor(length(c)/2);

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/deltaT);         

fband_var = zeros(length(spOmega),1);
species_pwr = zeros(length(spOmega),1);
species_SNR = zeros(length(spOmega),1);
species_modRC = zeros(length(spOmega),1);
fband_ci = zeros(length(spOmega),1);

df = f(2)-f(1);

for k=1:length(spOmega)
    ind = find(f > freq_bands(k,1) & f < freq_bands(k,2));
    [~,freq_ind] = min(abs(f(ind)-spOmega(k)*((1)/(2*pi))));
    fband_pwr = pxx(ind);
    a = species.recon_pwr(ind);
    species_pwr(k) = a(freq_ind);
    fband_pwr(freq_ind) = 0;
    fband_var(k) = sum(df*fband_pwr);
    fband_ci(k) = sqrt(fband_var(k)) * 1.96;            % 95 % CI of normal dist 
    species_SNR(k) = a(freq_ind)*df/(fband_var(k)/3600/24);
    species_modRC(k) = (spFiltLength(k)*3600*sqrt(species_SNR(k)))^-1;
end


species.fband.bands = freq_bands*3600*24;
species.fband.pwrs = pxx/3600/24;
species.fband.freqs = f * 3600 * 24;
species.fband.var = fband_var;
species.fband.ci = fband_ci;
species.resid = c;
species.fband.constit_pwr = species_pwr;
species.fband.SNR = species_SNR;
species.modRC = species_modRC;

species.reconHi(nanInd)     = NaN;
species.reconAll(nanInd)    = NaN;

constits.reconHi(nanInd)     = NaN;
constits.reconAll(nanInd)    = NaN;

%% Calculcate RMSE for reconstructed wl
if statWindowBool==0
    statStart = 1;
    statEnd = length(constits.resid);
else
    [~,statStart] = min(abs(statWindow(1)-dtimeIn));
    [~,statEnd]   = min(abs(statWindow(2)-dtimeIn));
end

if spOnlyBool==0
    constits.rmse = sqrt(mean(constits.resid(statStart:statEnd).^2));
    constits.bias = mean(constits.resid(statStart:statEnd));
end

species.rmse = sqrt(mean(species.resid(statStart:statEnd).^2));
species.bias = mean(species.resid(statStart:statEnd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function for bootstrap approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Moving Block Bootstrap for species
% using notation from Innocenti et al., 2022

if bootstrapFlag==1
    
    bootstrapActive = 1;
    
    L=30;           % length in days of residual blocks, as in paper referenced above
    
    % define reconstructed series and residuals of interest
    N = length(species.reconAll);

    lblocks=ceil((L*24*3600)/deltaT);
    B = ceil(N/lblocks);
    Nt = N-1;
    
    % data structure for alternate realizations of observed data
    species_alt = zeros(length(species.reconAll),nBoot);
    constits_alt = zeros(length(constits.reconAll),nBoot);
    
    % residuals only for respective realizations
    species_boot_resid = zeros(length(species.reconAll),nBoot);
    constits_boot_resid = zeros(length(constits.reconAll),nBoot);
    
    for k=1:nBoot        
        i_b = int32(1 + (Nt-lblocks-1).*rand(B,1));

        species_resid_blocks = zeros(B*lblocks,1);
        constits_resid_blocks = zeros(B*lblocks,1);

        for j=1:B
            start_ind = (j-1)*lblocks+1;
            end_ind = start_ind+lblocks-1;
            species_resid_blocks(start_ind:end_ind) = (species.resid(i_b(j):(i_b(j)+lblocks-1)));
            constits_resid_blocks(start_ind:end_ind) = (constits.resid(i_b(j):(i_b(j)+lblocks-1)));
        end
        
        species_resid_b = species_resid_blocks(1:N);
        constits_resid_b = constits_resid_blocks(1:N);
        
        species_alt(:,k) = species.reconAll + species_resid_b.';
        species_boot_resid(:,k) = species_resid_b;
        
        constits_alt(:,k) = constits.reconAll + constits_resid_b.';
        constits_boot_resid(:,k) = constits_resid_b;
        
    end
    
    % these are your nBoot noise realizations
    species.boots = species_boot_resid;    
    constits.boots = constits_boot_resid;
    
    for k=1:nBoot
        
        % update data
        coDataIn = constits.reconAll.'+constits.boots(:,k);
        spDataIn = species.reconAll.'+species.boots(:,k);
        
        % general analysis
        cwtConvolveFilters;
        cwtSolve;

        % OUTPUTS
        n = 1;
        for j=species.names
            species.(j).bootamps(:,k) = spSolnAmps(n,:);
            % species.(j).phases = spSolnPhases(n,:);
            n=n+1;
        end

        constits.dectimes = coTimes;

        n = 1;
        for j=constits.names
            constits.(j).bootamps(:,k) = coSolnAmps(n,:);
            % constits.(j).phases = coSolnPhases(n,:);
            n=n+1;
        end
        
        disp("Bootstrap error round "+num2str(k)+" at " + datestr(now, 'dd/mm/yy-HH:MM'))
        
%         if k==25
%             save boot25.mat species constits
%         elseif k==50
%             save boot50.mat species constits
%         elseif k==100
%             save boot100.mat species constits
%         elseif k==250
%             save boot250.mat species constits
%         elseif k==500
%             save boot500.mat species constits
%         elseif k==750
%             save boot750.mat species constits
%         end
        
    end
    
    % we now have all bootstrap amplitudes
    % we want to create confidence intervals for each species/constituent
    
%     for k=species.names
%         
%         amps = species.(k).bootamps - species.(k).amps.';
%         
%         boots = reshape(amps,nBoot*length(amps),1);
%         
%         species.(k).bootci = 1.96*std(boots);
%         
%         species.(k).boots = boots;
%     end
% 
%     for k=constits.names
%         
%         amps = constits.(k).bootamps - constits.(k).amps.';
%         
%         boots = reshape(amps,nBoot*length(amps),1);
%         
%         constits.(k).bootci = 1.96*std(boots);
%         
%         constits.(k).boots = boots;
%     end
%     
    bootstrapActive = 0;
    
end



