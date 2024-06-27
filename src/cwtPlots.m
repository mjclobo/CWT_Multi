%% #################  cwtPlots  ################# %%
% 
% This file plots results from the analysis, all contingent upon true (1)
% plot flags, as input by the user.
% By default, the routine doesn't plot anything.
% 
%%

if pfResid == 1
    
    x = species.alltimes;
    
    p = figure();
    plot(x,spDataInHi .* wtsIn.')
    hold on
    plot(x,species.reconHi)
    hold off
    grid on
    legend('original','recon (species)')
    ylabel('Amplitude')
    xlabel('Time [days]')
    xlim([min(x) max(x)])
    title('highpassed data vs recon (species)')
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])
    
    p = figure();
    plot(x,species.resid)
    grid on
    ylabel('Residuals (species)')
    xlabel('Time [days]')
    xlim([min(x) max(x)])
    title('Residuals (species)')
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])
    
    x = constits.alltimes;
    
    p = figure();
    plot(x,coDataInHi .* wtsIn.')
    hold on
    plot(x,constits.reconHi)
    hold off
    grid on
    legend('original','recon (constits)')
    ylabel('Amplitude')
    xlabel('Time [days]')
    xlim([min(x) max(x)])
    title('highpassed data vs recon (constits)') 
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])
    
    p = figure();
    plot(x,constits.resid)
    grid on
    ylabel('Residuals (constits)')
    xlabel('Time [days]')
    xlim([min(x) max(x)])
    title('Residuals (constits)')
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])

end

if pfAmps==1
    
    x = species.decTimesAll;
    p = figure();
    cols = hsv(length(species.names));
    leg_ind = zeros(length(species.names),1);
    childs = zeros(length(species.names)*3,1);
    l=1;
    i=1;
    hold on
    for k=species.names
        leg_ind(i) = l;
        childs(l) = plot(x,species.(k).amps,'color',cols(i,:),'LineWidth',2);
        childs(l+1) = plot(x,species.(k).amps-species.fband.ci(i),'--','color',cols(i,:),'LineWidth',0.5);
        childs(l+2) = plot(x,species.(k).amps+species.fband.ci(i),'--','color',cols(i,:),'LineWidth',0.5);
        i=i+1;
        l=l+3;
    end
%     if bootstrapFlag==1
%         leg_ind2 = zeros(length(species.names),1);
%         childs2 = zeros(length(species.names)*3,1);
%         l=1;
%         i=1;
%         for k=species.names
%             leg_ind2(i) = l;
%             childs2(l+1) = plot(x,species.(k).amps-species.(k).bootci,':','color',cols(i,:),'LineWidth',0.5);
%             childs2(l+2) = plot(x,species.(k).amps+species.(k).bootci,':','color',cols(i,:),'LineWidth',0.5);
%             i=i+1;
%             l=l+3;
%         end
%     end
    hold off
    grid on
    legend([childs(leg_ind)],species.names)
    grid off
    ylabel('Amplitude')
    xlabel('Time')    
    title('Species amplitudes')
    xlim([min(x) max(x)])
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])
    
    if spOnlyBool==0
        x = constits.decTimesAll;

        cols = hsv(length(constits.names));
        p = figure();
        leg_ind = zeros(length(constits.names),1);
        childs = zeros(length(constits.names)*3,1);
        l=1;
        i=1;
        hold on
        for k=constits.names
            leg_ind(i) = l;
            childs(l) = plot(x,constits.(k).amps,'color',cols(i,:),'LineWidth',2);
            childs(l+1) = plot(x,constits.(k).amps-constits.fband.ci(i),'--','color',cols(i,:),'LineWidth',0.5);
            childs(l+2) = plot(x,constits.(k).amps+constits.fband.ci(i),'--','color',cols(i,:),'LineWidth',0.5);
            i=i+1;
            l=l+3;
        end
    %     if bootstrapFlag==1
    %         leg_ind2 = zeros(length(constits.names),1);
    %         childs2 = zeros(length(constits.names)*3,1);
    %         l=1;
    %         i=1;
    %         for k=constits.names
    %             leg_ind2(i) = l;
    %             childs2(l+1) = plot(x,constits.(k).amps-constits.(k).bootci,':','color',cols(i,:),'LineWidth',0.5);
    %             childs2(l+2) = plot(x,constits.(k).amps+constits.(k).bootci,':','color',cols(i,:),'LineWidth',0.5);
    %             i=i+1;
    %             l=l+3;
    %         end
    %     end
        hold off
        grid on
        legend([childs(leg_ind)],constits.names)
        grid off
        ylabel('Amplitude')
        xlabel('Time')    
        title('Constits amplitudes')
        xlim([min(x) max(x)])
        ax = gca;
        ax.FontSize = 14;
        set(gcf,'Position',[100 100 1500 600])
    end
    
end


if pfFilt == 1
    
    p = figure();
    filt = spNFilt{4};
    t_half = (length(filt)-1)/2;
    x = (-t_half:t_half)*deltaT;
    plot(x,real(filt),'linewidth',2)
    hold on
    plot(x,imag(filt),'linewidth',2)
    hold off
    grid on
    title('Normalized D_2 filter (species)')
    xlabel('Time [s]')
    xlim([min(x) max(x)])
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])
    
    if spOnlyBool==0
        p = figure();
        filt = coNFiltD2{1};
        t_half = (length(filt)-1)/2;
        x = (-t_half:t_half)*deltaT;
        plot(x,real(filt),'linewidth',2)
        hold on
        plot(x,imag(filt),'linewidth',2)
        hold off
        grid on
        title('Normalized M_2 filter (constits)')
        xlabel('Time [s]')
        xlim([min(x) max(x)])   
        ax = gca;
        ax.FontSize = 14;
        set(gcf,'Position',[100 100 1500 600])
    end
    
end


if pfResp == 1
    
    N = length(coOmegaFlat);
    
    N_om = 501;                             % must be odd
    
    responses = zeros(N,N_om);
    freqs = zeros(N,N_om);
    
    for k=1:N
        [responses(k,:),freqs(k,:)] = cwt_freq_response(coFilt{k},coOmegaFlat(k),abs(coOmegaFlat(k)/500),N_om,deltaT);
    end
    
    p = figure();
    hold on
    for k=1:N
        plot(((freqs(k,:)*24*3600)/(2*pi)),responses(k,:),'linewidth',2)
    end

    nfft = 2^12;
    
%     for k=constits.names
%         [constits.(k).respMag,constits.(k).respFreq] = freq_resp(constits.(k).filter,1/deltaT,nfft);
%         plot(constits.(k).respFreq*3600*24,constits.(k).respMag,'linewidth',2)
%     end
    
    % set (gca,'xdir','reverse')
    hold off
    grid on
    xlabel('Freq [cy day^{-1}]')
    ylabel('Response')
    ylim([0 1])
    grid on
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])
    
end

if spOnlyBool==0
    if pfM4M2 == 1

        x = constits.decTimesAll;

        p = figure();
        semilogy(x,(constits.M4.amps ./ constits.M2.amps).^2)
        hold on
        plot(x,constits.M2.amps)
        plot(x,constits.M4.amps)
        hold off
        grid on
        legend('(M_4 / M_2)^2','M_2','M_4')
        ylabel('Amplitude')
        xlabel('Time [days]')
        xlim([min(x) max(x)])
        title('semilogy of (M_4/M_2)^2')
        ax = gca;
        ax.FontSize = 14;
        set(gcf,'Position',[100 100 1500 600])

    end
end

if pfResidSpectra == 1
    
    p1 = figure();
    loglog(species.fband.freqs.^-1,species.fband.pwrs,'linewidth',2)
    hold on
    plot(input.freqs.^-1,input.pwr,'r-','linewidth',1)
    plot(species.recon_freqs.^-1,species.recon_pwr,'b-','linewidth',0.5)
    hold off
    grid on
    yma = max(species.fband.pwrs)*1.2;
    ymi = min(species.fband.pwrs);
    set(gca,'xdir','reverse')
    xlabel('Period [days]')
    ylabel('Residual PSD [m^2 cpd^{-1}]')
    legend('residuals','input data','reconstructed data','Location','southwest')
    axis tight
    title('Species residual spectrum')
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])

    if spOnlyBool==0
        p2 = figure();
        loglog(constits.fband.freqs.^-1,constits.fband.pwrs,'linewidth',2)
        hold on
        plot(input.freqs.^-1,input.pwr,'r-','linewidth',1)
        plot(constits.recon_freqs.^-1,constits.recon_pwr,'b-','linewidth',0.5)
        hold off
        grid on
        yma = max(constits.fband.pwrs)*1.2;
        ymi = min(constits.fband.pwrs);
        set(gca,'xdir','reverse')
        xlabel('Period [days]')
        ylabel('Residual PSD [m^2 cpd^{-1}]')
        legend('residuals','input data','reconstructed data','Location','southwest')
        axis tight
        title('Constits residual spectrum')
        ax = gca;
        ax.FontSize = 14;
        set(gcf,'Position',[100 100 1500 600])
    end
end

if pfLoD2 == 1
    
    x2 = species.alltimes;
    
    p = figure();
    plot(species.dectimes,species.d2.amps,'LineWidth',2)
    hold on
    yyaxis right
    plot(x2,spDataInLo,'LineWidth',2)
    hold off
    legend('species D_2 (L axis)','low-passed data (R axis)')
    grid on
    ylabel('Amplitude')
    xlabel('Time [days]')    
    title('Species D_2 Amplitude vs. Low-passed Data')
    xlim auto
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])
    
    if spOnlyBool==0
        p = figure();
        cols = hsv(length(constits.names(coND1+1:coND1+coND2)));
        childs = zeros(length(species.names(coND1+1:coND1+coND2))+1,1);
        i=1;
        hold on
        for k=constits.names(coND1+1:coND1+coND2)
            childs(i) = plot(constits.dectimes,constits.(k).amps,'color',cols(i,:),'LineWidth',2);
            i=i+1;
        end
        ylabel('Amplitude')
        yyaxis right
        childs(length(species.names(coND1+1:coND1+coND2))+1) = plot(constits.alltimes,coDataInLo,'LineWidth',2);
        hold off    
        legend(childs,[constits.names(coND1+1:coND1+coND2)+' (L axis)','Low-passed data (R axis)'])
        grid on
        ylabel('Water level')
        xlabel('Time [days]')    
        title('Semidiurnal Constits Amps vs. Low-passed Data')
        xlim auto
        ax = gca;
        ax.FontSize = 14;
        set(gcf,'Position',[100 100 1500 600])
    end
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotvert(xin,N,ymin,ymax,col,i)
    if xin>0
        x = ones(N+1,1)*xin;
    else 
        x = ones(N+1,1)*10^-4;
    end
    dy = (ymax-ymin)/N;
    y = (ymin:dy:ymax)';
    plot(x,y,'--','color',col(i,:),'linewidth',2)
end

function [oneside_pwr,oneside_freq] = freq_resp(sig,fs,nfft)

Nzeros = nfft-length(sig);

sig_padded = horzcat(sig,zeros(1,Nzeros));

pwr_all = flip(fft(sig_padded));

oneside_pwr = 2*abs(pwr_all(1:(nfft/2)));

oneside_freq = fs*(0:(nfft/2-1))/nfft;

end

function [resp,om] = cwt_freq_response(filt,cent_freq,delta_om,N_om,deltaT)
% find frequency reponse of a filter using sine-wave analysis
%

N = length(filt);
N2 = (length(filt)-1)/2;

t = (-N2:N2) * deltaT;

N_om2 = (N_om - 1)/2;
om = (cent_freq - delta_om*N_om2:delta_om:cent_freq+delta_om*N_om2);

resp = zeros(length(om),1);

for k=1:N_om
    c = cos(om(k) * t);
    
    % for NConvolved() data must be same dim
    resp(k) = real(NConvolved([],c,filt,ones(size(filt)),1,0.5));
end

end
