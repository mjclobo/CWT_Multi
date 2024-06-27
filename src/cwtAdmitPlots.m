% various plots, called using various args

if plot_flag_amps ==1
    
    x = species.dectimes;
    p = figure();
    cols = hsv(length(species.names));
    leg_ind = zeros(length(species.names),1);
    childs = zeros(length(species.names)*3,1);
    l=1;
    i=1;
    hold on
    for k=species.names
        leg_ind(i) = l;
        childs(l) = plot(x,admit.species.(k).amps,'color',cols(i,:),'LineWidth',2);
        childs(l+1) = plot(x,admit.species.(k).amps-admit.species.fband.ci(i),'--','color',cols(i,:),'LineWidth',2);
        childs(l+2) = plot(x,admit.species.(k).amps+admit.species.fband.ci(i),'--','color',cols(i,:),'LineWidth',2);
        i=i+1;
        l=l+3;
    end
    if bootstrap_flag==1
        leg_ind2 = zeros(length(species.names),1);
        childs2 = zeros(length(species.names)*3,1);
        l=1;
        i=1;
        for k=species.names
            leg_ind2(i) = l;
            childs2(l+1) = plot(x,admit.species.(k).amps-admit.species.(k).bootci,':','color',cols(i,:),'LineWidth',2);
            childs2(l+2) = plot(x,admit.species.(k).amps+admit.species.(k).bootci,':','color',cols(i,:),'LineWidth',2);
            i=i+1;
            l=l+3;
        end
    end
    hold off    
    legend([childs(leg_ind)],species.names)
    grid off
    ylabel('Amplitude')
    xlabel('Time')    
    title('Species amplitudes')
    xlim([min(x) max(x)])
    ax = gca;
    ax.FontSize = 14;
    set(gcf,'Position',[100 100 1500 600])
    
    
    x = constits.dectimes;
    
    cols = hsv(length(constits.names));
    p = figure();
    leg_ind = zeros(length(constits.names),1);
    childs = zeros(length(constits.names)*3,1);
    l=1;
    i=1;
    hold on
    for k=constits.names
        leg_ind(i) = l;
        childs(l) = plot(x,admit.constits.(k).amps,'color',cols(i,:),'LineWidth',2);
        childs(l+1) = plot(x,admit.constits.(k).amps-admit.constits.fband.ci(i),'--','color',cols(i,:),'LineWidth',2);
        childs(l+2) = plot(x,admit.constits.(k).amps+admit.constits.fband.ci(i),'--','color',cols(i,:),'LineWidth',2);
        i=i+1;
        l=l+3;
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

function [resp,om] = cwt_freq_response(filt,cent_freq,delta_om,N_om,deltat)
% find frequency reponse of a filter using sine-wave analysis
%

N = length(filt);
N2 = (length(filt)-1)/2;

t = (-N2:N2) * deltat;

N_om2 = (N_om - 1)/2;
om = (cent_freq - delta_om*N_om2:delta_om:cent_freq+delta_om*N_om2);

resp = zeros(length(om),1);

for k=1:N_om
    c = cos(om(k) * t);
    
    % for NConvolved() data must be same dim
    resp(k) = real(NConvolved(c,filt,ones(size(filt)),1,0.5));
end

end
