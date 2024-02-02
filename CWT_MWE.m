%% Basic example of CWT_Multi usage

%% load data
load('./data/Astoria/astoria_wl.mat')
load('./data/Vancouver/vancouver_wl.mat')

% align station data
startDate = datetime(2010,06,03);
endDate = datetime(2012,06,03);

Van.start = find(Van.dates==startDate);
Van.end = find(Van.dates==endDate);

Van.dates = Van.dates(Van.start:Van.end);
Van.wl = Van.wl(Van.start:Van.end);

Ast.start = find(Ast.dates==startDate);
Ast.end = find(Ast.dates==endDate);

Ast.dates = Ast.dates(Ast.start:Ast.end);
Ast.wl = Ast.wl(Ast.start:Ast.end);

% Overlap all missing data between two stations
nan_ind = unique(sort([find(isnan(Van.wl));find(isnan(Ast.wl))]));

Van.wl(nan_ind) = NaN;
Ast.wl(nan_ind) = NaN;

% lon, lat for Astoria
lon = 123+46.1/60;
lat = 46+12.4/60;

% set time zone
Van.dates.TimeZone = 'America/Los_Angeles';
Ast.dates.TimeZone = 'America/Los_Angeles';

%% run cwt routine on Vancouver/Astoria data, using Astoria as reference station
[constits,species,ref,admit] = cwtMulti(Van.dates,Van.wl,'performAdmittance',lon,lat,'refStation',Ast.wl); 

%% PLOT SEMIDIURNAL AMPLITUDES AND ADMITTANCES
figure(); p=tiledlayout(3,1);

ax0=nexttile;
plot(constits.decTimesAll,constits.M2.amps,'linewidth',2,'DisplayName','M_2')
hold(ax0,'on')
plot(constits.decTimesAll,constits.S2.amps,'linewidth',2,'DisplayName','S_2')
plot(constits.decTimesAll,constits.N2.amps,'linewidth',2,'DisplayName','N_2')
hold(ax0,'off')
grid on
legend('location','northeast')
title('Vancouver D_2', 'Units', 'normalized', 'Position', [0.5, 0.85, 0]) 
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;

ax1=nexttile;
plot(constits.decTimesAll,ref.constits.M2.amps,'linewidth',2,'DisplayName','M_2')
hold(ax1,'on')
plot(constits.decTimesAll,ref.constits.S2.amps,'linewidth',2,'DisplayName','S_2')
plot(constits.decTimesAll,ref.constits.N2.amps,'linewidth',2,'DisplayName','N_2')
hold(ax1,'off')
grid on
legend('location','northeast')
title('Astoria D_2', 'Units', 'normalized', 'Position', [0.5, 0.85, 0]) 
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;

ax2=nexttile;
plot(constits.decTimesAll,admit.constits.M2.amps,'linewidth',2,'DisplayName','M_2')
hold(ax2,'on')
plot(constits.decTimesAll,admit.constits.S2.amps,'linewidth',2,'DisplayName','S_2')
plot(constits.decTimesAll,admit.constits.N2.amps,'linewidth',2,'DisplayName','N_2')
hold(ax2,'off')
grid on
legend('location','northeast')
title('Vancouver/Astoria D_2admittances', 'Units', 'normalized', 'Position', [0.5, 0.85, 0]) 
set(gca,'xtick',[])
xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
ax = gca;
ax.FontSize = 14;

% xlabel(p,'Time','FontSize',18)
ylabel(p,'Amplitude (m)','FontSize',18)
set(gcf,'Position',[100 100 1600 500])

