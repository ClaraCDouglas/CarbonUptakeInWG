%%
set(0,'units','pixels')
Pix_SS=get(0,'screensize')
% 1           1        1536         864


shelfcolor=[0.9290 0.6940 0.1250];
opencolor=[0.5 0.3 0.8];

%% NPP

% Monthly timeseries
aix=4;
rix=1;
TotNPP_monthly=figure(1);
NPPmonth=plot(timedec,OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC);
ylabel('NPP (Tg C)');
title('Total NPP per month, Weddell Gyre');
set(gca,'XLim',[2002.4 2020.4], 'FontSize', 12);
xticks(2003:1:2020);

% add other algorithms
fig.TotNPPalgo=figure(2);
for aix=1:length(algorithm)
rix=1;
plot(timedec,OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC,'LineWidth',1.5);
hold on
end
fig.NPPlgd=legend(algorithm);
ylabel('NPP (Tg C)');
title('Total NPP per month, Weddell Gyre');
set(gca,'XLim',[2002.4 2020.4], 'FontSize', 12);
xticks(2003:1:2020);

% Annual total
aix=4;
x=2003:1:2019
y=OceanProd.(algorithm{aix}).(region_sublist{1}).NPP_tot_TgC_annual(2:18,2);

shelf_open_annual(:,1)=OceanProd.(algorithm{aix}).(region_sublist{2}).NPP_tot_TgC_annual(2:18,2);
shelf_open_annual(:,2)=OceanProd.(algorithm{aix}).(region_sublist{3}).NPP_tot_TgC_annual(2:18,2);
y2=shelf_open_annual;

totalNPPplot=figure;
subplot(2,1,1);
totalNPPplot.Position=[150 80 800 700];
 bar(x,y,'FaceColor',[0.1 0.6 0.5])
  hold on
  xtickangle(45)
  set(gca, 'XTick', [2003:1:2018])
 ylabel('Annual NPP (Tg C)')
title('Annual NPP in Weddell Gyre')

subplot(2,1,2);
b=bar(x,y2,'stacked');
  hold on
b(1).FaceColor=[0.9290 0.6940 0.1250];
b(2).FaceColor=[0.5 0.3 0.8];
  xtickangle(45)
  set(gca, 'XTick', [2003:1:2018])
  ylim([0 150])
 ylabel('Annual NPP (Tg C)')
title('Annual NPP on Shelf and in Open Ocean')
legend('Shelf Region','Open Ocean', 'Position',[0.5 0.33 0.1 0.1])


%% Rates
% daily rates per month - timeseries
aix=4;
rix=1;
fig.NPP_mgm2day=figure(3);
fig.mgm2day=plot(timedec,OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2);
ylabel('NPP (mg m^-^2 day^-^1)');
title('Average daily NPP per month, Weddell Gyre');
set(gca,'XLim',[2002.4 2020.4], 'FontSize', 12);
xticks(2003:1:2020);

% add other algorithms
fig.rateNPPalgo=figure(4);
for aix=1:length(algorithm)
rix=1;
plot(timedec,OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2,'LineWidth',1.5);
hold on
end
fig.NPPlgd=legend(algorithm);
xtickangle(45)
ylabel('NPP (mg m^-^2 day^-^1)');
title('Average daily NPP per month, Weddell Gyre');
set(gca,'XLim',[2002.4 2020.4], 'FontSize', 12);
xticks(2003:1:2020);



graphical.shelfcolor=[0.8 0.4 0];
graphical.opencolor=[0.4 0.2 0.8];


    %rates in regions
    aix=4;
    ratesplot=figure;
    ratesplot.Position=[150 100 1000 400];
    p(1)=plot(timedec,OceanProd.(algorithm{aix}).(region_sublist{2}).NPP_av_mgm2);
    hold on
    p(2)=plot(timedec,OceanProd.(algorithm{aix}).(region_sublist{3}).NPP_av_mgm2);
    p(1).Color= graphical.shelfcolor;%[0.9290 0.6940 0.1250];
    p(2).Color= graphical.opencolor%[0.5 0.3 0.8];
    p(1).LineWidth=1.55;
    p(2).LineWidth=1.25;
    ylabel('NPP (mg m^-^2 day^-^1)');
    regionlgd=legend('Shelf', 'Open Ocean');
    title('Average daily NPP per month');
    set(gca,'XLim',[2002.4 2020.4], 'FontSize', 12);
    xticks(2003:2:2020);
    ratesplot.Color='none';
    
% annual rates
clear andailyplot
years=[2002:1:2020]';
aix=4;
rix=1;
andailyplot=figure;
andailyplot.Position=[150 100 1000 400];
yyaxis left
% d(1)=plot(years,OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(:,2));
% hold on
d(2)=plot(years,OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(:,2));
hold on
d(4)=plot(years,OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(:,2));
yyaxis right
d(3)=plot(years,OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(:,2));
hold on
d(1).Color=[0.9290 0.6940 0.1250];
d(2).Color=[0.5 0.3 0.8];
d(1).LineWidth=1.55;
d(2).LineWidth=1.25;
ylabel('NPP (mg m^-^2 day^-^1)');
regionlgd=legend('Shelf', 'Open Ocean');
title('Average daily NPP per month');
set(gca,'XLim',[2002.4 2020.4], 'FontSize', 12);
xticks(2003:2:2020);
ratesplot.Color='none';


%% Open Water Area
% area timeseries (technically the max open water area per month)
aix=4;
rix=1;
OW_monthly=figure(10);
OWmonth=plot(timedec,OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_km2);
ylabel('NPP (km^2)');
title('Open (Ice-Free) Water Area per month, Weddell Gyre');
set(gca,'XLim',[2002.4 2020.4], 'FontSize', 12);
xticks(2003:1:2020);

% average annual open water area
clear avannualIceFreeplot
years=[2002:1:2020]';
aix=4;
rix=1;
avannualIceFreeplot=figure;
avannualIceFreeplot.Position=[150 100 1000 400];
d(1)=plot(years,OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(:,2));
hold on
ylabel('Open Water Area (km^2)');
set(gca,'YLim', [0 1.9e6], 'FontSize',12);
title('Mean open (ice-free) water area per year');
set(gca,'XLim',[2002.4 2020.4], 'FontSize', 12);
xticks(2003:2:2020);
ratesplot.Color='none';

% with NPP
weddellcolor=[0.1 0.6 0.5]
icecolor=[0 0 1]

NPP_icefree=figure;
NPP_icefree.Position=[150 80 800 700];
set(NPP_icefree,'defaultAxesColorOrder',[weddellcolor; icecolor]);

yyaxis left
 bar(x,y,'FaceColor',[0.1 0.6 0.5])
  hold on
 ylabel('Annual NPP (Tg C)')
  yyaxis right
d(1)=plot(years,OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(:,2),'LineWidth',1.5);
 ylabel('Open Water Area (km^2)')
set(gca,'YLim', [1e6 1.9e6], 'FontSize',12);
  xtickangle(45)
  set(gca, 'XTick', [2003:1:2019])
title('Annual NPP and mean open water area in Weddell Gyre')


%% Relationship

