%% data ready
load('ProcessedDataMonthly_Jan22.mat')
SeaIce.Regression.MeanSIE=mean(SeaIce.Weddell.SIExtent_aus,1,'omitnan');
SeaIce.Regression.MeanSIE=SeaIce.Regression.MeanSIE';
SeaIce.Regression.NPP_AnTot=OceanProd.cafe.Weddell.NPP_tot_TgC_annual(2:end,2);
SeaIce.Regression.Year=yearrange0320;

%% Regression Jan 22
SeaIce.Regression.tbl=table(SeaIce.Regression.Year,SeaIce.Regression.MeanSIE,SeaIce.Regression.NPP_AnTot,'VariableNames',{'Year','SIE','NPP'});

SeaIce.Regression.lm=fitlm(SeaIce.Regression.tbl,'NPP~SIE')

regres=fitlm(SeaIce.Regression.tbl,'NPP~SIE')

%% plot
figure;
t = tiledlayout(1,1);
ax1 = axes(t);
plot(regres)
% ax1.XColor = 'r';
% ax1.YColor = 'r';
hold on
ax2 = axes(t);
plot(regres.lm.(region_sublist{3}))
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';

