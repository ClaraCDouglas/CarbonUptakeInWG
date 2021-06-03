%% Regression analysis
clearvars regres lmtest
aix = 4;
for rix = 1:length(region_sublist)

regres.year=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,1);
regres.icesmall.(region_sublist{1})=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,2)/1e6;
regres.icesmall.(region_sublist{3})=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,2)/1e6;
regres.icesmall.(region_sublist{2})=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,2)/1e4;
regres.ice.(region_sublist{rix})=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,2);
regres.NPP.(region_sublist{rix})=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(2:18,2);

regres.tbl.(region_sublist{rix})=table(regres.year,regres.ice.(region_sublist{rix}),regres.NPP.(region_sublist{rix}),'VariableNames',{'Year','icefree','NPP'});
regres.tbl.(region_sublist{rix})(1:5,:)
regres.tblsmall.(region_sublist{rix})=table(regres.year,regres.icesmall.(region_sublist{rix}),regres.NPP.(region_sublist{rix}),'VariableNames',{'Year','icefree','NPP'});
regres.tblsmall.(region_sublist{rix})(1:5,:)

% regres.tbl.(region_sublist{rix})=table(OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,1),...
%     OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,2),...
%     OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(2:18,2),'VariableNames',{'Year','NPP','icefree'});
% regres.tbl.(region_sublist{rix})(1:5,:)

regres.lm.(region_sublist{rix})=fitlm(regres.tbl.(region_sublist{rix}),'NPP~icefree')
regres.lmsmall.(region_sublist{rix})=fitlm(regres.tblsmall.(region_sublist{rix}),'NPP~icefree')

% lmtest=fitlm(regres.tbl.(region_sublist{rix}),'NPP~icefree')

end

% figure - need 2 x and y axes
clear NPPice_regression
NPPice_regression=figure;
NPPice_regression.Position=[150 80 800 700];
interesting=plot(regres.lm.(region_sublist{2}))
hold on
plot(regres.lm.(region_sublist{3}))

% 2 x and y axes, but doesn't work?
figure;
t = tiledlayout(1,1);
ax1 = axes(t);
plot(regres.lm.(region_sublist{2}))
% ax1.XColor = 'r';
% ax1.YColor = 'r';
hold on
ax2 = axes(t);
plot(regres.lm.(region_sublist{3}))
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';

% attempt to use predict function for confidence bounds...
[ypred,yci] = predict(regres.lm.(region_sublist{2}),regres.ice.(region_sublist{2}));
figure;
plot(regres.ice.(region_sublist{2}),regres.NPP.(region_sublist{2});
hold on
plot(regres.ice.(region_sublist{2}),yci);


% just scatter plot with regression line included
%VARIABLES
graphical.shelfcolor=[0.9290 0.6940 0.1250];
graphical.opencolor=[0.5 0.3 0.8];
graphical.SRsize=52;
graphical.OOsize=65;

clear scatterlm
scatterlm=figure(1);
scatterlm.Position=[150 80 700 700];
ax1=axes(scatterlm);
scatterSR=scatter((regres.ice.(region_sublist{2})/1e4),regres.NPP.(region_sublist{2}),graphical.SRsize,graphical.shelfcolor,'filled')
ylim([0 12.2])
xlabel('Annual mean area of open ice-free water on shelf (x10^4 km^2)','Interpreter','tex','FontSize',12)
ylabel('Total Annual NPP on shelf (Tg C)','Interpreter','tex','FontSize',12)
hold on
ax1.XColor=graphical.shelfcolor;
ax1.YColor=graphical.shelfcolor;
    % add regression line
        % REGRESSION LINE VARIABLES
graphical.SRIcemax=max(regres.ice.(region_sublist{2})/1e4);
graphical.SRIcemin=min(regres.ice.(region_sublist{2})/1e4);
graphical.xLSR = graphical.SRIcemin:0.02:graphical.SRIcemax;
graphical.mLSR = table2array(regres.lmsmall.(region_sublist{2}).Coefficients(2,1));
graphical.cLSR = table2array(regres.lmsmall.(region_sublist{2}).Coefficients(1,1));
graphical.yLSR = graphical.mLSR * graphical.xLSR + graphical.cLSR;

plot(graphical.xLSR,graphical.yLSR,'Color',graphical.shelfcolor);
graphical.text1=text(3.6,6, 'R^2=0.83, p<0.001','Color',graphical.shelfcolor,'Interpreter','tex','FontSize',12);

% add 2nd axis with open ocean data
ax2=axes(scatterlm);
scatterOO=scatter(regres.icesmall.(region_sublist{3}),regres.NPP.(region_sublist{3}),graphical.OOsize,graphical.opencolor,'filled')
hold on
%'s',
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color='none';
ax2.XColor=graphical.opencolor;
ax2.YColor=graphical.opencolor;
xlabel('Annual mean area of open ice-free water in open ocean (x10^6 km^2)','Interpreter','tex','FontSize',12)
ylabel('Total Annual NPP in open ocean(Tg C)','Interpreter','tex','FontSize',12)
graphical.OOIcemax=max(regres.icesmall.(region_sublist{3}))
graphical.OOIcemin=min(regres.icesmall.(region_sublist{3}))
graphical.xLOO = graphical.OOIcemin:0.02:graphical.OOIcemax;
graphical.mLOO = table2array(regres.lmsmall.(region_sublist{3}).Coefficients(2,1));
graphical.cLOO = table2array(regres.lmsmall.(region_sublist{3}).Coefficients(1,1));
graphical.yLOO = graphical.mLOO * graphical.xLOO + graphical.cLOO;
plot(graphical.xLOO,graphical.yLOO,'Color',graphical.opencolor);
lg2=legend([scatterSR, scatterOO],'Shelf','Open Ocean','Location','northwest')

graphical.text2=text(1.4,80, 'R^2=0.54, p<0.001','Color',graphical.opencolor,'Interpreter','tex','FontSize',12);


%% model 2 needed???

aix = 4;
for rix = 1:length(region_sublist)

regres.year=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,1);
regres.icesmall.(region_sublist{1})=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,2)/1e6;
regres.icesmall.(region_sublist{3})=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,2)/1e6;
regres.icesmall.(region_sublist{2})=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,2)/1e4;
regres.ice.(region_sublist{rix})=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,2);
regres.NPP.(region_sublist{rix})=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(2:18,2);

regres.tbl2.(region_sublist{rix})=table(regres.year,regres.NPP.(region_sublist{rix}),regres.ice.(region_sublist{rix}),'VariableNames',{'Year','NPP','icefree'});
regres.tbl2.(region_sublist{rix})(1:5,:)
regres.tblsmall2.(region_sublist{rix})=table(regres.year,regres.NPP.(region_sublist{rix}),regres.ice.(region_sublist{rix}),'VariableNames',{'Year','NPP','icefree'});
regres.tblsmall2.(region_sublist{rix})(1:5,:)

% regres.tbl.(region_sublist{rix})=table(OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,1),...
%     OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(2:18,2),...
%     OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(2:18,2),'VariableNames',{'Year','NPP','icefree'});
% regres.tbl.(region_sublist{rix})(1:5,:)

regres.lm2.(region_sublist{rix})=fitlm(regres.tbl.(region_sublist{rix}),'icefree~NPP')
regres.lmsmall2.(region_sublist{rix})=fitlm(regres.tblsmall.(region_sublist{rix}),'icefree~NPP')

% lmtest=fitlm(regres.tbl.(region_sublist{rix}),'NPP~icefree')

end

