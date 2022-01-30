%% data ready
% clearvars
cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
load('ProcessedData_8day_Jan22_wWAP.mat') %ProcessedData_8day_Jan22
load('SeaIceDaily_Jan22_wWAP.mat') %SeaIceDaily_Jan22

% Make max-min SIE
for rix = 1:3%length(region_sublist)
    SeaIce.(region_sublist{rix}).SIE_ausMAX=max(SeaIce.(region_sublist{rix}).SIExtent_aus,[],1,'omitnan');
    SeaIce.(region_sublist{rix}).SIE_ausMIN=min(SeaIce.(region_sublist{rix}).SIExtent_aus,[],1,'omitnan');
    SeaIce.(region_sublist{rix}).SIE_ausDIF=(SeaIce.(region_sublist{rix}).SIE_ausMAX-SeaIce.(region_sublist{rix}).SIE_ausMIN)';
figure; bar(SeaIce.(region_sublist{rix}).SIE_ausMAX); hold on; bar(SeaIce.(region_sublist{rix}).SIE_ausMIN); plot(SeaIce.(region_sublist{rix}).SIE_ausDIF); 
% ylim([3.5e6 6e6])
end

% make structure of all variables for corr testing
for rix = 1:3 %length(region_sublist)

Regression.MeanSIE.(region_sublist{rix})=mean(SeaIce.(region_sublist{rix}).SIExtent_aus,1,'omitnan');
Regression.MeanSIE.(region_sublist{rix})=(Regression.MeanSIE.(region_sublist{rix})/1e6)';
Regression.MeanSIA.(region_sublist{rix})=mean(SeaIce.(region_sublist{rix}).SIArea_aus,1,'omitnan');
Regression.MeanSIA.(region_sublist{rix})=(Regression.MeanSIA.(region_sublist{rix})/1e6)';

Regression.MeanIFE.(region_sublist{rix})=(SeaIce.g_area.(region_sublist{rix})/1e6)-Regression.MeanSIE.(region_sublist{rix});
Regression.MeanIFA.(region_sublist{rix})=(SeaIce.g_area.(region_sublist{rix})/1e6)-Regression.MeanSIA.(region_sublist{rix});

Regression.SIE_dif.(region_sublist{rix})=SeaIce.(region_sublist{rix}).SIE_ausDIF;

Regression.NPP_AnTot.(region_sublist{rix})=OceanProd_8day.cafe.(region_sublist{rix}).NPP_tot_TgC_annual(2:end,2);
Regression.NPP_AnRate.(region_sublist{rix})=OceanProd_8day.cafe.(region_sublist{rix}).annual_rate; %_chunk (prev version) % annual rate of NPP (gC m^-2 a^-1)
Regression.NPP_AvGSRate.(region_sublist{rix})=OceanProd_8day.cafe.(region_sublist{rix}).GSav_dayrate; % average daily rate of NPP (mgC m^-2 d^-1)

Regression.Year=yearrange0320;
end
Regression.NCP=NCP(:,2);
rix=2
%% Spearman Correlation (for non-normally distributed data)
for rix = 1:3%length(region_sublist)
lev=kstest(Regression.NPP_AnRate.(region_sublist{rix})) % etc - all not normal
end

for rix = 1:3%length(region_sublist)

[rho.IFEvNPPTot.(region_sublist{rix}),pval.IFEvNPPTot.(region_sublist{rix})] = corr(Regression.MeanIFE.(region_sublist{rix}),Regression.NPP_AnTot.(region_sublist{rix}), 'type', 'Spearman');
[rho.IFAvNPPTot.(region_sublist{rix}),pval.IFAvNPPTot.(region_sublist{rix})] = corr(Regression.MeanIFA.(region_sublist{rix}),Regression.NPP_AnTot.(region_sublist{rix}), 'type', 'Spearman');
[rho.IFEvNPPRate.(region_sublist{rix}),pval.IFEvNPPRate.(region_sublist{rix})] = corr(Regression.MeanIFE.(region_sublist{rix}),Regression.NPP_AnRate.(region_sublist{rix}), 'type', 'Spearman');
[rho.IFAvNPPRate.(region_sublist{rix}),pval.IFAvNPPRate.(region_sublist{rix})] = corr(Regression.MeanIFA.(region_sublist{rix}),Regression.NPP_AnRate.(region_sublist{rix}), 'type', 'Spearman');
[rho.IFEvNPPGSRate.(region_sublist{rix}),pval.IFEvNPPGSRate.(region_sublist{rix})] = corr(Regression.MeanIFE.(region_sublist{rix}),Regression.NPP_AvGSRate.(region_sublist{rix}), 'type', 'Spearman');
[rho.IFAvNPPGSRate.(region_sublist{rix}),pval.IFAvNPPGSRate.(region_sublist{rix})] = corr(Regression.MeanIFA.(region_sublist{rix}),Regression.NPP_AvGSRate.(region_sublist{rix}), 'type', 'Spearman');

[rho.IFAvIFE.(region_sublist{rix}),pval.IFAvIFE.(region_sublist{rix})] = corr(Regression.MeanIFA.(region_sublist{rix}),Regression.MeanIFE.(region_sublist{rix}), 'type', 'Spearman');

[rho.DIFvNPPTot.(region_sublist{rix}),pval.DIFvNPPTot.(region_sublist{rix})] = corr(Regression.SIE_dif.(region_sublist{rix}),Regression.NPP_AnTot.(region_sublist{rix}), 'type', 'Spearman');
[rho.DIFvNPPRate.(region_sublist{rix}),pval.DIFvNPPRate.(region_sublist{rix})] = corr(Regression.SIE_dif.(region_sublist{rix}),Regression.NPP_AnRate.(region_sublist{rix}), 'type', 'Spearman');
[rho.DIFvNPPGSRate.(region_sublist{rix}),pval.DIFvNPPGSRate.(region_sublist{rix})] = corr(Regression.SIE_dif.(region_sublist{rix}),Regression.NPP_AvGSRate.(region_sublist{rix}), 'type', 'Spearman');

end

[rho.NCPvNPPTot,pval.NCPvNPPTot] = corr(Regression.NCP(13:end),Regression.NPP_AnTot.Weddell(13:end), 'type', 'Spearman');
[rho.IFEvNCP,pval.IFEvNCP] = corr(Regression.MeanIFE.Open(13:end),Regression.NCP(13:end), 'type', 'Spearman');
[rho.IFAvNCP,pval.IFAvNCP] = corr(Regression.MeanIFA.Weddell(13:end),Regression.NCP(13:end), 'type', 'Spearman');
[rho.DIFvNCP,pval.DIFvNCP] = corr(Regression.SIE_dif.Weddell(13:end),Regression.NCP(13:end), 'type', 'Spearman');

figure; scatter(Regression.SIE_dif.Weddell(13:end),Regression.NCP(13:end))
%% Regression Jan 22
rix = 1
Regression.tbl=table(Regression.Year,Regression.MeanSIE.(region_sublist{rix}),Regression.MeanIFE.(region_sublist{rix}),...
    Regression.MeanSIA.(region_sublist{rix}),Regression.MeanIFA.(region_sublist{rix}),Regression.NPP_AnTot.(region_sublist{rix}),Regression.NPP_AnRate.(region_sublist{rix}),...
    'VariableNames',{'Year','SIE','IFE','SIA','IFA','NPP','NPPrate'});

Regression.lmNPPSIE=fitlm(Regression.tbl,'NPP~SIE');
Regression.lmNPPIFE=fitlm(Regression.tbl,'NPP~IFE');
Regression.lmNPPSIA=fitlm(Regression.tbl,'NPP~SIA');
Regression.lmNPPIFA=fitlm(Regression.tbl,'NPP~IFA');

disp(Regression.lmNPPSIE)
disp(Regression.lmNPPIFE)
disp(Regression.lmNPPSIA)
disp(Regression.lmNPPIFA)

%% plot
figure;
t = tiledlayout('flow')
% nexttile
% plot(Regression.lmNPPSIE)
nexttile
str={'R2=0.488, p<0.001';'NPP=102.28*IFE-11.534'}
plot(Regression.lmNPPIFE); annotation('textbox', [0.2, 0.8, 0.1, 0.1], 'String', str)
% nexttile
% plot(Regression.lmNPPSIA); annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "R2=0.488, p<0.001")
% nexttile
% plot(Regression.lmNPPIFA); annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "hi")

R=residuals(Regression.lmNPPIFA);
figure;
plotResiduals(Regression.lmNPPIFA)
figure; plotDiagnostics(Regression.lmNPPIFA,'cookd')
legend('show') % Show the legend
resid=table2array(Regression.lmNPPSIE.Residuals(:,1));
lev = leverage(resid);
fitted=Regression.lmNPPSIE.Fitted;
figure; tiledlayout('flow'); 
nexttile; scatter(fitted,resid); title('Residuals vs Fitted')
nexttile; qqplot(resid); title('Normal Q-Q') 
nexttile; histogram(resid); title('Distribution of residuals')
nexttile; scatter(lev,resid); title('Residuals vs Leverage')
%% SIA vs SIE
Regression.lmSIESIA=fitlm(Regression.tbl,'SIE~SIA');
disp(Regression.lmSIESIA)
figure;
t = tiledlayout('flow')
nexttile
plot(Regression.lmSIESIA)

%% rates vs ice

Regression.lmRateSIE=fitlm(Regression.tbl,'NPPrate~SIE');
Regression.lmRateIFE=fitlm(Regression.tbl,'NPPrate~IFE');
Regression.lmRateSIA=fitlm(Regression.tbl,'NPPrate~SIA');
Regression.lmRateIFA=fitlm(Regression.tbl,'NPPrate~IFA');

disp(Regression.lmRateSIE)
disp(Regression.lmRateIFE)
disp(Regression.lmRateSIA)
disp(Regression.lmRateIFA)

% plot
figure;
t = tiledlayout('flow')
nexttile
plot(Regression.lmRateSIE)
nexttile
plot(Regression.lmRateIFE)
nexttile
plot(Regression.lmRateSIA)
nexttile
plot(Regression.lmRateIFA)



%% pixel by pixel

%IceFree_pixels_years
%AnnualNPPRate_pixels_years
clearvars find*
% remove NaN NPP entries from all
for rix = 1:length(region_sublist)
    %     IceFree_pixels_years.(region_sublist{rix})(isnan(IceFree_pixels_years.(region_sublist{rix})))=0;
    %     AnnualNPPRate_pixels_years.(region_sublist{rix})(isnan(AnnualNPPRate_pixels_years.(region_sublist{rix})))=0;
    
    IceFree_pixels_years_COL.(region_sublist{rix})=reshape(IceFree_pixels_years.(region_sublist{rix}),[],1);
    AnnualNPPRate_pixels_years_COL.(region_sublist{rix})=reshape(AnnualNPPRate_pixels_years.(region_sublist{rix}),[],1);
    
    findNaN_NPP=isnan(AnnualNPPRate_pixels_years_COL.(region_sublist{rix}));
    findzero_ice=(IceFree_pixels_years_COL.(region_sublist{rix})==0);
    
        if findNaN_NPP==findzero_ice
            IceFree_pixels_years_COL.(region_sublist{rix})(findNaN_NPP)=[];
            AnnualNPPRate_pixels_years_COL.(region_sublist{rix})(findNaN_NPP)=[];
        else
            disp('No match')
        end
end

%
for rix = 1:length(region_sublist)
clearvars find*
IceFree_pixels_years.(region_sublist{rix})(IceFree_pixels_years.(region_sublist{rix})==0)=NaN;
   
    IceFree_pixels_years_mean.(region_sublist{rix})=mean(IceFree_pixels_years.(region_sublist{rix}),2,'omitnan');
%     IceFree_pixels_years_median.(region_sublist{rix})=median(IceFree_pixels_years.(region_sublist{rix}),2,'omitnan');
    AnnualNPPRate_pixels_years_mean.(region_sublist{rix})=mean(AnnualNPPRate_pixels_years.(region_sublist{rix}),2,'omitnan');
%     AnnualNPPRate_pixels_years_median.(region_sublist{rix})=median(AnnualNPPRate_pixels_years.(region_sublist{rix}),2,'omitnan');
    
%     if rix == 2
%         findWAPtip=find(AnnualNPPRate_pixels_years_mean.(region_sublist{rix})>200);
%         IceFree_pixels_years_mean.(region_sublist{rix})(findWAPtip)=[];
%         AnnualNPPRate_pixels_years_mean.(region_sublist{rix})(findWAPtip)=[];
%     end
        findNaN_NPP=isnan(AnnualNPPRate_pixels_years_mean.(region_sublist{rix}));
    findNaN_ice=isnan(IceFree_pixels_years_mean.(region_sublist{rix}));
    
        if findNaN_NPP==findNaN_ice
            IceFree_pixels_years_mean.(region_sublist{rix})(findNaN_NPP)=[];
            AnnualNPPRate_pixels_years_mean.(region_sublist{rix})(findNaN_NPP)=[];
        else
            disp('No match')
        end

    Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_years_COL.(region_sublist{rix}),AnnualNPPRate_pixels_years_COL.(region_sublist{rix}),'VariableNames',{'IceFree','AnNPPrate'});
    Regression.tbl_pixels_meanmed.(region_sublist{rix})=table(IceFree_pixels_years_mean.(region_sublist{rix}),...
        AnnualNPPRate_pixels_years_mean.(region_sublist{rix}),...
        'VariableNames',{'IceFreeMean','AnNPPrateMean'}); %AnnualNPPRate_pixels_years_median.(region_sublist{rix}),,'AnNPPrateMedian',IceFree_pixels_years_median.(region_sublist{rix}),'IceFreeMedian'
end

%% One-sample Kolmogorov-Smirnov test and Spearman
    % all data is not from a normal distribution
for rix = 1:length(region_sublist)
lev=kstest(IceFree_pixels_years_COL.(region_sublist{rix}))
lev=kstest(AnnualNPPRate_pixels_years_COL.(region_sublist{rix}))
lev=kstest(IceFree_pixels_years_mean.(region_sublist{rix}))
lev=kstest(AnnualNPPRate_pixels_years_mean.(region_sublist{rix}))
end

[rho,pval] = corr(IceFree_pixels_years_COL.(region_sublist{1}),AnnualNPPRate_pixels_years_COL.(region_sublist{1}), 'type', 'Spearman')
[rho,pval] = corr(IceFree_pixels_years_COL.(region_sublist{2}),AnnualNPPRate_pixels_years_COL.(region_sublist{2}), 'type', 'Spearman')
[rho,pval] = corr(IceFree_pixels_years_COL.(region_sublist{3}),AnnualNPPRate_pixels_years_COL.(region_sublist{3}), 'type', 'Spearman')

[rho,pval] = corr(IceFree_pixels_years_mean.(region_sublist{1}),AnnualNPPRate_pixels_years_mean.(region_sublist{1}), 'type', 'Spearman')
[rho,pval] = corr(IceFree_pixels_years_mean.(region_sublist{2}),AnnualNPPRate_pixels_years_mean.(region_sublist{2}), 'type', 'Spearman')
[rho,pval] = corr(IceFree_pixels_years_mean.(region_sublist{3}),AnnualNPPRate_pixels_years_mean.(region_sublist{3}), 'type', 'Spearman')

%% Regression
for rix=1:length(region_sublist)
Regression.lm_Pixel.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'AnNPPrate~IceFree');
 disp(Regression.lm_Pixel.(region_sublist{rix}))
Regression.lm_Pixelmean.(region_sublist{rix})=fitlm(Regression.tbl_pixels_meanmed.(region_sublist{rix}),'AnNPPrateMean~IceFreeMean');
 disp(Regression.lm_Pixelmean.(region_sublist{rix}))
% Regression.lm_Pixelmedian.(region_sublist{rix})=fitlm(Regression.tbl_pixels_meanmed.(region_sublist{rix}),'AnNPPrateMedian~IceFreeMedian');
%  disp(Regression.lm_Pixelmedian.(region_sublist{rix}))
end

figure;
t = tiledlayout(4,2)
for rix=1:length(region_sublist)
nexttile
plot(Regression.lm_Pixel.(region_sublist{rix}))
nexttile
plot(Regression.lm_Pixelmean.(region_sublist{rix}))
% nexttile
% plot(Regression.lm_Pixelmedian.(region_sublist{2}))
end

for rix = 1:length(region_sublist)
    [rho.(region_sublist{rix}),p.(region_sublist{rix})]=corr(AnnualNPPRate_pixels_years_mean.(region_sublist{rix}),IceFree_pixels_years_mean.(region_sublist{rix}),'Type','Pearson')
end