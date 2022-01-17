%% data ready
% clearvars
cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
load('ProcessedData_8day_Jan22.mat')
load('SeaIceDaily_Jan22.mat')
Regression.MeanSIE=mean(SeaIce.Weddell.SIExtent_aus,1,'omitnan');
Regression.MeanSIE=Regression.MeanSIE';
Regression.MeanSIA=mean(SeaIce.Weddell.SIArea_aus,1,'omitnan');
Regression.MeanSIA=Regression.MeanSIA';

Regression.MeanIFE=SeaIce.g_area.Weddell-Regression.MeanSIE;
Regression.MeanIFA=SeaIce.g_area.Weddell-Regression.MeanSIA;

Regression.NPP_AnTot=OceanProd_8day.cafe.Weddell.NPP_tot_TgC_annual(2:end,2);
Regression.NPP_AnRate=OceanProd_8day.cafe.Weddell.annual_rate_chunk; % annual rate of NPP (gC m^-2 a^-1)

Regression.Year=yearrange0320;

%% Regression Jan 22
Regression.tbl=table(Regression.Year,Regression.MeanSIE,Regression.MeanIFE,...
    Regression.MeanSIA,Regression.MeanIFA,Regression.NPP_AnTot,Regression.NPP_AnRate,...
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
nexttile
plot(Regression.lmNPPSIE)
nexttile
plot(Regression.lmNPPIFE)
nexttile
plot(Regression.lmNPPSIA)
nexttile
plot(Regression.lmNPPIFA)


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
IceFree_pixels_years.(region_sublist{rix})(IceFree_pixels_years.(region_sublist{rix})==0)=NaN;
for rix = 1:length(region_sublist)
   
    IceFree_pixels_years_mean.(region_sublist{rix})=mean(IceFree_pixels_years.(region_sublist{rix}),2,'omitnan');
    IceFree_pixels_years_median.(region_sublist{rix})=median(IceFree_pixels_years.(region_sublist{rix}),2,'omitnan');
    AnnualNPPRate_pixels_years_mean.(region_sublist{rix})=mean(AnnualNPPRate_pixels_years.(region_sublist{rix}),2,'omitnan');
    AnnualNPPRate_pixels_years_median.(region_sublist{rix})=median(AnnualNPPRate_pixels_years.(region_sublist{rix}),2,'omitnan');
    
    if rix == 2
        findWAPtip=find(AnnualNPPRate_pixels_years_mean.(region_sublist{rix})>200);
        IceFree_pixels_years_mean.(region_sublist{rix})(findWAPtip)=[];
        AnnualNPPRate_pixels_years_mean.(region_sublist{rix})(findWAPtip)=[];
    end
    
    Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_years_COL.(region_sublist{rix}),AnnualNPPRate_pixels_years_COL.(region_sublist{rix}),'VariableNames',{'IceFree','AnNPPrate'});
    Regression.tbl_pixels_meanmed.(region_sublist{rix})=table(IceFree_pixels_years_mean.(region_sublist{rix}),IceFree_pixels_years_median.(region_sublist{rix}),...
        AnnualNPPRate_pixels_years_mean.(region_sublist{rix}),AnnualNPPRate_pixels_years_median.(region_sublist{rix}),...
        'VariableNames',{'IceFreeMean','IceFreeMedian','AnNPPrateMean','AnNPPrateMedian'});
end
rix=1
for rix=2%1:length(region_sublist)
Regression.lm_Pixel.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'AnNPPrate~IceFree');
 disp(Regression.lm_Pixel.(region_sublist{rix}))
Regression.lm_Pixelmean.(region_sublist{rix})=fitlm(Regression.tbl_pixels_meanmed.(region_sublist{rix}),'AnNPPrateMean~IceFreeMean');
 disp(Regression.lm_Pixelmean.(region_sublist{rix}))
% Regression.lm_Pixelmedian.(region_sublist{rix})=fitlm(Regression.tbl_pixels_meanmed.(region_sublist{rix}),'AnNPPrateMedian~IceFreeMedian');
%  disp(Regression.lm_Pixelmedian.(region_sublist{rix}))
end

figure;
t = tiledlayout(2,2)
nexttile
plot(Regression.lm_Pixel.(region_sublist{2}))
nexttile
plot(Regression.lm_Pixelmean.(region_sublist{2}))
% nexttile
% plot(Regression.lm_Pixelmedian.(region_sublist{2}))
nexttile
plot(Regression.lm_Pixel.(region_sublist{3}))
nexttile
plot(Regression.lm_Pixelmean.(region_sublist{3}))
% nexttile
% plot(Regression.lm_Pixelmedian.(region_sublist{3}))

for rix = 1:length(region_sublist)
    [rho.(region_sublist{rix}),p.(region_sublist{rix})]=corr(AnnualNPPRate_pixels_years_mean.(region_sublist{rix}),IceFree_pixels_years_mean.(region_sublist{rix}),'Type','Pearson')
end