%% data ready
% clearvar
desktop = 0;
laptop=1;
if desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg'))
    cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Matlab Add-ins'));
elseif laptop
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG'))
    cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\m_map'));
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\SetCMap'));
end
load('ProcessedData_8day_May22.mat') %ProcessedData_8day_Jan22
load('SeaIcev4_8day_May22.mat') %SeaIceDaily_Jan22
close all

% Make max-min SIE
figure; tiledlayout('flow')
for rix = 1:length(region_sublist)
    SeaIce.(region_sublist{rix}).SIE_ausMAX=max(SeaIce.(region_sublist{rix}).SIExtent_aus,[],1,'omitnan');
    SeaIce.(region_sublist{rix}).SIE_ausMIN=min(SeaIce.(region_sublist{rix}).SIExtent_aus,[],1,'omitnan');
    SeaIce.(region_sublist{rix}).SIE_ausDIF=(SeaIce.(region_sublist{rix}).SIE_ausMAX-SeaIce.(region_sublist{rix}).SIE_ausMIN)';
    nexttile; bar(SeaIce.(region_sublist{rix}).SIE_ausMAX,'FaceColor',[0.2 0.55 0.75]); hold on; bar(SeaIce.(region_sublist{rix}).SIE_ausMIN,'FaceColor',[0.7 0.7 0.7]); plot(SeaIce.(region_sublist{rix}).SIE_ausDIF,'Color','k','LineWidth',2);
    if rix==4
        title('AP')
    else
        title((region_sublist{rix}))
    end
    % ylim([3.5e6 6e6])
    xticks([1:1:19])
    xticklabels(2003:1:2021)
    xtickangle(45)
end
sgtitle('Maximum and minimum sea ice extent and difference between winter max and summer min')
yearrange0321=permute(2003:1:2021,[2,1]);

% make structure of all variables for corr testing
for rix = 1:length(region_sublist)
Regression.MeanSIE.(region_sublist{rix})=mean(SeaIce.(region_sublist{rix}).SIExtent_aus,1,'omitnan');
Regression.MeanSIE.(region_sublist{rix})=(Regression.MeanSIE.(region_sublist{rix})/1e6)';
Regression.MeanSIA.(region_sublist{rix})=mean(SeaIce.(region_sublist{rix}).SIArea_aus,1,'omitnan');
Regression.MeanSIA.(region_sublist{rix})=(Regression.MeanSIA.(region_sublist{rix})/1e6)';

Regression.MeanIFE.(region_sublist{rix})=(SeaIce.g_area.(region_sublist{rix})/1e6)-Regression.MeanSIE.(region_sublist{rix});
Regression.MeanIFA.(region_sublist{rix})=(SeaIce.g_area.(region_sublist{rix})/1e6)-Regression.MeanSIA.(region_sublist{rix});

Regression.SIE_dif.(region_sublist{rix})=SeaIce.(region_sublist{rix}).SIE_ausDIF;

Regression.NPP_AnTot.(region_sublist{rix})=OceanProd_8day.cafe.(region_sublist{rix}).NPP_tot_TgC_annual(2:end,2);
Regression.NPP_AnRate.(region_sublist{rix})=OceanProd_8day.cafe.(region_sublist{rix}).AnnualNPPRate_gm2peryear; % % annual rate of NPP (gC m^-2 a^-1)
Regression.NPP_AvGSRate.(region_sublist{rix})=OceanProd_8day.cafe.(region_sublist{rix}).AnAvDayRate_mgm2d1; % average daily rate of NPP (mgC m^-2 d^-1)

Regression.IFDav.(region_sublist{rix})=OceanProd_8day.cafe.(region_sublist{rix}).IFD_av; % average daily rate of NPP (mgC m^-2 d^-1)
Regression.IFDmax.(region_sublist{rix})=OceanProd_8day.cafe.(region_sublist{rix}).IFD_max; % average daily rate of NPP (mgC m^-2 d^-1)

Regression.Year=yearrange0321;
Regression.NPP_AnTot_sqrt.(region_sublist{rix})=sqrt(OceanProd_8day.cafe.(region_sublist{rix}).NPP_tot_TgC_annual(2:end,2));

end
% Regression.NCP=NCP(:,2);
% rix=2

%% Spearman Correlation (for non-normally distributed data)
for rix = 1:3%length(region_sublist)
lev=kstest(Regression.NPP_AnRate.(region_sublist{rix})) % etc - all not normal
end

for rix = 1:length(region_sublist)
    [rho.IFAvIFD.(region_sublist{rix}),pval.IFAvIFD.(region_sublist{rix})] = corr(Regression.MeanIFA.(region_sublist{rix}),Regression.IFDav.(region_sublist{rix}), 'type', 'Spearman');
    [rho.IFEvIFD.(region_sublist{rix}),pval.IFEvIFD.(region_sublist{rix})] = corr(Regression.MeanIFE.(region_sublist{rix}),Regression.IFDav.(region_sublist{rix}), 'type', 'Spearman');
    [rho.IFAvIFE.(region_sublist{rix}),pval.IFAvIFE.(region_sublist{rix})] = corr(Regression.MeanIFA.(region_sublist{rix}),Regression.MeanIFE.(region_sublist{rix}), 'type', 'Spearman');

    
    [rho.IFAvNPPtot.(region_sublist{rix}),pval.IFAvNPPtot.(region_sublist{rix})] = corr(Regression.MeanIFA.(region_sublist{rix}),Regression.NPP_AnTot.(region_sublist{rix}), 'type', 'Spearman');
    [rho.IFAvNPPnorm.(region_sublist{rix}),pval.IFAvNPPnorm.(region_sublist{rix})] = corr(Regression.MeanIFA.(region_sublist{rix}),Regression.NPP_AnRate.(region_sublist{rix}), 'type', 'Spearman');
    [rho.NPPtotvIFD.(region_sublist{rix}),pval.NPPtotvIFD.(region_sublist{rix})] = corr(Regression.NPP_AnTot.(region_sublist{rix}),Regression.IFDav.(region_sublist{rix}), 'type', 'Spearman');
    [rho.NPPnormvIFD.(region_sublist{rix}),pval.NPPnormvIFD.(region_sublist{rix})] = corr(Regression.NPP_AnRate.(region_sublist{rix}),Regression.IFDav.(region_sublist{rix}), 'type', 'Spearman');
    [rho.IFAvNPPday.(region_sublist{rix}),pval.IFAvNPPday.(region_sublist{rix})] = corr(Regression.MeanIFA.(region_sublist{rix}),Regression.NPP_AvGSRate.(region_sublist{rix}), 'type', 'Spearman');
    [rho.NPPdayvIFD.(region_sublist{rix}),pval.NPPdayvIFD.(region_sublist{rix})] = corr(Regression.NPP_AvGSRate.(region_sublist{rix}),Regression.IFDav.(region_sublist{rix}), 'type', 'Spearman');

    [rho.NPPdayvtot.(region_sublist{rix}),pval.NPPdayvtot.(region_sublist{rix})] = corr(Regression.NPP_AvGSRate.(region_sublist{rix}),Regression.NPP_AnTot.(region_sublist{rix}), 'type', 'Spearman');
    [rho.NPPdayvrate.(region_sublist{rix}),pval.NPPdayvrate.(region_sublist{rix})] = corr(Regression.NPP_AvGSRate.(region_sublist{rix}),Regression.NPP_AnRate.(region_sublist{rix}), 'type', 'Spearman');
    [rho.NPPtotvrate.(region_sublist{rix}),pval.NPPtotvrate.(region_sublist{rix})] = corr(Regression.NPP_AnTot.(region_sublist{rix}),Regression.NPP_AnRate.(region_sublist{rix}), 'type', 'Spearman');
end
%           OR Just Use Correlation Matrix!
for rix = 1:length(region_sublist)
    all=[Regression.Year,Regression.MeanSIE.(region_sublist{rix}),Regression.MeanIFE.(region_sublist{rix}),...
        Regression.MeanSIA.(region_sublist{rix}),Regression.MeanIFA.(region_sublist{rix}),Regression.SIE_dif.(region_sublist{rix}),...
        Regression.NPP_AnTot.(region_sublist{rix}),Regression.NPP_AnRate.(region_sublist{rix}),Regression.NPP_AvGSRate.(region_sublist{rix}),...
        Regression.IFDav.(region_sublist{rix}),Regression.IFDmax.(region_sublist{rix})];
    %'VariableNames',{'Year','SIE','IFE','SIA','IFA','SIE_dif','NPP','NPPrate','NPPGSrate','IFDav','IFDmax'}); %
    
    [R.(region_sublist{rix}),P.(region_sublist{rix})] = corrcoef(all);
end
writetable(struct2table(R), 'someexcelfile.xlsx')
writetable(struct2table(P), 'someexcelfileP.xlsx')


%% Region means
for rix = 1:length(region_sublist)

Regression.tblmeans.(region_sublist{rix})=table(mean(Regression.MeanIFE.(region_sublist{rix})),...
    mean(Regression.MeanIFA.(region_sublist{rix})),...
    mean(Regression.NPP_AnTot.(region_sublist{rix})),mean(Regression.NPP_AnRate.(region_sublist{rix})),mean(Regression.NPP_AvGSRate.(region_sublist{rix})),...
    mean(Regression.IFDav.(region_sublist{rix})),mean(Regression.IFDmax.(region_sublist{rix})),...
    'VariableNames',{'IFE','IFA','NPP','NPPrate','NPPGSrate','IFDav','IFDmax'}); %
end
Regression.tblmeans.all=cat(1,Regression.tblmeans.Weddell,Regression.tblmeans.Shelf,Regression.tblmeans.Open,Regression.tblmeans.WAP);

%% Regression
for rix = 1:length(region_sublist)
Regression.tbl.(region_sublist{rix})=table(Regression.Year,Regression.MeanSIE.(region_sublist{rix}),Regression.MeanIFE.(region_sublist{rix}),...
    Regression.MeanSIA.(region_sublist{rix}),Regression.MeanIFA.(region_sublist{rix}),Regression.SIE_dif.(region_sublist{rix}),...
    Regression.NPP_AnTot.(region_sublist{rix}),Regression.NPP_AnRate.(region_sublist{rix}),Regression.NPP_AvGSRate.(region_sublist{rix}),...
    Regression.IFDav.(region_sublist{rix}),Regression.IFDmax.(region_sublist{rix}),...
    Regression.NPP_AnTot_sqrt.(region_sublist{rix}),...
    'VariableNames',{'Year','SIE','IFE','SIA','IFA','SIE_dif','NPP','NPPrate','NPPGSrate','IFDav','IFDmax','NPP_sqrt'}); %

% Regression.lmNPPSIE.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~SIE');
% Regression.lmNPPSIA.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~SIA');

%       Total NPP vs IFE/IFA
Regression.lmNPPIFE.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFE');
Regression.lmNPPIFA.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFA');
%       Area-normalised NPP vs IFE
Regression.lmRateIFE.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPPrate~IFE');
%       Mean daily rate vs IFA
Regression.lmGSRateIFA.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPPGSrate~IFA');
%       Total NPP vs mean IFD
Regression.lmNPPIFD.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFDav');
%       Total NPP vs max IFD
Regression.lmNPPIFDmax.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFDmax');
%       Total NPP vs IFA&meanIFD
Regression.lmNPP_IFAIFD.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFA+IFDav');
Regression.lmNPP_IFAIFDtest.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFA*IFDav');
%       Area-normalised NPP vs IFD
Regression.lmRateIFD.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPPrate~IFDav');

Regression.lmNPPdif.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~SIE_dif');
Regression.lmRatedif.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPPrate~SIE_dif');

end
model_list={'lmNPPIFE','lmNPPIFA','lmRateIFE','lmGSRateIFA','lmNPPIFD','lmNPPIFDmax','lmRateIFD','lmNPPdif'};

%       DISPLAY MODEL OUTPUTS
for rix = 1:length(region_sublist)
% disp(Regression.lmNPPSIE.(region_sublist{rix}))
% disp(Regression.lmNPPIFE.(region_sublist{rix}))
% disp(Regression.lmNPPSIA.(region_sublist{rix}))
% disp(Regression.lmNPPIFA.(region_sublist{rix}))
% disp(Regression.lmRateIFA.(region_sublist{rix}))
% disp(Regression.lmNPPdif.(region_sublist{rix}))
% disp(Regression.lmRatedif.(region_sublist{rix}))
% disp(Regression.lmGSRateIFA.(region_sublist{rix}))
% disp(Regression.lmNPPIFD.(region_sublist{rix}))
% disp(Regression.lmNPPIFDmax.(region_sublist{rix}))
% disp(Regression.lmNPP_IFAIFD.(region_sublist{rix}))
disp(Regression.lmNPP_IFAIFDtest.(region_sublist{rix}))
% disp(Regression.lmRateIFD.(region_sublist{rix}))
end

%       AIC
for rix = 1:length(region_sublist)
    for mix=1:length(model_list)
        disp([(region_sublist{rix}) (model_list{mix})])
        Regression.(model_list{mix}).(region_sublist{rix}).ModelCriterion.AIC
    end
end




%% plot
close(figure(5))
figure(5);
tiledlayout(2,2)
anpos=[0.1, 0.8, 0.1, 0.1;0.55, 0.8, 0.1, 0.1;0.1, 0.3, 0.1, 0.1;0.55, 0.3, 0.1, 0.1];
for rix = 1:length(region_sublist)
    nexttile
    plot(Regression.lmNPPIFE.(region_sublist{rix})) %lmNPPIFA %lmRateIFE %lmGSRateIFA %lmRateIFD
    title((region_sublist{rix}))
    xtxt={'Mean Ice Free Extent (10^6 km^2)'}; %Mean Ice Free Area (10^6 km^2) %Average # Ice Free Days
    xlabel(xtxt,'Interpreter','tex')
    ylabel('Annual NPP (TgC)','Interpreter','tex') %Annual NPP (TgC) %Annual NPP (gC m^2 a^{-1}) %Growing Season NPP (mg m^{-2} d^{-1})
legend('off')
    %     l=legend;
%     l.Location='southeast';
    if rix==2
        ylim([-0.2 inf])
        yline(0,':k')
    end
    c=Regression.lmNPPIFE.(region_sublist{rix}).Coefficients{1,1};
    m=Regression.lmNPPIFE.(region_sublist{rix}).Coefficients{2,1};
    R2=Regression.lmNPPIFE.(region_sublist{rix}).Rsquared.Adjusted;
    p=Regression.lmNPPIFE.(region_sublist{rix}).Coefficients{2,4};
    str={'R^2=' num2str(R2),' p=' num2str(p);'NPP=' num2str(m) '*IFE+' num2str(c)};
    str2={strjoin(str(1:2:8));strjoin(str(2:2:8))};
    annotation('textbox', anpos(rix,:),'String',str2,'FitBoxToText','on')
end

close(figure(10))
figure(10);
t=tiledlayout(2,2)
anpos=[0.1, 0.82, 0.1, 0.1;0.55, 0.82, 0.1, 0.1;0.1, 0.34, 0.1, 0.1;0.55, 0.34, 0.1, 0.1];
for rix = 1:length(region_sublist)
    nexttile
    plt=plot(Regression.lmRateIFE.(region_sublist{rix})) %lmNPPIFA %lmRateIFE %lmGSRateIFA
    hold on
    scatter(Regression.MeanIFE.(region_sublist{rix}),Regression.NPP_AnRate.(region_sublist{rix}),40,Regression.IFDav.(region_sublist{rix}),'filled') %lmNPPIFA %lmRateIFA %lmGSRateIFA
    title((region_sublist{rix}))
    cmocean('haline')
    cb=colorbar
    if rix==2 || rix==4
    ylabel(cb,'Mean Ice Free Days','FontSize',11);
    end
    xtxt={'Mean Ice Free Extent (10^6 km^2)'}; %Area
    xlabel(xtxt,'Interpreter','tex')
    ylabel('Annual NPP (gC m^2 a^{-1})','Interpreter','tex') %Annual NPP (TgC) %Annual NPP (gC m^2 a^{-1}) %Growing Season NPP (mg m^{-2} d^{-1})
    l=legend([plt]);
    l.Location='southeast';
%     if rix==2
%         ylim([-0.2 inf])
%         yline(0,':k')
%     end
    if rix>1
        legend('off')
    end
%     if rix==1
%         ylim([100 240])
%     end
    c=Regression.lmRateIFE.(region_sublist{rix}).Coefficients{1,1};
    m=Regression.lmRateIFE.(region_sublist{rix}).Coefficients{2,1};
    R2=Regression.lmRateIFE.(region_sublist{rix}).Rsquared.Adjusted;
    p=Regression.lmRateIFE.(region_sublist{rix}).Coefficients{2,4};
    str={'R^2=' num2str(R2),' p=' num2str(p);'NPP=' num2str(m) '*IFE+' num2str(c)};
    str2={strjoin(str(1:2:8));strjoin(str(2:2:8))};
    annotation('textbox', anpos(rix,:),'String',str2,'FitBoxToText','on')
    if rix==4
        title('Antarctic Peninsula')
    end
    set(gca,'FontSize', 11);
end
% t.Padding = 'compact';

%% Curvefitting tool
XX=Regression.MeanIFE.Shelf;
XXa=Regression.MeanIFA.Shelf;
YY=Regression.NPP_AnTot.Shelf;
ZZ=Regression.IFDav.Shelf;
YYr=Regression.NPP_AnRate.Shelf;
%[fitresult, gof] = createFit(XX, YY)
yysqrt=sqrt(YY);
yyRsqrt=sqrt(YYr);
XXo=Regression.MeanIFE.Open;
YYo=Regression.NPP_AnTot.Open;
YYosqrt=sqrt(YYo);

[fitresult, gof] = createFit(Regression.MeanIFE.Shelf, Regression.NPP_AnTot.Shelf)
[fitresultA, gofA] = createFit(Regression.MeanIFA.Shelf, Regression.NPP_AnTot.Shelf)

[fitresult, gof] = createFit(Regression.MeanIFE.WAP, Regression.NPP_AnTot.WAP)

lev=kstest(YYr) % etc - all not normal
levsqrt=kstest(yyRsqrt) % etc - all not normal
figure; tiledlayout('flow'); nexttile;
histogram(YYr,10)
nexttile; histogram(yyRsqrt,10)
%% Spearman Correlation (for non-normally distributed data)
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

%% residuals
R=residuals(Regression.lmRateIFE);
figure;
plotResiduals(Regression.lmNPPIFE.Open)
figure; plotDiagnostics(Regression.lmNPPIFE.Shelf,'cookd')
legend('show') % Show the legend
resid=table2array(Regression.lmNPPIFE.Open.Residuals(:,1));
lev = leverage(resid);
fitted=Regression.lmNPPIFE.Open.Fitted;
figure; tiledlayout('flow'); 
nexttile; scatter(fitted,resid); title('Residuals vs Fitted')
nexttile; qqplot(resid); title('Normal Q-Q') 
nexttile; histogram(resid,10); title('Distribution of residuals')
nexttile; scatter(lev,resid); title('Residuals vs Leverage')


%% SIA vs SIE
figure;
t = tiledlayout('flow')
for rix = 1:length(region_sublist)
Regression.lmIFEIFA.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'IFE~IFA');
disp(Regression.lmIFEIFA.(region_sublist{rix}))
nexttile
plot(Regression.lmIFEIFA.(region_sublist{rix}))
legend('off')
hold on 
x=0:3
plot(x,x)
xlim([(min(Regression.MeanIFA.(region_sublist{rix}))) (max(Regression.MeanIFA.(region_sublist{rix})))])
end
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
data_daily=true;
data_8day=false;
load data
cd 'D:\Data\SeaIceNIMBUS';
if data_daily
    load('SeaIce_daily_20022020.mat')
    load('SeaIce_8day_20022020.mat','g_area','g_lat','g_lon')
elseif data_8day
    load('SeaIce_8day_20022020.mat')
end
cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop

load('ProcessedCAFEArrays_May22.mat')

% Selecting values within region boxes only
for aix = 1:length(algorithm)
    for rix = 1:length(region_sublist)
        for yix=2003:2021
            temp_icefreedays=IceFreeDays_peryear(:,:,yix-2002);
            temp_icefreedays=round(temp_icefreedays);
            IceFree_pixels_years.(region_sublist{rix})(:,yix-2002)=temp_icefreedays(temp.(region_sublist{rix}).box_logic);
            
            temp_AnnualNPPRate=AnnualNPPRate_gperyear(:,:,yix-2002);
            AnnualNPPRate_pixels_years.(region_sublist{rix})(:,yix-2002)=temp_AnnualNPPRate(temp.(region_sublist{rix}).box_logic);
            Year_cols.(region_sublist{rix})(:,yix-2002)=AnnualNPPRate_pixels_years.(region_sublist{rix})(:,yix-2002);
            Year_cols.(region_sublist{rix})(:,yix-2002)=yix;
%             temp_AnAvDayRate=AnAvDayRate_mgm2d1(:,:,yix-2002);
%             AnAvDayRate_pixels_years.(region_sublist{rix})(:,yix-2002)=temp_AnAvDayRate(temp.(region_sublist{rix}).box_logic);
        end
    end
end
clearvars temp_*

% remove NaN NPP entries from all
for rix = 1:length(region_sublist)
    %     IceFree_pixels_years.(region_sublist{rix})(isnan(IceFree_pixels_years.(region_sublist{rix})))=0;
    %     AnnualNPPRate_pixels_years.(region_sublist{rix})(isnan(AnnualNPPRate_pixels_years.(region_sublist{rix})))=0;
    
    IceFree_pixels_years_COL.(region_sublist{rix})=reshape(IceFree_pixels_years.(region_sublist{rix}),[],1);
    AnnualNPPRate_pixels_years_COL.(region_sublist{rix})=reshape(AnnualNPPRate_pixels_years.(region_sublist{rix}),[],1);
    YEARS_COL.(region_sublist{rix})=reshape(Year_cols.(region_sublist{rix}),[],1);
%     AnAvDayRate_pixels_years_COL.(region_sublist{rix})=reshape(AnAvDayRate_pixels_years.(region_sublist{rix}),[],1);
    
    findNaN_NPP=isnan(AnnualNPPRate_pixels_years_COL.(region_sublist{rix}));
    findzero_ice=(IceFree_pixels_years_COL.(region_sublist{rix})==0);
%     findNaN_NPPd=isnan(AnAvDayRate_pixels_years_COL.(region_sublist{rix}));
    
        if findNaN_NPP==findzero_ice %& findzero_ice==findNaN_NPPd
            IceFree_pixels_years_COL.(region_sublist{rix})(findNaN_NPP)=[];
            AnnualNPPRate_pixels_years_COL.(region_sublist{rix})(findNaN_NPP)=[];
            YEARS_COL.(region_sublist{rix})(findNaN_NPP)=[];
%             AnAvDayRate_pixels_years_COL.(region_sublist{rix})(findNaN_NPPd)=[];
        else
            disp('No match')
        end
end


%% Pixel Regression

for rix=1:length(region_sublist)
    
    Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_years_COL.(region_sublist{rix}),...
        AnnualNPPRate_pixels_years_COL.(region_sublist{rix}),YEARS_COL.(region_sublist{rix}),...
        'VariableNames',{'IceFree','AnNPPrate','YEAR'}); %AnAvDayRate_pixels_years_COL.(region_sublist{rix}),...,'AnAvDayRate'
    
    Regression.lm_Pixel.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'AnNPPrate~IceFree');
    disp(Regression.lm_Pixel.(region_sublist{rix}))
    Regression.lm_Pixely.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'AnNPPrate~YEAR');
    disp(Regression.lm_Pixely.(region_sublist{rix}))

    Regression.lm_PixelIFy.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'AnNPPrate~IceFree*YEAR');
    disp(Regression.lm_PixelIFy.(region_sublist{rix}))

%     Regression.lm_Pixeld.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'AnAvDayRate~IceFree');
%     disp(Regression.lm_Pixeld.(region_sublist{rix}))
    % Regression.lm_Pixelmedian.(region_sublist{rix})=fitlm(Regression.tbl_pixels_meanmed.(region_sublist{rix}),'AnNPPrateMedian~IceFreeMedian');
    %  disp(Regression.lm_Pixelmedian.(region_sublist{rix}))
end

figure;
t = tiledlayout(2,4)
for rix=1:3%length(region_sublist)
    if rix == 1
        nexttile(2,[1 2])
    elseif rix==2
        nexttile(5,[1 2])
    elseif rix==3
        nexttile(7,[1 2])
    end
plot(Regression.lm_Pixel.(region_sublist{rix}))
    l=legend;
    l.Location='southeast';
    title((region_sublist{rix}))
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
% nexttile
% plot(Regression.lm_Pixelmean.(region_sublist{rix}))
% nexttile
% plot(Regression.lm_Pixelmedian.(region_sublist{2}))
end

figure;
t = tiledlayout(2,2)
for rix=1%:length(region_sublist)
    nexttile
    plot(Regression.lm_Pixel.(region_sublist{rix}))
    hold on
    scatter(IceFree_pixels_years_COL.(region_sublist{rix}),AnnualNPPRate_pixels_years_COL.(region_sublist{rix}),...
        20,YEARS_COL.(region_sublist{rix}),'filled') %lmNPPIFA %lmRateIFA %lmGSRateIFA
    title((region_sublist{rix}))
    cmocean('haline')
    cb=colorbar
    l=legend;
    l.Location='northeast';
    title((region_sublist{rix}))
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
    xlim([0 300])
    % nexttile
    % plot(Regression.lm_Pixelmean.(region_sublist{rix}))
    % nexttile
    % plot(Regression.lm_Pixelmedian.(region_sublist{2}))
end


figure;
t = tiledlayout(2,2)
for rix=1:length(region_sublist)
nexttile
boxplot(AnAvDayRate_pixels_years_COL.(region_sublist{rix}),IceFree_pixels_years_COL.(region_sublist{rix}))
    title((region_sublist{rix}))
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Mean daily NPP (g m^{-2} d^{-1})','Interpreter','tex')
end


for rix = 1:length(region_sublist)
    [rho.(region_sublist{rix}),p.(region_sublist{rix})]=corr(AnnualNPPRate_pixels_years_mean.(region_sublist{rix}),IceFree_pixels_years_mean.(region_sublist{rix}),'Type','Pearson')
end

%% means regressions by location
for rix = 1:length(region_sublist)
    IceFree_pixels_yearsN.(region_sublist{rix})=IceFree_pixels_years.(region_sublist{rix});
    IceFree_pixels_yearsN.(region_sublist{rix})(IceFree_pixels_years.(region_sublist{rix})==0)=NaN;
    IceFree_pixels_MEANS.(region_sublist{rix})=mean(IceFree_pixels_yearsN.(region_sublist{rix}),2,'omitnan');
    AnnualNPPRate_pixels_MEAN.(region_sublist{rix})=mean(AnnualNPPRate_pixels_years.(region_sublist{rix}),2,'omitnan');
    AnAvDayRate_pixels_MEAN.(region_sublist{rix})=mean(AnAvDayRate_pixels_years.(region_sublist{rix}),2,'omitnan');
end
for rix=1:length(region_sublist)
    
    Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_MEANS.(region_sublist{rix}),...
        AnnualNPPRate_pixels_MEAN.(region_sublist{rix}),AnAvDayRate_pixels_MEAN.(region_sublist{rix}),...
        'VariableNames',{'IceFreeM','AnNPPrateM','AnAvDayRateM'});
    
    Regression.lm_PixelM.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'AnNPPrateM~IceFreeM');
    disp(Regression.lm_PixelM.(region_sublist{rix}))
    Regression.lm_PixeldM.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'AnAvDayRateM~IceFreeM');
    disp(Regression.lm_PixeldM.(region_sublist{rix}))
    % Regression.lm_Pixelmedian.(region_sublist{rix})=fitlm(Regression.tbl_pixels_meanmed.(region_sublist{rix}),'AnNPPrateMedian~IceFreeMedian');
    %  disp(Regression.lm_Pixelmedian.(region_sublist{rix}))
end

figure;
t = tiledlayout(2,2)
for rix=1:length(region_sublist)
nexttile
plot(Regression.lm_PixelM.(region_sublist{rix}))
    l=legend;
    l.Location='southeast';
    title((region_sublist{rix}))
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
% nexttile
% plot(Regression.lm_Pixelmean.(region_sublist{rix}))
% nexttile
% plot(Regression.lm_Pixelmedian.(region_sublist{2}))
end

figure;
t = tiledlayout(2,2)
for rix=1:length(region_sublist)
nexttile
plot(Regression.lm_PixeldM.(region_sublist{rix}))
    l=legend;
    l.Location='northeast';
    title((region_sublist{rix}))
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Mean daily NPP (g m^{-2} d^{-1})','Interpreter','tex')
    xlim([0 300])
% nexttile
% plot(Regression.lm_Pixelmean.(region_sublist{rix}))
% nexttile
% plot(Regression.lm_Pixelmedian.(region_sublist{2}))
end


%% means regressions by year
for rix = 1:length(region_sublist)
    IceFree_pixels_MEANSy.(region_sublist{rix})=permute(mean(IceFree_pixels_yearsN.(region_sublist{rix}),1,'omitnan'),[2 1]);
    AnnualNPPRate_pixels_MEANy.(region_sublist{rix})=permute(mean(AnnualNPPRate_pixels_years.(region_sublist{rix}),1,'omitnan'),[2 1]);
    AnAvDayRate_pixels_MEANy.(region_sublist{rix})=permute(mean(AnAvDayRate_pixels_years.(region_sublist{rix}),1,'omitnan'),[2 1]);
end
for rix=1:length(region_sublist)
    Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_MEANSy.(region_sublist{rix}),...
        AnnualNPPRate_pixels_MEANy.(region_sublist{rix}),AnAvDayRate_pixels_MEANy.(region_sublist{rix}),...
        'VariableNames',{'IceFreeM','AnNPPrateM','AnAvDayRateM'});
    
    Regression.lm_PixelMy.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'AnNPPrateM~IceFreeM');
    disp(Regression.lm_PixelMy.(region_sublist{rix}))
    Regression.lm_PixeldMy.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'AnAvDayRateM~IceFreeM');
    disp(Regression.lm_PixeldMy.(region_sublist{rix}))
    % Regression.lm_Pixelmedian.(region_sublist{rix})=fitlm(Regression.tbl_pixels_meanmed.(region_sublist{rix}),'AnNPPrateMedian~IceFreeMedian');
    %  disp(Regression.lm_Pixelmedian.(region_sublist{rix}))
end

figure;
t = tiledlayout(2,2)
for rix=1:length(region_sublist)
nexttile
plot(Regression.lm_PixelMy.(region_sublist{rix}))
    l=legend;
    l.Location='southeast';
    title((region_sublist{rix}))
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
% nexttile
% plot(Regression.lm_Pixelmean.(region_sublist{rix}))
% nexttile
% plot(Regression.lm_Pixelmedian.(region_sublist{2}))
end

figure;
t = tiledlayout(2,2)
for rix=1:length(region_sublist)
nexttile
plot(Regression.lm_PixeldMy.(region_sublist{rix}))
    l=legend;
    l.Location='northeast';
    title((region_sublist{rix}))
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Mean daily NPP (g m^{-2} d^{-1})','Interpreter','tex')
% nexttile
% plot(Regression.lm_Pixelmean.(region_sublist{rix}))
% nexttile
% plot(Regression.lm_Pixelmedian.(region_sublist{2}))
end

%% Mean IFE vs #IFd
for rix=1:length(region_sublist)
    Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_MEANSy.(region_sublist{rix}),...
        Regression.MeanIFE.(region_sublist{rix}),...
        'VariableNames',{'IceFreeM','MeanIFE'});
    
    Regression.lm_.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'MeanIFE~IceFreeM');
    disp(Regression.lm_.(region_sublist{rix}))
end
figure;
t = tiledlayout(2,2)
for rix=1:length(region_sublist)
nexttile
plot(Regression.lm_.(region_sublist{rix}))
    l=legend;
    l.Location='southeast';
    title((region_sublist{rix}))
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Mean Ice Free Extent (10^6 km^2)','Interpreter','tex')
% nexttile
% plot(Regression.lm_Pixelmean.(region_sublist{rix}))
% nexttile
% plot(Regression.lm_Pixelmedian.(region_sublist{2}))
end
for rix=1:length(region_sublist)
    [rho,pval] = corr(IceFree_pixels_MEANSy.(region_sublist{rix}),...
        Regression.MeanIFE.(region_sublist{rix}), 'type', 'Spearman')
end

%% Boxplots

figure;
tiledlayout(2,6)
nexttile(5,[1 2])
rix=1;
boxplot(IceFree_pixels_years.(region_sublist{rix}))
    title((region_sublist{rix}))
nexttile(1,[1 1])
boxplot(IceFree_pixels_years_COL.(region_sublist{rix}))
ylim([0 250])
    title({'All' (region_sublist{rix})})
    
nexttile(7,[1 2])
rix=2;
boxplot(IceFree_pixels_years.(region_sublist{rix}))
    title((region_sublist{rix}))
nexttile(2,[1 1])
boxplot(IceFree_pixels_years_COL.(region_sublist{rix}))
ylim([0 250])
    title({'All' (region_sublist{rix})})

    nexttile(9,[1 2])
    rix=3;
boxplot(IceFree_pixels_years.(region_sublist{rix}))
    title((region_sublist{rix}))
nexttile(3,[1 1])
boxplot(IceFree_pixels_years_COL.(region_sublist{rix}))
ylim([0 250])
    title({'All' (region_sublist{rix})})

    nexttile(11,[1 2])
    rix=4;
boxplot(IceFree_pixels_years.(region_sublist{rix}))
    title((region_sublist{rix}))
nexttile(4,[1 1])
boxplot(IceFree_pixels_years_COL.(region_sublist{rix}))
ylim([0 250])
    title({'All' (region_sublist{rix})})

%% MEANS... 
for rix = 1:length(region_sublist)
clearvars find*
IceFree_pixels_years.(region_sublist{rix})(IceFree_pixels_years.(region_sublist{rix})==0)=NaN;
   
    IceFree_pixels_years_mean.(region_sublist{rix})=mean(IceFree_pixels_years.(region_sublist{rix}),2,'omitnan');
%     IceFree_pixels_years_median.(region_sublist{rix})=median(IceFree_pixels_years.(region_sublist{rix}),2,'omitnan');
    AnnualNPPRate_pixels_years_mean.(region_sublist{rix})=mean(AnnualNPPRate_pixels_years.(region_sublist{rix}),2,'omitnan');
%     AnnualNPPRate_pixels_years_median.(region_sublist{rix})=median(AnnualNPPRate_pixels_years.(region_sublist{rix}),2,'omitnan');
        AnnualNPPRate_pixels_years_mean.(region_sublist{rix})=mean(AnnualNPPRate_pixels_years.(region_sublist{rix}),2,'omitnan');

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





%% Shelf Regressions
rix=2;
Regression.lmShelfall.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFE+SIE_dif+IFDav+Year+IFE*Year+SIE_dif*Year+IFDav*Year');
disp(Regression.lmShelfall.(region_sublist{rix}))

Regression.lmShelf_lessIFE.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~SIE_dif+IFDav+Year+SIE_dif*Year+IFDav*Year');
disp(Regression.lmShelf_lessIFE.(region_sublist{rix}))

Regression.lmShelfall_noint.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFE+SIE_dif+IFDav+Year');
disp(Regression.lmShelfall_noint.(region_sublist{rix}))

Regression.lmShelf_YIFD.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFDav+Year+IFDav*Year');
disp(Regression.lmShelf_YIFD.(region_sublist{rix}))
Regression.lmShelf_Ydif.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~SIE_dif+Year+SIE_dif*Year');
disp(Regression.lmShelf_Ydif.(region_sublist{rix}))

Regression.lmShelf_noyr.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFE+SIE_dif+IFDav');
disp(Regression.lmShelf_noyr.(region_sublist{rix}))
Regression.lmShelf_noyrd.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFE+IFDav');
disp(Regression.lmShelf_noyrd.(region_sublist{rix}))

Regression.lmShelf_year.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~Year');
disp(Regression.lmShelf_year.(region_sublist{rix}))
Regression.lmShelf_IFE.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFE');
disp(Regression.lmShelf_IFE.(region_sublist{rix}))
Regression.lmShelf_IFD.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~IFDav');
disp(Regression.lmShelf_IFD.(region_sublist{rix}))
Regression.lmShelf_SIE_dif.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~SIE_dif');
disp(Regression.lmShelf_SIE_dif.(region_sublist{rix}))

Regression.lmShelf_IFE_sqrt.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP_sqrt~IFE');
disp(Regression.lmShelf_IFE_sqrt.(region_sublist{rix}))

model_list={'lmShelfall','lmShelf_lessIFE','lmShelfall_noint','lmShelf_YIFD','lmShelf_Ydif','lmShelf_noyr','lmShelf_noyrd','lmShelf_year','lmShelf_IFE','lmShelf_IFD','lmShelf_SIE_dif'};

%       AIC
for rix = 2%1:length(region_sublist)
    for mix=1:length(model_list)
        disp([(region_sublist{rix}) (model_list{mix})])
        AIC(mix,1)=Regression.(model_list{mix}).(region_sublist{rix}).ModelCriterion.AIC;
    end
end


close(figure(22))
figure(22);
tiledlayout(2,2)
anpos=[0.1, 0.8, 0.1, 0.1;0.55, 0.8, 0.1, 0.1;0.1, 0.3, 0.1, 0.1;0.55, 0.3, 0.1, 0.1];
for rix = 2%1:length(region_sublist)
    nexttile
    plot(Regression.lmShelf_Ydif.(region_sublist{rix})) %lmNPPIFA %lmRateIFE %lmGSRateIFA %lmRateIFD
    nexttile
    scatter(Regression.Year,Regression.NPP_AnTot.(region_sublist{rix}))
    nexttile
    scatter(Regression.SIE_dif.(region_sublist{rix}),Regression.NPP_AnTot.(region_sublist{rix}))
    nexttile
    scatter(Regression.Year,Regression.SIE_dif.(region_sublist{rix}))
%     title((region_sublist{rix}))
%     xtxt={'Mean Ice Free Extent (10^6 km^2)'}; %Mean Ice Free Area (10^6 km^2) %Average # Ice Free Days
%     xlabel(xtxt,'Interpreter','tex')
%     ylabel('Annual NPP (TgC)','Interpreter','tex') %Annual NPP (TgC) %Annual NPP (gC m^2 a^{-1}) %Growing Season NPP (mg m^{-2} d^{-1})
% legend('off')
%     %     l=legend;
% %     l.Location='southeast';
%     if rix==2
%         ylim([-0.2 inf])
%         yline(0,':k')
%     end
%     c=Regression.lmNPPIFE.(region_sublist{rix}).Coefficients{1,1};
%     m=Regression.lmNPPIFE.(region_sublist{rix}).Coefficients{2,1};
%     R2=Regression.lmNPPIFE.(region_sublist{rix}).Rsquared.Adjusted;
%     p=Regression.lmNPPIFE.(region_sublist{rix}).Coefficients{2,4};
%     str={'R^2=' num2str(R2),' p=' num2str(p);'NPP=' num2str(m) '*IFE+' num2str(c)};
%     str2={strjoin(str(1:2:8));strjoin(str(2:2:8))};
%     annotation('textbox', anpos(rix,:),'String',str2,'FitBoxToText','on')
end

% residuals
figure;
plotResiduals(Regression.lmShelf_noyr.(region_sublist{rix}))
figure; plotDiagnostics(Regression.lmShelf_noyr.(region_sublist{rix}),'cookd')
legend('show') % Show the legend
resid=table2array(Regression.lmShelf_noyr.(region_sublist{rix}).Residuals(:,1));
lev = leverage(resid);
fitted=Regression.lmShelf_noyr.(region_sublist{rix}).Fitted;
figure; tiledlayout('flow'); 
nexttile; scatter(fitted,resid); title('Residuals vs Fitted')
nexttile; qqplot(resid); title('Normal Q-Q') 
nexttile; histogram(resid,10); title('Distribution of residuals')
nexttile; scatter(lev,resid); title('Residuals vs Leverage')

figure;
plotResiduals(Regression.lmShelf_noyrd.(region_sublist{rix}))
figure; plotDiagnostics(Regression.lmShelf_noyrd.(region_sublist{rix}),'cookd')
legend('show') % Show the legend
resid=table2array(Regression.lmShelf_noyrd.(region_sublist{rix}).Residuals(:,1));
lev = leverage(resid);
fitted=Regression.lmShelf_noyrd.(region_sublist{rix}).Fitted;
figure; tiledlayout('flow'); 
nexttile; scatter(fitted,resid); title('Residuals vs Fitted')
nexttile; qqplot(resid); title('Normal Q-Q') 
nexttile; histogram(resid,10); title('Distribution of residuals')
nexttile; scatter(lev,resid); title('Residuals vs Leverage')


figure;
plotResiduals(Regression.lmShelf_IFE.(region_sublist{rix}))
figure; plotDiagnostics(Regression.lmShelf_IFE.(region_sublist{rix}),'cookd')
legend('show') % Show the legend
resid=table2array(Regression.lmShelf_IFE.(region_sublist{rix}).Residuals(:,1));
lev = leverage(resid);
fitted=Regression.lmShelf_IFE.(region_sublist{rix}).Fitted;
figure; tiledlayout('flow'); 
nexttile; scatter(fitted,resid); title('Residuals vs Fitted')
nexttile; qqplot(resid); title('Normal Q-Q') 
nexttile; histogram(resid,10); title('Distribution of residuals')
nexttile; scatter(lev,resid); title('Residuals vs Leverage')


figure;
plotResiduals(Regression.lmShelf_IFE_sqrt.(region_sublist{rix}))
figure; plotDiagnostics(Regression.lmShelf_IFE_sqrt.(region_sublist{rix}),'cookd')
legend('show') % Show the legend
resid=table2array(Regression.lmShelf_IFE_sqrt.(region_sublist{rix}).Residuals(:,1));
lev = leverage(resid);
fitted=Regression.lmShelf_IFE_sqrt.(region_sublist{rix}).Fitted;
figure; tiledlayout('flow'); 
nexttile; scatter(fitted,resid); title('Residuals vs Fitted')
nexttile; qqplot(resid); title('Normal Q-Q') 
nexttile; histogram(resid,5); title('Distribution of residuals')
nexttile; scatter(lev,resid); title('Residuals vs Leverage')

%% Shelf 2006 onwards
for rix = 1:length(region_sublist)
Regression.tbl_2007on.(region_sublist{rix})=Regression.tbl.(region_sublist{rix})(5:end,:);
end

figure; plot(Regression.Year(5:end),Regression.NPP_AnTot.Shelf(5:end))
Regression.lmSR_YvTot.Shelf=fitlm(Regression.tbl_2007on.Shelf,'NPP~Year');
disp(Regression.lmSR_YvTot.Shelf)

figure(5);
    plot(Regression.lmSR_YvTot.Shelf) %lmNPPIFA %lmRateIFE %lmGSRateIFA %lmRateIFD
    title('Shelf 2006 onwards')
    xtxt={'Year'}; %Mean Ice Free Area (10^6 km^2) %Average # Ice Free Days
    xlabel(xtxt,'Interpreter','tex')
    ylabel('Annual NPP (TgC)','Interpreter','tex') %Annual NPP (TgC) %Annual NPP (gC m^2 a^{-1}) %Growing Season NPP (mg m^{-2} d^{-1})
legend('off')
    %     l=legend;
%     l.Location='southeast';
ylim([-0.2 inf])
yline(0,':k')

figure;
plotResiduals(Regression.lmSR_YvTot.Shelf)
figure; plotDiagnostics(Regression.lmSR_YvTot.Shelf,'cookd')
legend('show') % Show the legend
resid=table2array(Regression.lmSR_YvTot.Shelf.Residuals(:,1));
lev = leverage(resid);
fitted=Regression.lmSR_YvTot.Shelf.Fitted;
figure; tiledlayout('flow'); 
nexttile; scatter(fitted,resid); title('Residuals vs Fitted')
nexttile; qqplot(resid); title('Normal Q-Q') 
nexttile; histogram(resid,10); title('Distribution of residuals')
nexttile; scatter(lev,resid); title('Residuals vs Leverage')

