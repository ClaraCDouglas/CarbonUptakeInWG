%Import annual SAM values from E:\Data\Wind\SAM

SAM0220=Annual20022020(:,2);
SAM0220=table2array(SAM0220);

% Since annual SAM values are provided for calendar year:
    % could do regressions of total annual NPP vs SAM of preceding and
    % current year
        %i.e. 2003 austral total NPP ~ 2002 SAM
        % and 2003 austral total NPP ~ 2003 SAM
        % my predition is that preceding year's SAM will be more important
    % maybe want to get mean of SON-JJA or JJA-MAM SAM values too, to
    % also try using the dodge way to calc austral SAM
    
    
    


for rix = 1:length(region_sublist)
Regression.tblSAM.(region_sublist{rix})=table(Regression.MeanSIE.(region_sublist{rix})(1:18),...
    Regression.MeanIFE.(region_sublist{rix})(1:18),...
    Regression.SIE_dif.(region_sublist{rix})(1:18),...
    Regression.NPP_AnTot.(region_sublist{rix})(1:18),...
    Regression.NPP_AnRate.(region_sublist{rix})(1:18),...
    Regression.NPP_AvGSRate.(region_sublist{rix})(1:18),...
    SAM0220(2:19),...
    'VariableNames',{'SIE','IFE','SIE_dif','NPP','NPPrate','NPPGSrate','SAM'}); %

% Regression.lmNPPSIE.(region_sublist{rix})=fitlm(Regression.tblSAM.(region_sublist{rix}),'NPP~SIE');
% disp(Regression.lmNPPSIE.(region_sublist{rix}))
Regression.lmNPPrSIE.(region_sublist{rix})=fitlm(Regression.tblSAM.(region_sublist{rix}),'NPPrate~SIE');
% disp(Regression.lmNPPrSIE.(region_sublist{rix}))

% Regression.lmNPPSAM.(region_sublist{rix})=fitlm(Regression.tblSAM.(region_sublist{rix}),'NPP~SAM');
% disp(Regression.lmNPPSAM.(region_sublist{rix}))

% Regression.lmSIESAM.(region_sublist{rix})=fitlm(Regression.tblSAM.(region_sublist{rix}),'SIE~SAM');
% disp(Regression.lmSIESAM.(region_sublist{rix}))

% Regression.lmNPPrSAM.(region_sublist{rix})=fitlm(Regression.tblSAM.(region_sublist{rix}),'NPPrate~SAM');
% disp(Regression.lmNPPrSAM.(region_sublist{rix}))
% 
% Regression.lmNPPSIESAM.(region_sublist{rix})=fitlm(Regression.tblSAM.(region_sublist{rix}),'NPP~SIE+SAM');
% disp(Regression.lmNPPSIESAM.(region_sublist{rix}))
% Regression.lmNPPrSIESAM.(region_sublist{rix})=fitlm(Regression.tblSAM.(region_sublist{rix}),'NPPrate~SIE+SAM');
% disp(Regression.lmNPPrSIESAM.(region_sublist{rix}))

%% SIdif
% Regression.lmNPPSIEd.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPP~SIE_dif');
% disp(Regression.lmNPPSIEd.(region_sublist{rix}))
Regression.lmNPPrSIEd.(region_sublist{rix})=fitlm(Regression.tbl.(region_sublist{rix}),'NPPrate~SIE_dif');
% disp(Regression.lmNPPrSIEd.(region_sublist{rix}))

% Regression.lmSIEdSAM.(region_sublist{rix})=fitlm(Regression.tblSAM.(region_sublist{rix}),'SIE_dif~SAM');
% disp(Regression.lmSIEdSAM.(region_sublist{rix}))
 
% Regression.lmNPPSIEdSAM.(region_sublist{rix})=fitlm(Regression.tblSAM.(region_sublist{rix}),'NPP~SIE_dif+SAM');
% disp(Regression.lmNPPSIEdSAM.(region_sublist{rix}))
Regression.lmNPPrSIEdSAM.(region_sublist{rix})=fitlm(Regression.tblSAM.(region_sublist{rix}),'NPPrate~SIE_dif+SAM');
disp(Regression.lmNPPrSIEdSAM.(region_sublist{rix}))

end
%%
figure;
tiledlayout(2,2)
anpos=[0.1, 0.8, 0.1, 0.1;0.55, 0.8, 0.1, 0.1;0.1, 0.3, 0.1, 0.1;0.55, 0.3, 0.1, 0.1];
for rix = 1:length(region_sublist)
    nexttile
    plot(Regression.lmNPPSIEd.(region_sublist{rix}))
    title((region_sublist{rix}))
    xtxt={'SIE difference (Winter Max SIE - Summer Min SIE)'}; %Calendar Year SAM Index
    xlabel(xtxt,'Interpreter','tex')
    ylabel('Annual NPP (TgC)','Interpreter','tex') %Annual NPP (TgC) %Sea Ice Extent (x10^6 km^2) % Annual NPP (gC m^2 a^{-1})
    legend('off')
    %     l=legend;
%     l.Location='southeast';
    if rix==2 || rix==4
        ylim([-0.2 inf])
        yline(0,':k')
     end
    c=Regression.lmNPPSIEd.(region_sublist{rix}).Coefficients{1,1};
    m=Regression.lmNPPSIEd.(region_sublist{rix}).Coefficients{2,1};
    R2=Regression.lmNPPSIEd.(region_sublist{rix}).Rsquared.Adjusted;
    p=Regression.lmNPPSIEd.(region_sublist{rix}).Coefficients{2,4};
    str={'R^2=' num2str(R2),' p=' num2str(p);'NPP=' num2str(m) '*IFE+' num2str(c)};
    str2={strjoin(str(1:2:8));strjoin(str(2:2:8))};
    annotation('textbox', anpos(rix,:),'String',str2,'FitBoxToText','on')
end

%% AIC
model_list={'lmNPPSIE','lmNPPrSIE','lmNPPSAM',...
    'lmNPPrSAM','lmNPPSIESAM','lmNPPrSIESAM','lmNPPSIEd','lmSIEdSAM',...
    'lmNPPSIEdSAM','lmNPPrSIEdSAM'};
NPPmodel_list={'lmNPPSIE','lmNPPSAM',...
    'lmNPPSIESAM','lmNPPSIEd',...
    'lmNPPSIEdSAM'};
Ratemodel_list={'lmNPPrSIE',...
    'lmNPPrSAM','lmNPPrSIESAM',...
    'lmNPPrSIEdSAM'};

%       AIC
for rix = 1:length(region_sublist)
    for mix=1:length(Ratemodel_list)
        disp([(region_sublist{rix}) (Ratemodel_list{mix})])
        Regression.(Ratemodel_list{mix}).(region_sublist{rix}).ModelCriterion.AIC
    end
end

