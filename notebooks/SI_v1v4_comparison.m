%% Compare SIv1 and SIv4

for rix = 1:length(region_sublist)
    
    SeaIcet1.(region_sublist{rix}).SIE_ausMAX=max(SeaIce.(region_sublist{rix}).SIExtent_aus,[],1,'omitnan');
    SeaIcet1.(region_sublist{rix}).SIE_ausMIN=min(SeaIce.(region_sublist{rix}).SIExtent_aus,[],1,'omitnan');
    SeaIcet1.(region_sublist{rix}).SIE_ausDIF=(SeaIcet1.(region_sublist{rix}).SIE_ausMAX-SeaIcet1.(region_sublist{rix}).SIE_ausMIN)';
    SeaIcet4.(region_sublist{rix}).SIE_ausMAX=max(SeaIcev4.(region_sublist{rix}).SIExtent_aus,[],1,'omitnan');
    SeaIcet4.(region_sublist{rix}).SIE_ausMIN=min(SeaIcev4.(region_sublist{rix}).SIExtent_aus,[],1,'omitnan');
    SeaIcet4.(region_sublist{rix}).SIE_ausDIF=(SeaIcet4.(region_sublist{rix}).SIE_ausMAX-SeaIcet4.(region_sublist{rix}).SIE_ausMIN)';

Regression.MeanSIEv1.(region_sublist{rix})=mean(SeaIce.(region_sublist{rix}).SIExtent_aus,1,'omitnan');
Regression.MeanSIEv1.(region_sublist{rix})=(Regression.MeanSIEv1.(region_sublist{rix})/1e6)';
Regression.MeanSIAv1.(region_sublist{rix})=mean(SeaIce.(region_sublist{rix}).SIArea_aus,1,'omitnan');
Regression.MeanSIAv1.(region_sublist{rix})=(Regression.MeanSIAv1.(region_sublist{rix})/1e6)';

Regression.SIE_difv1.(region_sublist{rix})=SeaIcet1.(region_sublist{rix}).SIE_ausDIF;

Regression.MeanSIEv4.(region_sublist{rix})=mean(SeaIcev4.(region_sublist{rix}).SIExtent_aus,1,'omitnan');
Regression.MeanSIEv4.(region_sublist{rix})=(Regression.MeanSIEv4.(region_sublist{rix})/1e6)';
Regression.MeanSIAv4.(region_sublist{rix})=mean(SeaIcev4.(region_sublist{rix}).SIArea_aus,1,'omitnan');
Regression.MeanSIAv4.(region_sublist{rix})=(Regression.MeanSIAv4.(region_sublist{rix})/1e6)';

Regression.SIE_difv4.(region_sublist{rix})=SeaIcet4.(region_sublist{rix}).SIE_ausDIF;
end



for rix = 1:length(region_sublist)
Regression.tblSIE.(region_sublist{rix})=table(Regression.MeanSIEv1.(region_sublist{rix})(1:18),...
    Regression.MeanSIEv4.(region_sublist{rix})(1:18),...
    Regression.MeanSIAv1.(region_sublist{rix})(1:18),...
    Regression.MeanSIAv4.(region_sublist{rix})(1:18),...
    Regression.SIE_difv1.(region_sublist{rix})(1:18),...
    Regression.SIE_difv4.(region_sublist{rix})(1:18),...
    'VariableNames',{'SIEv1','SIEv4','SIAv1','SIAv4','SIE_difv1','SIE_difv4'}); %
% Regression.lmSIE.(region_sublist{rix})=fitlm(Regression.tblSIE.(region_sublist{rix}),'SIEv1~SIEv4');
% disp(Regression.lmSIE.(region_sublist{rix}))
% Regression.lmSIA.(region_sublist{rix})=fitlm(Regression.tblSIE.(region_sublist{rix}),'SIAv1~SIAv4');
% disp(Regression.lmSIA.(region_sublist{rix}))
Regression.lmSIE_dif.(region_sublist{rix})=fitlm(Regression.tblSIE.(region_sublist{rix}),'SIE_difv1~SIE_difv4');
disp(Regression.lmSIE_dif.(region_sublist{rix}))

end


figure;
tiledlayout(2,2)
anpos=[0.1, 0.8, 0.1, 0.1;0.55, 0.8, 0.1, 0.1;0.1, 0.3, 0.1, 0.1;0.55, 0.3, 0.1, 0.1];
for rix = 1:length(region_sublist)
    nexttile
    plot(Regression.lmSIE_dif.(region_sublist{rix}))
    title((region_sublist{rix}))
    xtxt={'SIE_difv1'}; %Calendar Year SAM Index
    xlabel(xtxt,'Interpreter','tex')
    ylabel('SIE_difv4','Interpreter','tex') %Annual NPP (TgC) %Sea Ice Extent (x10^6 km^2) % Annual NPP (gC m^2 a^{-1})
    legend('off')
    %     l=legend;
%     l.Location='southeast';
    c=Regression.lmSIE_dif.(region_sublist{rix}).Coefficients{1,1};
    m=Regression.lmSIE_dif.(region_sublist{rix}).Coefficients{2,1};
    R2=Regression.lmSIE_dif.(region_sublist{rix}).Rsquared.Adjusted;
    p=Regression.lmSIE_dif.(region_sublist{rix}).Coefficients{2,4};
    str={'R^2=' num2str(R2),' p=' num2str(p);'NPP=' num2str(m) '*IFE+' num2str(c)};
    str2={strjoin(str(1:2:8));strjoin(str(2:2:8))};
    annotation('textbox', anpos(rix,:),'String',str2,'FitBoxToText','on')
end
sgtitle('SIE difference v1 vs v4')


%% SIE vs SIA v4
for rix = 1:length(region_sublist)
Regression.tblSIE.(region_sublist{rix})=table(Regression.MeanSIEv1.(region_sublist{rix})(1:18),...
    Regression.MeanSIEv4.(region_sublist{rix})(1:18),...
    Regression.MeanSIAv1.(region_sublist{rix})(1:18),...
    Regression.MeanSIAv4.(region_sublist{rix})(1:18),...
    Regression.SIE_difv1.(region_sublist{rix})(1:18),...
    Regression.SIE_difv4.(region_sublist{rix})(1:18),...
    'VariableNames',{'SIEv1','SIEv4','SIAv1','SIAv4','SIE_difv1','SIE_difv4'}); %
Regression.lmSIESIA.(region_sublist{rix})=fitlm(Regression.tblSIE.(region_sublist{rix}),'SIEv4~SIAv4');
disp(Regression.lmSIESIA.(region_sublist{rix}))
end


figure;
tiledlayout(2,2)
anpos=[0.1, 0.8, 0.1, 0.1;0.55, 0.8, 0.1, 0.1;0.1, 0.3, 0.1, 0.1;0.55, 0.3, 0.1, 0.1];
for rix = 1:length(region_sublist)
    nexttile
    plot(Regression.lmSIESIA.(region_sublist{rix}))
    title((region_sublist{rix}))
    xtxt={'SIE'}; %Calendar Year SAM Index
    xlabel(xtxt,'Interpreter','tex')
    ylabel('SIA','Interpreter','tex') %Annual NPP (TgC) %Sea Ice Extent (x10^6 km^2) % Annual NPP (gC m^2 a^{-1})
    legend('off')
    %     l=legend;
%     l.Location='southeast';
    c=Regression.lmSIESIA.(region_sublist{rix}).Coefficients{1,1};
    m=Regression.lmSIESIA.(region_sublist{rix}).Coefficients{2,1};
    R2=Regression.lmSIESIA.(region_sublist{rix}).Rsquared.Adjusted;
    p=Regression.lmSIESIA.(region_sublist{rix}).Coefficients{2,4};
    str={'R^2=' num2str(R2),' p=' num2str(p);'SIA=' num2str(m) '*SIE+' num2str(c)};
    str2={strjoin(str(1:2:8));strjoin(str(2:2:8))};
    annotation('textbox', anpos(rix,:),'String',str2,'FitBoxToText','on')
end
sgtitle('SIE v4 vs SIA v4')