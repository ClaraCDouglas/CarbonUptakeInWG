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
algorithm='cafe'
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
glmFIT='polynomial3'%'quadratic'%'linear'
for rix=1:length(region_sublist)
    switch glmFIT
        case 'linear'
            Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_years_COL.(region_sublist{rix}),...
                AnnualNPPRate_pixels_years_COL.(region_sublist{rix}),categorical(YEARS_COL.(region_sublist{rix})),...
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
        case 'quadratic'
            Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_years_COL.(region_sublist{rix}),...
                AnnualNPPRate_pixels_years_COL.(region_sublist{rix}),...
                'VariableNames',{'IceFree','AnNPPrate'}); %AnAvDayRate_pixels_years_COL.(region_sublist{rix}),...,'AnAvDayRate'
            Regression.lm_Pixel.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'purequadratic');
            disp(Regression.lm_Pixel.(region_sublist{rix}))
        case 'polynomial3'
            Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_years_COL.(region_sublist{rix}),...
                AnnualNPPRate_pixels_years_COL.(region_sublist{rix}),categorical(YEARS_COL.(region_sublist{rix})),...
                'VariableNames',{'IceFree','AnNPPrate','YEAR'}); %AnAvDayRate_pixels_years_COL.(region_sublist{rix}),...,'AnAvDayRate'           
            Regression.lm_Pixel.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'poly12');
            disp(Regression.lm_Pixel.(region_sublist{rix}))
    end
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

%% hist2d
figure;
tiledlayout('flow')
for rix=1:length(region_sublist)
nexttile
Yedges = [0:5:(max(AnnualNPPRate_pixels_years_COL.(region_sublist{rix}))+10)];
Xedges = [0:10:(max(IceFree_pixels_years_COL.(region_sublist{rix}))+10)];
h=histogram2(IceFree_pixels_years_COL.(region_sublist{rix}),AnnualNPPRate_pixels_years_COL.(region_sublist{rix}),...
    Xedges,Yedges,'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','off') %lmNPPIFA %lmRateIFA %lmGSRateIFA
h.EdgeColor='None';
colormap(cmocean('matter'));
colorbar
xlim([0 Inf])
ylim([0 Inf])
xlabel('Number of ice free days (number of days NPP data is available)')
ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
title((region_sublist{rix}))
end

%calc mean for each IFD bin
for rix=1:length(region_sublist)
z=IceFree_pixels_years_COL.(region_sublist{rix});
setedges=[0:10:(max(IceFree_pixels_years_COL.(region_sublist{rix}))+10)];
[count_in_bin,edges,bin_IFD]=histcounts(z,setedges);

w = AnnualNPPRate_pixels_years_COL.(region_sublist{rix});
IFDbinMeanNPP.(region_sublist{rix})=nan(size(count_in_bin));
midbin.(region_sublist{rix})=nan(size(count_in_bin));
for ix=1:length(count_in_bin)
    IFDbinMeanNPP.(region_sublist{rix})(ix)=(nanmean(w(bin_IFD == ix)));
    midbin.(region_sublist{rix})(ix)=mean([edges(ix),edges(ix+1)]);
end

figure;
plot(midbin.(region_sublist{rix}),IFDbinMeanNPP.(region_sublist{rix}),'-o')
end

% together
figure;
tiledlayout('flow')
for rix=1:length(region_sublist)
nexttile
Yedges = [0:5:(max(AnnualNPPRate_pixels_years_COL.(region_sublist{rix}))+10)];
Xedges = [0:10:(max(IceFree_pixels_years_COL.(region_sublist{rix}))+10)];
h=histogram2(IceFree_pixels_years_COL.(region_sublist{rix}),AnnualNPPRate_pixels_years_COL.(region_sublist{rix}),...
    Xedges,Yedges,'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','off') %lmNPPIFA %lmRateIFA %lmGSRateIFA
h.EdgeColor='None';
colormap(cmocean('matter'));
colorbar
xlim([0 Inf])
ylim([0 Inf])
xlabel('Number of ice free days (number of days NPP data is available)')
ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
hold on
plot(midbin.(region_sublist{rix}),IFDbinMeanNPP.(region_sublist{rix}),'k-o')
title((region_sublist{rix}))
end
txt={'Density distribution of IFD and Area Normalised Annual NPP','Black circles represent mean NPP in each IFD bin'};
sgtitle(txt)

%% categorical years
for rix=1:length(region_sublist)
    figure;
    tiledlayout('flow')
    for yix=2003:2021
        Regression.tbl_pixelsYEARsel.(region_sublist{rix})=table(IceFree_pixels_years_COL.(region_sublist{rix})(YEARS_COL.(region_sublist{rix})==yix),...
            AnnualNPPRate_pixels_years_COL.(region_sublist{rix})(YEARS_COL.(region_sublist{rix})==yix),...
            ... %categorical(YEARS_COL.(region_sublist{rix})(YEARS_COL.(region_sublist{rix})==yix)),
            'VariableNames',{'IceFree','AnNPPrate'}); %,'YEAR',AnAvDayRate_pixels_years_COL.(region_sublist{rix}),...,'AnAvDayRate'
        Regression.lm_PixelCAT.(region_sublist{rix})=fitlm(Regression.tbl_pixelsYEARsel.(region_sublist{rix}),'purequadratic'); %,'AnNPPrate~IceFree'
        disp(Regression.lm_PixelCAT.(region_sublist{rix}))
        
        nexttile
        plot(Regression.lm_PixelCAT.(region_sublist{rix}))
        title(yix)
        legend('off')
        sgtitle(region_sublist{rix})
        ylim([0 Inf])
        xlim([0 Inf])
    end
end

%curvefitting tool says polynomial has best R2 and RMSE:
CurveFits = struct();
makeplot=0;
for rix=1:length(region_sublist)
figure; tiledlayout('flow')
    for yix=2003:2021
        x_IFD=IceFree_pixels_years_COL.(region_sublist{rix})(YEARS_COL.(region_sublist{rix})==yix);
        y_NPP=AnnualNPPRate_pixels_years_COL.(region_sublist{rix})(YEARS_COL.(region_sublist{rix})==yix);
        
        fn = sprintf('n%d', yix );
        % A.(fn) = 2; % use the struct
        [CurveFits.(region_sublist{rix}).(fn).p2fitobject,CurveFits.(region_sublist{rix}).(fn).p2gof] = fit(x_IFD,y_NPP,'poly2');
        if makeplot
        nexttile; plot(CurveFits.(region_sublist{rix}).(fn).p2fitobject,x_IFD,y_NPP); 
        title({['poly2 - ',num2str(yix)];...
            ['R^2=' num2str(round(CurveFits.(region_sublist{rix}).(fn).p2gof.adjrsquare,2)),...
            ' RMSE=' num2str(round(CurveFits.(region_sublist{rix}).(fn).p2gof.rmse,2))]}); legend('off');
        end
        [CurveFits.(region_sublist{rix}).(fn).p1fitobject,CurveFits.(region_sublist{rix}).(fn).p1gof] = fit(x_IFD,y_NPP,'poly1');
        if makeplot
        nexttile; plot(CurveFits.(region_sublist{rix}).(fn).p1fitobject,x_IFD,y_NPP); 
        title({['poly1 - ',num2str(yix)];...
            ['R^2=' num2str(round(CurveFits.(region_sublist{rix}).(fn).p1gof.adjrsquare,2)),...
            ' RMSE=' num2str(round(CurveFits.(region_sublist{rix}).(fn).p1gof.rmse,2))]}); legend('off');
        end
        % for the plots:
        CurveFits.(region_sublist{rix}).(fn).p = polyfit(x_IFD,y_NPP,2);
    end
    sgtitle(region_sublist{rix})
end

% together again
figure;
tiledlayout('flow')
for rix=1:length(region_sublist)
    nexttile
    Yedges = [0:5:(max(AnnualNPPRate_pixels_years_COL.(region_sublist{rix}))+10)];
    Xedges = [0:10:(max(IceFree_pixels_years_COL.(region_sublist{rix}))+10)];
    h=histogram2(IceFree_pixels_years_COL.(region_sublist{rix}),AnnualNPPRate_pixels_years_COL.(region_sublist{rix}),...
        Xedges,Yedges,'Normalization','probability',...
        'DisplayStyle','tile','ShowEmptyBins','off') %lmNPPIFA %lmRateIFA %lmGSRateIFA
    h.EdgeColor='None';
    colormap(cmocean('matter'));
    colorbar
    xlim([0 Inf])
    ylim([0 Inf])
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
    hold on
    title((region_sublist{rix}))
    
    for yix=2003:2021
        x_IFD=IceFree_pixels_years_COL.(region_sublist{rix})(YEARS_COL.(region_sublist{rix})==yix);
        x1 = linspace(0,max(x_IFD));
        fn = sprintf('n%d', yix );
        y1 = polyval(CurveFits.(region_sublist{rix}).(fn).p,x1);
        plot(x1,y1,'k','LineWidth',1)
    end
    plot(midbin.(region_sublist{rix}),IFDbinMeanNPP.(region_sublist{rix}),'b-o','MarkerFaceColor','b')
    %axis=gca; uistack(gca,'top');
    set(gca,'Layer','top');
end
txt={'Density distribution of IFD and Area Normalised Annual NPP','Blue circles represent mean NPP in each IFD bin','Black lines are the annual regressions'};
sgtitle(txt)
    
    
%% means regressions by pixel
for rix = 1:length(region_sublist)
    IceFree_pixels_yearsN.(region_sublist{rix})=IceFree_pixels_years.(region_sublist{rix});
    IceFree_pixels_yearsN.(region_sublist{rix})(IceFree_pixels_years.(region_sublist{rix})==0)=NaN;
    IceFree_pixels_MEANS.(region_sublist{rix})=mean(IceFree_pixels_yearsN.(region_sublist{rix}),2,'omitnan');
    AnnualNPPRate_pixels_MEAN.(region_sublist{rix})=mean(AnnualNPPRate_pixels_years.(region_sublist{rix}),2,'omitnan');

    IceFree_pixels_STD.(region_sublist{rix})=std(IceFree_pixels_yearsN.(region_sublist{rix}),0,2,'omitnan');
    AnnualNPPRate_pixels_STD.(region_sublist{rix})=std(AnnualNPPRate_pixels_years.(region_sublist{rix}),0,2,'omitnan');

end
for rix=1:length(region_sublist)
    
    Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_MEANS.(region_sublist{rix}),...
        AnnualNPPRate_pixels_MEAN.(region_sublist{rix}),...
        'VariableNames',{'IceFreeM','AnNPPrateM'});
    
    Regression.lm_PixelM.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'purequadratic');
    disp(Regression.lm_PixelM.(region_sublist{rix}))
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

for rix=1:length(region_sublist)
    figure;
    t = tiledlayout(2,2)
    yneg=AnnualNPPRate_pixels_STD.(region_sublist{rix});
    ypos=AnnualNPPRate_pixels_STD.(region_sublist{rix});
    xneg=IceFree_pixels_STD.(region_sublist{rix});
    xpos=IceFree_pixels_STD.(region_sublist{rix});
    nexttile
    errorbar(IceFree_pixels_MEANS.(region_sublist{rix}),...
        AnnualNPPRate_pixels_MEAN.(region_sublist{rix}),...
        yneg,'o','MarkerSize',5,...
        'MarkerEdgeColor','red','MarkerFaceColor','red');
    l=legend;
    l.Location='northeast';
    title('NPP STD')
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')

    nexttile
    errorbar(IceFree_pixels_MEANS.(region_sublist{rix}),...
        AnnualNPPRate_pixels_MEAN.(region_sublist{rix}),...
        xneg,'horizontal','o','MarkerSize',5,...
        'MarkerEdgeColor','red','MarkerFaceColor','red');
    l=legend;
    l.Location='northeast';
    title('IFD STD')
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
    
    nexttile(3,[1 2])
    errorbar(IceFree_pixels_MEANS.(region_sublist{rix}),...
        AnnualNPPRate_pixels_MEAN.(region_sublist{rix}),...
        yneg,ypos,xneg,xpos,'o','MarkerSize',5,...
        'MarkerEdgeColor','red','MarkerFaceColor','red');
    l=legend;
    l.Location='northeast';
    title('Both STD')
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
    
    sgtitle({'Pixel mean IFD vs pixel mean NPP with STD',(region_sublist{rix})})
end


%% means regressions by year
for rix = 1:length(region_sublist)
    IceFree_pixels_yearsN.(region_sublist{rix})=IceFree_pixels_years.(region_sublist{rix});
    IceFree_pixels_yearsN.(region_sublist{rix})(IceFree_pixels_years.(region_sublist{rix})==0)=NaN;
    IceFree_pixels_MEANSy.(region_sublist{rix})=permute(mean(IceFree_pixels_yearsN.(region_sublist{rix}),1,'omitnan'),[2 1]);
    AnnualNPPRate_pixels_MEANy.(region_sublist{rix})=permute(mean(AnnualNPPRate_pixels_years.(region_sublist{rix}),1,'omitnan'),[2 1]);
    IceFree_pixels_STDy.(region_sublist{rix})=permute(std(IceFree_pixels_yearsN.(region_sublist{rix}),0,1,'omitnan'),[2 1]);
    AnnualNPPRate_pixels_STDy.(region_sublist{rix})=permute(std(AnnualNPPRate_pixels_years.(region_sublist{rix}),0,1,'omitnan'),[2 1]);
end
for rix=1:length(region_sublist)
    Regression.tbl_pixels.(region_sublist{rix})=table(IceFree_pixels_MEANSy.(region_sublist{rix}),...
        AnnualNPPRate_pixels_MEANy.(region_sublist{rix}),...
        'VariableNames',{'IceFreeM','AnNPPrateM'});
    
    Regression.lm_PixelMy.(region_sublist{rix})=fitlm(Regression.tbl_pixels.(region_sublist{rix}),'purequadratic');
    disp(Regression.lm_PixelMy.(region_sublist{rix}))
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

for rix=1:length(region_sublist)
    figure;
    t = tiledlayout(2,2)
    yneg=AnnualNPPRate_pixels_STDy.(region_sublist{rix});
    ypos=AnnualNPPRate_pixels_STDy.(region_sublist{rix});
    xneg=IceFree_pixels_STDy.(region_sublist{rix});
    xpos=IceFree_pixels_STDy.(region_sublist{rix});
    nexttile
    errorbar(IceFree_pixels_MEANSy.(region_sublist{rix}),...
        AnnualNPPRate_pixels_MEANy.(region_sublist{rix}),...
        yneg,'o','MarkerSize',5,...
        'MarkerEdgeColor','red','MarkerFaceColor','red');
    l=legend;
    l.Location='northeast';
    title('NPP STD')
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')

    nexttile
    errorbar(IceFree_pixels_MEANSy.(region_sublist{rix}),...
        AnnualNPPRate_pixels_MEANy.(region_sublist{rix}),...
        xneg,'horizontal','o','MarkerSize',5,...
        'MarkerEdgeColor','red','MarkerFaceColor','red');
    l=legend;
    l.Location='northeast';
    title('IFD STD')
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
    
    nexttile(3,[1 2])
    errorbar(IceFree_pixels_MEANSy.(region_sublist{rix}),...
        AnnualNPPRate_pixels_MEANy.(region_sublist{rix}),...
        yneg,ypos,xneg,xpos,'o','MarkerSize',5,...
        'MarkerEdgeColor','red','MarkerFaceColor','red');
    l=legend;
    l.Location='northeast';
    title('Both STD')
    xlabel('Number of ice free days (number of days NPP data is available)')
    ylabel('Annual NPP (g m^{-2})','Interpreter','tex')
    
    sgtitle({'Year mean IFD vs Year mean NPP with STD',(region_sublist{rix})})
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

