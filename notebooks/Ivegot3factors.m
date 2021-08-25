clearvars
desktop = 1; 
if desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg'))
    cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Matlab Add-ins'));
else
    cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
end
% load('vgpm_imported.mat', 'vgpm_npp_all') % vgpm_npp_all as from OceanProductivity site - average (A-W) daily rates per month
% load('cafe_imported.mat', 'cafe_npp_all') % cafe_npp_all as from OceanProductivity site - average (A-W) daily rates per month
% load('cbpm_imported.mat', 'cbpm_npp_all') % cbpm_npp_all as from OceanProductivity site - average (A-W) daily rates per month
% load('eppley_imported.mat', 'eppley_npp_all') % eppley_npp_all as from OceanProductivity site - average (A-W) daily rates per month
% load('vgpm_imported.mat', 'time_start_all')
load('ProcessedData.mat')

algorithm={'cafe','cbpm','eppley','vgpm'};
% temp.cafe.rates=cafe_npp_all;
% temp.cbpm.rates=cafe_npp_all;
% temp.eppley.rates=eppley_npp_all;
% temp.vgpm.rates=vgpm_npp_all;
% clearvars cafe_npp_all cbpm_npp_all eppley_npp_all vgpm_npp_all

%% box logic
%
load('latlon_m.mat')
load('openshelf_coord.mat')
load('box_lat_lons.mat', 'andrex_box')
IN_and=inpolygon(lon_m,lat_m,andrex_box(:,1),andrex_box(:,2));
findweddell=find(IN_and==1);
IN_shelf=inpolygon(lon_m,lat_m,shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2));
findshelf=find(IN_shelf==1);
IN_open=inpolygon(lon_m,lat_m,open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2));
findopen=find(IN_open==1);

region_sublist={'Weddell','Shelf','Open'};
regionfindlist= {'findweddell','findshelf','findopen'};

clearvars temp
for rix = 1:length(region_sublist)
    %box,box logic
    temp.(region_sublist{rix}).box = lat_m;
    temp.(region_sublist{rix}).box(:) = 0;
    eval(['temp.',region_sublist{rix},'.box(',regionfindlist{rix},')=1;']);
    temp.(region_sublist{rix}).box_logic=logical(temp.(region_sublist{rix}).box);
end

%% accounting for the area of regions
for aix = 4%1:length(algorithm)
    for rix = 1:length(region_sublist)
        temp.(region_sublist{rix}).monthly_AW_rate_gmpermonth=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC./OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_m2;
    end
end
        

for aix = 4%1:length(algorithm)
    for rix = 1:length(region_sublist)
        for yix = 2003:2019
            findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
            temp.(region_sublist{rix}).NPP_annual_AWrate(yix-2002,1)=yix;

            if size(findyear)<12
                temp.(region_sublist{rix}).NPP_annual_AWrate(yix-2002,2)=NaN;
            else
                temp.(region_sublist{rix}).NPP_annual_AWrate(yix-2002,2)=nansum(temp.(region_sublist{rix}).monthly_AW_rate_gmpermonth(findyear,1));
            end
        end
    end
end

x=temp.Weddell.NPP_annual_AWrate(:,1)
y(:,1)=temp.Shelf.NPP_annual_AWrate(:,2);
y(:,2)=temp.Open.NPP_annual_AWrate(:,2);
figure;
subplot(2,1,1);
b=bar(x,y);
shelfcolor=[0.9290 0.6940 0.1250];
opencolor=[0.5 0.3 0.8];
b(1).FaceColor=[0.9290 0.6940 0.1250];
b(2).FaceColor=[0.5 0.3 0.8];
ylim([0 65])
xtickangle(45)
set(gca, 'XTick', [2003:1:2019])
legend('Shelf', 'Open Ocean')
title('Annual area-weighted PP')
xlabel('Year')
ylabel('NPP (g C m^-^2 d^-^1)')


x=temp.Weddell.NPP_annual_AWrate(:,1)
y2(:,1)=OceanProd.vgpm.Shelf.NPP_tot_TgC_annual(2:18,2);
y2(:,2)=OceanProd.vgpm.Open.NPP_tot_TgC_annual(2:18,2);
subplot(2,1,2);
b=bar(x,y2);
shelfcolor=[0.9290 0.6940 0.1250];
opencolor=[0.5 0.3 0.8];
b(1).FaceColor=[0.9290 0.6940 0.1250];
b(2).FaceColor=[0.5 0.3 0.8];
xtickangle(45)
set(gca, 'XTick', [2003:1:2019])
legend('Shelf', 'Open Ocean')
title('Annual Total NPP')
xlabel('Year')
ylabel('NPP (Tg C a^-^1)')

% annual average area weighted NPP
for aix = 4%1:length(algorithm)
    for rix = 1:length(region_sublist)
        temp.(region_sublist{rix}).av_ann_tot=nanmean(OceanProd.vgpm.(region_sublist{rix}).NPP_tot_TgC_annual(2:18,2));
    end
end

% stats - are the annual area-weighted rates on the shelf different to
% those in the open ocean? 
    % (null hypothesis is that they are not different)
[p,tbl]=anova1(y)
years=x;
shelfvopen=y;
mdl=fitlm(shelfvopen,years)
figure;
plotResiduals(mdl)
figure;
Res = table2array(mdl.Residuals);
 boxplot(Res)
 figure;
 plotResiduals(mdl,'probability')
 
figure;
histogram(y(:,2));
figure; 
boxplot(y);

[h,p,ci,stats] = ttest2(y(:,1),y(:,2))
[h,p,ci,stats] = ttest2(y2(:,1),y2(:,2))

% regression and correlation between annual mean area of ice-free water vs A-W NPP
% run some in AnalyseData.m

lm_AW_Shelf=fitlm(regres.tblAW.(region_sublist{2}),'NPP_AW~icefree')
lm_AW_Open=fitlm(regres.tblAW.(region_sublist{3}),'NPP_AW~icefree')
lmsmall_AW_Shelf=fitlm(regres.tblAWsmall.(region_sublist{2}),'NPP_AW~icefree')
lmsmall_AW_Open=fitlm(regres.tblAWsmall.(region_sublist{3}),'NPP_AW~icefree')

[rhow,pw]=corr(regres.NPP_AW.(region_sublist{1}),regres.ice.(region_sublist{1}),'Type','Pearson')
[rhoo,po]=corr(regres.NPP_AW.(region_sublist{3}),regres.ice.(region_sublist{3}),'Type','Pearson')
[rhos,ps]=corr(regres.NPP_AW.(region_sublist{2}),regres.ice.(region_sublist{2}),'Type','Pearson')

figure; 
subplot(2,1,1)
plot(lmsmall_AW_Shelf); title('Shelf')
txtS = {'y=0.59x+40.726','rho=0.2172','p=0.402'};
text(2.4,54,txtS)
xlabel('Annual mean area of ice-free water (km^2 x10^4)','Interpreter','tex')
ylabel('Area-weighted NPP (g C m^-^2 a^-^1)','Interpreter','tex')
subplot(2,1,2)
plot(lmsmall_AW_Open); title('Open Ocean')
txtO = {'y=-2.3933x+42.069','rho=-0.0987','p=0.706'};
text(1.1,30,txtO)
xlabel('Annual mean area of ice-free water (km ^2 x10 ^6)','Interpreter','tex')
ylabel('Area-weighted NPP (g C m ^-^2 a^-^1)','Interpreter','tex')
sgtitle('Regression between area normalised/weighted annual NPP and mean ice-free area')


%% length of growing season
NPP_nans=vgpm_npp_all;
NPP_nans(NPP_nans<0)=NaN;
icecover=isnan(NPP_nans);

% for permanantly ice covered - should change any ==12( or 7) to NaN before
% summing
% orrrrr value will = length of time ice covered
    % so will want growing season length
        % so 12 - months ice covered = # months ice free
        % so numbers == 0 can be set to NaN
%
for yix = 2003:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    ice_covered(:,:,yix-2002)=nansum(icecover(:,:,findyear),3);
    growing_season(:,:,yix-2002)=12.-ice_covered(:,:,yix-2002);
end

growing_season(growing_season==0)=NaN;

for rix = 1:length(region_sublist)
    for yix = 2003:2019
        gstemp=growing_season(:,:,yix-2002);
        temp.(region_sublist{rix}).growseason(yix-2002,1)=nanmean(nanmean(gstemp(temp.(region_sublist{rix}).box_logic)));
%         temp.(region_sublist{rix}).growseason(yix-2002,1)=nanmean(nanmean(growing_season(:,:,yix-2002).*(temp.(region_sublist{rix}).box)));
    end
end

figure;
plot(x,temp.Shelf.growseason)
hold on
plot(x,temp.Open.growseason)
legend('Shelf','Open')
title('Number of months that are considered ice-free ')

figure;
pcolor(lon_m,lat_m,growing_season(:,:,1)); shading flat
cm = get(gca,'Colormap');

plot_folder=['C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\figures\'];
load('openshelf_coord.mat')
for yix = 2003:2019
figure('units','normalized','outerposition',[0 0 1 1])
pcolor(lon_m,lat_m,growing_season(:,:,yix-2002)); shading flat
colormap(cm)
colorbar
geoshow('landareas.shp','facecolor','k')
title(num2str(yix))
xlim([-65,47]);
ylim([-80,-50])
hold on
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2)%'#80471C'
OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2)%,'LineStyle','--')
print('-dpng',[plot_folder ['IceFreeMonths',num2str(yix), '.png']])
close
end

for ix = 1:12
figure(ix);
pcolor(lon_m,lat_m,NPP_nans(:,:,ix)); shading flat
colormap(cm)
colorbar
geoshow('landareas.shp','facecolor','k')
title(num2str(time_start_all(ix,2)))
end

% Only including months where there is full coverage (so excludes April-August)
for yix = 2003:2019
    findyear=find(timedec>yix-0.3 & timedec<yix+0.25); % removes April-August (partial/no coverage)
    
    ice_covered(:,:,yix-2002)=nansum(icecover(:,:,findyear),3);
    growing_season(:,:,yix-2002)=7.-ice_covered(:,:,yix-2002);
end

growing_season(growing_season==0)=NaN;
growing_season_meanmap=mean(growing_season,3,'omitnan');
gs_temp=growing_season(:,:,1);

figure;
% subplot(2,1,1)
histogram(gs_temp(temp.(region_sublist{2}).box_logic),'Normalization','cdf')
hold on
% subplot(2,1,2)
histogram(gs_temp(temp.(region_sublist{3}).box_logic),'Normalization','cdf')

figure;
edges=[0.5:0.5:7];
yyaxis left
histogram(growing_season_meanmap(temp.(region_sublist{3}).box_logic),edges)%,'Normalization','pdf')
hold on
yyaxis right
histogram(growing_season_meanmap(temp.(region_sublist{2}).box_logic),edges)%,'Normalization','pdf')
legend('Open Ocean','Shelf')

edges=[0:1:7];
figure;
for yix = 2003:2019
    subplot(6,3,yix-2002)
    gstemp=growing_season(:,:,yix-2002);

    yyaxis left
    histogram(gstemp(temp.(region_sublist{3}).box_logic),edges)%,'Normalization','pdf')
    hold on
    
    yyaxis right
    histogram(gstemp(temp.(region_sublist{2}).box_logic),edges)%,'Normalization','pdf')
    hold on
    title(num2str(yix))
    if yix==2003
    legend('Open Ocean','Shelf')
    end
end



figure;
pcolor(lon_m,lat_m,growing_season(:,:,1)); shading flat
% cmblueice = get(gca,'Colormap');
for yix = 2003:2019
figure(yix-2002);
pcolor(lon_m,lat_m,growing_season(:,:,yix-2002)); shading flat
colormap(cm)
colorbar
geoshow('landareas.shp','facecolor','k')
title(num2str(yix))
end

for rix = 1:length(region_sublist)
    for yix = 2003:2019
        gstemp=growing_season(:,:,yix-2002);
        temp.(region_sublist{rix}).growseason(yix-2002,1)=mean(gstemp(temp.(region_sublist{rix}).box_logic),'omitnan');
        temp.(region_sublist{rix}).growseason_med(yix-2002,1)=median(gstemp(temp.(region_sublist{rix}).box_logic),'omitnan');
        temp.(region_sublist{rix}).growseason_mode(yix-2002,1)=mode(gstemp(temp.(region_sublist{rix}).box_logic));

        temp.(region_sublist{rix}).growseason_meanstd(yix-2002,1)=std(gstemp(temp.(region_sublist{rix}).box_logic),'omitnan');
    end
end

figure;
p1=plot(x,temp.Shelf.growseason_mode,'LineWidth',1.5,'Color',[0.8 0.4 0])%,'LineStyle','--')
hold on
p2=plot(x,temp.Open.growseason_mode,'LineWidth',1.5,'Color',[0.6 0.2 0.8])%,'LineStyle','--')
legend('Shelf','Open')
title('Mode number of months that are considered ice-free ')

x=2003:2019;
figure;
bar(x,[temp.Shelf.growseason_med temp.Open.growseason_med])
legend('Shelf','Open')
title('Median number of months that are considered ice-free ')
ylim([0 5])

figure;
p1=plot(x,temp.Shelf.growseason,'LineWidth',1.5,'Color',[0.8 0.4 0])%,'LineStyle','--')
hold on
% plot(x,temp.Shelf.growseason_med,'LineWidth',2,'Color',[0.7 0.3 0.1])
% hold on
p2=plot(x,temp.Open.growseason,'LineWidth',1.5,'Color',[0.6 0.2 0.8])%,'LineStyle','--')
% plot(x,temp.Open.growseason_med,'LineWidth',1.8,'Color',[0.5 0.1 0.7])

p3=plot(x,(temp.Shelf.growseason+temp.Shelf.growseason_meanstd),'LineWidth',1,'Color',[0.6 0.2 0],'LineStyle','--')
p4=plot(x,(temp.Shelf.growseason-temp.Shelf.growseason_meanstd),'LineWidth',1,'Color',[0.6 0.2 0],'LineStyle','--')

p5=plot(x,(temp.Open.growseason+temp.Open.growseason_meanstd),'LineWidth',1,'Color',[0.7 0.4 0.8],'LineStyle','--')
p6=plot(x,(temp.Open.growseason-temp.Open.growseason_meanstd),'LineWidth',1,'Color',[0.7 0.4 0.8],'LineStyle','--')

title('Number of months that are considered ice-free')
% legend('Shelf Mean','Shelf Median','Open Ocean Mean','Open Ocean Median','Shelf STD','Shelf STD','Open STD','Open STD')
legend([p1 p2 p3 p5],'Shelf Mean','Open Ocean Mean','Shelf STD','Open STD')
ylim([0 7])
ylabel('Months Ice-Free')
xlabel('Year')


% boxplot...
for rix = 1:length(region_sublist)
    for yix = 2003:2019
        gstemp=growing_season(:,:,yix-2002);
        temp.(region_sublist{rix}).growseason_map(:,:,yix-2002)=gstemp(temp.(region_sublist{rix}).box_logic);
    end
%     temp.(region_sublist{rix}).growseason_map(isnan(temp.(region_sublist{rix}).growseason_map))=0;
    temp.(region_sublist{rix}).growseason_map(temp.(region_sublist{rix}).growseason_map==0)=NaN;

end


figure;
boxplot([temp.Shelf.growseason_map(:,:,1) temp.Shelf.growseason_map(:,:,2) temp.Shelf.growseason_map(:,:,6)])

test=reshape(temp.Shelf.growseason_map,[9317,17]);
test2=reshape(temp.Open.growseason_map,[35282,17]);
figure;
subplot(1,2,1)
boxplot(test,x)
ylabel('Months ice-free')
title('Shelf Region')

subplot(1,2,2)
boxplot(test2,x)
ylabel('Months ice-free')
title('Open Ocean')

sgtitle('Number of months that are ice-free between September-March')


figure;
boxplot([temp.Open.growseason_map(:,:,1) temp.Open.growseason_map(:,:,8) temp.Open.growseason_map(:,:,14)])
