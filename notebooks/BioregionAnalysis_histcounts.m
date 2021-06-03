clearvars

home = 'C:\Users\ccd1n18\Documents\Projects\Carbon-Uptake-in-WG_Manuscript\'; % Clara laptop Doc folder
plot_folder = [home 'Figures\'];

LONG_MIN=-180; LONG_MAX=180; LAT_MIN=-90; LAT_MAX=90;
[Z,LONG,LAT]=m_tbase([LONG_MIN LONG_MAX LAT_MIN LAT_MAX]);

load('latlon_m.mat')
load('vgpm_imported.mat')
load('ProcessedData.mat')

% looking at documentation for interp2
% X = LONG; Y = LAT; V = Z;
% Xq = lon_m; Yq = lat_m
% Vq = interpolated elevation

Z_interp = interp2(LONG,LAT,Z,lon_m,lat_m);

figure
surf(LONG,LAT,Z) % or pcolor
shading flat

figure
surf(lon_m,lat_m,Z_interp) % or pcolor
shading flat

cd C:\Users\ccd1n18\Documents\Projects\Carbon-Uptake-in-WG_Manuscript\Workspace_Variables
load('openshelf_coord.mat')
load('box_lat_lons.mat', 'andrex_box')
%% Import regions: shelf and open ocean then andrex box
% load('latlon_m.mat')

IN_and=inpolygon(lon_m,lat_m,andrex_box(:,1),andrex_box(:,2));
findweddell=find(IN_and==1);

region_sublist={'Weddell'};
regionfindlist= {'findweddell'};

for rix = 1:length(region_sublist)
    %box,box logic
    temp.(region_sublist{rix}).box = lat_m;
    temp.(region_sublist{rix}).box(:) = 0;
    eval(['temp.',region_sublist{rix},'.box(',regionfindlist{rix},')=1;']);
    
    temp.(region_sublist{rix}).box_logic=logical(temp.(region_sublist{rix}).box);
    
    OceanProd.vgpm.(region_sublist{rix}).bathymetry=Z_interp.*(temp.(region_sublist{rix}).box);
end

%% plot distribution of depths in WG box
figure(1)
histogram(OceanProd.vgpm.Weddell.bathymetry(findweddell));
title('Count of Bathymetric Depths in the Weddell Gyre box');

xline(nanmean(OceanProd.vgpm.Weddell.bathymetry(findweddell)),'Color','r','LineWidth',2);
xline(median(OceanProd.vgpm.Weddell.bathymetry(findweddell),'omitnan'),'Color','b','LineWidth',2);
legend('','mean','median')
ylabel('Counts')
xlabel('Depth (m BSL)')
print('-dpng',[plot_folder  'Distr_WGboxDepths' '.png'])

% could even just do:
% figure
% histogram(Z_interp(findweddell));

%% NPP area-weighted rate 
% Average NPP during ice free conditions
temp.vgpm=vgpm_npp_all;
temp.findneg=find(temp.vgpm<0);
temp.vgpm(temp.findneg)=NaN;
vgpm_av_day_nan=nanmean(temp.vgpm,3);

%print('-dpng',[plot_folder  'Distr_NPPvBathy_WGbox' '.png'])

%% Area of regions
OceanProd.vgpm.Weddell.bathymetry(findweddell)
% area of WG
WeddellBoxArea=sum(sum(IN_and.*area_MODISVGPM_km2))
    % 6.079e6 km2
    
% area of shelf

% area of bathy bins beyond shelf

%% Total NPP 
% VGPM_npp_tot_gC_all
% VGPM_npp_tot_gC_nans

% temp.totalNPP_nans = VGPM_npp_tot_gC_nans(:,:,1);
temp.vgpmtot=VGPM_npp_tot_gC_nans;
% vgpm_av_tot_nan=nanmean(temp.vgpmtot,3); % this has calculated the average total NPP in each pixel

for yix = 2003:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    if size (findyear)<12
    DepthFunction.AnnIntegNPP(:,:,yix-2002)=NaN(size(temp.vgpmtot(:,:,1)));
    else
    DepthFunction.AnnIntegNPP(:,:,yix-2002)=nansum(temp.vgpmtot(:,:,findyear),3);    
    end
end

DepthFunction.AvAnnIntegNPP=nanmean(DepthFunction.AnnIntegNPP,3);

%% Total NPP in each bathy bin using histcounts
% need to find where in weddell box is <=0m so as to not include islands in
    % binning if specifying number of bins
% Instead, specify bin edges
    % -0.001 1000 2000 3000 4000 5000 6000
    % or if using bathy(findweddell), then bin sizes can be 0:1000:6000
x=OceanProd.vgpm.Weddell.bathymetry(findweddell);
setedges=[min(x) -5000 -4000 -3000 -2000 -1000 0 max(x)]
[N,edges,bin]=histcounts(OceanProd.vgpm.Weddell.bathymetry(findweddell),setedges)
        % N = number of occurrences within each bin
        % edges = in ascending order from -6000m to 0m
        % bin = bin number allocation (where 1 = -6000 - -5000m, and 6 is -1000m to 0m

% sum of NPP within each bin for average total NPP
temp.testfind = DepthFunction.AvAnnIntegNPP(findweddell);
for ix=1:1:6
    temp.histcountsSUM(:,ix)=(nansum(temp.testfind(bin == ix)))/1e12;
end
temp.histcountsSUM
temp.plotX = categorical({'0-1000m','1000-2000m','2000-3000m','3000-4000m','4000-5000m','5000m+'}) %'Above Sea Level',
temp.plotX = reordercats(temp.plotX,{'0-1000m','1000-2000m','2000-3000m','3000-4000m','4000-5000m','5000m+'}) %'Above Sea Level',
temp.plotY=temp.histcountsSUM'
temp.plotY=flipud(temp.plotY)

figure
bar(temp.plotX,temp.plotY)

% for all the years
for yix = 2003:2019
    temp_totNPParray=DepthFunction.AnnIntegNPP(:,:,yix-2002);
    temp_testfind=temp_totNPParray(findweddell);
    for ix=1:1:6
        temp.histcounts_IntNPP(ix,yix-2002)=(nansum(temp_testfind(bin == ix)))/1e12;
    end
    clear temp_totNPParray temp_testfind
end

temp.plotY=flipud(temp.histcounts_IntNPP);
temp.NamesS=repmat({'0-1000m'},17,1)
temp.NamesS1=repmat({'1000-2000m'},17,1)
temp.Names2=repmat({'2000-3000m'},17,1)
temp.Names3=repmat({'3000-4000m'},17,1)
temp.Names4=repmat({'4000-5000m'},17,1)
temp.Names5=repmat({'5000m+'},17,1)
temp.Names=[temp.NamesS ; temp.NamesS1 ; temp.Names2 ; temp.Names3 ; temp.Names4 ; temp.Names5];
figure
boxplot(temp.plotY',temp.Names)
xlabel('Bathymetric bin')
ylabel('Total Annual NPP (Tg)')
xtickangle(45)


% anova?...
% [p,tbl]=anova1(temp.plotY(1:2,:))



%% Daily NPP 
% Average daily NPP during ice free conditions
temp.vgpm=vgpm_npp_all;
temp.findneg=find(temp.vgpm<0);
temp.vgpm(temp.findneg)=NaN;
% vgpm_av_day_nan=nanmean(temp.vgpm,3);

for yix = 2003:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    if size (findyear)<12
    DepthFunction.DailyNPP(:,:,yix-2002)=NaN(size(temp.vgpm(:,:,1)));
    else
    DepthFunction.DailyNPP(:,:,yix-2002)=nanmean(temp.vgpm(:,:,findyear),3);    
    end
end

DepthFunction.AvDailyNPP=nanmean(DepthFunction.DailyNPP,3);


%% Average NPP (mg/m2/day) in each bathy bin using histcounts
temp.testfind_rate = DepthFunction.AvDailyNPP(findweddell);
for ix=1:1:6
    temp.histcountsRATES(:,ix)=(nanmean(temp.testfind_rate(bin == ix)));
    temp.histcountsRATES_max(:,ix)=max(temp.testfind_rate(bin == ix));
    temp.histcountsRATES_min(:,ix)=min(temp.testfind_rate(bin == ix));
end
temp.histcountsRATES
temp.histcountsRATES_max

%plot
temp.plotX = categorical({'0-1000m','1000-2000m','2000-3000m','3000-4000m','4000-5000m','5000m+'}) %'Above Sea Level',
temp.plotX = reordercats(temp.plotX,{'0-1000m','1000-2000m','2000-3000m','3000-4000m','4000-5000m','5000m+'}) %'Above Sea Level',
temp.plotY=temp.histcountsRATES'
temp.plotY=flipud(temp.plotY)

figure
bar(temp.plotX,temp.plotY) % average rates per bathy bin (doesn't include the spread of rates in each bin)

temp.plotY=temp.histcountsRATES_max'
temp.plotY=flipud(temp.plotY)
figure
bar(temp.plotX,temp.plotY) % max rates per bathy bin 

temp.plotY=temp.histcountsRATES_min'
temp.plotY=flipud(temp.plotY)
figure
bar(temp.plotX,temp.plotY) % min rates per bathy bin (doesn't go down to zero.. because it is an average?)


% for all the years
for yix = 2003:2019
    temp_dailyNPParray=DepthFunction.DailyNPP(:,:,yix-2002);
    temp_testfind=temp_dailyNPParray(findweddell);
    for ix=1:1:6
        temp.histcounts_dailyNPP(ix,yix-2002)=(nanmean(temp_testfind(bin == ix)));
        temp.histcounts_dailyNPP(ix,yix-2002)=(max(temp_testfind(bin == ix)));
    end
    %clear temp_dailyNPParray temp_testfind
end

temp.plotY=flipud(temp.histcounts_dailyNPP);
temp.NamesS=repmat({'0-1000m'},17,1)
temp.NamesS1=repmat({'1000-2000m'},17,1)
temp.Names2=repmat({'2000-3000m'},17,1)
temp.Names3=repmat({'3000-4000m'},17,1)
temp.Names4=repmat({'4000-5000m'},17,1)
temp.Names5=repmat({'5000m+'},17,1)
temp.Names=[temp.NamesS ; temp.NamesS1 ; temp.Names2 ; temp.Names3 ; temp.Names4 ; temp.Names5];
figure
boxplot(temp.plotY',temp.Names)
xlabel('Bathymetric bin')
ylabel('Average NPP (mg m^-^2 day^-^1) for each austral year')
xtickangle(45)



%% Bathymetry distribution within high NPP regions
x_NPP=vgpm_av_day_nan(findweddell);
x_NPP(isnan(x_NPP))=-9999;
setedgesNPP=[-9999 0:100:500 800:300:1800]
[N_NPP,edges_NPP,bin_NPP]=histcounts(x_NPP,setedgesNPP)


% boxplot using bin# vs bathy
x_bathy=OceanProd.vgpm.Weddell.bathymetry(findweddell);

figure
boxplot(OceanProd.vgpm.Weddell.bathymetry(findweddell),bin_NPP)
xlim([1.5 10.5])
xticklabels(edges_NPP(2:end))
xlabel('NPP (mg m^-^2 day^-^1) upper bin edge')
ylabel('Bathymetric depth (m)')
%% histcounts 2
% using bins 

x=OceanProd.vgpm.Weddell.bathymetry(findweddell); % bathymetry within the Weddell Box
y=temp.testfind_rate; % average NPP when ice-free (<-9999 was treated as NaNs, not 0)
    temp.nan=find(isnan(y));
    y(temp.nan)=-9999; % permanent ice covered reverted back to -9999
Xedges=[-6030 -5000 -4000 -3000 -2000 -1000 0 max(x)];
Yedges=[-9999 0 200 max(y)]
[N2,Xedgesout,Yedgesout,binX,binY] = histcounts2(x,y,Xedges,Yedges)
        % N2 comes out with sig less counts that there should be
        sum(sum(N2))
        
% this creates bins for:
    % 5000m BSL and deeper
    % 4000-5000m BSL
    % 3000-4000m BSL
    % 2000-3000m BSL
    % 1000-2000m BSL
    % 0-1000m BSL
    % above sea level (44 pixel counts)
    
    % permanently ice covered
    % 0-200 mg/m2/day NPP (rates)
    % greater than 200 mg/m2/day NPP 
    

for ix=1:6 % 1:6 to remove above sea level or 1:length(N2(:,1))
    for nix = 2:3 % 2:3 and nix-1 for column input to remove perm ice covered or 1:length(N2(1,:))
        temp.histcounts2SUM(ix,nix-1)=(nansum(temp.testfind(binX == ix & binY ==nix)))/1e12;
    end
end

temp.plotY=temp.histcounts2SUM;
temp.plotY=flipud(temp.plotY);
% temp.plotY=temp.plotY';
temp.plotX = categorical({'0-1000m','1000-2000m','2000-3000m','3000-4000m','4000-5000m','5000m+'}) %'Above Sea Level',
temp.plotX = reordercats(temp.plotX,{'0-1000m','1000-2000m','2000-3000m','3000-4000m','4000-5000m','5000m+'}) %'Above Sea Level',

figure
b=bar(temp.plotX,temp.plotY);
b(1).FaceColor=[0.5 0.75 1];
b(2).FaceColor=[0 0.5 0];
legend('Low productivity (<200 mg m^-^2 day^-^1)','High productivity (>200 mg m^-^2 day^-^1)','Location','northwest') % 'Permanently Ice Covered',
xlabel('Bathymetry bin')
ylabel('Annual NPP (Tg)')
        
temp.histcounts2SUM_total=sum(sum(temp.histcounts2SUM))
for ix=1:6 % 1:6 to remove above sea level or 1:length(N2(:,1))
    for nix = 1:2 % 2:3 and nix-1 for column input to remove perm ice covered or 1:length(N2(1,:))
        temp.histcounts2SUM_contribution(ix,nix)=(temp.histcounts2SUM(ix,nix)/temp.histcounts2SUM_total)*100;        
    end
end

figure
histogram2(x,y,Xedges,Yedges(2:end))