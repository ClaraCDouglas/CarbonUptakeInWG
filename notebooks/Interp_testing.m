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

%% NPP rate vs depth
% Average NPP during ice free conditions
temp.vgpm=vgpm_npp_all;
temp.findneg=find(temp.vgpm<0);
temp.vgpm(temp.findneg)=NaN;
vgpm_av_day_nan=nanmean(temp.vgpm,3);

% find NPP in bathy groups
DepthFunction.findShelf=find(OceanProd.vgpm.Weddell.bathymetry<=-0.001 & OceanProd.vgpm.Weddell.bathymetry>=-1999.99); % CAREFUL WITH THIS BECAUSE IT WILL SELECT SOME OF THE SHALLOW AREAS ALONG THE NORTHERN EDGE AND MAUD RISE
DepthFunction.NPPShelf=vgpm_av_day_nan(DepthFunction.findShelf);
DepthFunction.NPPShelf_mean=nanmean(DepthFunction.NPPShelf)
DepthFunction.NPPShelf_median=median(DepthFunction.NPPShelf,'omitnan')
figure
histogram(DepthFunction.NPPShelf);

DepthFunction.findOFFS_2to3=find(OceanProd.vgpm.Weddell.bathymetry<=-2000 & OceanProd.vgpm.Weddell.bathymetry>=-2999.99);
DepthFunction.NPPOFFS_2to3=vgpm_av_day_nan(DepthFunction.findOFFS_2to3);
DepthFunction.NPPOFFS_2to3_mean=nanmean(DepthFunction.NPPOFFS_2to3)
DepthFunction.NPPOFFS_2to3_median=median(DepthFunction.NPPOFFS_2to3,'omitnan')
figure
histogram(DepthFunction.NPPOFFS_2to3);

DepthFunction.findOFFS_3to4=find(OceanProd.vgpm.Weddell.bathymetry<=-3000 & OceanProd.vgpm.Weddell.bathymetry>=-3999.99);
DepthFunction.NPPOFFS_3to4=vgpm_av_day_nan(DepthFunction.findOFFS_3to4);
DepthFunction.NPPOFFS_3to4_mean=nanmean(DepthFunction.NPPOFFS_3to4)
DepthFunction.NPPOFFS_3to4_median=median(DepthFunction.NPPOFFS_3to4,'omitnan')
figure
histogram(DepthFunction.NPPOFFS_3to4);

DepthFunction.findOFFS_4to5=find(OceanProd.vgpm.Weddell.bathymetry<=-4000 & OceanProd.vgpm.Weddell.bathymetry>=-4999.99);
DepthFunction.NPPOFFS_4to5=vgpm_av_day_nan(DepthFunction.findOFFS_4to5);
DepthFunction.NPPOFFS_4to5_mean=nanmean(DepthFunction.NPPOFFS_4to5)
DepthFunction.NPPOFFS_4to5_median=median(DepthFunction.NPPOFFS_4to5,'omitnan')
figure
histogram(DepthFunction.NPPOFFS_4to5);

DepthFunction.findOFFS_5to6=find(OceanProd.vgpm.Weddell.bathymetry<=-5000 & OceanProd.vgpm.Weddell.bathymetry>=-6000);
DepthFunction.NPPOFFS_5to6=vgpm_av_day_nan(DepthFunction.findOFFS_5to6);
DepthFunction.NPPOFFS_5to6_mean=nanmean(DepthFunction.NPPOFFS_5to6)
DepthFunction.NPPOFFS_5to6_median=median(DepthFunction.NPPOFFS_5to6,'omitnan')
figure
histogram(DepthFunction.NPPOFFS_5to6);


%% plotting histograms of NPP rates in bathy bins
graphical.shelfcolor=[0.8 0.4 0];
graphical.opencolor=[0.4 0.2 0.8];


figure
subplot(2,3,1)
histogram(DepthFunction.NPPShelf,0:20:1600,'FaceAlpha',0.8,'FaceColor',graphical.shelfcolor);
hold on
histogram(DepthFunction.NPPOFFS_2to3,'FaceAlpha',0.8,'FaceColor',[0.50 0.08 0.5]);
histogram(DepthFunction.NPPOFFS_3to4,'FaceAlpha',0.8,'FaceColor',[0.50 0.08 0]);
histogram(DepthFunction.NPPOFFS_4to5,'FaceAlpha',0.8,'FaceColor',[0.05 0.3 0]);
histogram(DepthFunction.NPPOFFS_5to6,'FaceAlpha',0.8,'FaceColor',[0.2 0 0.4]);
legend('Shelf','2000','3000','4000','5000')

subplot(2,3,2)
histogram(DepthFunction.NPPShelf,0:20:1600,'FaceAlpha',0.8,'FaceColor',graphical.shelfcolor);
xline(DepthFunction.NPPShelf_mean,'Color','r','LineWidth',2);
xline(DepthFunction.NPPShelf_median,'Color','b','LineWidth',2);
legend('Shelf','mean','median')
subplot(2,3,3)
histogram(DepthFunction.NPPOFFS_2to3,'FaceAlpha',0.8,'FaceColor',[0.50 0.08 0.5]);
xline(DepthFunction.NPPOFFS_2to3_mean,'Color','r','LineWidth',2);
xline(DepthFunction.NPPOFFS_2to3_median,'Color','b','LineWidth',2);
legend('2000-2999m','mean','median')

subplot(2,3,4)
histogram(DepthFunction.NPPOFFS_3to4,'FaceAlpha',0.8,'FaceColor',[0.50 0.08 0]);
xline(DepthFunction.NPPOFFS_3to4_mean,'Color','r','LineWidth',2);
xline(DepthFunction.NPPOFFS_3to4_median,'Color','b','LineWidth',2);
legend('3000-3999m','mean','median')
subplot(2,3,5)
histogram(DepthFunction.NPPOFFS_4to5,'FaceAlpha',0.8,'FaceColor',[0.05 0.3 0]);
xline(DepthFunction.NPPOFFS_4to5_mean,'Color','r','LineWidth',2);
xline(DepthFunction.NPPOFFS_4to5_median,'Color','b','LineWidth',2);
legend('4000-4999m','mean','median')
subplot(2,3,6)
histogram(DepthFunction.NPPOFFS_5to6,'FaceAlpha',0.8,'FaceColor',[0.2 0 0.4]);
xline(DepthFunction.NPPOFFS_5to6_mean,'Color','r','LineWidth',2);
xline(DepthFunction.NPPOFFS_5to6_median,'Color','b','LineWidth',2);
legend('5000-6000m','mean','median')
title('')
suptitle('Distribution of Rates of NPP in bathymetry groups')

print('-dpng',[plot_folder  'Distr_NPPvBathy_WGbox' '.png'])

%% ANOVA
tab4stats=[DepthFunction.NPPShelf_mean DepthFunction.NPPOFFS_3to4_mean DepthFunction.NPPOFFS_4to5_mean DepthFunction.NPPOFFS_5to6_mean]
[p,tbl,stats] = anova1(tab4stats)

%% Area of regions
OceanProd.vgpm.Weddell.bathymetry(findweddell)
% area of WG
WeddellBoxArea=sum(sum(IN_and.*area_MODISVGPM_km2))
    % 6.079e6 km2
    
% area of shelf

% area of bathy bins beyond shelf

%% Total NPP in each bathy bin 
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


DepthFunction.NPPtotShelf=DepthFunction.AvAnnIntegNPP(DepthFunction.findShelf); % this has listed the average totals for each pixel within the bathy bin
                                                                    % need to sum these cells to get the total NPP in the bathy bin
                                                                    % then to compare to each of the other bins and determine percentage contribution of total
DepthFunction.NPPtotShelf_total=nansum(DepthFunction.NPPtotShelf)
DepthFunction.NPPtotShelf_totalTg=DepthFunction.NPPtotShelf_total/1e12

DepthFunction.NPPtotOFFS_2to3=DepthFunction.AvAnnIntegNPP(DepthFunction.findOFFS_2to3);
DepthFunction.NPPtotOFFS_2to3_total=nansum(DepthFunction.NPPtotOFFS_2to3)
DepthFunction.NPPtotOFFS_2to3_totalTg=DepthFunction.NPPtotOFFS_2to3_total/1e12

DepthFunction.NPPtotOFFS_3to4=DepthFunction.AvAnnIntegNPP(DepthFunction.findOFFS_3to4);
DepthFunction.NPPtotOFFS_3to4_total=nansum(DepthFunction.NPPtotOFFS_3to4)
DepthFunction.NPPtotOFFS_3to4_totalTg=DepthFunction.NPPtotOFFS_3to4_total/1e12

DepthFunction.NPPtotOFFS_4to5=DepthFunction.AvAnnIntegNPP(DepthFunction.findOFFS_4to5);
DepthFunction.NPPtotOFFS_4to5_total=nansum(DepthFunction.NPPtotOFFS_4to5)
DepthFunction.NPPtotOFFS_4to5_totalTg=DepthFunction.NPPtotOFFS_4to5_total/1e12

DepthFunction.NPPtotOFFS_5to6=DepthFunction.AvAnnIntegNPP(DepthFunction.findOFFS_5to6);
DepthFunction.NPPtotOFFS_5to6_total=nansum(DepthFunction.NPPtotOFFS_5to6)
DepthFunction.NPPtotOFFS_5to6_totalTg=DepthFunction.NPPtotOFFS_5to6_total/1e12

temp.checktotal=DepthFunction.NPPtotOFFS_5to6_totalTg+DepthFunction.NPPtotOFFS_4to5_totalTg+DepthFunction.NPPtotOFFS_3to4_totalTg+ ...
    DepthFunction.NPPtotOFFS_2to3_totalTg+DepthFunction.NPPtotShelf_totalTg


X = categorical({'Shelf','2000','3000','4000','5000'});
X = reordercats(X,{'Shelf','2000','3000','4000','5000'});
Y = [DepthFunction.NPPtotShelf_totalTg DepthFunction.NPPtotOFFS_2to3_totalTg DepthFunction.NPPtotOFFS_3to4_totalTg ...
     DepthFunction.NPPtotOFFS_4to5_totalTg DepthFunction.NPPtotOFFS_5to6_totalTg];
figure
bar(X,Y)

% for all the years

for ix=1:1:17
    temp.totalNPPall=DepthFunction.AnnIntegNPP(:,:,ix);
    tempNPPtotShelf_all=temp.totalNPPall(DepthFunction.findShelf);
    DepthFunction.NPPtotShelf_all(:,ix)=sum(tempNPPtotShelf_all);
    clear tempNPPtotShelf_all
end

for ix=1:1:17
    temp.totalNPPall=DepthFunction.AnnIntegNPP(:,:,ix);
    tempNPPtotOFFS_2to3_all=temp.totalNPPall(DepthFunction.findOFFS_2to3);
    DepthFunction.NPPtotOFFS_2to3_all(:,ix)=sum(tempNPPtotOFFS_2to3_all);
    clear tempNPPtotOFFS_2to3_all
end
for ix=1:1:17
    temp.totalNPPall=DepthFunction.AnnIntegNPP(:,:,ix);
    tempNPPtotOFFS_3to4_all=temp.totalNPPall(DepthFunction.findOFFS_3to4);
    DepthFunction.NPPtotOFFS_3to4_all(:,ix)=sum(tempNPPtotOFFS_3to4_all);
    clear tempNPPtotOFFS_3to4_all
end
for ix=1:1:17
    temp.totalNPPall=DepthFunction.AnnIntegNPP(:,:,ix);
    tempNPPtotOFFS_4to5_all=temp.totalNPPall(DepthFunction.findOFFS_4to5);
    DepthFunction.NPPtotOFFS_4to5_all(:,ix)=sum(tempNPPtotOFFS_4to5_all);
    clear tempNPPtotOFFS_4to5_all
end
for ix=1:1:17
    temp.totalNPPall=DepthFunction.AnnIntegNPP(:,:,ix);
    tempNPPtotOFFS_5to6_all=temp.totalNPPall(DepthFunction.findOFFS_5to6);
    DepthFunction.NPPtotOFFS_5to6_all(:,ix)=sum(tempNPPtotOFFS_5to6_all);
    clear tempNPPtotOFFS_5to6_all
end

temp.Numbers=[DepthFunction.NPPtotShelf_all' ; DepthFunction.NPPtotOFFS_2to3_all' ; DepthFunction.NPPtotOFFS_3to4_all' ; ...
    DepthFunction.NPPtotOFFS_4to5_all' ; DepthFunction.NPPtotOFFS_5to6_all']
temp.NumbersTg=temp.Numbers/1e12;
temp.NamesS=repmat({'Shelf'},17,1)
temp.Names2=repmat({'2000'},17,1)
temp.Names3=repmat({'3000'},17,1)
temp.Names4=repmat({'4000'},17,1)
temp.Names5=repmat({'5000'},17,1)
temp.Names=[temp.NamesS ; temp.Names2 ; temp.Names3 ; temp.Names4 ; temp.Names5]

figure
boxplot(temp.NumbersTg,temp.Names)
xlabel('Bathy bin')
ylabel('Total Annual NPP (Tg)')

temp.NumbersCols=[DepthFunction.NPPtotShelf_all' DepthFunction.NPPtotOFFS_2to3_all' DepthFunction.NPPtotOFFS_3to4_all' ...
    DepthFunction.NPPtotOFFS_4to5_all' DepthFunction.NPPtotOFFS_5to6_all']
temp.NumbersColsTg=temp.NumbersCols/1e12;
temp.NamesCols=[temp.NamesS temp.Names2 temp.Names3 temp.Names4 temp.Names5]

[p,tbl]=anova1(temp.NumbersColsTg(:,4:5))


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

% sum of NPP within each bin
temp.testfind = DepthFunction.AvAnnIntegNPP(findweddell);
for ix=1:1:6
    temp.histcountsSUM(:,ix)=(nansum(temp.testfind(bin == ix)))/1e12;
end
temp.histcountsSUM

%% Average NPP (mg/m2/day) in each bathy bin using histcounts
% Average NPP during ice free conditions
temp.vgpm=vgpm_npp_all;
temp.findneg=find(temp.vgpm<0);
temp.vgpm(temp.findneg)=NaN;
vgpm_av_day_nan=nanmean(temp.vgpm,3);

temp.testfind_rate = vgpm_av_day_nan(findweddell);
for ix=1:1:6
    temp.histcountsRATES(:,ix)=(nanmean(temp.testfind_rate(bin == ix)));
    temp.histcountsRATES_max(:,ix)=max(temp.testfind_rate(bin == ix));
end
temp.histcountsRATES
temp.histcountsRATES_max

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

temp.pieX = {'0-1000m L','1000-2000m L','2000-3000m L','3000-4000m L','4000-5000m L','5000m+ L';...
    '0-1000m H','1000-2000m H','2000-3000m H','3000-4000m H','4000-5000m H','5000m+ H'} %'Above Sea Level',
temp.
figure
p=pie(temp.plotY,temp.pieX);

