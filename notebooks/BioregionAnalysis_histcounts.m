% clearvars

home = 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg'; % desktop computer

plot_folder = 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\figures'; % desktop computer
make_plots =0;

LONG_MIN=-180; LONG_MAX=180; LAT_MIN=-90; LAT_MAX=90;
[Z,LONG,LAT]=m_tbase([LONG_MIN LONG_MAX LAT_MIN LAT_MAX]);

cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop computer
load('latlon_m.mat')
load('box_lat_lons.mat', 'andrex_box')
load('vgpm_imported.mat')
load('ProcessedData.mat')

% looking at documentation for interp2
% X = LONG; Y = LAT; V = Z;
% Xq = lon_m; Yq = lat_m
% Vq = interpolated elevation

Z_interp = interp2(LONG,LAT,Z,lon_m,lat_m);

if make_plots
figure
surf(LONG,LAT,Z) % or pcolor
shading flat

figure
surf(lon_m,lat_m,Z_interp) % or pcolor
shading flat
end

load('openshelf_coord.mat')
load('box_lat_lons.mat', 'andrex_box')
%% Allocate bathymetry to ANDREX box
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
if make_plots
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
end
%% NPP area-weighted rate 
% Average NPP during ice free conditions
temp.vgpm=vgpm_npp_all;
temp.findneg=find(temp.vgpm<0);
temp.vgpm(temp.findneg)=NaN;
	% average growing season daily rate across the whole timeseries. 
        vgpm_av_day_nan=nanmean(temp.vgpm,3); 
        
% average growing season daily rate per year
NPP.vgpm_annual_day_nan(:,:,length([2003:2019]))=NaN(size(temp.vgpm(:,:,1)));    
for yix = 2003:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    if size (findyear)<12
    NPP_years.vgpm_annual_day_nan(:,:,yix-2002)=NaN(size(temp.vgpm(:,:,1)));
    else
    NPP_years.vgpm_annual_day_nan(:,:,yix-2002)=nanmean(temp.vgpm(:,:,findyear),3);    
    end
end
% climatology of average growing season daily rate per year
NPP.vgpm_annual_av_day_nan=nanmean(NPP_years.vgpm_annual_day_nan,3);
    
%% Area of regions
OceanProd.vgpm.Weddell.bathymetry(findweddell)
% area of WG
WeddellBoxArea=sum(sum(IN_and.*area_MODISVGPM_km2))
    % 6.079e6 km2
    
% area of shelf

% area of bathy bins beyond shelf

%% Total NPP 
% VGPM_npp_tot_gC_all - land/permanent ice are set as zeros
% VGPM_npp_tot_gC_nans - land/permanent ice are set as NaNs
load('vgpm_imported.mat', 'VGPM_npp_tot_gC_all', 'VGPM_npp_tot_gC_nans')
    % monthly climatology of total NPP per pixel for whole time series 
        % this has calculated the average total NPP in each pixel
        NPP.vgpm_tot_monthclim=nanmean(VGPM_npp_tot_gC_all,3);
        NPP.vgpm_tot_monthclim(find(NPP.vgpm_tot_monthclim==0))=NaN;

% annual climatology of total NPP per year
temp.vgpmtot=VGPM_npp_tot_gC_nans;
% vgpm_av_tot_nan=nanmean(temp.vgpmtot,3); % this has calculated the average total NPP in each pixel
NPP_years.vgpm_tot_years(:,:,length([2003:2019]))=NaN(size(temp.vgpmtot(:,:,1)));    
for yix = 2003:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    if size (findyear)<12
    NPP_years.vgpm_tot_years(:,:,yix-2002)=NaN(size(temp.vgpmtot(:,:,1)));
    else
    NPP_years.vgpm_tot_years(:,:,yix-2002)=nansum(temp.vgpmtot(:,:,findyear),3);    
    end
end
NPP.vgpm_tot_av_years=nanmean(NPP_years.vgpm_tot_years,3);
NPP.vgpm_tot_av_years(find(NPP.vgpm_tot_av_years==0))=NaN;

%% Total NPP in each bathy bin using histcounts
% need to find where in weddell box is <=0m so as to not include islands in
    % binning if specifying number of bins
% Instead, specify bin edges
    % -0.001 1000 2000 3000 4000 5000 6000
    % or if using bathy(findweddell), then bin sizes can be 0:1000:6000
z=OceanProd.vgpm.Weddell.bathymetry(findweddell);
setedges=[min(z) -5000 -4000 -3000 -2000 -1000 0 max(z)];
[count_in_bin,edges,bin_bathy]=histcounts(OceanProd.vgpm.Weddell.bathymetry(findweddell),setedges);
        % N = number of occurrences within each bin
        % edges = in ascending order from -6000m to 0m
        % bin = bin number allocation (where 1 = -6000 - -5000m, and 6 is -1000m to 0m

% sum of NPP within each bin for average total NPP (TgC)
temp.testfind = NPP.vgpm_tot_av_years(findweddell);
for ix=1:1:6 %length(setedges) if using all the bins (the last 2 are above sea level)
    temp.histcountsSUM(:,ix)=(nansum(temp.testfind(bin_bathy == ix)))/1e12;
end
temp.histcountsSUM

if make_plots
    temp.plotX = categorical({'0-1000m','1000-2000m','2000-3000m','3000-4000m','4000-5000m','5000m+'}) %'Above Sea Level',
    temp.plotX = reordercats(temp.plotX,{'0-1000m','1000-2000m','2000-3000m','3000-4000m','4000-5000m','5000m+'}) %'Above Sea Level',
    temp.plotY=temp.histcountsSUM'
    temp.plotY=flipud(temp.plotY)
    
    figure
    bar(temp.plotX,temp.plotY)
end

% for all the years
for yix = 2003:2019
    temp_totNPParray=NPP_years.vgpm_tot_years(:,:,yix-2002);
    temp_testfind=temp_totNPParray(findweddell);
    for ix=1:1:6
        temp.histcounts_IntNPP(ix,yix-2002)=(nansum(temp_testfind(bin_bathy == ix)))/1e12;
    end
    clear temp_totNPParray temp_testfind
end

if make_plots
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
end

%% Average NPP (mg/m2/day) in each bathy bin using histcounts
temp.testfind_rate = NPP.vgpm_annual_av_day_nan(findweddell);
for ix=1:1:6
    temp.histcountsRATES(:,ix)=(nanmean(temp.testfind_rate(bin_bathy == ix)));
    temp.histcountsRATES_max(:,ix)=max(temp.testfind_rate(bin_bathy == ix));
    temp.histcountsRATES_min(:,ix)=min(temp.testfind_rate(bin_bathy == ix));
end
temp.histcountsRATES
temp.histcountsRATES_max

%plot
if make_plot
    temp.plotYratemean=temp.histcountsRATES'
    temp.plotYratemean=flipud(temp.plotYratemean)
    
    figure
    bar(temp.plotX,temp.plotYratemean) % average rates per bathy bin (doesn't include the spread of rates in each bin)
    
    temp.plotYratemax=temp.histcountsRATES_max'
    temp.plotYratemax=flipud(temp.plotYratemax)
    figure
    bar(temp.plotX,temp.plotYratemax) % max rates per bathy bin
    
%     temp.plotYratetog(:,1)=temp.plotYratemean;
%     temp.plotYratetog(:,2)=temp.plotYratemax;
%     figure
%         bar(temp.plotX,temp.plotYratetog) % max rates per bathy bin

%     temp.plotY=temp.histcountsRATES_min'
%     temp.plotY=flipud(temp.plotY)
%     figure
%     bar(temp.plotX,temp.plotY) % min rates per bathy bin (doesn't go down to zero.. because it is an average?)
end

% for all the years
for yix = 2003:2019
    temp_dailyNPParray=NPP_years.vgpm_annual_day_nan(:,:,yix-2002);
    temp_testfind=temp_dailyNPParray(findweddell);
    for ix=1:1:6
        temp.histcounts_dailyNPP(ix,yix-2002)=(nanmean(temp_testfind(bin_bathy == ix)));
        temp.histcounts_dailyNPP(ix,yix-2002)=(max(temp_testfind(bin_bathy == ix)));
    end
    clear temp_dailyNPParray temp_testfind
end

if make_plots
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
    title('make a useful title')
end

%% Bathymetry distribution within high NPP regions - ignore just now
x_NPP=vgpm_av_day_nan(findweddell);
x_NPP(isnan(x_NPP))=-9999;
setedgesNPP=[-9999 0:100:500 800:300:1800]
[N_NPP,edges_NPP,bin_NPP]=histcounts(x_NPP,setedgesNPP)
        % N_NPP = number of occurrences within each bin
        % edges_NPP = in ascending order from -9999 (NaN) to 1800
        % bin = bin number allocation (where 1 = -9999 - 0, 2 is 0-100 and so on


% boxplot using bin# vs bathy
x_bathy=OceanProd.vgpm.Weddell.bathymetry(findweddell);

figure
boxplot(OceanProd.vgpm.Weddell.bathymetry(findweddell),bin_NPP)
xlim([1.5 10.5])
xticklabels(edges_NPP(2:end))
xlabel('NPP (mg m^-^2 day^-^1) upper bin edge')
ylabel('Bathymetric depth (m)')
%% histcounts 2 - binning by bathy and daily NPP rates
% using bins 

z=OceanProd.vgpm.Weddell.bathymetry(findweddell); % bathymetry within the Weddell Box
y=NPP.vgpm_annual_av_day_nan(findweddell); % average NPP during growing season = when ice-free (<-9999 was treated as NaNs, not 0)
temp.nan=find(isnan(y));
y(temp.nan)=-9999; % permanent ice covered reverted back to -9999
bathyedges=[min(z) -5000 -4000 -3000 -2000 -1000 0 max(z)];
NPPsetedges=[-9999 0 200 max(y)];
[count_in_hist2bins,bathyedgesout,NPPedgesout,binX_bathy,binY_NPP] = histcounts2(z,y,bathyedges,NPPsetedges);
        % N2 comes out with sig less counts that there should be
        sum(sum(count_in_hist2bins))
        
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
        temp.histcounts2SUM(ix,nix-1)=(nansum(temp.testfind(binX_bathy == ix & binY_NPP ==nix)))/1e12;
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


%% Making a-w growing season productivity bins for all years
%we want:
bin_NPP_allyears=zeros(length(y),1,length([2003:2019]));
count_in_NPPbin_allyears=zeros(1,3,length([2003:2019]));

for yix = 1:length([2003:2019])
    NPPin=NPP_years.vgpm_annual_day_nan(:,:,yix);
    NPPin=NPPin(findweddell);
    NPPin(isnan(NPPin))=-9999;
    NPPsetedges=[-9999 0 200 max(NPPin)];
    [count_in_NPPbin,edgesNPP,bin_NPP]=histcounts(NPPin,NPPsetedges);
    bin_NPP_allyears(:,1,yix)=bin_NPP;
    count_in_NPPbin_allyears(:,:,yix)=count_in_NPPbin;
    if yix<2019
        clear NPPin NPPsetedges count_in_NPPbin bin_NPP
    elseif yix==2019
        clear NPPin count_in_NPPbin bin_NPP
    end
end

% map NPP rates area-weighted for each year
temp.mapNPPbins=zeros(1080,2160);
mapNPPbins=zeros(1080,2160,length([2003:2019]));
for yix = 1:length([2003:2019])
for ix=1:length(edgesNPP)
    prodbins=findweddell(bin_NPP_allyears(:,:,yix)==ix);
    temp.mapNPPbins(prodbins)=ix;
end
    mapNPPbins(:,:,yix)=temp.mapNPPbins;
end
mapNPPbins(find(mapNPPbins==0))=NaN;

% map NPP rate bins for climatology
mapprodbins=zeros(1080,2160);
for ix=1:length(NPPedgesout)
    curious=findweddell(binY_NPP==ix);
    mapprodbins(curious)=ix;
    clear curious
end
temp.nanan=find(mapprodbins==0);
mapprodbins(temp.nanan)=NaN;

% make the maps
map2 = [0 0 0
    0.5 0.75 1
    0 0.5 0];

years=2003:2019;

P=figure('units','normalized','outerposition',[0 0 1 1]);
for pix = 1:length([2003:2019])
    subplot(3,6,pix)
    pcolor(lon_m,lat_m,mapNPPbins(:,:,pix)); shading flat; hold on
    colormap(map2)
    xlim([-65,47]);
    ylim([-80,-50]);
    title(['',num2str(years(pix)), ''])
end 
% add mean and color bar to last plot
    subplot(3,6,18)
    pcolor(lon_m,lat_m,mapprodbins); shading flat; hold on
    colormap(map2)
    xlim([-65,47]);
    ylim([-80,-50]);
    caxis([1 3]);
    title(['Climatology'])

pp=colorbar
pp.Ticks=[1.45:0.65:3];
pp.TickLabels={'None', '0-200','>200'};
pp.Label.String = 'Area-weighted productivity bins (mg m^-^2 d^-^1)';
pp.Label.Rotation = 90;

sgtitle('Productivity bins based on average area-weighted growing season NPP (mg C m^-^2 d^-^1)')

% add floats
load('WG_floats.mat')
t = datetime(biofloat.SD5904397.M_TIME,'ConvertFrom','datenum');

P(2)=plot(coast.long,coast.lat,'k'); hold on % coast variable from floats script

for fix = 1:length(gfixs)
    float_name = float_names{fix};
    if ~isempty(biofloat.(float_name).LONGITUDE)
        P(fix+2) = plot(biofloat.(float_name).LONGITUDE,biofloat.(float_name).LATITUDE,'linewidth',3,'color',colors(fix,:));
        hold on
    end
    
end
legend(P(3:end),float_names(gfixs))
Pax=gca
Pax.FontSize=14
title('Productivity bins')

