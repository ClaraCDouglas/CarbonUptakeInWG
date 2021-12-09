%% Processing data brought in by ImportData_concise.m
clearvars
desktop = 0;
if desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg'))
    cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Matlab Add-ins'));
else
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG'))
    cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\m_map'));
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\SetCMap'));
end

setup.checkregions=false;
setup.plotfigures=false;
setup.startyear=2002;
setup.endyear=2020;

%% load data
% load sea ice data
cd 'E:\Data\SeaIceNIMBUS';
load('SeaIce_8day_20022020.mat')
% cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
% load('cafe_8day_imported_recalc.mat', 'time_start_all');
load('openshelfisobath_clean21.mat')
load('box_lat_lons.mat', 'andrex_box')
% time_start_ice8
% time_end_ice8

% time_start_all=time_start_all(1:848,:); % NPP time_start for each 8 day slice ending in 2020
time_start_ice8=time_start_ice8(24:end,:); %there are more time slices for the sea ice data because the 2nd-4th weeks in Aug 2020 are missing from the NPP dataset
time_end_ice8=time_end_ice8(24:end,:); %there are more time slices for the sea ice data because the 2nd-4th weeks in Aug 2020 are missing from the NPP dataset
ice_conc_8day=ice_conc_8day(:,:,24:end); % remove ice data for Jan-June 2002 (no NPP data for that)

% if removing the data that is also missing from NPP data:
time_start_ice8(834:836,:)=[];
time_end_ice8(834:836,:)=[];
ice_conc_8day(:,:,834:836)=[];

datenum8day=datenum(time_start_ice8);
clearvars time_start_all

%% Import regions: shelf and open ocean then andrex box

IN_and=inpolygon(g_lon,g_lat,andrex_box(:,1),andrex_box(:,2));
findweddell=find(IN_and==1);
IN_shelf=inpolygon(g_lon,g_lat,shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2));
findshelf=find(IN_shelf==1);
IN_open=inpolygon(g_lon,g_lat,open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2));
findopen=find(IN_open==1);

region_sublist={'Weddell','Shelf','Open'};
regionfindlist= {'findweddell','findshelf','findopen'};

if setup.checkregions
    % % To check regions are within the same place
    for rix = 1:length(region_sublist)
        eval(['box_check.',region_sublist{rix},'=g_lat;']);
        eval(['box_check.',region_sublist{rix},'(:)=0;']);
        eval(['box_check.',region_sublist{rix},'(',regionfindlist{rix},')=1;']);
        eval(['box_logic_check.',region_sublist{rix},'=logical(box_check.',region_sublist{rix},');']);
    end
    
    % simple plot
    %     figure;
    %     t=tiledlayout(1,2)
    %     ax1 = nexttile;
    %     pcolor(ice_conc_8day(:,:,1)'); shading flat
    %     title('July 4th-11th 2002 Sea ice concentration')
    %     ax2 = nexttile;
    %     pcolor(IN_and'); shading flat
    %     title('ANDREX box')
    
    % whole S.Hemi view of data
    figure;
    t=tiledlayout(1,2);
    ax1 = nexttile;
    latlim = [-90 -50];
    lonlim = [-100 170];
    worldmap(latlim,lonlim)
    geoshow(antarctica)
    hold on
    pcolorm(g_lat,g_lon,ice_conc_8day(:,:,1))
    title('July 4th-11th 2002 Sea ice concentration')
    ax2=nexttile;
    worldmap(latlim,lonlim)
    hold on
    pcolorm(g_lat,g_lon,box_check.Weddell)
    %pcolorm(g_lat,g_lon,box_check.Shelf);
    %pcolorm(g_lat,g_lon,box_check.Open);
    geoshow(antarctica)
    title('ANDREX box')
end

%% Extract sea ice data for study region
    % % % % NSIDC defines ice extent/ ice cover by <15% sea ice
        % concentration (by looking at plots on
        % https://nsidc.org/data/seaice_index/)
        
        % in the sea ice extent time series plots, a 5-day trailing mean is
        % used = value for a plotted day is the average of that and the 4
        % prev days.

% plot S. Pole ice cover vs selected region
icewed=testice.*box_logic_check.Weddell;
icenotwg=find(box_logic_check.Weddell==0);
icewed(icenotwg)=NaN;

    figure;
        t=tiledlayout(1,2)
    ax1 = nexttile;

    latlim = [-90 -50];
    lonlim = [-100 170];
    worldmap(latlim,lonlim)
    geoshow(antarctica)
    hold on
    pcolorm(g_lat,g_lon,ice_conc_8day(:,:,10))

    ax2 = nexttile;

    latlim = [-90 -50];
    lonlim = [-100 170];
    worldmap(latlim,lonlim)
    geoshow(antarctica)
    hold on
    pcolorm(g_lat,g_lon,icewed)

%% Area of WG, SR, OO based on sea ice grid
SeaIce.g_area.Weddell=sum(g_area(findweddell));
% SeaIce.g_area.Weddell2=sum(sum(g_area.*box_logic_check.Weddell)); %same same
SeaIce.g_area.Shelf=sum(g_area(findshelf));
SeaIce.g_area.Open=sum(g_area(findopen));
SeaIce.g_area.OOSRtot=SeaIce.g_area.Open+SeaIce.g_area.Shelf; % same same

%% Test extent calculation
% test run to get area of ice free (<15% ice concentration) in regions
testice=ice_conc_8day(:,:,1);
testicewed=testice(findweddell);
icefree=find(testicewed<=0.15);
icefreearea=sum(g_area(icefree))
clearvars testice* icefree*

% for rix = 1:length(region_sublist)
%     for tix=1:length(time_start_ice8)
%        tempice=ice_conc_8day(:,:,tix);
%        tempice_reg=tempice(find(regionfindlist{rix}));
%        icefreefind=find(tempice_reg<=0.1);
%        SeaIce.region_sublist{rix}.IceFreeArea(tix,1)=sum(g_area(icefreefind));
%     end
% end

%% Sea Ice Extent 
    %(area of ice-covered i.e. total area of pixels that have >15% sea ice concentration)
for tix=1:length(time_start_ice8)
    tempice=ice_conc_8day(:,:,tix);
    tempice_reg=tempice(findweddell);
    g_area_reg=g_area(findweddell);
    icefind=find(tempice_reg>0.15);
    SeaIce.Weddell.SIExtent(tix,1)=sum(g_area_reg(icefind));
end
for tix=1:length(time_start_ice8)
    tempice=ice_conc_8day(:,:,tix);
    tempice_reg=tempice(findshelf);
    g_area_reg=g_area(findshelf);
    icefind=find(tempice_reg>0.15);
    SeaIce.Shelf.SIExtent(tix,1)=sum(g_area_reg(icefind));
end
for tix=1:length(time_start_ice8)
    tempice=ice_conc_8day(:,:,tix);
    tempice_reg=tempice(findopen);
    g_area_reg=g_area(findopen);
    icefind=find(tempice_reg>0.15);
    SeaIce.Open.SIExtent(tix,1)=sum(g_area_reg(icefind));
end
SeaIce.Weddell.IceFree_Extent=SeaIce.g_area.Weddell-SeaIce.Weddell.SIExtent;
SeaIce.Shelf.IceFree_Extent=SeaIce.g_area.Open-SeaIce.Shelf.SIExtent;
SeaIce.Open.IceFree_Extent=SeaIce.g_area.Open-SeaIce.Open.SIExtent;
clearvars tempice icefind tix

% timeseries figures in longform
figure;
t=tiledlayout(3,1);
ax1 = nexttile;
plot(datetime(time_start_ice8),SeaIce.Weddell.SIExtent,'LineWidth',2)
title('Weddell')
hold on
ax2=nexttile;
plot(datetime(time_start_ice8),SeaIce.Shelf.SIExtent,'LineWidth',2)
title('Shelf')
ax3=nexttile;
plot(datetime(time_start_ice8),SeaIce.Open.SIExtent,'LineWidth',2)
title('Open')
title(t,'Extent of ice free water (km^2)')

% timeseries superimposed over JASONDJFMAMJ
    % get data arranged in austral years
yearrange0320=yearrange(2:end,:);
x_months= {'J','A','S','O','N','D','J','F','M','A','M','J'};
temptime=datenum(time_start_ice8);
for yix = 1:length(yearrange0320)
    year=yearrange0320(yix)
    findyear1=find(time_start_ice8(:,1)==year-1 & time_start_ice8(:,2)>=7);
    findyear2=find(time_start_ice8(:,1)==year & time_start_ice8(:,2)<=6);
    findyearaus=cat(1,findyear1,findyear2);
    
    SeaIce.Weddell.SIExtent_aus(:,year-2002)=SeaIce.Weddell.SIExtent(findyearaus);
    SeaIce.Shelf.SIExtent_aus(:,year-2002)=SeaIce.Shelf.SIExtent(findyearaus);
    SeaIce.Open.SIExtent_aus(:,year-2002)=SeaIce.Open.SIExtent(findyearaus);
    
end 
SeaIce.Weddell.austral_meanmonth=mean(SeaIce.Weddell.SIExtent_aus,2);
    % get time for x-axis
time_start_ice8_datetime=datenum(time_start_ice8);
time_end_ice8_datetime=datenum(time_end_ice8);
time_av_ice8=mean([time_start_ice8_datetime, time_end_ice8_datetime],2)
time_av_ice8_datetime=datetime(time_av_ice8,'ConvertFrom','datenum');

linecolors = jet(length(yearrange0320));

    % plot figure - jet colourful
figure;
for cix=1:length(yearrange0320)
    plot(SeaIce.Weddell.SIExtent_aus(:,cix),'Color',linecolors(cix,:),'LineWidth',1.5);
    hold on;
end
set(gca,'xtick',[1 5 8 12 16 20 24 28 31 35 39 43],'xticklabel',x_months,'FontSize',14)
plot(SeaIce.Weddell.austral_meanmonth,'Color','k','LineWidth',3,'LineStyle',':');
xlim([0.8 46.2])
ylabel('Sea Ice Extent (km^2)','FontSize',14)
title('Weddell Gyre Sea Ice Extent','FontSize',14)
cmap = colormap(jet(length(yearrange0320)));
ticks=[0.03125:0.125:0.96875];
colorbar('Ticks',[0.03125:0.055:0.96875],'TickLabels',{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'})

    % with most lines grey
linecolors_grey=linecolors;
jetty=parula(10);
for yix = 1:length(yearrange0320)
    year=yearrange0320(yix);
    if year==2014
        linecolors_grey(yix,:)=jetty(2,:);
    elseif year ==2015
        linecolors_grey(yix,:)=jetty(4,:);
    elseif year ==2016
        linecolors_grey(yix,:)=jetty(6,:);
    elseif year ==2017
        linecolors_grey(yix,:)=jetty(8,:);
    elseif year ==2018
        linecolors_grey(yix,:)=jetty(10,:);
    else
        linecolors_grey(yix,:)=0.8;
    end
end
figure;
for cix=1:length(yearrange0320)
    plot(SeaIce.Weddell.SIExtent_aus(:,cix),'Color',linecolors_grey(cix,:),'LineWidth',2);
    hold on;
end
set(gca,'xtick',[1 5 8 12 16 20 24 28 31 35 39 43],'xticklabel',x_months,'FontSize',14)
plot(SeaIce.Weddell.austral_meanmonth,'Color','k','LineWidth',3,'LineStyle',':');
xlim([0.8 46.2])
ylabel('Sea Ice Extent (km^2)','FontSize',14)
title('Weddell Gyre Sea Ice Extent from 8-day means','FontSize',14)
cmap = colormap(linecolors_grey);
ticks=[0.03125:0.125:0.96875];
colorbar('Ticks',[0.03125:0.055:0.96875],'TickLabels',{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'})

%% Sea ice area
SeaIce.Weddell.ice_covered_area=NaN(length(time_start_ice8),1);
SeaIce.Weddell.ice_free_area=NaN(length(time_start_ice8),1);
SeaIce.Weddell.ice_covered_prop=NaN(length(time_start_ice8),1);

for ix=1:length(time_start_ice8)
    testice=ice_conc_8day(:,:,ix); % select time slice
    icefind=find(testice<=0.15);
    testice(icefind)=0;
    weightice=testice.*g_area; % get the area per pixel covered by ice
    SeaIce.Weddell.ice_covered_area(ix,1)=sum(weightice(findweddell),'omitnan'); % total area covered by ice in region
    SeaIce.Weddell.ice_free_area(ix,1)=SeaIce.g_area.Weddell-SeaIce.Weddell.ice_covered_area(ix,1);
    SeaIce.Weddell.ice_covered_prop(ix,1)=SeaIce.Weddell.ice_covered_area(ix,1)/SeaIce.g_area.Weddell; % divided by area of region to get proportion of region that is ice covered
clearvars testice weightice
end

% quick compare of SIA and SIE
figure;
plot(datetime(time_start_ice8),SeaIce.Weddell.SIExtent,'LineWidth',2)
title('Weddell Sea Ice Extent vs Sea Ice Area','FontSize',14)
hold on
plot(datetime(time_start_ice8),SeaIce.Weddell.ice_covered_area,'LineWidth',2)
legend('SIE','SIA')
ylabel('Ice covered waters (km^2)','FontSize',14)
set(gca,'FontSize',14,'XMinorGrid','on')

% timeseries superimposed over JASONDJFMAMJ
    % get data arranged in austral years
for yix = 1:length(yearrange0320)
    year=yearrange0320(yix)
    findyear1=find(time_start_ice8(:,1)==year-1 & time_start_ice8(:,2)>=7);
    findyear2=find(time_start_ice8(:,1)==year & time_start_ice8(:,2)<=6);
    findyearaus=cat(1,findyear1,findyear2);
    
    SeaIce.Weddell.SIArea_aus(:,year-2002)=SeaIce.Weddell.ice_covered_area(findyearaus);
%     SeaIce.Shelf.SIArea_aus(:,year-2002)=SeaIce.Shelf.ice_covered_area(findyearaus);
%     SeaIce.Open.SIArea_aus(:,year-2002)=SeaIce.Open.ice_covered_area(findyearaus);
end 
SeaIce.Weddell.SIArea_aus_meanmonth=mean(SeaIce.Weddell.SIArea_aus,2);

    % with most lines grey
linecolors_grey=linecolors;
jetty=parula(10);
for yix = 1:length(yearrange0320)
    year=yearrange0320(yix);
    if year==2014
        linecolors_grey(yix,:)=jetty(2,:);
    elseif year ==2015
        linecolors_grey(yix,:)=jetty(4,:);
    elseif year ==2016
        linecolors_grey(yix,:)=jetty(6,:);
    elseif year ==2017
        linecolors_grey(yix,:)=jetty(8,:);
    elseif year ==2018
        linecolors_grey(yix,:)=jetty(10,:);
    else
        linecolors_grey(yix,:)=0.8;
    end
end
figure;
for cix=1:length(yearrange0320)
    plot(SeaIce.Weddell.SIArea_aus(:,cix),'Color',linecolors_grey(cix,:),'LineWidth',2);
    hold on;
end
set(gca,'xtick',[1 5 8 12 16 20 24 28 31 35 39 43],'xticklabel',x_months,'FontSize',14)
plot(SeaIce.Weddell.SIArea_aus_meanmonth,'Color','k','LineWidth',3,'LineStyle',':');
xlim([0.8 46.2])
ylabel('Sea Ice Extent (km^2)','FontSize',14)
title('Weddell Gyre Sea Ice Area from 8-day means','FontSize',14)
cmap = colormap(linecolors_grey);
ticks=[0.03125:0.125:0.96875];
colorbar('Ticks',[0.03125:0.055:0.96875],'TickLabels',{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'})

%% Mean ice-free water per austral year
%SeaIce.Weddell.ice_free_area
%SeaIce.Weddell.IceFree_Extent
for rix = 1%:length(region_sublist)
    for yix = 2003:2020
        findyear1=find(time_start_ice8(:,1)==yix-1 & time_start_ice8(:,2)>=7);
        findyear2=find(time_start_ice8(:,1)==yix & time_start_ice8(:,2)<=6);
        findyear=cat(1,findyear1,findyear2);
        
        SeaIce.(region_sublist{rix}).meanIceFreeArea(yix-2002,1)=yix;
        SeaIce.(region_sublist{rix}).meanIceFreeExtent(yix-2002,1)=yix;
        if size(findyear)<46
            SeaIce.(region_sublist{rix}).meanIceFreeArea(yix-2002,2)=NaN;
            SeaIce.(region_sublist{rix}).meanIceFreeExtent(yix-2002,2)=NaN;
        else
            SeaIce.(region_sublist{rix}).meanIceFreeArea(yix-2002,2)=mean(SeaIce.(region_sublist{rix}).ice_free_area(findyear,1));
            SeaIce.(region_sublist{rix}).meanIceFreeExtent(yix-2002,2)=mean(SeaIce.(region_sublist{rix}).IceFree_Extent(findyear,1));
        end
    end
end



figure;
% t=tiledlayout(3,1)
% ax1 = nexttile;
bar(SeaIce.Weddell.meanIceFreeArea(:,1),SeaIce.Weddell.meanIceFreeArea(:,2))
txt=[{'Weddell: Annual mean ice-free water - calculated from SIA and SIE'};{'note, flipped colors from timeseries'}];
title(txt)
hold on
bar(SeaIce.Weddell.meanIceFreeExtent(:,1),SeaIce.Weddell.meanIceFreeExtent(:,2))
legend('SIA','SIE')
% ax2=nexttile;
% bar(SeaIce.Shelf.meanIceFreeArea(:,1),SeaIce.Shelf.meanIceFreeArea(:,2))
% title('Shelf')
% ax3=nexttile;
% bar(SeaIce.Open.meanIceFreeArea(:,1),SeaIce.Open.meanIceFreeArea(:,2))
% title('Open')
% title(t,'Average ice free water per year (km^2)')%t,
ylabel('Ice free water area (km^2)')%t,
set(gca,'FontSize',12);

%% Average sea ice area
SeaIce.ice_covered_area=NaN(length(time_start_ice8),1);
SeaIce.ice_free_area=NaN(length(time_start_ice8),1);
SeaIce.ice_covered_prop=NaN(length(time_start_ice8),1);
for ix=1:length(time_start_ice8)
    testice=ice_conc_8day(:,:,ix); % select time slice
    weightice=testice.*g_area; % get the area per pixel covered by ice
    SeaIce.ice_covered_area(ix,1)=sum(weightice(findweddell),'omitnan'); % total area covered by ice in region
    SeaIce.ice_free_area(ix,1)=SeaIce.g_area.Weddell-SeaIce.ice_covered_area(ix,1);
    SeaIce.ice_covered_prop(ix,1)=SeaIce.ice_covered_area(ix,1)/SeaIce.g_area.Weddell; % divided by area of region to get proportion of region that is ice covered
%     clearvars testice weigtice
end

figure; plot(datetime(time_start_ice8),SeaIce.ice_covered_prop);









