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

setup.checkregions=true;
setup.plotfigures=true;
setup.startyear=2002;
setup.endyear=2021;
data_daily=false;
data_8day=true;
setup.continuerun=false;
%% load data
% load sea ice data
% SeaIce_daily_20022020.mat, SeaIce_8day_20022020.mat
cd E:\Data\SeaIceV4\sidads.colorado.edu\AllFiles
if data_daily
    load('SICv4_daily.mat')
elseif data_8day
    load('SICv4_8day.mat')
end

if data_daily
    % make 3 column time start variable that matches daily seaice data
    time_start_allice=[year month day];
    yearrange=unique(year);
    clearvars year month day filenames
    % cut to start of NPP dataset - 8days starts on 4th July 2002
    % = 185 in seaice daily data
    time_start_allice=time_start_allice(185:end,:);
    ice_conc=ice_conc(:,:,185:end);
    %make dec time %(and datenum)
    timedec_allice=decyear(time_start_allice);
end

if data_8day
    % time_start_all=time_start_all(1:848,:); % NPP time_start for each 8 day slice ending in 2020
%     time_start_8day_SICv4=time_start_8day_SICv4(24:end,:); %there are more time slices for the sea ice data because the 2nd-4th weeks in Aug 2020 are missing from the NPP dataset
%     time_end_ice8=time_end_ice8(24:end,:); %there are more time slices for the sea ice data because the 2nd-4th weeks in Aug 2020 are missing from the NPP dataset
%     SICv4_8day=SICv4_8day(:,:,24:end); % remove ice data for Jan-June 2002 (no NPP data for that)
%     
%     % if removing the data that is also missing from NPP data:
%     time_start_8day_SICv4(834:836,:)=[];
%     time_end_ice8(834:836,:)=[];
%     SICv4_8day(:,:,834:836)=[];
%     
    %make dec time %(and datenum)
    timedec8dayice_start=decyear(time_start_8day_SICv4);
%     timedec8dayice_end=decyear(time_end_ice8);
%     timedec8dayice=mean([timedec8dayice_start,timedec8dayice_end],2);
    % datenum8day=datenum(time_start_8day_SICv4);
    clearvars time_start_all
end


if desktop
    cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
else
    cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
end
% load('cafe_8day_imported_recalc.mat', 'time_start_all');

%load('openshelfisobath_clean21.mat')
load('WAPSHelfOpenJan22.mat')
load('box_lat_lons.mat', 'andrex_box')
% time_start_8day_SICv4
% time_end_ice8

%% Import regions: shelf and open ocean then andrex box

BoxIn=andrex_box;
%BoxIn=[-55,-55,-30,-30;-45,-38,-38,-45]';
ShelfBox=ShelfMinusWAPJan22;%shelf_region_ANDbox%ShelfMinusWAPJan22
OOBox=OpenOceanMinusWAPJan22;%open_ocean_ANDbox%OpenOceanMinusWAPJan22
WAPBox=WAPJan22;
IN_and=inpolygon(g_lon,g_lat,BoxIn(:,1),BoxIn(:,2));
findweddell=find(IN_and==1);
IN_shelf=inpolygon(g_lon,g_lat,ShelfBox(:,1),ShelfBox(:,2));
findshelf=find(IN_shelf==1);
IN_open=inpolygon(g_lon,g_lat,OOBox(:,1),OOBox(:,2));
findopen=find(IN_open==1);
IN_wap=inpolygon(g_lon,g_lat,WAPBox(:,1),WAPBox(:,2));
findWAP=find(IN_wap==1);

region_sublist={'Weddell','Shelf','Open','AP'}; %
regionfindlist= {'findweddell','findshelf','findopen','findWAP'}; %

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
    %     pcolor(SICv4_8day(:,:,1)'); shading flat
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
    antarctica = shaperead('landareas.shp', 'UseGeoCoords', true,...
        'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});
    worldmap(latlim,lonlim)
    %     patchm(antarctica.Lat, antarctica.Lon, [0.5 1 0.5])
    %     geoshow(antarctica)
    hold on
    if data_8day
        pcolorm(g_lat,g_lon,SICv4_8day(:,:,1))
    elseif data_daily
        pcolorm(g_lat,g_lon,SICv4_0221(:,:,1))
    end
    title('July 4th-11th 2002 Sea ice concentration')
    ax2=nexttile;
    worldmap(latlim,lonlim)
    hold on
    pcolorm(g_lat,g_lon,box_check.Weddell)
    pcolorm(g_lat,g_lon,box_check.Shelf);
    pcolorm(g_lat,g_lon,box_check.Open);
    pcolorm(g_lat,g_lon,box_check.AP);
    geoshow(antarctica)
    title('ANDREX box')
end

%% Test extent calculation
% test run to get area of ice free (<15% ice concentration) in regions
if data_8day
    testice=SICv4_8day(:,:,1);
elseif data_daily
    testice=SICv4_0221(:,:,1);
end
testicewed=testice(findweddell);
icefree=find(testicewed<=0.15);
icefreearea=sum(g_area(icefree))

% for rix = 1:length(region_sublist)
%     for tix=1:length(time_start_8day_SICv4)
%        tempice=SICv4_8day(:,:,tix);
%        tempice_reg=tempice(find(regionfindlist{rix}));
%        icefreefind=find(tempice_reg<=0.1);
%        SeaIce.region_sublist{rix}.IceFreeArea(tix,1)=sum(g_area(icefreefind));
%     end
% end
%% Extract sea ice data for study region
% % % % NSIDC defines ice extent/ ice cover by <15% sea ice
% concentration (by looking at plots on
% https://nsidc.org/data/seaice_index/)

% in the sea ice extent time series plots, a 5-day trailing mean is
% used = value for a plotted day is the average of that and the 4
% prev days.

% plot S. Pole ice cover vs selected region
if setup.checkregions
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
    if data_8day
        pcolorm(g_lat,g_lon,SICv4_8day(:,:,1))
    elseif data_daily
        pcolorm(g_lat,g_lon,ice_conc(:,:,1))
    end
    ax2 = nexttile;
    latlim = [-90 -50];
    lonlim = [-100 170];
    worldmap(latlim,lonlim)
    geoshow(antarctica)
    hold on
    pcolorm(g_lat,g_lon,icewed)
end

clearvars testice* icefree*

%% Area of WG, SR, OO based on sea ice grid
% SeaIce.g_area.Weddell2=sum(sum(g_area.*box_logic_check.Weddell)); %same same
SeaIce.g_area.Weddell=sum(g_area(findweddell));
SeaIce.g_area.Shelf=sum(g_area(findshelf));
SeaIce.g_area.Open=sum(g_area(findopen));
SeaIce.g_area.AP=sum(g_area(findWAP));
%SeaIce.g_area.OOSRtot=SeaIce.g_area.Open+SeaIce.g_area.Shelf; % same same
SeaIce.g_area.OOSRWAPtot=SeaIce.g_area.Open+SeaIce.g_area.Shelf+SeaIce.g_area.AP; % same same

SeaIce.g_area.Weddell-SeaIce.g_area.OOSRWAPtot

%% 8day analysis
if data_8day
    for tix=1:length(time_start_8day_SICv4)
        tempice=SICv4_8day(:,:,tix);
        tempice_reg=tempice(findweddell);
        g_area_reg=g_area(findweddell);
        icefind=find(tempice_reg>0.15);
        SeaIce.Weddell.SIExtent(tix,1)=sum(g_area_reg(icefind));
    end
    for tix=1:length(time_start_8day_SICv4)
        tempice=SICv4_8day(:,:,tix);
        tempice_reg=tempice(findshelf);
        g_area_reg=g_area(findshelf);
        icefind=find(tempice_reg>0.15);
        SeaIce.Shelf.SIExtent(tix,1)=sum(g_area_reg(icefind));
    end
    for tix=1:length(time_start_8day_SICv4)
        tempice=SICv4_8day(:,:,tix);
        tempice_reg=tempice(findopen);
        g_area_reg=g_area(findopen);
        icefind=find(tempice_reg>0.15);
        SeaIce.Open.SIExtent(tix,1)=sum(g_area_reg(icefind));
    end
    for tix=1:length(time_start_8day_SICv4)
        tempice=SICv4_8day(:,:,tix);
        tempice_reg=tempice(findWAP);
        g_area_reg=g_area(findWAP);
        icefind=find(tempice_reg>0.15);
        SeaIce.AP.SIExtent(tix,1)=sum(g_area_reg(icefind));
    end
    for rix = 1:length(region_sublist)
        SeaIce.(region_sublist{rix}).IceFree_Extent=SeaIce.g_area.(region_sublist{rix})-SeaIce.(region_sublist{rix}).SIExtent;
    end
    clearvars tempice icefind tix
    
    if setup.plotfigures
        % timeseries figures in longform
        figure;
        t=tiledlayout(4,1);
        ax1 = nexttile;
        plot(datetime(time_start_8day_SICv4),SeaIce.Weddell.SIExtent,'LineWidth',2)
        title('Weddell')
        hold on
        ax2=nexttile;
        plot(datetime(time_start_8day_SICv4),SeaIce.Shelf.SIExtent,'LineWidth',2)
        title('Shelf')
        ax3=nexttile;
        plot(datetime(time_start_8day_SICv4),SeaIce.Open.SIExtent,'LineWidth',2)
        title('Open')
        ax4=nexttile;
        plot(datetime(time_start_8day_SICv4),SeaIce.AP.SIExtent,'LineWidth',2)
        title('AP')
        title(t,'Sea Ice Extent (km^2)')
    end
    
    % timeseries superimposed over JASONDJFMAMJ
    % get data arranged in austral years
    yearrange=unique(2002:2021);
    yearrange=yearrange';
    yearrange0321=yearrange(2:end,:);
    temptime=datenum(time_start_8day_SICv4);
    for rix = 1:length(region_sublist)
        SeaIce.(region_sublist{rix}).SIExtent_aus=NaN(46,length(yearrange0321));
        for yix = 2003:2021
            findyearaus=find(timedec8dayice_start>yix-0.5 & timedec8dayice_start<yix+0.5);
            year = yix;
            SeaIce.(region_sublist{rix}).SIExtent_aus(1:length(findyearaus),year-2002)=SeaIce.(region_sublist{rix}).SIExtent(findyearaus);
        end
        SeaIce.(region_sublist{rix}).austral_meanmonth=mean(SeaIce.(region_sublist{rix}).SIExtent_aus,2);
    end
    time_start_8day_SICv4_datetime=datenum(time_start_8day_SICv4);
%     time_end_ice8_datetime=datenum(time_end_ice8);
%     time_av_ice8=mean([time_start_8day_SICv4_datetime, time_end_ice8_datetime],2);
%     time_av_ice8_datetime=datetime(time_av_ice8,'ConvertFrom','datenum');

    % get time for x-axis
    if setup.plotfigures
        linecolors = jet(length(yearrange0321));
        % plot figure - jet colourful
        x_months= {'J','A','S','O','N','D','J','F','M','A','M','J'};
        figure;
        for cix=1:length(yearrange0321)
            plot(SeaIce.Weddell.SIExtent_aus(:,cix),'Color',linecolors(cix,:),'LineWidth',1.5);
            hold on;
        end
        set(gca,'xtick',[1 5 8 12 16 20 24 28 31 35 39 43],'xticklabel',x_months,'FontSize',14)
        plot(SeaIce.Weddell.austral_meanmonth,'Color','k','LineWidth',3,'LineStyle',':');
        %xlim([0.8 366.2])
        ylabel('Sea Ice Extent (km^2)','FontSize',14)
        title('Weddell Gyre Sea Ice Extent','FontSize',14)
        cmap = colormap(jet(length(yearrange0321)));
        ticks=[0.03125:0.125:0.96875];
        colorbar('Ticks',[0.03125:0.055:0.96875],'TickLabels',{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'})
        
        % with most lines grey
        linecolors_grey=linecolors;
        jetty=parula(10);
        for yix = 1:length(yearrange0321)
            year=yearrange0321(yix);
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
        for cix=1:length(yearrange0321)
            plot(SeaIce.Weddell.SIExtent_aus(:,cix),'Color',linecolors_grey(cix,:),'LineWidth',2);
            hold on;
        end
        set(gca,'xtick',[1 5 8 12 16 20 24 28 31 35 39 43],'xticklabel',x_months,'FontSize',14)
        plot(SeaIce.Weddell.austral_meanmonth,'Color','k','LineWidth',3,'LineStyle',':');
        %xlim([0.8 366])
        ylabel('Sea Ice Extent (km^2)','FontSize',14)
        title('Weddell Gyre Sea Ice Extent from 8-day means','FontSize',14)
        cmap = colormap(linecolors_grey);
        ticks=[0.03125:0.125:0.96875];
        colorbar('Ticks',[0.03125:0.055:0.96875],'TickLabels',{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'})
    end
    %% Sea ice area
    for rix = 1:length(region_sublist)
        SeaIce.(region_sublist{rix}).ice_covered_area=NaN(length(time_start_8day_SICv4),1);
        SeaIce.(region_sublist{rix}).ice_free_area=NaN(length(time_start_8day_SICv4),1);
        SeaIce.(region_sublist{rix}).ice_covered_prop=NaN(length(time_start_8day_SICv4),1);
    end

    for rix = 1:length(region_sublist)
        for ix=1:length(time_start_8day_SICv4)
            testice=SICv4_8day(:,:,ix); % select time slice
            icefind=find(testice<=0.15);
            testice(icefind)=0;
            weightice=testice.*g_area; % get the area per pixel covered by ice
            SeaIce.Weddell.ice_covered_area(ix,1)=sum(weightice(findweddell),'omitnan'); % total area covered by ice in region
            SeaIce.Shelf.ice_covered_area(ix,1)=sum(weightice(findshelf),'omitnan'); % total area covered by ice in region
            SeaIce.Open.ice_covered_area(ix,1)=sum(weightice(findopen),'omitnan'); % total area covered by ice in region
            SeaIce.AP.ice_covered_area(ix,1)=sum(weightice(findWAP),'omitnan'); % total area covered by ice in region
            
            SeaIce.(region_sublist{rix}).ice_free_area(ix,1)=SeaIce.g_area.(region_sublist{rix})-SeaIce.(region_sublist{rix}).ice_covered_area(ix,1);
            SeaIce.(region_sublist{rix}).ice_covered_prop(ix,1)=SeaIce.(region_sublist{rix}).ice_covered_area(ix,1)/SeaIce.g_area.(region_sublist{rix}); % divided by area of region to get proportion of region that is ice covered
            clearvars testice weightice
        end
    end
    
    if setup.plotfigures
        % quick compare of SIA and SIE
        figure;
        plot(datetime(time_start_8day_SICv4),SeaIce.Weddell.SIExtent,'LineWidth',2)
        title('Weddell Sea Ice Extent vs Sea Ice Area','FontSize',14)
        hold on
        plot(datetime(time_start_8day_SICv4),SeaIce.Weddell.ice_covered_area,'LineWidth',2)
        legend('SIE','SIA')
        ylabel('Ice covered waters (km^2)','FontSize',14)
        set(gca,'FontSize',14,'XMinorGrid','on')
    end
    % timeseries superimposed over JASONDJFMAMJ
    % get data arranged in austral years
    
    for rix = 1:length(region_sublist)
        SeaIce.(region_sublist{rix}).SIArea_aus=NaN(46,length(yearrange0321));
        for yix = 2003:2021
            findyearaus=find(timedec8dayice_start>yix-0.5 & timedec8dayice_start<yix+0.5);
            year = yix;
            SeaIce.(region_sublist{rix}).SIArea_aus(1:length(findyearaus),year-2002)=SeaIce.(region_sublist{rix}).ice_covered_area(findyearaus);
        end
        SeaIce.(region_sublist{rix}).SIArea_aus_meanmonth=mean(SeaIce.(region_sublist{rix}).SIArea_aus,2);
    end
    
    if setup.plotfigures
        % with most lines grey
        linecolors_grey=linecolors;
        jetty=parula(10);
        for yix = 1:length(yearrange0321)
            year=yearrange0321(yix);
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
        for cix=1:length(yearrange0321)
            cix;
            plot(SeaIce.Weddell.SIArea_aus(:,cix),'Color',linecolors_grey(cix,:),'LineWidth',2);
            hold on;
        end
        set(gca,'xtick',[1 5 8 12 16 20 24 28 31 35 39 43],'xticklabel',x_months,'FontSize',14)
        plot(SeaIce.Weddell.SIArea_aus_meanmonth,'Color','k','LineWidth',3,'LineStyle',':');
%         xlim([0.8 366])
        ylabel('Sea Ice Extent (km^2)','FontSize',14)
        title('Weddell Gyre Sea Ice Area from 8-day means','FontSize',14)
        cmap = colormap(linecolors_grey);
        ticks=[0.03125:0.125:0.96875];
        colorbar('Ticks',[0.03125:0.055:0.96875],'TickLabels',{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'})
    end
    % clearvars -except data* SeaIce time* yearrange*
    
    desktop=0;
    laptop=1;
    if desktop
        cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
    elseif laptop
        cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
    end
    clearvars desktop laptop filenames BoxIN *ix *Box g* find* temp* ticks
    save('SeaIcev4_8day_May22.mat','SeaIce','time_start_8day_SICv4','timedec8dayice_start','data_8day','data_daily'); %SeaIceDaily_Jan22_wWAP
end

%% daily analysis
if data_daily
    for tix=1:length(time_start_allice)
        tempice=ice_conc(:,:,tix);
        tempice_reg=tempice(findweddell);
        g_area_reg=g_area(findweddell);
        icefind=find(tempice_reg>0.15);
        SeaIce.Weddell.SIExtent(tix,1)=sum(g_area_reg(icefind));
    end
    for tix=1:length(time_start_allice)
        tempice=ice_conc(:,:,tix);
        tempice_reg=tempice(findshelf);
        g_area_reg=g_area(findshelf);
        icefind=find(tempice_reg>0.15);
        SeaIce.Shelf.SIExtent(tix,1)=sum(g_area_reg(icefind));
    end
    for tix=1:length(time_start_allice)
        tempice=ice_conc(:,:,tix);
        tempice_reg=tempice(findopen);
        g_area_reg=g_area(findopen);
        icefind=find(tempice_reg>0.15);
        SeaIce.Open.SIExtent(tix,1)=sum(g_area_reg(icefind));
    end
    for tix=1:length(time_start_allice)
        tempice=ice_conc(:,:,tix);
        tempice_reg=tempice(findWAP);
        g_area_reg=g_area(findWAP);
        icefind=find(tempice_reg>0.15);
        SeaIce.WAP.SIExtent(tix,1)=sum(g_area_reg(icefind));
    end
    for rix = 1:length(region_sublist)
        SeaIce.(region_sublist{rix}).IceFree_Extent=SeaIce.g_area.(region_sublist{rix})-SeaIce.(region_sublist{rix}).SIExtent;
    end
    clearvars tempice icefind tix
    
    if setup.plotfigures
        % timeseries figures in longform
        figure;
        t=tiledlayout(4,1);
        ax1 = nexttile;
        plot(datetime(time_start_allice),SeaIce.Weddell.SIExtent,'LineWidth',2)
        title('Weddell')
        hold on
        ax2=nexttile;
        plot(datetime(time_start_allice),SeaIce.Shelf.SIExtent,'LineWidth',2)
        title('Shelf')
        ax3=nexttile;
        plot(datetime(time_start_allice),SeaIce.Open.SIExtent,'LineWidth',2)
        title('Open')
        ax4=nexttile;
        plot(datetime(time_start_allice),SeaIce.WAP.SIExtent,'LineWidth',2)
        title('WAP')
        title(t,'Sea Ice Extent (km^2)')
    end
    
    % timeseries superimposed over JASONDJFMAMJ
    % get data arranged in austral years
    yearrange0321=yearrange(2:end,:);
    temptime=datenum(time_start_allice);
    for rix = 1:length(region_sublist)
        SeaIce.(region_sublist{rix}).SIExtent_aus=NaN(366,length(yearrange0321));
        for yix = 2003:2021
            findyearaus=find(timedec_allice>yix-0.5 & timedec_allice<yix+0.5);
            year = yix;
            SeaIce.(region_sublist{rix}).SIExtent_aus(1:length(findyearaus),year-2002)=SeaIce.(region_sublist{rix}).SIExtent(findyearaus);
        end
        SeaIce.(region_sublist{rix}).austral_meanmonth=mean(SeaIce.(region_sublist{rix}).SIExtent_aus,2);
    end
    
    % get time for x-axis
    if setup.plotfigures
        linecolors = jet(length(yearrange0321));
        % plot figure - jet colourful
        x_months= {'J','A','S','O','N','D','J','F','M','A','M','J'};
        figure;
        for cix=1:length(yearrange0321)
            plot(SeaIce.Weddell.SIExtent_aus(:,cix),'Color',linecolors(cix,:),'LineWidth',1.5);
            hold on;
        end
        %     set(gca,'xtick',[1 5 8 12 16 20 24 28 31 35 39 43],'xticklabel',x_months,'FontSize',14)
        plot(SeaIce.Weddell.austral_meanmonth,'Color','k','LineWidth',3,'LineStyle',':');
        xlim([0.8 366.2])
        ylabel('Sea Ice Extent (km^2)','FontSize',14)
        title('Weddell Gyre Sea Ice Extent','FontSize',14)
        cmap = colormap(jet(length(yearrange0321)));
        ticks=[0.03125:0.125:0.96875];
        colorbar('Ticks',[0.03125:0.055:0.96875],'TickLabels',{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'})
        
        % with most lines grey
        linecolors_grey=linecolors;
        jetty=parula(10);
        for yix = 1:length(yearrange0321)
            year=yearrange0321(yix);
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
        for cix=1:length(yearrange0321)
            plot(SeaIce.Weddell.SIExtent_aus(:,cix),'Color',linecolors_grey(cix,:),'LineWidth',2);
            hold on;
        end
        set(gca,'xtick',[1:30:365],'FontSize',14)
        plot(SeaIce.Weddell.austral_meanmonth,'Color','k','LineWidth',3,'LineStyle',':');
        xlim([0.8 366])
        ylabel('Sea Ice Extent (km^2)','FontSize',14)
        title('Weddell Gyre Sea Ice Extent from 8-day means','FontSize',14)
        cmap = colormap(linecolors_grey);
        ticks=[0.03125:0.125:0.96875];
        colorbar('Ticks',[0.03125:0.055:0.96875],'TickLabels',{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'})
    end
    %% Sea ice area
    for rix = 1:length(region_sublist)
        SeaIce.(region_sublist{rix}).ice_covered_area=NaN(length(time_start_allice),1);
        SeaIce.(region_sublist{rix}).ice_free_area=NaN(length(time_start_allice),1);
        SeaIce.(region_sublist{rix}).ice_covered_prop=NaN(length(time_start_allice),1);
    end

    for rix = 1:length(region_sublist)
        for ix=1:length(time_start_allice)
            testice=ice_conc(:,:,ix); % select time slice
            icefind=find(testice<=0.15);
            testice(icefind)=0;
            weightice=testice.*g_area; % get the area per pixel covered by ice
            SeaIce.Weddell.ice_covered_area(ix,1)=sum(weightice(findweddell),'omitnan'); % total area covered by ice in region
            SeaIce.Shelf.ice_covered_area(ix,1)=sum(weightice(findshelf),'omitnan'); % total area covered by ice in region
            SeaIce.Open.ice_covered_area(ix,1)=sum(weightice(findopen),'omitnan'); % total area covered by ice in region
            SeaIce.WAP.ice_covered_area(ix,1)=sum(weightice(findWAP),'omitnan'); % total area covered by ice in region
            
            SeaIce.(region_sublist{rix}).ice_free_area(ix,1)=SeaIce.g_area.(region_sublist{rix})-SeaIce.(region_sublist{rix}).ice_covered_area(ix,1);
            SeaIce.(region_sublist{rix}).ice_covered_prop(ix,1)=SeaIce.(region_sublist{rix}).ice_covered_area(ix,1)/SeaIce.g_area.(region_sublist{rix}); % divided by area of region to get proportion of region that is ice covered
            clearvars testice weightice
        end
    end
    
    if setup.plotfigures
        % quick compare of SIA and SIE
        figure;
        plot(datetime(time_start_allice),SeaIce.Weddell.SIExtent,'LineWidth',2)
        title('Weddell Sea Ice Extent vs Sea Ice Area','FontSize',14)
        hold on
        plot(datetime(time_start_allice),SeaIce.Weddell.ice_covered_area,'LineWidth',2)
        legend('SIE','SIA')
        ylabel('Ice covered waters (km^2)','FontSize',14)
        set(gca,'FontSize',14,'XMinorGrid','on')
    end
    % timeseries superimposed over JASONDJFMAMJ
    % get data arranged in austral years
    
    for rix = 1:length(region_sublist)
        SeaIce.(region_sublist{rix}).SIArea_aus=NaN(366,length(yearrange0321));
        for yix = 2003:2021
            findyearaus=find(timedec_allice>yix-0.5 & timedec_allice<yix+0.5);
            year = yix;
            SeaIce.(region_sublist{rix}).SIArea_aus(1:length(findyearaus),year-2002)=SeaIce.(region_sublist{rix}).ice_covered_area(findyearaus);
        end
        SeaIce.(region_sublist{rix}).SIArea_aus_meanmonth=mean(SeaIce.(region_sublist{rix}).SIArea_aus,2);
    end
    
    if setup.plotfigures
        % with most lines grey
        linecolors_grey=linecolors;
        jetty=parula(10);
        for yix = 1:length(yearrange0321)
            year=yearrange0321(yix);
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
        for cix=1:length(yearrange0321)
            cix;
            plot(SeaIce.Weddell.SIArea_aus(:,cix),'Color',linecolors_grey(cix,:),'LineWidth',2);
            hold on;
        end
        set(gca,'xtick',[1:30:365],'FontSize',14)
        plot(SeaIce.Weddell.SIArea_aus_meanmonth,'Color','k','LineWidth',3,'LineStyle',':');
        xlim([0.8 366])
        ylabel('Sea Ice Extent (km^2)','FontSize',14)
        title('Weddell Gyre Sea Ice Area from 8-day means','FontSize',14)
        cmap = colormap(linecolors_grey);
        ticks=[0.03125:0.125:0.96875];
        colorbar('Ticks',[0.03125:0.055:0.96875],'TickLabels',{'2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'})
    end
    % clearvars -except data* SeaIce time* yearrange*
    
    desktop=0;
    laptop=1;
    if desktop
        cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
    elseif laptop
        cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
    end
    clearvars desktop laptop
    save('SeaIceDaily_Jan22_wWAPauto.mat','SeaIce','time_start_allice','timedec_allice','yearrange','yearrange0320','data_8day','data_daily'); %SeaIceDaily_Jan22_wWAP
    
end
clearvars ans antarctica ax* box* *ix filenames *find* t temp* IN* icenotwg icewed
if setup.continuerun
%% Ice free waters
% quick compare of SIA and SIE
figure;
plot(timedec8dayice,SeaIce.Weddell.ice_free_area,'LineWidth',2)
title('Weddell: Annual mean ice-free water - calculated from SIA and SIE and NPP data','FontSize',14)
hold on
plot(timedec8dayice,SeaIce.Weddell.IceFree_Extent,'LineWidth',2)
plot(timedec8day,OceanProd_8day.cafe.Weddell.area_8day_km2,'LineWidth',2)
legend('SIA','SIE','NPP ice-free')
xlim([2002.8 2021.2])
ylabel('Ice free waters (km^2)','FontSize',14)
set(gca,'FontSize',14,'XMinorGrid','on')

SIEvNPP=SeaIce.Weddell.IceFree_Extent(1:848)-OceanProd_8day.cafe.Weddell.area_8day_km2(1:848);
SIAvNPP=SeaIce.Weddell.ice_free_area(1:848)-OceanProd_8day.cafe.Weddell.area_8day_km2(1:848);
SIAvSIE=SeaIce.Weddell.ice_free_area(1:848)-SeaIce.Weddell.IceFree_Extent(1:848);

figure;
plot(timedec8dayice,SIAvSIE,'LineWidth',2)
title('Weddell: Annual mean ice-free water - calculated from SIA and SIE and NPP data','FontSize',14)
hold on
plot(timedec8dayice,SIAvNPP,'LineWidth',2)
plot(timedec8dayice,SIEvNPP,'LineWidth',2)
legend('SIA','SIE','NPP ice-free')
xlim([2002.8 2021.2])
ylabel('Ice free waters (km^2)','FontSize',14)
set(gca,'FontSize',14,'XMinorGrid','on')
ylim([0 inf])

%% Mean ice-free water per austral year
%SeaIce.Weddell.ice_free_area
%SeaIce.Weddell.IceFree_Extent
for rix = 1%:length(region_sublist)
    for yix = 2003:2021
        findyear1=find(time_start_8day_SICv4(:,1)==yix-1 & time_start_8day_SICv4(:,2)>=7);
        findyear2=find(time_start_8day_SICv4(:,1)==yix & time_start_8day_SICv4(:,2)<=6);
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
txt=[{'Weddell: Annual mean ice-free water - calculated from SIA and SIE and NPP ice-free'};{'note, flipped colors from timeseries'}];
title(txt)
hold on
bar(SeaIce.Weddell.meanIceFreeExtent(:,1),SeaIce.Weddell.meanIceFreeExtent(:,2))
bar(OceanProd_8day.cafe.Weddell.IceFree_annualMEAN(:,1),OceanProd_8day.cafe.Weddell.IceFree_annualMEAN(:,2))
legend('SIA','SIE','NPP Ice-Free')
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
SeaIce.ice_covered_area=NaN(length(time_start_8day_SICv4),1);
SeaIce.ice_free_area=NaN(length(time_start_8day_SICv4),1);
SeaIce.ice_covered_prop=NaN(length(time_start_8day_SICv4),1);
for ix=1:length(time_start_8day_SICv4)
    testice=SICv4_8day(:,:,ix); % select time slice
    weightice=testice.*g_area; % get the area per pixel covered by ice
    SeaIce.ice_covered_area(ix,1)=sum(weightice(findweddell),'omitnan'); % total area covered by ice in region
    SeaIce.ice_free_area(ix,1)=SeaIce.g_area.Weddell-SeaIce.ice_covered_area(ix,1);
    SeaIce.ice_covered_prop(ix,1)=SeaIce.ice_covered_area(ix,1)/SeaIce.g_area.Weddell; % divided by area of region to get proportion of region that is ice covered
    %     clearvars testice weigtice
end

figure; plot(datetime(time_start_8day_SICv4),SeaIce.ice_covered_prop);

end







