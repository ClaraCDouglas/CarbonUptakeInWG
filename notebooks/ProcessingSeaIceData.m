%% Processing data brought in by ImportData_concise.m
clearvars
desktop = 1;
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
    figure;
    pcolor(g_lon,g_lat,box_check.Weddell);
    shading flat
    figure;
    pcolor(g_lon,g_lat,box_check.Shelf);
    shading flat
    figure;
    pcolor(g_lon,g_lat,box_check.Open);
    shading flat
end

%% Extract sea ice data for study region
% Area of WG, SR, OO based on sea ice grid
SeaIce.g_area.Weddell=sum(g_area(findweddell));
SeaIce.g_area.Shelf=sum(g_area(findshelf));
SeaIce.g_area.Open=sum(g_area(findopen));
SeaIce.g_area.OOSRtot=SeaIce.g_area.Open+SeaIce.g_area.Shelf;

% test run to get area of ice free (<10% ice concentration) in regions
testice=ice_conc_8day(:,:,1);
testicewed=testice(findweddell);
icefree=find(testicewed<=0.1);
icefreearea=sum(g_area(icefree))

% for rix = 1:length(region_sublist)
%     for tix=1:length(time_start_ice8)
%        tempice=ice_conc_8day(:,:,tix);
%        tempice_reg=tempice(find(regionfindlist{rix}));
%        icefreefind=find(tempice_reg<=0.1);
%        SeaIce.region_sublist{rix}.IceFreeArea(tix,1)=sum(g_area(icefreefind));
%     end
% end

for tix=1:length(time_start_ice8)
    tempice=ice_conc_8day(:,:,tix);
    tempice_reg=tempice(findweddell);
    icefreefind=find(tempice_reg<=0.1);
    SeaIce.Weddell.IceFreeArea(tix,1)=sum(g_area(icefreefind));
end
for tix=1:length(time_start_ice8)
    tempice=ice_conc_8day(:,:,tix);
    tempice_reg=tempice(findshelf);
    icefreefind=find(tempice_reg<=0.1);
    SeaIce.Shelf.IceFreeArea(tix,1)=sum(g_area(icefreefind));
end
for tix=1:length(time_start_ice8)
    tempice=ice_conc_8day(:,:,tix);
    tempice_reg=tempice(findopen);
    icefreefind=find(tempice_reg<=0.1);
    SeaIce.Open.IceFreeArea(tix,1)=sum(g_area(icefreefind));
end

figure;
t=tiledlayout(3,1)
ax1 = nexttile;
plot(datetime(time_start_ice8),SeaIce.Weddell.IceFreeArea,'LineWidth',2)
title('Weddell')
hold on
ax2=nexttile;
plot(datetime(time_start_ice8),SeaIce.Shelf.IceFreeArea,'LineWidth',2)
title('Shelf')
ax3=nexttile;
plot(datetime(time_start_ice8),SeaIce.Open.IceFreeArea,'LineWidth',2)
title('Open')
title(t,'Area of ice free water (km^2)')
% legend('Weddell','Shelf','Open');

%% Mean ice-free water per austral year
for rix = 1:length(region_sublist)
    for yix = 2002:2020
        findyear1=find(time_start_ice8(:,1)==yix-1 & time_start_ice8(:,2)>=7);
        findyear2=find(time_start_ice8(:,1)==yix & time_start_ice8(:,2)<=6);
        findyear=cat(1,findyear1,findyear2);

        SeaIce.(region_sublist{rix}).meanIceFreeArea(yix-2001,1)=yix;
        if size(findyear)<46
            SeaIce.(region_sublist{rix}).meanIceFreeArea(yix-2001,2)=NaN;
        else
            SeaIce.(region_sublist{rix}).meanIceFreeArea(yix-2001,2)=mean(SeaIce.(region_sublist{rix}).IceFreeArea(findyear,1));
        end
    end
end

figure;
t=tiledlayout(3,1)
ax1 = nexttile;
bar(SeaIce.Weddell.meanIceFreeArea(:,1),SeaIce.Weddell.meanIceFreeArea(:,2))
title('Weddell')
hold on
ax2=nexttile;
bar(SeaIce.Shelf.meanIceFreeArea(:,1),SeaIce.Shelf.meanIceFreeArea(:,2))
title('Shelf')
ax3=nexttile;
bar(SeaIce.Open.meanIceFreeArea(:,1),SeaIce.Open.meanIceFreeArea(:,2))
title('Open')
title(t,'Average ice free water per year (km^2)')
ylabel(t, 'Ice free water area (km^2)')

