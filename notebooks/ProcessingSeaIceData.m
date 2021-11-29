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

setup.checkregions=true;
setup.plotfigures=false;
setup.startyear=2002;
setup.endyear=2020;

%% load data
% load sea ice data

timedec=time_start_all(:,1)+(time_start_all(:,2)/12)-1/24;
datenum8day=datenum(time_start_all);
load('openshelfisobath_clean21.mat')
load('box_lat_lons.mat', 'andrex_box')


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
        eval(['box_check.',region_sublist{rix},'=lat_m;']);
        eval(['box_check.',region_sublist{rix},'(:)=0;']);
        eval(['box_check.',region_sublist{rix},'(',regionfindlist{rix},')=1;']);
        eval(['box_logic_check.',region_sublist{rix},'=logical(box_check.',region_sublist{rix},');']);
    end
    figure;
    pcolor(lon_m,lat_m,box_check.Weddell);
    shading flat
    figure;
    pcolor(lon_m,lat_m,box_check.Shelf);
    shading flat
    figure;
    pcolor(lon_m,lat_m,box_check.Open);
    shading flat
end

%% Extract sea ice data for study region 

testing=sum(g_area(findweddell));
testice=ice_conc(:,:,1);
testicewed=testice(findweddell);
icefree=find(testicewed<=0.1);
icefreearea=sum(g_area(icefree))
