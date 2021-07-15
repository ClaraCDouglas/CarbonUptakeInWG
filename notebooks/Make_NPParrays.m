clearvars
desktop =0; 
laptop =0;
remote=1;
if desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg'))
    cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Matlab Add-ins'));
elseif laptop
    cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
elseif remote
    addpath(genpath('/noc/users/ccd1n18/Documents/Projects/CarbonUptakeInWG'))
    cd '/noc/users/ccd1n18/Documents/Projects/CarbonUptakeInWG/data/processed' % desktop
end
load('vgpm_npp_all.mat') % vgpm_npp_all as from OceanProductivity site - average (A-W) daily rates per month
load('cafe_npp_all.mat') % cafe_npp_all as from OceanProductivity site - average (A-W) daily rates per month
load('cbpm_npp_all.mat') % cbpm_npp_all as from OceanProductivity site - average (A-W) daily rates per month
load('eppley_npp_all.mat') % eppley_npp_all as from OceanProductivity site - average (A-W) daily rates per month
load('time_start_all.mat')
load('ProcessedData.mat', 'timedec')

algorithm={'cafe','cbpm','eppley','vgpm'};
temp.cafe.rates=cafe_npp_all;
temp.cbpm.rates=cafe_npp_all;
temp.eppley.rates=eppley_npp_all;
temp.vgpm.rates=vgpm_npp_all;
clearvars cafe_npp_all cbpm_npp_all eppley_npp_all vgpm_npp_all

%% NPP area-weighted rate 
for aix = 1:length(algorithm)
    temp.findneg=find(temp.(algorithm{aix}).rates<0);
    temp.(algorithm{aix}).rates(temp.findneg)=NaN;
    % annual climatology of average daily rates per year
    NPP.(algorithm{aix}).annual_day_nan(:,:,length([2003:2019]))=NaN(size(temp.(algorithm{aix}).rates(:,:,1)));
    for yix = 2003:2019
        findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
        if size (findyear)<12
            NPP_years.(algorithm{aix}).annual_day_nan(:,:,yix-2002)=NaN(size(temp.(algorithm{aix}).rates(:,:,1)));
        else
            NPP_years.(algorithm{aix}).annual_day_nan(:,:,yix-2002)=nanmean(temp.(algorithm{aix}).rates(:,:,findyear),3);
        end
    end
    NPP.(algorithm{aix}).annual_av_day_nan=nanmean(NPP_years.(algorithm{aix}).annual_day_nan,3);
temp = rmfield(temp, 'findneg');
end
%% Januarys/monthly rates

% January maps
for aix = 1:length(algorithm)
    for mix = 1%:12
        findmonth=find(time_start_all(:,2)==mix);
        NPP_years.(algorithm{aix}).jan_rates=temp.(algorithm{aix}).rates(:,:,findmonth);
    end
end
clearvars findmonth

month={'jan','feb','march','apr','may','june','july','aug','sept','oct','nov','dec'};
% Monthly maps
for aix = 1:length(algorithm)
    for mix = 1:12%[1:1:4 8:1:12]
        findmonth=find(time_start_all(:,2)==mix);
        NPP_years.(algorithm{aix}).month_rates.(month{mix})=temp.(algorithm{aix})(:,:,findmonth);
    end
end
clearvars findmonth

%% Total NPP 
% VGPM_npp_tot_gC_nans - land/permanent ice are set as NaNs
load('vgpm_imported.mat', 'VGPM_npp_tot_gC_nans')
    % annual climatology of total NPP per year
temp.vgpmtot=VGPM_npp_tot_gC_nans;

temp.cafe.totals=cafe_npp_all;
temp.cbpm.totals=cafe_npp_all;
temp.eppley.totals=eppley_npp_all;
temp.vgpm2.totals=load('vgpm_imported.mat', 'VGPM_npp_tot_gC_nans');

% vgpm_av_tot_nan=nanmean(temp.vgpmtot,3); % this has calculated the average total NPP in each pixel
for aix = 1:length(algorithm)
NPP_years.(algorithm{aix}).tot_years(:,:,length([2003:2019]))=NaN(size(temp.(algorithm{aix}).totals(:,:,1)));    
for yix = 2003:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    if size (findyear)<12
    NPP_years.vgpm_tot_years(:,:,yix-2002)=NaN(size(temp.vgpmtot(:,:,1)));
    else
    NPP_years.vgpm_tot_years(:,:,yix-2002)=nansum(temp.vgpmtot(:,:,findyear),3);    
    end
end
end 
NPP.vgpm_tot_av_years=nanmean(NPP_years.vgpm_tot_years,3);
NPP.vgpm_tot_av_years(find(NPP.vgpm_tot_av_years==0))=NaN;

%% make WG box
load('vgpm_imported.mat', 'area_MODISVGPM_m2')
load('latlon_m.mat')
load('box_lat_lons.mat', 'andrex_box')
    % make box
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
end
