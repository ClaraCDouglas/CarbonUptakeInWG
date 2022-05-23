clearvars
close all

desktop = 1;
laptop=0;
if desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg'))
    cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Matlab Add-ins'));
elseif laptop
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG'))
    cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\m_map'));
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\SetCMap'));
end

%algorithm={'cafe','cbpm','eppley','vgpm'};
algorithm={'cafe'};
setup.checkregions=true;
setup.plotfigures=false;
setup.startyear=2002;
setup.endyear=2021;
setup.monthly=false;
setup.eightday=true;
setup.yearrange=setup.startyear:1:setup.endyear; setup.yearrange=setup.yearrange';
setup.yearrange0321=setup.yearrange(2:end,:);

%% load data

for aix = 1:length(algorithm)
    %     load([algorithm{aix} '_imported.mat'], 'area*','time*',[algorithm{aix} '*']);
    %     if aix==4
    %         load([algorithm{aix} '_imported.mat'], 'VGPM*');
    %     end
    if setup.monthly
        if desktop
            cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
            load([algorithm{aix} '_imported.mat']);
        elseif laptop
            cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
            load([algorithm{aix} '_imported.mat']);
        end
    elseif setup.eightday
        if desktop
            cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
            load('cafe_8day_imported_May22.mat') % data has been corrected and saved in processed folder. ignored for commits
        elseif laptop
            cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
            load('cafe_8day_imported_May22.mat') 
        end
    end
end
clearvars algo_choice b filebase filedir


%% Allocate to ANDREX box
load('box_lat_lons.mat', 'andrex_box')
IN_and=inpolygon(lon_wg,lat_wg,andrex_box(:,1),andrex_box(:,2));
findweddell=find(IN_and==1);

region_sublist={'Weddell'};
regionfindlist= {'findweddell'};

for rix = 1:length(region_sublist)
    temp.(region_sublist{rix}).box = lat_m;
    temp.(region_sublist{rix}).box(:) = 0;
    eval(['temp.',region_sublist{rix},'.box(',regionfindlist{rix},')=1;']);
    temp.(region_sublist{rix}).box_logic=logical(temp.(region_sublist{rix}).box);
    
    AND_box_selection.NPPRate_gperyear=Z_interp.*(temp.(region_sublist{rix}).box);
end

        for yix=2003:2020
            temp_icefreedays=IceFreeDays_peryear(:,:,yix-2002);
            temp_icefreedays=round(temp_icefreedays);
            IceFree_pixels_years.(region_sublist{rix})(:,yix-2002)=temp_icefreedays(temp.(region_sublist{rix}).box_logic);
        end
%% to get high/low prod

% AnnualDayRate_mgperdmean=NaN(1080,1380,length([2003:1:2021]));
AnnualDayRate_mgperdmax=NaN(1080,1380,length([2003:1:2021]));
for yix = 2003:2021
    clearvars temp1 temp2
    temp_anntot=zeros(1080,1380);
    findyear=find(timedec8day>yix-0.5 & timedec8day<yix+0.5);
    tempNPP_yeararray=cafe_npp_all_8day(:,:,findyear); %austral year time slice of NPP av daily rates
    tempNPP_yeararray(tempNPP_yeararray<0)=NaN;
    
%     AnnualDayRate_mgperdmean(:,:,yix-2002)=mean(tempNPP_yeararray,3,'omitnan');
    AnnualDayRate_mgperdmax(:,:,yix-2002)=max(tempNPP_yeararray,[],3,'omitnan');

end



NPP_200bin=cafe_npp_all_8day;
NPP_200bin(NPP_200bin<0)=NaN; % or 0 if binning doesnt work
NPP_200bin(NPP_200bin<=200)=1;
NPP_200bin(NPP_200bin>200)=2;


setedges=[1 2];
[count_in_bin,edges,bin_bathy]=histcounts(OceanProd.vgpm.Weddell.bathymetry(findweddell),setedges);

