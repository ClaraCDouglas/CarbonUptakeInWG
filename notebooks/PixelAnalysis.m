%% pixel by pixel analysis

%% load data
clearvars
close all

desktop = 0;
laptop=1;
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
setup.endyear=2020;
setup.monthly=false;
setup.eightday=true;
setup.yearrange=setup.startyear:1:setup.endyear; setup.yearrange=setup.yearrange';
setup.yearrange0320=setup.yearrange(2:end,:);

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
            load('cafe_8day_imported_eqWG_withNaNs.mat') % data has been corrected and saved in processed folder. ignored for commits
        elseif laptop
            cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
            load('cafe_8day_imported_eqWG_withNaNs.mat') % data has been corrected and saved in processed folder. ignored for commits
        end
    end
end
% vgpm_npp_tot_gC_all=VGPM_npp_tot_gC_all;
% vgpm_npp_tot_gC_nans=VGPM_npp_tot_gC_nans;
clearvars VGPM_npp_tot_gC_nans VGPM_npp_tot_gC_all algo_choice b filebase filedir *vgpm*
if setup.monthly
    timedec=time_start_all(:,1)+(time_start_all(:,2)/12)-1/24;
end
if setup.eightday
    %     datenum8day=datenum(time_start_all);
    timedec8day=decyear(time_start_all);
    timedec8day_end=decyear(time_end_all);
end


%% day chunks
leapyear1=time_start_all(:,1);
leapyear2=rem(leapyear1,4)==0;
leapyear3=zeros(length(leapyear1),1);
leapyear3(leapyear2==1)=366;
leapyear3(leapyear3==0)=365;
days_test(:,1)=(timedec8day_end-timedec8day);
days_test(:,2)=days_test(:,1)*365;
days_test(:,3)=days_test(:,1).*leapyear3;
day_chunk=round(days_test(:,1).*leapyear3);

clearvars leap* *test
%% Ice Free Days per pixel (map)
timedec8day_mid=mean([timedec8day,timedec8day_end],2);
% temp4=zeros(1080,1380);
IceFreeDays_peryear=NaN(1080,1380,length([2003:1:2020]));
for yix = 2003:2020
    temp4=zeros(1080,1380);
    findyear=find(timedec8day_mid>yix-0.5 & timedec8day_mid<yix+0.5);
    tempNPP_yeararray=cafe_npp_all_8day(:,:,findyear); %austral year time slice of NPP av daily rates
    %     tempNPP_yeararray(tempNPP_yeararray<0)=NaN;
    temp_daychunk=day_chunk(findyear);
    daystot=sum(temp_daychunk);
    if daystot<365
        error('days tot less than a year')
    end
    for qix=1:1:length(findyear)
        temparray=tempNPP_yeararray(:,:,qix); % one by one, go through 8 day slices for each aus year
        temp2=ones(size(temparray)); % make blank array
        temp2(find(temparray<0))=0; % put in 1s where there is data
        temp3=temp2*temp_daychunk(qix); % multiply 1s (data present) by the number of days in the time slice = number of days ice free
        temp4=temp4+temp3;
    end
    IceFreeDays_peryear(:,:,yix-2002)=temp4;
end

clearvars temp*
%% Annual rate
AnnualNPPRate_mgperyear=NaN(1080,1380,length([2003:1:2020]));
for yix = 2003:2020
    clearvars temp1 temp2
    temp_anntot=zeros(1080,1380);
    findyear=find(timedec8day_mid>yix-0.5 & timedec8day_mid<yix+0.5);
    tempNPP_yeararray=cafe_npp_all_8day(:,:,findyear); %austral year time slice of NPP av daily rates
    tempNPP_yeararray(tempNPP_yeararray<0)=NaN;
    temp_daychunk=day_chunk(findyear);
    
    for qix=1:1:length(findyear)
        temparray=tempNPP_yeararray(:,:,qix); % one by one, go through 8 day slices for each aus year
        temp2=temparray*temp_daychunk(qix);
        temp2(isnan(temp2))=0;
        temp_anntot=temp_anntot+temp2;
    end
    AnnualNPPRate_mgperyear(:,:,yix-2002)=temp_anntot;
end

AnnualNPPRate_mgperyear(AnnualNPPRate_mgperyear==0)=NaN;
AnnualNPPRate_gperyear=AnnualNPPRate_mgperyear./1000;

%% Average daily rate per year when NPP was recorded
% temp_anntot=zeros(1080,1380);
% AvDailyNPPRate_mgperday=NaN(1080,1380,length([2003:1:2020]));
% for yix = 2003:2020
%     findyear=find(timedec8day_mid>yix-0.5 & timedec8day_mid<yix+0.5);
%     tempNPP_yeararray=cafe_npp_all_8day(:,:,findyear); %austral year time slice of NPP av daily rates
%     tempNPP_yeararray(tempNPP_yeararray<0)=NaN;
%
%     AvDailyNPPRate_mgperday(:,:,yix-2002)=mean(tempNPP_yeararray,3,'omitnan');
% end

%%
load('openshelfisobath_clean21.mat')
load('box_lat_lons.mat', 'andrex_box')


%% 8-DAY DATA
%% Check data region
if setup.checkregions
    figure;
    pcolor(lon_wg,lat_wg,cafe_npp_all_8day(:,:,1)); shading flat; caxis([0 400]); colorbar;
    figure;
    pcolor(lon_wg,lat_wg,cafe_npp_all_8day(:,:,10)); shading flat; caxis([0 400]); colorbar;
    figure;
    pcolor(lon_wg,lat_wg,cafe_npp_all_8day(:,:,30)); shading flat; caxis([0 600]); colorbar;
    hold on
    SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2)%'#80471C'
    O_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2)%,'LineStyle','--')
end
clearvars SR_line O_line
%% Import regions: shelf and open ocean then andrex box
BoxIn=andrex_box;
%BoxIn=[-55,-55,-30,-30;-45,-38,-38,-45]';
IN_and=inpolygon(lon_wg,lat_wg,BoxIn(:,1),BoxIn(:,2));
findweddell=find(IN_and==1);
IN_shelf=inpolygon(lon_wg,lat_wg,shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2));
findshelf=find(IN_shelf==1);
IN_open=inpolygon(lon_wg,lat_wg,open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2));
findopen=find(IN_open==1);

region_sublist={'Weddell','Shelf','Open'};
regionfindlist= {'findweddell','findshelf','findopen'};

if setup.checkregions
    % % To check regions are within the same place
    for rix = 1:length(region_sublist)
        eval(['box_check.',region_sublist{rix},'=lat_wg;']);
        eval(['box_check.',region_sublist{rix},'(:)=0;']);
        eval(['box_check.',region_sublist{rix},'(',regionfindlist{rix},')=1;']);
        eval(['box_logic_check.',region_sublist{rix},'=logical(box_check.',region_sublist{rix},');']);
    end
    figure;
    pcolor(lon_wg,lat_wg,box_check.Weddell);
    shading flat
    figure;
    pcolor(lon_wg,lat_wg,box_check.Shelf);
    shading flat
    figure;
    pcolor(lon_wg,lat_wg,box_check.Open);
    shading flat
    hold on
    SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',5)%'#80471C'
    O_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',3.5)%,'LineStyle','--')
    
    % they look correct now! Phew!!
end
% close all

%%  loop for box,box logic, area calcs
% (varible area of open ice-free water calculated)
for aix = 1:length(algorithm)
    for tix=1:length(timedec8day)
        % make mask for ocean/area where NPP is recorded
        tempNPP=cafe_npp_all_8day(:,:,tix); %using mg m^-2 array for this
        findNPP=find(tempNPP>fill); %>=0
        NPPmask=zeros(1080,1380);
        NPPmask(findNPP)=1;
        for rix = 1:length(region_sublist)
            %box,box logic
            temp.(region_sublist{rix}).box = lat_wg;
            temp.(region_sublist{rix}).box(:) = 0;
            eval(['temp.',region_sublist{rix},'.box(',regionfindlist{rix},')=1;']);
            temp.(region_sublist{rix}).box_logic=logical(temp.(region_sublist{rix}).box);
        end
    end
end


%% Selecting values within region boxes only
for aix = 1:length(algorithm)
    for rix = 1:length(region_sublist)
        for yix=2003:2020
            temp_icefreedays=IceFreeDays_peryear(:,:,yix-2002);
            IceFree_pixels_years.(region_sublist{rix})(:,yix-2002)=temp_icefreedays(temp.(region_sublist{rix}).box_logic);
            
            temp_AnnualNPPRate=AnnualNPPRate_gperyear(:,:,yix-2002);
            AnnualNPPRate_pixels_years.(region_sublist{rix})(:,yix-2002)=temp_AnnualNPPRate(temp.(region_sublist{rix}).box_logic);
        end
    end
end


%%

% maps version

for rix = 1%:length(region_sublist)
    for yix=2003:2020
        temp_icefreedays=IceFreeDays_peryear(:,:,yix-2002);
        IceFree_pixels_yearsMAP(:,:,yix-2002)=temp_icefreedays.*(temp.(region_sublist{rix}).box_logic);
        
        temp_AnnualNPPRate=AnnualNPPRate_gperyear(:,:,yix-2002);
        AnnualNPPRate_pixels_yearsMAP(:,:,yix-2002)=temp_AnnualNPPRate.*(temp.(region_sublist{rix}).box_logic);
    end
end


figure; tiledlayout('flow'); nexttile; pcolor(IceFree_pixels_yearsMAP(:,:,1)); shading flat; nexttile; pcolor(AnnualNPPRate_pixels_yearsMAP(:,:,1)); shading flat


figure;
tiledlayout('flow'); nexttile;
pcolor(lon_wg,lat_wg,IceFree_pixels_yearsMAP(:,:,4)); shading flat; %caxis([0 600]); colorbar;
hold on
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2)%'#80471C'
O_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2)%,'LineStyle','--')
nexttile;
pcolor(lon_wg,lat_wg,AnnualNPPRate_pixels_yearsMAP(:,:,4)); shading flat; %caxis([0 600]); colorbar;
hold on
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2)%'#80471C'
O_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2)%,'LineStyle','--')
