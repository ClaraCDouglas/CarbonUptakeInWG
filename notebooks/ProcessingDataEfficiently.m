%% Processing data brought in by ImportData_concise.m
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
load('openshelfisobath_clean21.mat')
load('box_lat_lons.mat', 'andrex_box')
desktop = 1;
laptop=0;

%% MONTHLY DATA
if setup.monthly
    %% Check data region
if setup.checkregions
    figure;
    pcolor(lon_m,lat_m,cafe_npp_all(:,:,1)); shading flat; caxis([0 400]); colorbar;
    figure;
    pcolor(lon_m,lat_m,cafe_npp_all(:,:,10)); shading flat; caxis([0 400]); colorbar;
    figure;
    pcolor(lon_m,lat_m,cafe_npp_all(:,:,30)); shading flat; caxis([0 600]); colorbar;
    hold on
    SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2)%'#80471C'
    O_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2)%,'LineStyle','--')
end
clearvars SR_line O_line
%% Import regions: shelf and open ocean then andrex box
% load('latlon_m.mat')
%BoxIn=andrex_box;
BoxIn=[-55,-55,-30,-30;-45,-38,-38,-45]';
IN_and=inpolygon(lon_m,lat_m,BoxIn(:,1),BoxIn(:,2));
findweddell=find(IN_and==1);
IN_shelf=inpolygon(lon_m,lat_m,shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2));
findshelf=find(IN_shelf==1);
IN_open=inpolygon(lon_m,lat_m,open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2));
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
    hold on
    SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',5)%'#80471C'
    O_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',3.5)%,'LineStyle','--')

    % they look correct now! Phew!!
end

    %%  loop for box,box logic, area calcs
    % (varible area of open ice-free water calculated)
    for aix = 1:length(algorithm)
        for tix=1:length(timedec)
            % make mask for ocean/area where NPP is recorded
            %         eval(['temp.NPP=',algorithm{aix},'_npp_all(:,:,tix);']); %using mg m^-2 array for this
            %         temp.findNPP=find(temp.NPP>=0);
            %         temp.NPPmask=zeros(1080,2160);
            %         temp.NPPmask(temp.findNPP)=1;
            tempNPP=cafe_npp_all(:,:,tix); %using mg m^-2 array for this
            findNPP=find(tempNPP>=0);
            NPPmask=zeros(1080,2160);
            NPPmask(findNPP)=1;
            for rix = 1:length(region_sublist)
                %box,box logic
                temp.(region_sublist{rix}).box = lat_m;
                temp.(region_sublist{rix}).box(:) = 0;
                eval(['temp.',region_sublist{rix},'.box(',regionfindlist{rix},')=1;']);

                temp.(region_sublist{rix}).box_logic=logical(temp.(region_sublist{rix}).box);

                % Area of Open Ice-Free water (per month)
                %             eval(['temp.area_MODIS=area_MODIS',algorithm{aix},'_m2;']);
                %             OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_m2(tix,1)=nansum(nansum((temp.NPPmask).*(temp.area_MODIS).*(temp.(region_sublist{rix}).box)));
                OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_m2(tix,1)=sum(sum(NPPmask.*area_MODIScafe_m2.*(temp.(region_sublist{rix}).box),'omitnan'),'omitnan');
                OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_km2=OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_m2/1e6;
            end
        end
    end
    clear IN_and IN_shelf IN_open tempNPP findNPP NPPmask regionfindlist findweddell findshelf findopen

    %% Calculating variables

    % NPP monthly timeseries
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            for tix=1:length(timedec)
                % Total NPP per month
                temp.NPP_tot=ones(1080,2160);
                eval(['temp.NPP(:,:)=',algorithm{aix},'_npp_tot_gC_nans(:,:,tix);']);
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC(tix,1)=sum(sum(temp.NPP(temp.(region_sublist{rix}).box_logic),'omitnan'),'omitnan');
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC/1e12;

                % Mean daily rates for each month
                temp.NPP_daily=ones(1080,2160);
                eval(['temp.NPP_daily(:,:)=',algorithm{aix},'_npp_all(:,:,tix);']);
                temp.findneg=find(temp.NPP_daily<0);
                temp.NPP_daily(temp.findneg)=NaN;
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans(tix,1)=mean(mean(temp.NPP_daily(temp.(region_sublist{rix}).box_logic),'omitnan'),'omitnan');
                % replacing NaNs with zeros for the NaN months
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans;
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2(find(isnan(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2)))=0;

                % mean daily rates for region for each month where ice covered is 0
                %temp.findneg=find(temp.NPP_daily<0);
                temp.NPP_daily(temp.findneg)=0;
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_zeros(tix,1)=mean(mean(temp.NPP_daily(temp.(region_sublist{rix}).box_logic),'omitnan'),'omitnan');

                % Monthly NPP rates from daily rates (daily*number of days in month)
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_mgm2_month(tix,1)=(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2(tix,1)).*time_end_all(tix,3);

                % Monthly NPP rates calculated by total NPP per month/(max)area of open water per month
                OceanProd.(algorithm{aix}).(region_sublist{rix}).monthly_NPP_gm2(tix,1)=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC(tix,1)./OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_m2(tix,1);
            end
        end
    end

    % Monthly means and anomalies
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            for mix = 1:12
                % mean/anom NPP
                temp.tot=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC;
                temp.totTG=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC;
                findmonth=find(time_start_all(:,2)==mix);

                temp.mean=mean(temp.tot(findmonth));
                OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_gC(mix,1)=temp.mean;

                temp.meanTG=mean(temp.totTG(findmonth));
                OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_TgC(mix,1)=temp.meanTG;

                OceanProd.(algorithm{aix}).(region_sublist{rix}).anommonth_NPP_tot_gC(findmonth,1)=temp.tot(findmonth)-OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_gC(mix);
                OceanProd.(algorithm{aix}).(region_sublist{rix}).anommonth_NPP_tot_TgC(findmonth,1)=temp.totTG(findmonth)-OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_TgC(mix);

                % mean/anom for mean daily rates per month
                temp.daily=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2;
                temp.meandaily=mean(temp.daily(findmonth));

                OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_avday_mgm2(mix,1)=temp.meandaily;

                OceanProd.(algorithm{aix}).(region_sublist{rix}).anommonth_NPP_avday_mgm2(findmonth,1)=temp.daily(findmonth)-OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_avday_mgm2(mix);


                % mean/anom ice free
                temp.totkm2=OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_km2;
                temp.meankm2=mean(temp.totkm2(findmonth));
                OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_icefree(mix,1)=temp.meankm2;

                OceanProd.(algorithm{aix}).(region_sublist{rix}).anommonth_icefree(findmonth,1)=temp.totkm2(findmonth)-OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_icefree(mix);
            end
        end
    end
    clearvars findmonth

    % Annual total NPP, annual open water variables, annual rates of NPP
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            for yix = 2002:2020
                findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
                % Integrated NPP over austral year (June to June)
                % and total and mean annual open ice-free water area and max monthly
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,1)=yix;
                OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,1)=yix;
                OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,1)=yix;
                OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,1)=yix;

                % and annual rates of NPP
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(yix-2001,1)=yix;
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(yix-2001,1)=yix;
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(yix-2001,1)=yix;
                OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(yix-2001,1)=yix;

                if size(findyear)<12
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)=NaN;
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,2)=NaN;
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,2)=NaN;
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,2)=NaN;

                    OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(yix-2001,2)=NaN;
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(yix-2001,2)=NaN;
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(yix-2001,2)=NaN;
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(yix-2001,2)=NaN;
                else
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)=sum(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC(findyear,1));
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,2)=sum(OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_km2(findyear));
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,2)=OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,2)/12;
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,2)=max(OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_m2(findyear,1),[],'omitnan');

                    OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(yix-2001,2)=sum(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_mgm2_month(findyear),'omitnan');
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(yix-2001,2)=sum(OceanProd.(algorithm{aix}).(region_sublist{rix}).monthly_NPP_gm2(findyear,1),'omitnan');
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(yix-2001,2)=(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)*1e12)/(OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,2)*1e6);
                    OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(yix-2001,2)=(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)*1e12)/OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,2);
                end
            end
        end
    end

    % mean annual total NPP and ice-free water area
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_MEANtot_TgC_annual=mean(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(:,2),'omitnan');
            OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_MEANMEAN=mean(OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(:,2),'omitnan');
        end
    end

    for aix = 1:length(algorithm)
        OceanProd.(algorithm{aix}).(region_sublist{2}).IceFree_proportion(:,1)=(OceanProd.(algorithm{aix}).(region_sublist{2}).IceFree_annualMEAN(:,2)./OceanProd.(algorithm{aix}).(region_sublist{1}).IceFree_annualMEAN(:,2)).*100;
        OceanProd.(algorithm{aix}).(region_sublist{2}).IceFree_AVproportion=nanmean(OceanProd.(algorithm{aix}).(region_sublist{2}).IceFree_proportion(:,1));
        OceanProd.(algorithm{aix}).(region_sublist{3}).IceFree_proportion(:,1)=(OceanProd.(algorithm{aix}).(region_sublist{3}).IceFree_annualMEAN(:,2)./OceanProd.(algorithm{aix}).(region_sublist{1}).IceFree_annualMEAN(:,2)).*100;
        OceanProd.(algorithm{aix}).(region_sublist{3}).IceFree_AVproportion=nanmean(OceanProd.(algorithm{aix}).(region_sublist{3}).IceFree_proportion(:,1));
    end

    clearvars regionfindlist findweddell findshelf findopen temp setup

%     save('ProcessedDataMonthly_Jan22.mat','OceanProd','algorithm','region_sublist','timedec');

    clearvars -except OceanProd algorithm region_sublist timedec vgpm_npp_all time* lat_m lon_m andrex_box
    %% calculate NPP per m2 per day per year...
    % either: average daily rate per month (NPP_av_mgm2) for each year
    % or NPP_rates_fromdaily/# days in year
    % findyear, if findyear >12, then sum column 3 of time_end_all where column 1 = yix

end

%% 8-DAY DATA
if setup.eightday
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
%BoxIn=andrex_box;
BoxIn=[-55,-55,-30,-30;-45,-38,-38,-45]';
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

                % Area of Open Ice-Free water (per month)
                %             eval(['temp.area_MODIS=area_MODIS',algorithm{aix},'_m2;']);
                %             OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_m2(tix,1)=nansum(nansum((temp.NPPmask).*(temp.area_MODIS).*(temp.(region_sublist{rix}).box)));

                % where there is ocean * area in m2 * region box. then summed = total area that there is ocean for each 8 day period
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_8day_m2(tix,1)=sum(sum(NPPmask.*area_MODIScafe_m2_wg.*(temp.(region_sublist{rix}).box),'omitnan'),'omitnan');
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_8day_km2=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_8day_m2/1e6;
            end
        end
    end
    %% area of region and sub-region

    for rix = 1:length(region_sublist)
        OceanProd_8day.total_area.(region_sublist{rix})=sum(sum(area_MODIScafe_km2_wg.*(temp.(region_sublist{rix}).box),'omitnan'),'omitnan');
    end
    OceanProd_8day.total_area.SROOtotal=OceanProd_8day.total_area.Shelf+OceanProd_8day.total_area.Open;
    OceanProd_8day.total_area.SROOtot_WGdiff=OceanProd_8day.total_area.SROOtotal-OceanProd_8day.total_area.Weddell;
    %     clear IN_and IN_shelf IN_open tempNPP findNPP* NPPmask regionfindlist findweddell findshelf findopen

    %% Calculating variables

    % NPP 8-day timeseries
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            for tix=1:length(timedec8day)
                % Total NPP per 8-day period
                temp.NPP_tot=NaN(1080,1380);
                eval(['temp.NPP(:,:)=',algorithm{aix},'_npp_tot_gC_nans_8day(:,:,tix);']);
                %                 temp.NPP=cafe_npp_all_8day(:,:,tix); %using mg m^-2 array for this

                % temp.NPP is total NPP in 8 day period (land/ice/clouds=Nan)
                %summed for pixels within region
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC(tix,1)=sum(sum(temp.NPP(temp.(region_sublist{rix}).box_logic),'omitnan'),'omitnan');
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC/1e12;

                % Mean daily rates for each 8-day period
                temp.NPP_daily=ones(1080,1380);
                eval(['temp.NPP_daily(:,:)=',algorithm{aix},'_npp_all_8day(:,:,tix);']);
                temp.findneg=find(temp.NPP_daily<0);
                temp.NPP_daily(temp.findneg)=NaN;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans(tix,1)=mean(mean(temp.NPP_daily(temp.(region_sublist{rix}).box_logic),'omitnan'),'omitnan');
                % replace NaN estimates for 8 day period with zeros
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2(isnan(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2))=0;

                % mean daily rates for region for each 8-day period where ice covered is 0
                %temp.findneg=find(temp.NPP_daily<0);
                temp.NPP_daily(temp.findneg)=0;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_zeros(tix,1)=mean(mean(temp.NPP_daily(temp.(region_sublist{rix}).box_logic)));
            end
        end
    end
    clearvars cafe_npp_tot_gC_nans_8day

    % Monthly means and anomalies
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            for yix = 2003:2020
                findyear=find(timedec8day>yix-0.5 & timedec8day<yix+0.5);
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_zeros_cols(:,yix-2002)=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_zeros(findyear);
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_cols(:,yix-2002)=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC(findyear);
            end
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_zeros_mean=mean(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_zeros_cols,2);
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_zeros_anom=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_zeros_cols-OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_zeros_mean;
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_mean=mean(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_cols,2);
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_anom=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_cols-OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_mean;
        end
    end

    if setup.plotfigures
        linecolors = jet(length(setup.yearrange0320));
        figure;
        tiledlayout('flow')
        nexttile
        for cix=1:length(setup.yearrange0320)
            plot(OceanProd_8day.cafe.Shelf.NPP_av_mgm2_zeros_anom(:,cix),'Color',linecolors(cix,:),'LineWidth',1.5);
            hold on
        end
        legend(num2str(setup.yearrange0320))
        yline(0,'--')
        title('Daily rate anomalies')
        nexttile
        for cix=1:length(setup.yearrange0320)
            plot(OceanProd_8day.cafe.Shelf.NPP_av_mgm2_zeros_cols(:,cix),'Color',linecolors(cix,:),'LineWidth',1.5);
            hold on
            plot(OceanProd_8day.cafe.Shelf.NPP_av_mgm2_zeros_mean,'k');
        end
        title('Daily rate and mean')
        nexttile
        for cix=1:length(setup.yearrange0320)
            plot(OceanProd_8day.cafe.Weddell.NPP_tot_TgC_anom(:,cix),'Color',linecolors(cix,:),'LineWidth',1.5);
            hold on
        end
        legend(num2str(setup.yearrange0320))
        yline(0,'--')
        title('8-day total anomalies')
        nexttile
        for cix=1:length(setup.yearrange0320)
            plot(OceanProd_8day.cafe.Weddell.NPP_tot_TgC_cols(:,cix),'Color',linecolors(cix,:),'LineWidth',1.5);
            hold on
            plot(OceanProd_8day.cafe.Weddell.NPP_tot_TgC_mean,'k');
        end
        title('8-day total NPP and mean')
    end

    % Annual total NPP, annual open water variables, annual rates of NPP
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            for yix = 2002:2020
                findyear=find(timedec8day>yix-0.5 & timedec8day<yix+0.5);
                %findyear1=find(time_start_all(:,1)==yix-1 & time_start_all(:,2)>=7);
                %findyear2=find(time_start_all(:,1)==yix & time_start_all(:,2)<=6);
                %findyear=cat(1,findyear1,findyear2);

                % Integrated NPP over austral year (July to June)
                % and total and mean annual open ice-free water area and max monthly
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,1)=yix;

                if size(findyear)<46
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,2)=NaN;
                else
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)=sum(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC(findyear,1));
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,2)=sum(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_8day_km2(findyear));
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,2)=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,2)/46;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,2)=max(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_8day_m2(findyear,1),[],'omitnan');
                end
            end
        end
    end

    % Annual rates of NPP (mg C m^-2 yr^-1)
    % Via Integrated daily
    % (attempt 3 - see AnnualNPPratesCalc.m for previous workings)
    leapyear1=time_start_all(:,1);
    leapyear2=rem(leapyear1,4)==0;
    leapyear3=zeros(length(leapyear1),1);
    leapyear3(leapyear2==1)=366;
    leapyear3(leapyear3==0)=365;
    days_test(:,1)=(timedec8day_end-timedec8day);
    days_test(:,2)=days_test(:,1)*365;
    days_test(:,3)=days_test(:,1).*leapyear3;
    day_chunk=days_test(:,1).*leapyear3;

    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).day_chunk_rates_gm2=(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans.*day_chunk)/1000;
            for yix = 2003:2020
                findyear=find(timedec8day>yix-0.5 & timedec8day<yix+0.5);
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).annual_rate_chunk(yix-2002,1)= ...
                    sum(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).day_chunk_rates_gm2(findyear),'omitnan');
                day_chunk_IF=day_chunk(findyear);
                day_chunk_IF(isnan(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).day_chunk_rates_gm2(findyear)))=NaN;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).days_ice_free(yix-2002,1)=sum(day_chunk_IF,'omitnan')
            end
        end
    end
    clearvars leap* days_test day_chunk_IF

    % mean annual total NPP and ice-free water area
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_MEANtot_TgC_annual=mean(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(:,2),'omitnan');
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_MEANMEAN=mean(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(:,2),'omitnan');
        end
    end

    % proportion of ice-free area seen in OO & SR
    for aix = 1:length(algorithm)
        OceanProd_8day.(algorithm{aix}).(region_sublist{2}).IceFree_proportion(:,1)=(OceanProd_8day.(algorithm{aix}).(region_sublist{2}).IceFree_annualMEAN(:,2)./OceanProd_8day.(algorithm{aix}).(region_sublist{1}).IceFree_annualMEAN(:,2)).*100;
        OceanProd_8day.(algorithm{aix}).(region_sublist{2}).IceFree_AVproportion=nanmean(OceanProd_8day.(algorithm{aix}).(region_sublist{2}).IceFree_proportion(:,1));
        OceanProd_8day.(algorithm{aix}).(region_sublist{3}).IceFree_proportion(:,1)=(OceanProd_8day.(algorithm{aix}).(region_sublist{3}).IceFree_annualMEAN(:,2)./OceanProd_8day.(algorithm{aix}).(region_sublist{1}).IceFree_annualMEAN(:,2)).*100;
        OceanProd_8day.(algorithm{aix}).(region_sublist{3}).IceFree_AVproportion=nanmean(OceanProd_8day.(algorithm{aix}).(region_sublist{3}).IceFree_proportion(:,1));
    end

    %     clearvars regionfindlist findweddell findshelf findopen temp setup
         save('ProcessedData_8day_Jan22.mat','OceanProd_8day','algorithm','region_sublist','timedec8day');
         clearvars -except OceanProd* algorithm region_sublist time* cafe_npp_all time* lat_* lon_* andrex_box
    %% calculate NPP per m2 per day per year...
    % either: average daily rate per month (NPP_av_mgm2) for each year
    % or NPP_rates_fromdaily/# days in year
    % findyear, if findyear >12, then sum column 3 of time_end_all where column 1 = yix

end

%% check data
load('ProcessedDataMonthly_Jan22.mat', 'OceanProd')
load('ProcessedDataMonthly_Jan22.mat', 'timedec')
time_start_all_8day=time_start_all;
time_end_all_8day=time_end_all;
load('cafe_imported.mat', 'time_start_all');
load('cafe_imported.mat', 'time_end_all');

timedecend=time_end_all(:,1)+(time_end_all(:,2)/12)-1/24;
timedec_momid=mean([timedec,timedecend],2);
timedec8day_mid=mean([timedec8day,timedec8day_end],2);

figure;
plot(timedec_momid,OceanProd.cafe.Open.NPP_tot_gC)
hold on
plot(timedec8day_mid,OceanProd_8day.cafe.Open.NPP_tot_gC)

%145seconds to run on laptop

figure;
tiledlayout(2,6)
nexttile(2,[1 2])
y=[OceanProd.cafe.Weddell.NPP_tot_TgC_annual(:,2) OceanProd_8day.cafe.Weddell.NPP_tot_TgC_annual(:,2)];
bar(OceanProd.cafe.Weddell.NPP_tot_TgC_annual(:,1),y)
txt1={'Total NPP per year (Tg C) calculated from 8-day averages vs monthly averages', 'NORTH OF ACC...'}
title(txt1)
ylabel('Integrated annual NPP (Tg C yr^-^1)')
legend('Monthly','8-day')
nexttile(4,[1 2])
y=[OceanProd.cafe.Open.NPP_tot_TgC_annual(:,2) OceanProd_8day.cafe.Open.NPP_tot_TgC_annual(:,2)];
bar(OceanProd.cafe.Open.NPP_tot_TgC_annual(:,1),y)
txt2={'Total NPP per year (Tg C) calculated from 8-day averages vs monthly averages', 'OPEN OCEAN Weddell Gyre REGION...'}
title(txt2)
ylabel('Integrated annual NPP (Tg C yr^-^1)')
legend('Monthly','8-day')

nexttile(7,[1 2])
plot(timedec_momid,OceanProd.cafe.Open.NPP_tot_gC)
hold on
plot(timedec8day_mid,OceanProd_8day.cafe.Open.NPP_tot_gC)
legend('Monthly','8-day')
title('Total NPP (Tg C) per time chunk (month or 8 day) OPEN OCEAN WG REGION...')

nexttile(9,[1 2])
plot(timedec_momid,OceanProd.cafe.Open.NPP_tot_gC)
hold on
plot(timedec8day_mid,OceanProd_8day.cafe.Open.NPP_tot_gC.*3.8)
legend('Monthly','8-day')
txt={'Total NPP (Tg C) per time chunk (month or 8 day) OPEN OCEAN WG REGION','but 8-day values *3.8 to pretend they are month totals'}
title(txt)
nexttile(11,[1 2])
plot(timedec_momid,OceanProd.cafe.Open.NPP_av_mgm2,'LineWidth',3)
hold on
plot(timedec8day_mid,OceanProd_8day.cafe.Open.NPP_av_mgm2,'LineWidth',1.75)
ylabel('NPP average daily rate mg m^-^2 day^-^1')
title('Differences in daily rate timeseries calculated from 8-day averages vs monthly averages...')
legend('Monthly','8-day')

% close all


figure;
y=[OceanProd.cafe.Weddell.IceFree_annualMAX(:,2) OceanProd_8day.cafe.Weddell.IceFree_annualMAX(:,2)];
bar(OceanProd.cafe.Weddell.NPP_tot_TgC_annual(:,1),y)
title('Maximum ice-free water area calculated from 8-day averages vs monthly averages...')
ylabel('km^2')
legend('Monthly','8-day')

%% thinking
think=(OceanProd_8day.cafe.Open.NPP_tot_TgC_annual(:,2)./OceanProd.cafe.Open.NPP_tot_TgC_annual(:,2)).*100;
think(:,2)=(OceanProd.cafe.Open.NPP_tot_TgC_annual(:,2)./OceanProd_8day.cafe.Open.NPP_tot_TgC_annual(:,2)).*100;



sum(OceanProd.cafe.Open.NPP_tot_TgC_annual(:,2),'omitnan')
sum(OceanProd_8day.cafe.Open.NPP_tot_TgC_annual(:,2),'omitnan')






