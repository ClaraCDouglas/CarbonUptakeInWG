%% Processing data brought in by ImportData_concise.m
% clearvars

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

%algorithm={'cafe','cbpm','eppley','vgpm'};
algorithm={'cafe'};
setup.checkregions=true;
setup.plotfigures=false;
setup.startyear=2002;
setup.endyear=2020;
setup.monthly=false;
setup.eightday=true;
%% load data

% for aix = 1:length(algorithm)
% %     load([algorithm{aix} '_imported.mat'], 'area*','time*',[algorithm{aix} '*']);
% %     if aix==4
% %         load([algorithm{aix} '_imported.mat'], 'VGPM*');
% %     end
% 
% %     load([algorithm{aix} '_imported.mat']);
% 
%     load([algorithm{aix} '_8day_imported.mat']);
% end
% % vgpm_npp_tot_gC_all=VGPM_npp_tot_gC_all;
% % vgpm_npp_tot_gC_nans=VGPM_npp_tot_gC_nans;
area_MODIScafe_m2=area_MODISvgpm_m2;
clearvars VGPM_npp_tot_gC_nans VGPM_npp_tot_gC_all D0 algo_choice b filebase filedir fill *vgpm*
if setup.monthly
    timedec=time_start_all(:,1)+(time_start_all(:,2)/12)-1/24;
elseif setup.eightday
    datenum8day=datenum(time_start_all);
end
load('openshelfisobath_clean21.mat')
load('box_lat_lons.mat', 'andrex_box')


%% Import regions: shelf and open ocean then andrex box
% load('latlon_m.mat') 

IN_and=inpolygon(lon_m,lat_m,andrex_box(:,1),andrex_box(:,2));
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
    % they look correct now! Phew!!
end

%% MONTHLY DATA
if setup.monthly
    %%  loop for box,box logic, area calcs
    % (varible area of open ice-free water calculated)
    for aix = 1:length(algorithm)
        for tix=1:length(timedec)
            % make mask for ocean/area where NPP is recorded
            %         eval(['temp.NPP=',algorithm{aix},'_npp_all(:,:,tix);']); %using mg m^-2 array for this
            %         temp.findNPP=find(temp.NPP>=0);
            %         temp.NPPmask=zeros(1080,2160);
            %         temp.NPPmask(temp.findNPP)=1;
            tempNPP=vgpm_npp_all(:,:,tix); %using mg m^-2 array for this
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
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_m2(tix,1)=nansum(nansum(NPPmask.*area_MODISVGPM_m2.*(temp.(region_sublist{rix}).box)));
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_km2=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_m2/1e6;
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
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC(tix,1)=nansum(nansum(temp.NPP(temp.(region_sublist{rix}).box_logic)));
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC/1e12;

                % Mean daily rates for each month
                temp.NPP_daily=ones(1080,2160);
                eval(['temp.NPP_daily(:,:)=',algorithm{aix},'_npp_all(:,:,tix);']);
                temp.findneg=find(temp.NPP_daily<0);
                temp.NPP_daily(temp.findneg)=NaN;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans(tix,1)=nanmean(nanmean(temp.NPP_daily(temp.(region_sublist{rix}).box_logic)));
                % with zeros for the NaN months
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2(find(isnan(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2)))=0;

                % Monthly NPP rates from daily rates (daily*number of days in month)
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_mgm2_month(tix,1)=(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2(tix,1)).*time_end_all(tix,3);

                % Monthly NPP rates calculated by total NPP per month/(max)area of open water per month
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).monthly_NPP_gm2(tix,1)=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC(tix,1)./OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_m2(tix,1);
            end
        end
    end

    % Monthly means and anomalies
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            for mix = 1:12
                % mean/anom NPP
                temp.tot=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC;
                temp.totTG=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC;
                findmonth=find(time_start_all(:,2)==mix);

                temp.mean=mean(temp.tot(findmonth));
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_gC(mix,1)=temp.mean;

                temp.meanTG=mean(temp.totTG(findmonth));
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_TgC(mix,1)=temp.meanTG;

                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).anommonth_NPP_tot_gC(findmonth,1)=temp.tot(findmonth)-OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_gC(mix);
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).anommonth_NPP_tot_TgC(findmonth,1)=temp.totTG(findmonth)-OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_TgC(mix);

                % mean/anom for mean daily rates per month
                temp.daily=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2;
                temp.meandaily=mean(temp.daily(findmonth));

                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_avday_mgm2(mix,1)=temp.meandaily;

                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).anommonth_NPP_avday_mgm2(findmonth,1)=temp.daily(findmonth)-OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_avday_mgm2(mix);


                % mean/anom ice free
                temp.totkm2=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_km2;
                temp.meankm2=mean(temp.totkm2(findmonth));
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).meanmonth_icefree(mix,1)=temp.meankm2;

                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).anommonth_icefree(findmonth,1)=temp.totkm2(findmonth)-OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).meanmonth_icefree(mix);
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
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,1)=yix;

                % and annual rates of NPP
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(yix-2001,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(yix-2001,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(yix-2001,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(yix-2001,1)=yix;

                if size(findyear)<12
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,2)=NaN;

                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(yix-2001,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(yix-2001,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(yix-2001,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(yix-2001,2)=NaN;
                else
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)=sum(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC(findyear,1));
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,2)=sum(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_km2(findyear));
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,2)=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2001,2)/12;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,2)=max(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_m2(findyear,1),[],'omitnan');

                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(yix-2001,2)=nansum(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_mgm2_month(findyear));
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(yix-2001,2)=nansum(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).monthly_NPP_mgm2(findyear,1));
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(yix-2001,2)=(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)*1e12)/(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,2)*1e6);
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(yix-2001,2)=(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)*1e12)/OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,2);
                end
            end
        end
    end

    % mean annual total NPP and ice-free water area
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_MEANtot_TgC_annual=nanmean(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(:,2));
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_MEANMEAN=nanmean(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(:,2));
        end
    end

    for aix = 1:length(algorithm)
        OceanProd_8day.(algorithm{aix}).(region_sublist{2}).IceFree_proportion(:,1)=(OceanProd_8day.(algorithm{aix}).(region_sublist{2}).IceFree_annualMEAN(:,2)./OceanProd_8day.(algorithm{aix}).(region_sublist{1}).IceFree_annualMEAN(:,2)).*100;
        OceanProd_8day.(algorithm{aix}).(region_sublist{2}).IceFree_AVproportion=nanmean(OceanProd_8day.(algorithm{aix}).(region_sublist{2}).IceFree_proportion(:,1));
        OceanProd_8day.(algorithm{aix}).(region_sublist{3}).IceFree_proportion(:,1)=(OceanProd_8day.(algorithm{aix}).(region_sublist{3}).IceFree_annualMEAN(:,2)./OceanProd_8day.(algorithm{aix}).(region_sublist{1}).IceFree_annualMEAN(:,2)).*100;
        OceanProd_8day.(algorithm{aix}).(region_sublist{3}).IceFree_AVproportion=nanmean(OceanProd_8day.(algorithm{aix}).(region_sublist{3}).IceFree_proportion(:,1));
    end

    clearvars regionfindlist findweddell findshelf findopen temp setup

    save('ProcessedData.mat','OceanProd_8day','algorithm','region_sublist','timedec');

    clearvars -except OceanProd algorithm region_sublist timedec vgpm_npp_all time* lat_m lon_m andrex_box
    %% calculate NPP per m2 per day per year...
    % either: average daily rate per month (NPP_av_mgm2) for each year
    % or NPP_rates_fromdaily/# days in year
    % findyear, if findyear >12, then sum column 3 of time_end_all where column 1 = yix

end

%% 8-DAY DATA
if setup.eightday
    %%  loop for box,box logic, area calcs
    % (varible area of open ice-free water calculated)
    for aix = 1:length(algorithm)
        for tix=1:length(datenum8day)
            % make mask for ocean/area where NPP is recorded
            tempNPP=cafe_npp_all_8day(:,:,tix); %using mg m^-2 array for this
            %findNPP=find(tempNPP>=0);
            findNPP2=find(~(tempNPP<0)); % would be even better coding if put ==FillValue
            NPPmask=zeros(540,2160);
            NPPmask(findNPP2)=1;
            for rix = 1:length(region_sublist)
                %box,box logic
                temp.(region_sublist{rix}).box = lat_m;
                temp.(region_sublist{rix}).box(:) = 0;
                eval(['temp.',region_sublist{rix},'.box(',regionfindlist{rix},')=1;']);

                temp.(region_sublist{rix}).box_logic=logical(temp.(region_sublist{rix}).box);

                % Area of Open Ice-Free water (per month)
                %             eval(['temp.area_MODIS=area_MODIS',algorithm{aix},'_m2;']);
                %             OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_m2(tix,1)=nansum(nansum((temp.NPPmask).*(temp.area_MODIS).*(temp.(region_sublist{rix}).box)));
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_m2(tix,1)=sum(sum(NPPmask.*area_MODIScafe_m2.*(temp.(region_sublist{rix}).box),'omitnan'),'omitnan');
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_km2=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_m2/1e6;
            end
        end
    end
    %% area of region and sub-region

    for rix = 1:length(region_sublist)
        OceanProd_8day.total_area.(region_sublist{rix})=sum(sum(area_MODIScafe_km2.*(temp.(region_sublist{rix}).box),'omitnan'),'omitnan');
    end
    OceanProd_8day.total_area.SROOtotal=OceanProd_8day.total_area.Shelf+OceanProd_8day.total_area.Open;
    OceanProd_8day.total_area.SROOtot_WGdiff=OceanProd_8day.total_area.SROOtotal-OceanProd_8day.total_area.Weddell;
    clear IN_and IN_shelf IN_open tempNPP findNPP* NPPmask regionfindlist findweddell findshelf findopen

    %% Calculating variables

    % NPP 8-day timeseries
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            for tix=1:length(datenum8day)
                % Total NPP per 8-day period
                temp.NPP_tot=ones(540,2160);
                eval(['temp.NPP(:,:)=',algorithm{aix},'_npp_tot_gC_nans_8day(:,:,tix);']);
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC(tix,1)=sum(sum(temp.NPP(temp.(region_sublist{rix}).box_logic),'omitnan'),'omitnan');
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC/1e12;

                % Mean daily rates for each 8-day period
                temp.NPP_daily=ones(540,2160);
                eval(['temp.NPP_daily(:,:)=',algorithm{aix},'_npp_all_8day(:,:,tix);']);
                temp.findneg=find(temp.NPP_daily<0);
                temp.NPP_daily(temp.findneg)=NaN;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans(tix,1)=mean(mean(temp.NPP_daily(temp.(region_sublist{rix}).box_logic),'omitnan'),'omitnan');
                % with zeros for the NaN months
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2(find(isnan(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2)))=0;

            end
        end
    end

    % Monthly means and anomalies - not run just now
    %  for aix = 1:length(algorithm)
    %      for rix = 1:length(region_sublist)
    %          for mix = 1:12
    %              % mean/anom NPP
    %              temp.tot=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC;
    %              temp.totTG=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC;
    %              findmonth=find(time_start_all(:,2)==mix);
    %
    %              temp.mean=mean(temp.tot(findmonth));
    %              OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_gC(mix,1)=temp.mean;
    %
    %              temp.meanTG=mean(temp.totTG(findmonth));
    %              OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_TgC(mix,1)=temp.meanTG;
    %
    %              OceanProd.(algorithm{aix}).(region_sublist{rix}).anommonth_NPP_tot_gC(findmonth,1)=temp.tot(findmonth)-OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_gC(mix);
    %              OceanProd.(algorithm{aix}).(region_sublist{rix}).anommonth_NPP_tot_TgC(findmonth,1)=temp.totTG(findmonth)-OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_tot_TgC(mix);
    %
    %              % mean/anom for mean daily rates per month
    %              temp.daily=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2;
    %              temp.meandaily=mean(temp.daily(findmonth));
    %
    %              OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_avday_mgm2(mix,1)=temp.meandaily;
    %
    %              OceanProd.(algorithm{aix}).(region_sublist{rix}).anommonth_NPP_avday_mgm2(findmonth,1)=temp.daily(findmonth)-OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_NPP_avday_mgm2(mix);
    %
    %
    %              % mean/anom ice free
    %              temp.totkm2=OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_km2;
    %              temp.meankm2=mean(temp.totkm2(findmonth));
    %              OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_icefree(mix,1)=temp.meankm2;
    %
    %              OceanProd.(algorithm{aix}).(region_sublist{rix}).anommonth_icefree(findmonth,1)=temp.totkm2(findmonth)-OceanProd.(algorithm{aix}).(region_sublist{rix}).meanmonth_icefree(mix);
    %          end
    %      end
    %  end
    %  clearvars findmonth

    % Annual total NPP, annual open water variables, annual rates of NPP
    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            for yix = 2003:2020
                findyear1=find(time_start_all(:,1)==yix-1 & time_start_all(:,2)>=7);
                findyear2=find(time_start_all(:,1)==yix & time_start_all(:,2)<=6);
                findyear=cat(1,findyear1,findyear2);
                % Integrated NPP over austral year (June to June)
                % and total and mean annual open ice-free water area and max monthly
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2002,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2002,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2002,1)=yix;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2002,1)=yix;

                % and annual rates of NPP
%                 OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(yix-2002,1)=yix;
%                 OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(yix-2002,1)=yix;
%                 OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(yix-2002,1)=yix;
%                 OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(yix-2002,1)=yix;

                if size(findyear)<46
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2002,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2002,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2002,2)=NaN;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2002,2)=NaN;

%                     OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(yix-2002,2)=NaN;
%                     OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(yix-2002,2)=NaN;
%                     OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(yix-2002,2)=NaN;
%                     OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(yix-2002,2)=NaN;
                else
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2002,2)=sum(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC(findyear,1));
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2002,2)=sum(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_km2(findyear));
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2002,2)=OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualTOT(yix-2002,2)/12;
                    OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2002,2)=max(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).area_month_m2(findyear,1),[],'omitnan');

                    %                  OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(yix-2002,2)=nansum(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_mgm2_month(findyear));
                    %                  OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(yix-2002,2)=nansum(OceanProd.(algorithm{aix}).(region_sublist{rix}).monthly_NPP_mgm2(findyear,1));
                    %                  OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(yix-2002,2)=(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2002,2)*1e12)/(OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2002,2)*1e6);
                    %                  OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(yix-2002,2)=(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2002,2)*1e12)/OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2002,2);
                end
            end
        end
    end

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

    clearvars regionfindlist findweddell findshelf findopen temp setup

%     save('ProcessedData_8day.mat','OceanProd','algorithm','region_sublist','timedec');

%     clearvars -except OceanProd algorithm region_sublist timedec cafe_npp_all time* lat_m lon_m andrex_box
    %% calculate NPP per m2 per day per year...
    % either: average daily rate per month (NPP_av_mgm2) for each year
    % or NPP_rates_fromdaily/# days in year
    % findyear, if findyear >12, then sum column 3 of time_end_all where column 1 = yix

end

%% check data
load('ProcessedData.mat', 'OceanProd')

figure;
plot(OceanProd.cafe.Open.NPP_tot_gC)

figure;
plot(OceanProd_8day.cafe.Open.NPP_tot_gC)
