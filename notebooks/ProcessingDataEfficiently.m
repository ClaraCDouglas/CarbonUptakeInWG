%% Processing data brought in by ImportData_concise.m
clearvars
cd C:\Users\ccd1n18\Documents\Projects\Carbon-Uptake-in-WG_Manuscript\Workspace_Variables
algorithm={'cafe','cbpm','eppley','vgpm'};
setup.checkregions=false;
setup.plotfigures=false;
setup.startyear=2002;
setup.endyear=2020;

%% load data

for aix = 1:length(algorithm)
%     load([algorithm{aix} '_imported.mat'], 'area*','time*',[algorithm{aix} '*']);
%     if aix==4
%         load([algorithm{aix} '_imported.mat'], 'VGPM*');
%     end
    load([algorithm{aix} '_imported.mat']);
end
vgpm_npp_tot_gC_all=VGPM_npp_tot_gC_all;
vgpm_npp_tot_gC_nans=VGPM_npp_tot_gC_nans;
clear VGPM_npp_tot_gC_nans VGPM_npp_tot_gC_all D0 algo_choice b filebase filedir fill
timedec=time_start_all(:,1)+(time_start_all(:,2)/12)-1/24;
load('openshelf_coord.mat')
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
            OceanProd.(algorithm{aix}).(region_sublist{rix}).area_month_m2(tix,1)=nansum(nansum(NPPmask.*area_MODISVGPM_m2.*(temp.(region_sublist{rix}).box)));
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
            OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC(tix,1)=nansum(nansum(temp.NPP(temp.(region_sublist{rix}).box_logic)));
            OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_gC/1e12;
            
% Mean daily rates for each month            
            temp.NPP_daily=ones(1080,2160);
            eval(['temp.NPP_daily(:,:)=',algorithm{aix},'_npp_all(:,:,tix);']);
            temp.findneg=find(temp.NPP_daily<0);
            temp.NPP_daily(temp.findneg)=NaN;
            OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans(tix,1)=nanmean(nanmean(temp.NPP_daily(temp.(region_sublist{rix}).box_logic)));
                % with zeros for the NaN months
            OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2=OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans;
            OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2(find(isnan(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2)))=0;
            
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
                 
                 OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromdaily(yix-2001,2)=nansum(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_mgm2_month(findyear));
                 OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_frommonth(yix-2001,2)=nansum(OceanProd.(algorithm{aix}).(region_sublist{rix}).monthly_NPP_mgm2(findyear,1));
                 OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_fromannual(yix-2001,2)=(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)*1e12)/(OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(yix-2001,2)*1e6);
                 OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_rates_yearmaxOW(yix-2001,2)=(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(yix-2001,2)*1e12)/OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMAX(yix-2001,2);
             end
         end
     end
 end

 % mean annual total NPP and ice-free water area
 for aix = 1:length(algorithm)
     for rix = 1:length(region_sublist)
         OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_MEANtot_TgC_annual=nanmean(OceanProd.(algorithm{aix}).(region_sublist{rix}).NPP_tot_TgC_annual(:,2));
         OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_MEANMEAN=nanmean(OceanProd.(algorithm{aix}).(region_sublist{rix}).IceFree_annualMEAN(:,2));
     end
 end
 
 
 
 for aix = 1:length(algorithm)
     OceanProd.(algorithm{aix}).(region_sublist{2}).IceFree_proportion(:,1)=(OceanProd.(algorithm{aix}).(region_sublist{2}).IceFree_annualMEAN(:,2)./OceanProd.(algorithm{aix}).(region_sublist{1}).IceFree_annualMEAN(:,2)).*100;
     OceanProd.(algorithm{aix}).(region_sublist{2}).IceFree_AVproportion=nanmean(OceanProd.(algorithm{aix}).(region_sublist{2}).IceFree_proportion(:,1));
     OceanProd.(algorithm{aix}).(region_sublist{3}).IceFree_proportion(:,1)=(OceanProd.(algorithm{aix}).(region_sublist{3}).IceFree_annualMEAN(:,2)./OceanProd.(algorithm{aix}).(region_sublist{1}).IceFree_annualMEAN(:,2)).*100;
     OceanProd.(algorithm{aix}).(region_sublist{3}).IceFree_AVproportion=nanmean(OceanProd.(algorithm{aix}).(region_sublist{3}).IceFree_proportion(:,1));
 end
 
 clearvars regionfindlist findweddell findshelf findopen temp setup

save('ProcessedData.mat','OceanProd','algorithm','region_sublist','timedec');
 
clearvars -except OceanProd algorithm region_sublist timedec vgpm_npp_all time* lat_m lon_m andrex_box
%% calculate NPP per m2 per day per year... 
    % either: average daily rate per month (NPP_av_mgm2) for each year
    % or NPP_rates_fromdaily/# days in year
            % findyear, if findyear >12, then sum column 3 of time_end_all where column 1 = yix




