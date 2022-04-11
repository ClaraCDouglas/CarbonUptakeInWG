load('WAPSHelfOpenJan22.mat')
load('box_lat_lons.mat', 'andrex_box')
ShelfBox=ShelfMinusWAPJan22;%shelf_region_ANDbox
OOBox=OpenOceanMinusWAPJan22;%open_ocean_ANDbox
WAPBox=WAPJan22;
load('cafe_8day_imported_eqWG_withNaNs.mat', 'cafe_npp_tot_gC_nans_8day')
load('cafe_8day_imported_eqWG_withNaNs.mat', 'cafe_npp_all_8day')
load('cafe_8day_imported_eqWG_withNaNs.mat', 'lat_wg')
load('cafe_8day_imported_eqWG_withNaNs.mat', 'lon_wg')
load('cafe_8day_imported_eqWG_withNaNs.mat', 'time_end_all')
load('cafe_8day_imported_eqWG_withNaNs.mat', 'time_start_all')
load('cafe_8day_imported_eqWG_withNaNs.mat', 'timedec8day')
load('cafe_8day_imported_eqWG_withNaNs.mat', 'timedec8day_end')
data_daily=false;
data_8day=true;
cd 'D:\Data\SeaIceNIMBUS'; % SeaIce_daily_20022020.mat, SeaIce_8day_20022020.mat
if data_daily
    load('SeaIce_daily_20022020.mat')
    load('SeaIce_8day_20022020.mat','g_area','g_lat','g_lon')
elseif data_8day
    load('SeaIce_8day_20022020.mat')
end
if data_8day
    % time_start_all=time_start_all(1:848,:); % NPP time_start for each 8 day slice ending in 2020
    time_start_ice8=time_start_ice8(24:end,:); %there are more time slices for the sea ice data because the 2nd-4th weeks in Aug 2020 are missing from the NPP dataset
    time_end_ice8=time_end_ice8(24:end,:); %there are more time slices for the sea ice data because the 2nd-4th weeks in Aug 2020 are missing from the NPP dataset
    ice_conc_8day=ice_conc_8day(:,:,24:end); % remove ice data for Jan-June 2002 (no NPP data for that)
    
    % if removing the data that is also missing from NPP data:
    time_start_ice8(834:836,:)=[];
    time_end_ice8(834:836,:)=[];
    ice_conc_8day(:,:,834:836)=[];
    
    %make dec time %(and datenum)
    timedec8dayice_start=decyear(time_start_ice8);
    timedec8dayice_end=decyear(time_end_ice8);
    timedec8dayice=mean([timedec8dayice_start,timedec8dayice_end],2);
    % datenum8day=datenum(time_start_ice8);
end

get(0,'ScreenSize')

% for 46 tiles
% anindex=[0.05 0.15 0.3:0.2:0.7];
% colors={'y','r','c','k','g','m','b'};
anindex=[0.05 0.15];
colors={'k','r'};
cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\figures\TimeSlice_perYear_Maps\'
NPP=1;

for yix=2003%:2020
    yyix=yix-2002;
    addno=(yyix-1)*46;
    f=figure(yyix+20); 
%     f = figure('visible','off');
    tlo=tiledlayout('flow');
    f.WindowState='maximized';
    for ix=26%addno+1:1:addno+46%addno+9:1:addno+35
        ax=nexttile(tlo);
        
        switch NPP
            case 1 %rate
        temptemp=cafe_npp_all_8day(:,:,ix); %cafe_npp_tot_gC_nans_8day %cafe_npp_all_8day
        temptemp(temptemp<0)=NaN;
            case 2 %total
        temptemp=cafe_npp_tot_gC_nans_8day(:,:,ix); %cafe_npp_tot_gC_nans_8day %cafe_npp_all_8day
        temptemp(temptemp<0)=NaN;
        end
        
        pcolor(lon_wg,lat_wg,temptemp); shading flat
        hold on
%         SR_line=plot(ShelfBox(:,1),ShelfBox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
%         O_line=plot(OOBox(:,1),OOBox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')
%         WAP_line=plot(WAPBox(:,1),WAPBox(:,2),'color',[0 0.4 0.2],'linewi',2);%,'LineStyle','--')
        
        for icex=1:length(anindex)
            v=[anindex(icex),anindex(icex)];
            [c,h]=contourm(g_lat,g_lon,ice_conc_8day(:,:,ix),v,'LineWidth',2,'Color',colors{icex}); %
        end
        xlim([-70 35]); ylim([-80 -50])
        switch NPP
            case 1 %rate
                caxis([0 1000])
            case 2 %total
                caxis([0 5e8])
        end
        geoshow('landareas.shp','facecolor','k')
        ti=datetime(time_start_all(ix,1:3)); % should put the mid-date instead
        title(datestr(ti))
    end
    c=colorbar;
    switch NPP
        case 1
            c.Label.String='Daily NPP (mg m^{-2} d^{-1}';%Total NPP (gC)%Daily NPP (mg m^{-2} d^{-1}
        case 2
            c.Label.String='Total NPP (gC)';%Total NPP (gC)%Daily NPP (mg m^{-2} d^{-1}
    end
    sgti=num2str(yearrange0320(yyix));
    sgtitle(sgti)
    nexttile;
    for icex=1:length(anindex)
        plot(NaN,'LineWidth',2,'Color',colors{icex});
        hold on
    end
%     lgd=legend({'5%','15%','30%','50%','70%'},'FontSize',11)
    lgd=legend({'5%','15%'},'FontSize',11)
    title(lgd,'SIC Contours')
    % name=['figure',num2str(yix),'.jpg']
    switch NPP
        case 1 %rate
            hgexport(gcf, ['Austral_DayRate_',num2str(yix),'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
        case 2
            hgexport(gcf, ['Austral_Total_',num2str(yix),'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
    end
    close all
end


NPP=1;
anindex=[0.05 0.15 0.4];
colors={'k','r','b'};
% for data tiles
for yix=2003:2020
    yyix=yix-2002;
    addno=(yyix-1)*46;
    f=figure(yyix+20); 
%     f = figure('visible','off');
    tlo=tiledlayout('flow');
    f.WindowState='maximized';
    for ix=addno+9:1:addno+35
        ax=nexttile(tlo);
        switch NPP
            case 1 %rate
        temptemp=cafe_npp_all_8day(:,:,ix); %cafe_npp_tot_gC_nans_8day %cafe_npp_all_8day
        temptemp(temptemp<0)=NaN;
            case 2 %total
        temptemp=cafe_npp_tot_gC_nans_8day(:,:,ix); %cafe_npp_tot_gC_nans_8day %cafe_npp_all_8day
        temptemp(temptemp<0)=NaN;
        end
        pcolor(lon_wg,lat_wg,temptemp); shading flat
        hold on
        c=[[.21 .44 .31];[0.8 0.4 0];[0.44 0.78 0.91]];
        SR_line=plot(ShelfBox(:,1),ShelfBox(:,2),'color',[.21 .44 .31],'linewi',4);%'#80471C'
        O_line=plot(OOBox(:,1),OOBox(:,2),'color',[0.44 0.4 0.91],'linewi',4);%,'LineStyle','--')
        WAP_line=plot(WAPBox(:,1),WAPBox(:,2),'color',[0.8 0.4 0],'linewi',4);%,'LineStyle','--')
        for icex=1:length(anindex)
            v=[anindex(icex),anindex(icex)];
            [c,h]=contourm(g_lat,g_lon,ice_conc_8day(:,:,ix),v,'LineWidth',2,'Color',colors{icex}); %
        end
        xlim([-70 40]); ylim([-80 -40])
        switch NPP
            case 1 %rate
                caxis([0 1000])
            case 2 %total
                caxis([0 5e8])
        end
        geoshow('landareas.shp','facecolor','k')
        ti=datetime(time_start_all(ix,1:3)); % should put the mid-date instead
        title(datestr(ti))
    end
    cb=colorbar;
    switch NPP
        case 1
            cb.Label.String='Daily NPP (mg m^{-2} d^{-1}';%Total NPP (gC)%Daily NPP (mg m^{-2} d^{-1}
        case 2
            cb.Label.String='Total NPP (gC)';%Total NPP (gC)%Daily NPP (mg m^{-2} d^{-1}
    end
    sgti=num2str(yearrange0320(yyix));
    sgtitle(sgti)
    nexttile;
    for icex=1:length(anindex)
        plot(NaN,'LineWidth',2,'Color',colors{icex});
        hold on
    end
    yticks([])
%     lgd=legend({'5%','15%','30%','50%','70%'},'FontSize',11)
    lgd=legend({'5%','15%','40%'},'FontSize',11);
    title(lgd,'SIC Contours')
    % name=['figure',num2str(yix),'.jpg']
    switch NPP
        case 1 %rate
            hgexport(gcf, ['Visible_DayRate_',num2str(yix),'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
        case 2
            hgexport(gcf, ['Visible_Total_',num2str(yix),'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
    end
    close all
end


%% ice concentrations, stereo proj
for yix=2017:2018%2003:2020
    yyix=yix-2002;
    addno=(yyix-1)*46;
    f=figure(yyix+20); tiledlayout('flow');
    f.WindowState='maximized';
    for ix=addno+1:1:addno+46%addno+9:1:addno+35
        nexttile
        latlim = [-90 -50];
        lonlim = [-100 170];
        antarctica = shaperead('landareas.shp', 'UseGeoCoords', true,...
            'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});
        worldmap(latlim,lonlim)
        %     patchm(antarctica.Lat, antarctica.Lon, [0.5 1 0.5])
        %     geoshow(antarctica)
        hold on
        pcolorm(g_lat,g_lon,ice_conc_8day(:,:,ix))
        ti=datetime(time_start_all(ix,1:3)); % should put the mid-date instead
        title(datestr(ti))
    end
    colorbar
end



