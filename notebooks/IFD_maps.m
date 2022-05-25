load('ProcessedCAFEArrays_May22.mat')
load('WAPSHelfOpenJan22.mat')
load('box_lat_lons.mat', 'andrex_box')

% ways to calc the mean in the most meaningful way... chosen to just stick
% with mean, and make sure colorbar is appropriate colour intervals to set
% 0-1 as white, meaning rest will show up.
%chosen way below 
    tempnans=IceFreeDays_peryear;
    tempnans(tempnans==0)=NaN;
    IFD_mean_nans=mean(tempnans(:,:,1:18),3,'omitnan');

    IFD_array=IceFreeDays_peryear;
    IFD_array(isnan(AnnualNPPRate_gperyear))=NaN;
    IFD_mean_nans2=mean(IFD_array(:,:,1:18),3,'omitnan');

    BB = reshape(AnnualNPPRate_gperyear(:,:,1:18),[],18);
    rowsrecord=[];
    for xix=1:length(BB)
       %disp(xix);
       if sum(isnan(BB(xix,:)))>0
           rowsrecord=cat(1,rowsrecord,xix);
       end
    end
    BB(rowsrecord,:)=-10;
    BBB= reshape(BB,1080,1380,18);
    IFD_arrayB=IceFreeDays_peryear;
    IFD_arrayB(BB==-10)=NaN;
    IFD_mean_nans2B=mean(IFD_arrayB(:,:,1:18),3,'omitnan');

% use this way:
IFD_mean=mean(IceFreeDays_peryear(:,:,1:18),3);
IFD_mean(IFD_mean==0)=NaN;
[IFD_max,I]=max(IceFreeDays_peryear(:,:,1:18),[],3);
IFD_max(isnan(IFD_mean))=NaN;
I(isnan(IFD_mean))=NaN;

% plot:
figure; tiledlayout(3,1)
nt1=nexttile;
pcolor(lon_wg,lat_wg,IFD_mean); shading flat
xlim([-70 35]); ylim([-80 -50])
caxis([0 250])
colorbar
geoshow('landareas.shp','facecolor',[0.8 0.8 0.8])
cmap=cmocean('balance',250,'pivot',120);
% cmap(1,:)=1;
colormap(nt1,cmap)
hold on
SR_line=plot(ShelfBox(:,1),ShelfBox(:,2),'color',[.21 .44 .31],'linewi',4);%'#80471C'
O_line=plot(OOBox(:,1),OOBox(:,2),'color',[0.44 0.4 0.91],'linewi',4);%,'LineStyle','--')
WAP_line=plot(WAPBox(:,1),WAPBox(:,2),'color',[0.8 0.4 0],'linewi',4);%,'LineStyle','--')
title('Mean IFD 2003-2020')

nt2=nexttile;
pcolor(lon_wg,lat_wg,IFD_max); shading flat
xlim([-70 35]); ylim([-80 -50])
caxis([0 250])
colorbar
geoshow('landareas.shp','facecolor',[0.8 0.8 0.8])
cmap=cmocean('balance',250,'pivot',120);
% cmap(1,:)=1;
colormap(nt2,cmap)
hold on
SR_line=plot(ShelfBox(:,1),ShelfBox(:,2),'color',[.21 .44 .31],'linewi',4);%'#80471C'
O_line=plot(OOBox(:,1),OOBox(:,2),'color',[0.44 0.4 0.91],'linewi',4);%,'LineStyle','--')
WAP_line=plot(WAPBox(:,1),WAPBox(:,2),'color',[0.8 0.4 0],'linewi',4);%,'LineStyle','--')
title('Max IFD (longest GS at each location) 2003-2020')

nt3=nexttile;
pcolor(lon_wg,lat_wg,I+2002); shading flat
xlim([-70 35]); ylim([-80 -50])
% caxis([0 250])
colorbar
geoshow('landareas.shp','facecolor',[0.8 0.8 0.8])
cmap=cmocean('balance',18);
% cmap(1,:)=1;
colormap(nt3,cmap)
hold on
SR_line=plot(ShelfBox(:,1),ShelfBox(:,2),'color',[.21 .44 .31],'linewi',4);%'#80471C'
O_line=plot(OOBox(:,1),OOBox(:,2),'color',[0.44 0.4 0.91],'linewi',4);%,'LineStyle','--')
WAP_line=plot(WAPBox(:,1),WAPBox(:,2),'color',[0.8 0.4 0],'linewi',4);%,'LineStyle','--')
title('Year when longest GS was seen at each location')


figure;
tiledlayout('flow');
for yix=2003:2020
    nexttile;
    pcolor(lon_wg,lat_wg,AnnualNPPRate_gperyear(:,:,yix-2002)); shading flat %
    xlim([-70 35]); ylim([-80 -50])
    caxis([0 260])
    colorbar
    geoshow('landareas.shp','facecolor',[0.8 0.8 0.8])
    cmap=cmocean('thermal',52);
    cmap(1,:)=1;
    colormap(cmap)
%     hold on
%     SR_line=plot(ShelfBox(:,1),ShelfBox(:,2),'color',[.21 .44 .31],'linewi',4);%'#80471C'
%     O_line=plot(OOBox(:,1),OOBox(:,2),'color',[0.44 0.4 0.91],'linewi',4);%,'LineStyle','--')
%     WAP_line=plot(WAPBox(:,1),WAPBox(:,2),'color',[0.8 0.4 0],'linewi',4);%,'LineStyle','--')
end


% Selecting values within region boxes only
for aix = 1:length(algorithm)
    for rix = 1:length(region_sublist)
        for yix=2003:2021
            temp_icefreedays=IceFreeDays_peryear(:,:,yix-2002);
            temp_icefreedays=round(temp_icefreedays);
            IceFree_pixels_years.(region_sublist{rix})(:,yix-2002)=temp_icefreedays(temp.(region_sublist{rix}).box_logic);
            
            temp_AnnualNPPRate=AnnualNPPRate_gperyear(:,:,yix-2002);
            AnnualNPPRate_pixels_years.(region_sublist{rix})(:,yix-2002)=temp_AnnualNPPRate(temp.(region_sublist{rix}).box_logic);
            Year_cols.(region_sublist{rix})(:,yix-2002)=AnnualNPPRate_pixels_years.(region_sublist{rix})(:,yix-2002);
            Year_cols.(region_sublist{rix})(:,yix-2002)=yix;
%             temp_AnAvDayRate=AnAvDayRate_mgm2d1(:,:,yix-2002);
%             AnAvDayRate_pixels_years.(region_sublist{rix})(:,yix-2002)=temp_AnAvDayRate(temp.(region_sublist{rix}).box_logic);
        end
    end
end
clearvars temp_*
