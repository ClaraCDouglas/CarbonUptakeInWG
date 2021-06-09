%% NPP area-weighted rate 
clearvars
cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
load('vgpm_imported.mat', 'vgpm_npp_all') % as from OceanProductivity site - average (A-W) daily rates per month
load('ProcessedData.mat', 'timedec')

% Average NPP during ice free conditions
    % monthly climatology of average daily rates for whole monthly time series
temp.vgpm=vgpm_npp_all;
temp.findneg=find(temp.vgpm<0);
temp.vgpm(temp.findneg)=NaN;
NPP.vgpm_av_day_nan=nanmean(temp.vgpm,3); 

    % annual climatology of average daily rates per year
NPP.vgpm_annual_day_nan(:,:,length([2003:2019]))=NaN(size(temp.vgpm(:,:,1)));    
for yix = 2003:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    if size (findyear)<12
    NPP_years.vgpm_annual_day_nan(:,:,yix-2002)=NaN(size(temp.vgpm(:,:,1)));
    else
    NPP_years.vgpm_annual_day_nan(:,:,yix-2002)=nanmean(temp.vgpm(:,:,findyear),3);    
    end
end
NPP.vgpm_annual_av_day_nan=nanmean(NPP_years.vgpm_annual_day_nan,3);

% Average NPP over whole year
    % monthly climatology of average daily rates for whole monthly time series
temp.vgpm_zeros=vgpm_npp_all;
temp.findneg=find(temp.vgpm_zeros<0);
temp.vgpm_zeros(temp.findneg)=0;
NPP.vgpm_av_day_zeros=nanmean(temp.vgpm_zeros,3); 

    % annual climatology of average daily rates per year
NPP.vgpm_annual_day_zeros(:,:,length([2003:2019]))=zeros(size(temp.vgpm(:,:,1)));    
for yix = 2003:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    if size (findyear)<12
    NPP_years.vgpm_annual_day_zeros(:,:,yix-2002)=NaN(size(temp.vgpm_zeros(:,:,1)));
    else
    NPP_years.vgpm_annual_day_zeros(:,:,yix-2002)=nanmean(temp.vgpm_zeros(:,:,findyear),3);    
    end
end
NPP.vgpm_annual_av_day_zeros=nanmean(NPP_years.vgpm_annual_day_zeros,3);

%% Total NPP 
% VGPM_npp_tot_gC_all - land/permanent ice are set as zeros
% VGPM_npp_tot_gC_nans - land/permanent ice are set as NaNs
load('vgpm_imported.mat', 'VGPM_npp_tot_gC_all', 'VGPM_npp_tot_gC_nans')
    % monthly climatology of total NPP per pixel for whole time series
NPP.vgpm_tot_monthclim=nanmean(VGPM_npp_tot_gC_all,3);
NPP.vgpm_tot_monthclim(find(NPP.vgpm_tot_monthclim==0))=NaN;

    % annual climatology of total NPP per year
temp.vgpmtot=VGPM_npp_tot_gC_nans;
% vgpm_av_tot_nan=nanmean(temp.vgpmtot,3); % this has calculated the average total NPP in each pixel
NPP_years.vgpm_tot_years(:,:,length([2003:2019]))=NaN(size(temp.vgpmtot(:,:,1)));    
for yix = 2003:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    if size (findyear)<12
    NPP_years.vgpm_tot_years(:,:,yix-2002)=NaN(size(temp.vgpmtot(:,:,1)));
    else
    NPP_years.vgpm_tot_years(:,:,yix-2002)=nansum(temp.vgpmtot(:,:,findyear),3);    
    end
end

NPP.vgpm_tot_av_years=nanmean(NPP_years.vgpm_tot_years,3);
NPP.vgpm_tot_av_years(find(NPP.vgpm_tot_av_years==0))=NaN;

% % change to Tg
% NPP.vgpm_tot_monthclim=NPP.vgpm_tot_monthclim;
% NPP.vgpm_tot_av_years=NPP.vgpm_tot_av_years;

% Get proportion - contribution of each pixel to total WG NPP
load('ProcessedData.mat', 'OceanProd')
WG_annual_tot=OceanProd.vgpm.Weddell.NPP_tot_TgC_annual(2:18,2);
WG_annual_tot=WG_annual_tot.*1e12 % back to g C

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

test4=NaN(size(NPP_years.vgpm_tot_years(:,:,:)));
for yix=2003:2019
    test=NPP_years.vgpm_tot_years(:,:,yix-2002).*temp.Weddell.box_logic;
    test(test==0)=NaN;
    test2=test./WG_annual_tot(yix-2002);
    test3=nansum(nansum(test2))
    test4(:,:,yix-2002)=test2;
    % then in polygon and set 0 to nan
    % then add the values to make sure they add to one as a sanity check
end

NPP.test5=nanmean(test4,3);
%% Maps
% load the files needed maps
load('box_lat_lons.mat', 'andrex_box')
load('Just2000mContour2.mat')
load('openshelf_coord.mat')
load('Colormap_Delta_NPP.mat')
load('latlon_m.mat')
load('SouthernBoundary_Orsi1995.mat')

LAT_MIN = -80; LAT_MAX = -50; LONG_MIN = -75; LONG_MAX = 40;
xin=lon_m; yin=lat_m;  
Proj_List = {'Stereographic';'Orthographic';'Azimuthal Equal-area';...
    'Azimuthal Equidistant';'Gnomonic';'Satellite';...
    'Albers Equal-Area Conic';'Lambert Conformal Conic';'Mercator';...
    'Miller Cylindrical';'Equidistant Cylindrical';'Oblique Mercator';...
    'Transverse Mercator';'Sinusoidal';'Gall-Peters';'Hammer-Aitoff';...
    'Mollweide';'Robinson';'UTM';};
proj=8;
Projection = Proj_List{proj};
% Projection inputs
if proj < 7
    %     m_proj(Projection, 'lon',(LONG_MIN+LONG_MAX)./2,'lat',(LAT_MIN+LAT_MAX)./2,'rad',90 | -30 -50, 'alt',1); %For projs 1-6\
    m_proj(Projection, 'lon',0,'lat',-90,'rad',-70, 'alt',45); %For projs 1-6\
elseif proj == 7 || proj == 8
    m_proj(Projection,'lon',[LONG_MIN LONG_MAX],'lat',[LAT_MIN LAT_MAX],'clo',-15,'par',[-80 -45],'ell','normal'); %For projs 7-8
elseif proj == 9 || proj == 10 || proj == 11
    m_proj(Projection,'lon',[LONG_MIN LONG_MAX],'lat',[LAT_MIN LAT_MAX]); %For projs 9-11
elseif proj == 12
    m_proj(Projection,'lon',[-200 LONG_MAX],'lat',[LAT_MIN LAT_MAX],'asp',100,'dir','horizontal') % For projs 12
elseif proj == 19
    m_proj(Projection,'lon',[LONG_MIN LONG_MAX],'lat',[LAT_MIN LAT_MAX],'hem',0,'ell','grs80'); %For proj 19
else
    m_proj(Projection,'lon',[LONG_MIN LONG_MAX],'lat',[LAT_MIN LAT_MAX],'clo',-30); %For projs 13-18 -30
end

NPP_names=fieldnames(NPP);

for zix=2:length(NPP_names)
    NPP_name=NPP_names{zix};
    zin=NPP.(NPP_name);
    
    figure(zix); clf; hold on

    m_pcolor(xin,yin,zin)
    shading interp
    set(0,'DefaultAxesColor','none')
    hh=colorbar;
    switch zix
        case 1 || 2
            caxis([0 400]);
            colormap(ColourmapDelta_NPP)
            title(NPP_name)
            ylabel(hh,'mg C m^-^2 day^-^1','FontSize',12);

        case 3 || 4
            caxis([0 200]);
            colormap(ColourmapDelta_NPP)
            title(NPP_name)
            ylabel(hh,'mg C m^-^2 day^-^1','FontSize',12);
        case 5 || 6
            caxis([0 9e9]);
            colormap(ColourmapDelta_NPP)
            title(NPP_name)
            ylabel(hh,'g C','FontSize',12);
        case 7
            caxis([0 0.000075]);
            colormap(flipud(hot))
            title('Proportion of total WG NPP')
            ylabel(hh,'Proportion','FontSize',12);
    end
% ylabel(hh,'mg C m^-^2 day^-^1','FontSize',12);

%  hh.Location='southoutside'; % want to make it a bit shorter still
% hh.Position = [0.25 0.13 0.4 0.05];
% hh.Position = [0.47 0.15 0.4 0.05];
%         pos2=[0.91 0.35 0.025 0.39]; % colorbar
%         hh.Position=pos2;
%         hh.Label.Rotation = 0;
%         hh.Label.Position = [0.7 10.65 0.7];

% plot 2000m isobath
v=[-2000,-2000]; % this sets level at -2000m isobath 
m_tbase('contour',v,'edgecolor','k','linewidth',0.8);
    % just to get a black line for the legend:
    contour_blank=[0 1 2;0 1 2]
    contour_blank_draw=m_plot(contour_blank(1,:),contour_blank(2,:),'color','k','linewidth',0.8)
% plot SB front
sb_line=m_plot(sbdy(:,1),sbdy(:,2),'color',[.7 .7 .7],'linewi',1.5,'LineStyle','-')
% plot region box and SR/OO line
% graphical.shelfcolor=[0.9290 0.6940 0.1250];
% graphical.opencolor=[0.6 0.4 0.8];

% OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',5)%,'LineStyle','--')
% SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',3.5)%'#80471C'

SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',5)%'#80471C'
OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',3.5)%,'LineStyle','--')

m_coast('patch',[0.7 0.7 0.7]);
m_grid('box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);
% caxis([0 500]);
set(gca,'FontSize',12)

get(0,'ScreenSize')
% set(gcf,'color','none','position',[150 80 850 620])
set(gcf,'color','white','position',[150 80 850 620])

set(gca,'color','none')
% title('Average Net Primary Production [MODIS VGPM] (2002 to 2020)')
set(get(gca,'title'),'Position',[-0.0169,0.25,-5000000000000000])

lgd_map=legend([contour_blank_draw sb_line SR_line OO_line],'2000m Isobath','Southern Boundary','Shelf Region','Open Ocean Region','Location',[0.5 0.17 0.4 0.05]);
legend('boxoff')

end

%% individually..

plot_folder=['C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\figures\'];

%% climatology of annual total 
zin=NPP.vgpm_tot_av_years;

figure(10); clf; hold on

m_pcolor(xin,yin,zin)
shading interp
set(0,'DefaultAxesColor','none')
hh=colorbar;
caxis([0 9e9]);
colormap(ColourmapDelta_NPP)
title('Climatology of Total Annual NPP (g C)')
ylabel(hh,'g C','FontSize',12);

%  hh.Location='southoutside'; % want to make it a bit shorter still
% hh.Position = [0.25 0.13 0.4 0.05];
% hh.Position = [0.47 0.15 0.4 0.05];
%         pos2=[0.91 0.35 0.025 0.39]; % colorbar
%         hh.Position=pos2;
%         hh.Label.Rotation = 0;
%         hh.Label.Position = [0.7 10.65 0.7];

% plot 2000m isobath
v=[-2000,-2000]; % this sets level at -2000m isobath
m_tbase('contour',v,'edgecolor','k','linewidth',0.8);
% just to get a black line for the legend:
contour_blank=[0 1 2;0 1 2]
contour_blank_draw=m_plot(contour_blank(1,:),contour_blank(2,:),'color','k','linewidth',0.8)
% plot SB front
sb_line=m_plot(sbdy(:,1),sbdy(:,2),'color',[.7 .7 .7],'linewi',1.5,'LineStyle','-')
% plot region box and SR/OO line
% graphical.shelfcolor=[0.9290 0.6940 0.1250];
% graphical.opencolor=[0.6 0.4 0.8];

% OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',5)%,'LineStyle','--')
% SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',3.5)%'#80471C'

SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',5)%'#80471C'
OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',3.5)%,'LineStyle','--')

m_coast('patch',[0.7 0.7 0.7]);
m_grid('box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);
% caxis([0 500]);
set(gca,'FontSize',12)

get(0,'ScreenSize')
% set(gcf,'color','none','position',[150 80 850 620])
set(gcf,'color','white','position',[150 80 850 620])

set(gca,'color','none')
% title('Average Net Primary Production [MODIS VGPM] (2002 to 2020)')
set(get(gca,'title'),'Position',[-0.0169,0.25,-5000000000000000])

lgd_map=legend([contour_blank_draw sb_line SR_line OO_line],'2000m Isobath','Southern Boundary','Shelf Region','Open Ocean Region','Location',[0.5 0.17 0.4 0.05]);
legend('boxoff')

print('-dpng',[plot_folder ['Map_TotNPP_Climatology' '.png']])

%% climatology of annual daily average in ice-free conditions 
zin=NPP.vgpm_annual_av_day_nan;

figure(11); clf; hold on

m_pcolor(xin,yin,zin)
shading interp
set(0,'DefaultAxesColor','none')
hh=colorbar;
caxis([0 400]);
colormap(ColourmapDelta_NPP)
title('Annual Climatology of Daily Average NPP in Ice-Free Conditions')
ylabel(hh,'mg C m^-^2 day^-^1','FontSize',12);
%  hh.Location='southoutside'; % want to make it a bit shorter still
% hh.Position = [0.25 0.13 0.4 0.05];
% hh.Position = [0.47 0.15 0.4 0.05];
%         pos2=[0.91 0.35 0.025 0.39]; % colorbar
%         hh.Position=pos2;
%         hh.Label.Rotation = 0;
%         hh.Label.Position = [0.7 10.65 0.7];

% plot 2000m isobath
v=[-2000,-2000]; % this sets level at -2000m isobath
m_tbase('contour',v,'edgecolor','k','linewidth',0.8);
% just to get a black line for the legend:
contour_blank=[0 1 2;0 1 2]
contour_blank_draw=m_plot(contour_blank(1,:),contour_blank(2,:),'color','k','linewidth',0.8)
% plot SB front
sb_line=m_plot(sbdy(:,1),sbdy(:,2),'color',[.7 .7 .7],'linewi',1.5,'LineStyle','-')
% plot region box and SR/OO line
% graphical.shelfcolor=[0.9290 0.6940 0.1250];
% graphical.opencolor=[0.6 0.4 0.8];

% OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',5)%,'LineStyle','--')
% SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',3.5)%'#80471C'

SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',5)%'#80471C'
OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',3.5)%,'LineStyle','--')

m_coast('patch',[0.7 0.7 0.7]);
m_grid('box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);
% caxis([0 500]);
set(gca,'FontSize',12)

get(0,'ScreenSize')
% set(gcf,'color','none','position',[150 80 850 620])
set(gcf,'color','white','position',[150 80 850 620])

set(gca,'color','none')
% title('Average Net Primary Production [MODIS VGPM] (2002 to 2020)')
set(get(gca,'title'),'Position',[-0.0169,0.25,-5000000000000000])

lgd_map=legend([contour_blank_draw sb_line SR_line OO_line],'2000m Isobath','Southern Boundary','Shelf Region','Open Ocean Region','Location',[0.5 0.17 0.4 0.05]);
legend('boxoff')

print('-dpng',[plot_folder ['Map_DailyNPP_IceFree_Climatology' '.png']])

%% climatology of annual daily average in ice-free conditions
NPP.vgpm_annual_av_day_zeros(NPP.vgpm_annual_av_day_zeros==0)=NaN;
zin=NPP.vgpm_annual_av_day_zeros;

figure(12); clf; hold on

m_pcolor(xin,yin,zin)
shading interp
set(0,'DefaultAxesColor','none')
hh=colorbar;
caxis([0 150]);
colormap(ColourmapDelta_NPP)
title('Annual Climatology of Daily Average NPP Over Whole Year')
ylabel(hh,'mg C m^-^2 day^-^1','FontSize',12);
%  hh.Location='southoutside'; % want to make it a bit shorter still
% hh.Position = [0.25 0.13 0.4 0.05];
% hh.Position = [0.47 0.15 0.4 0.05];
%         pos2=[0.91 0.35 0.025 0.39]; % colorbar
%         hh.Position=pos2;
%         hh.Label.Rotation = 0;
%         hh.Label.Position = [0.7 10.65 0.7];

% plot 2000m isobath
v=[-2000,-2000]; % this sets level at -2000m isobath
m_tbase('contour',v,'edgecolor','k','linewidth',0.8);
% just to get a black line for the legend:
contour_blank=[0 1 2;0 1 2]
contour_blank_draw=m_plot(contour_blank(1,:),contour_blank(2,:),'color','k','linewidth',0.8)
% plot SB front
sb_line=m_plot(sbdy(:,1),sbdy(:,2),'color',[.7 .7 .7],'linewi',1.5,'LineStyle','-')
% plot region box and SR/OO line
% graphical.shelfcolor=[0.9290 0.6940 0.1250];
% graphical.opencolor=[0.6 0.4 0.8];

% OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',5)%,'LineStyle','--')
% SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',3.5)%'#80471C'

SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',5)%'#80471C'
OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',3.5)%,'LineStyle','--')

m_coast('patch',[0.7 0.7 0.7]);
m_grid('box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);
% caxis([0 500]);
set(gca,'FontSize',12)

get(0,'ScreenSize')
% set(gcf,'color','none','position',[150 80 850 620])
set(gcf,'color','white','position',[150 80 850 620])

set(gca,'color','none')
% title('Average Net Primary Production [MODIS VGPM] (2002 to 2020)')
set(get(gca,'title'),'Position',[-0.0169,0.25,-5000000000000000])

lgd_map=legend([contour_blank_draw sb_line SR_line OO_line],'2000m Isobath','Southern Boundary','Shelf Region','Open Ocean Region','Location',[0.5 0.17 0.4 0.05]);
legend('boxoff')

print('-dpng',[plot_folder ['Map_DailyNPP_AllYear_Climatology' '.png']])

%% climatology of annual daily average in ice-free conditions
zin=NPP.test5;

figure(13); clf; hold on

m_pcolor(xin,yin,zin)
shading interp
set(0,'DefaultAxesColor','none')
hh=colorbar;
caxis([0 0.000075]);
colormap(flipud(hot))
title('Climatology of the Contribution/Proporion to Total WG NPP')
ylabel(hh,'Proportion','FontSize',12);
%  hh.Location='southoutside'; % want to make it a bit shorter still
% hh.Position = [0.25 0.13 0.4 0.05];
% hh.Position = [0.47 0.15 0.4 0.05];
%         pos2=[0.91 0.35 0.025 0.39]; % colorbar
%         hh.Position=pos2;
%         hh.Label.Rotation = 0;
%         hh.Label.Position = [0.7 10.65 0.7];

% plot 2000m isobath
v=[-2000,-2000]; % this sets level at -2000m isobath
m_tbase('contour',v,'edgecolor','k','linewidth',0.8);
% just to get a black line for the legend:
contour_blank=[0 1 2;0 1 2]
contour_blank_draw=m_plot(contour_blank(1,:),contour_blank(2,:),'color','k','linewidth',0.8)
% plot SB front
sb_line=m_plot(sbdy(:,1),sbdy(:,2),'color',[.7 .7 .7],'linewi',1.5,'LineStyle','-')
% plot region box and SR/OO line
% graphical.shelfcolor=[0.9290 0.6940 0.1250];
% graphical.opencolor=[0.6 0.4 0.8];

% OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',5)%,'LineStyle','--')
% SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',3.5)%'#80471C'

SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',5)%'#80471C'
OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',3.5)%,'LineStyle','--')

m_coast('patch',[0.7 0.7 0.7]);
m_grid('box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);
% caxis([0 500]);
set(gca,'FontSize',12)

get(0,'ScreenSize')
% set(gcf,'color','none','position',[150 80 850 620])
set(gcf,'color','white','position',[150 80 850 620])

set(gca,'color','none')
% title('Average Net Primary Production [MODIS VGPM] (2002 to 2020)')
set(get(gca,'title'),'Position',[-0.0169,0.25,-5000000000000000])

lgd_map=legend([contour_blank_draw sb_line SR_line OO_line],'2000m Isobath','Southern Boundary','Shelf Region','Open Ocean Region','Location',[0.5 0.17 0.4 0.05]);
legend('boxoff')

print('-dpng',[plot_folder ['Map_Contribution2TotNPP_Climatology' '.png']])
