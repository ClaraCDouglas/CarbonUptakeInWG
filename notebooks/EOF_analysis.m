NPP_years.anom_annual_day_nan_test = NPP_years.vgpm_annual_day_nan(:,:,1) - NPP.vgpm_annual_av_day_nan;
% anomalies ~ -500-500 
figure;
pcolor(lon_m, lat_m, NPP_years.anom_annual_day_nan_test(:,:,1));
shading flat


NPP_years.anom_annual_day_nan = NPP_years.vgpm_annual_day_nan - NPP.vgpm_annual_av_day_nan;
% anomalies ~ -500-500 


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

years=2003:2019;
plot_folder=['C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\figures\'];

for mix=1%:17
    zin=NPP_years.anom_annual_day_nan(:,:,mix);
    
    figure(mix+1); clf; hold on
    
    m_pcolor(xin,yin,zin)
    shading interp
    set(0,'DefaultAxesColor','none')
    hh=colorbar;
    caxis([-300 300]);
    colormap(ColourmapDelta_NPP)
    [cmap] = setcmappete('rest',28,'reg'); 
    title("Anomaly of Daily Average NPP in Ice-Free Conditions " +(years(1,mix))+ "")
    ylabel(hh,'mg C m^-^2 day^-^1','FontSize',12);

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
    
    SR_line=m_plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2)%'#80471C'
    OO_line=m_plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2)%,'LineStyle','--')
    
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
print('-dpng',[plot_folder ['Anom_avdailyNPP',num2str(years(mix)), '.png']])
close all
end
clearvars ans cmap ColormapDelta_NPP contour* hh Just* LAT_MAX LAT_MIN lgd* LONG* mix OO* open* proj* Proj* sb* shelf* SR* v *in

%% EOF analysis
% process data to only have within WG box
% NPP_years.anom_annual_day_nan = NPP_years.vgpm_annual_day_nan - NPP.vgpm_annual_av_day_nan;

NPP_years.vgpm_annual_day_nan_BOX=NPP_years.vgpm_annual_day_nan.*temp.Weddell.box_logic;
NPP_years.vgpm_annual_day_nan_BOX(NPP_years.vgpm_annual_day_nan_BOX==0)=NaN;

% cut down to box around WG box

NPP_years.vgpm_annual_day_nan_BOX=NPP_years.vgpm_annual_day_nan_BOX(861:1012,708:1279,:);
% to check region:
    % lat_m_box=lat_m(861:1012,708:1279);
    % lon_m_box=lon_m(861:1012,708:1279);
    % figure;
    % pcolor(lon_m_box,lat_m_box,NPP_years.vgpm_annual_day_nan_BOX(:,:,1)); shading flat

years=years';

yearsARRAY=zeros(size(NPP_years.vgpm_annual_day_nan_BOX));
for tix=1:17
    yearsARRAY(:,:,tix)=years(tix);    
end

time=1:length(years);
time=time';
[ev_index,tda,pev,trends] = calc_pigup_EOF2(NPP_years.vgpm_annual_day_nan_BOX,time)

MR_box=NPP_years.vgpm_annual_day_nan_BOX(62:83,372:404,:);
[ev_index,tda,pev,trends] = calc_pigup_EOF2(MR_box,time)
% to check region:
    lat_m_MR=lat_m_box(62:83,372:404);
    lon_m_MR=lon_m_box(62:83,372:404);
    figure;
    pcolor(lon_m_MR,lat_m_MR,MR_box(:,:,1)); shading flat

    % plot EOF output
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1);
plot(pev,'Color',[0,0.2,0.4],'LineWidth',2);
ylabel('PEV');
xlabel('Mode Number')
title('Percentage Explained Variance','fontsize',14);
subplot(2,2,2);
plot(tda(:,1),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 1 Time Dependent Amplidtude','fontsize',14);
subplot(2,2,[3,4]);
pcolor(lon_m_MR,lat_m_MR,ev_index(:,:,1)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 1 Pattern','fontsize',14);
sgt = sgtitle('VGPM growing season average NPP near Maud Rise SVD EOF','FontWeight','bold');
sgt.FontSize = 18;

figure;
plot(pev,'Color',[0,0.2,0.4],'LineWidth',2);
ylabel('PEV');
xlabel('Mode Number')
title('Percentage Explained Variance','fontsize',14);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1);
plot(tda(:,1),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 1 Time Dependent Amplidtude','fontsize',14);
subplot(2,4,2);
plot(tda(:,2),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 2 TDA','fontsize',14);
subplot(2,4,3);
plot(tda(:,3),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 3 TDA','fontsize',14);
subplot(2,4,4);
plot(tda(:,4),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 4 TDA','fontsize',14);

subplot(2,4,5);
pcolor(lon_m_MR,lat_m_MR,ev_index(:,:,1)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 1 Pattern','fontsize',14);
subplot(2,4,6);
pcolor(lon_m_MR,lat_m_MR,ev_index(:,:,2)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 1 Pattern','fontsize',14);
subplot(2,4,7);
pcolor(lon_m_MR,lat_m_MR,ev_index(:,:,3)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 1 Pattern','fontsize',14);
subplot(2,4,8);
pcolor(lon_m_MR,lat_m_MR,ev_index(:,:,4)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 1 Pattern','fontsize',14);





% .. testing the function
index3d = NPP_years.vgpm_annual_day_nan_BOX;
time=time;





%% smaller than whole WG, larger than MR box
% keeping data outwith the WG box, because they might be part of the patterns

% 12.5W
NPP_years.vgpm_annual_day_nan_east=NPP_years.vgpm_annual_day_nan(861:1012,1005:1279,:);

[ev_index,tda,pev,trends] = calc_pigup_EOF2(NPP_years.vgpm_annual_day_nan_east,time)

    lat_m_E=lat_m(861:1012,1005:1279);
    lon_m_E=lon_m(861:1012,1005:1279);
    figure;
    pcolor(lon_m_E,lat_m_E,NPP_years.vgpm_annual_day_nan_east(:,:,1)); shading flat

    % plot EOF output
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1);
plot(pev,'Color',[0,0.2,0.4],'LineWidth',2);
ylabel('PEV');
xlabel('Mode Number')
title('Percentage Explained Variance','fontsize',14);
subplot(2,2,2);
plot(tda(:,1),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 1 Time Dependent Amplidtude','fontsize',14);
subplot(2,2,[3]);
pcolor(lon_m_E,lat_m_E,ev_index(:,:,1)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 1 Pattern','fontsize',14);
sgt = sgtitle('VGPM growing season average NPP east WG SVD EOF','FontWeight','bold');
sgt.FontSize = 18;


figure;
plot(pev,'Color',[0,0.2,0.4],'LineWidth',2);
ylabel('PEV');
xlabel('Mode Number')
title('Percentage Explained Variance','fontsize',14);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1);
plot(tda(:,1),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 1 Time Dependent Amplidtude','fontsize',14);
subplot(2,4,2);
plot(tda(:,2),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 2 TDA','fontsize',14);
subplot(2,4,3);
plot(tda(:,3),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 3 TDA','fontsize',14);
subplot(2,4,4);
plot(tda(:,4),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 4 TDA','fontsize',14);

subplot(2,4,5);
pcolor(lon_m_E,lat_m_E,ev_index(:,:,1)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 1 Pattern','fontsize',14);
subplot(2,4,6);
pcolor(lon_m_E,lat_m_E,ev_index(:,:,2)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 2 Pattern','fontsize',14);
subplot(2,4,7);
pcolor(lon_m_E,lat_m_E,ev_index(:,:,3)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 3 Pattern','fontsize',14);
subplot(2,4,8);
pcolor(lon_m_E,lat_m_E,ev_index(:,:,4)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 4 Pattern','fontsize',14);



NPP_years.vgpm_annual_day_nan_weast=NPP_years.vgpm_annual_day_nan(861:1012,930:1279,:);

[ev_index,tda,pev,trends] = calc_pigup_EOF2(NPP_years.vgpm_annual_day_nan_weast,time)

    lat_m_Ew=lat_m(861:1012,930:1279);
    lon_m_Ew=lon_m(861:1012,930:1279);
    figure;
    pcolor(lon_m_Ew,lat_m_Ew,NPP_years.vgpm_annual_day_nan_weast(:,:,1)); shading flat

    % plot EOF output
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1);
plot(pev,'Color',[0,0.2,0.4],'LineWidth',2);
ylabel('PEV');
xlabel('Mode Number')
title('Percentage Explained Variance','fontsize',14);
subplot(2,2,2);
plot(tda(:,1),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 1 Time Dependent Amplidtude','fontsize',14);
subplot(2,2,[3]);
pcolor(lon_m_Ew,lat_m_Ew,ev_index(:,:,1)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 1 Pattern','fontsize',14);
sgt = sgtitle('VGPM growing season average NPP more east WG SVD EOF','FontWeight','bold');
sgt.FontSize = 18;



NPP_years.vgpm_annual_day_nan_mweast=NPP_years.vgpm_annual_day_nan(861:1012,770:1279,:); % too big %850:1279 also too big

[ev_index,tda,pev,trends] = calc_pigup_EOF2(NPP_years.vgpm_annual_day_nan_mweast,time) % too big


NPP_years.vgpm_annual_day_nan_mweast=NPP_years.vgpm_annual_day_nan(861:1012,880:1279,:); 

[ev_index,tda,pev,trends] = calc_pigup_EOF2(NPP_years.vgpm_annual_day_nan_mweast,time) 
    lat_m_Ew=lat_m(861:1012,880:1279);
    lon_m_Ew=lon_m(861:1012,880:1279);
    figure;
    pcolor(lon_m_Ew,lat_m_Ew,NPP_years.vgpm_annual_day_nan_mweast(:,:,1)); shading flat

    % plot EOF output
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1);
plot(pev,'Color',[0,0.2,0.4],'LineWidth',2);
ylabel('PEV');
xlabel('Mode Number')
title('Percentage Explained Variance','fontsize',14);
subplot(2,2,2);
plot(tda(:,1),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 1 Time Dependent Amplidtude','fontsize',14);
subplot(2,2,[3]);
pcolor(lon_m_Ew,lat_m_Ew,ev_index(:,:,1)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 1 Pattern','fontsize',14);
sgt = sgtitle('VGPM growing season average NPP more east WG SVD EOF','FontWeight','bold');
sgt.FontSize = 18;


figure;
plot(pev,'Color',[0,0.2,0.4],'LineWidth',2);
ylabel('PEV');
xlabel('Mode Number')
title('Percentage Explained Variance','fontsize',14);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1);
plot(tda(:,1),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 1 Time Dependent Amplidtude','fontsize',14);
subplot(2,4,2);
plot(tda(:,2),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 2 TDA','fontsize',14);
subplot(2,4,3);
plot(tda(:,3),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 3 TDA','fontsize',14);
subplot(2,4,4);
plot(tda(:,4),'r','LineWidth',2);ylabel('PEV');
ylabel('TDA');
xlabel('Year #')
title('Mode 4 TDA','fontsize',14);

subplot(2,4,5);
pcolor(lon_m_Ew,lat_m_Ew,ev_index(:,:,1)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 1 Pattern','fontsize',14);
subplot(2,4,6);
pcolor(lon_m_Ew,lat_m_Ew,ev_index(:,:,2)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 2 Pattern','fontsize',14);
subplot(2,4,7);
pcolor(lon_m_Ew,lat_m_Ew,ev_index(:,:,3)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 3 Pattern','fontsize',14);
subplot(2,4,8);
pcolor(lon_m_Ew,lat_m_Ew,ev_index(:,:,4)); 
shading flat
hold on 
plot(3,-66,'w*')
colorbar
colormap(jet);
geoshow('landareas.shp','facecolor','k')
title('Mode 4 Pattern','fontsize',14);






