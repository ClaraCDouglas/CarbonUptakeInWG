cd C:\Users\ccd1n18\Documents\Projects\Carbon-Uptake-in-WG_Manuscript\Workspace_Variables
%% Map
load('vgpm_imported.mat', 'vgpm_npp_all','lat_m','lon_m')
% where NPP in ice covered areas = NaN:
% Average NPP during ice free conditions
temp.vgpm=vgpm_npp_all;
temp.findneg=find(temp.vgpm<0);
temp.vgpm(temp.findneg)=NaN;
vgpm_av_day_nan=nanmean(temp.vgpm,3);
% figure; pcolor(lon_m,lat_m,vgpm_av_day_nan); shading flat;      % test image

% load the files needed for next step
load('box_lat_lons.mat', 'andrex_box')
load('Just2000mContour2.mat')
load('openshelf_coord.mat')
load('Colormap_Delta_NPP.mat')
load('SouthernBoundary_Orsi1995.mat')
% Using 'VGPM_av_day_nan' for plots for now 

% Plot the map
LAT_MIN = -80; LAT_MAX = -50; LONG_MIN = -75; LONG_MAX = 40;

xin=lon_m; yin=lat_m; zin=vgpm_av_day_nan; 
% zin2=VGPM_av_day;
% zin3=VGPM_av_day_nan;


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


figure(2)
m_pcolor(xin,yin,zin) 
shading interp
set(0,'DefaultAxesColor','none')
caxis([0 400]); %250 %380 %300
colormap(ColourmapDelta_NPP)

hold on

hh=colorbar;
ylabel(hh,'mg C m^-^2 day^-^1','FontSize',12);
%  hh.Location='southoutside'; % want to make it a bit shorter still
% hh.Position = [0.25 0.13 0.4 0.05];
% hh.Position = [0.47 0.15 0.4 0.05];
%         pos2=[0.91 0.35 0.025 0.39]; % colorbar
%         hh.Position=pos2;
%         hh.Label.Rotation = 0;
%         hh.Label.Position = [0.7 10.65 0.7];

% plot 2000m isobath
       % To draw the contours at one height (k), 
        % specify levels as a two-element row vector [k k]
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
set(gcf,'color','none','position',[150 80 850 620])
set(gcf,'color','white','position',[150 80 850 620])

set(gca,'color','none')
title('')
title('Average Net Primary Production [MODIS VGPM] (2002 to 2020)')
set(get(gca,'title'),'Position',[-0.0169,0.25,-5000000000000000])

lgd_map=legend([contour_blank_draw sb_line SR_line OO_line],'2000m Isobath','Southern Boundary','Shelf Region','Open Ocean Region','Location',[0.5 0.17 0.4 0.05]);
legend('boxoff')

%% Regression
% just scatter plot with regression line included
%VARIABLES
graphical.shelfcolor=[0.8 0.4 0];
graphical.opencolor=[0.4 0.2 0.8];
graphical.SRsize=65;
graphical.OOsize=52;

clear scatterlm
scatterlm=figure(2);
scatterlm.Position=[150 80 700 700];
ax1=axes(scatterlm);
scatterSR=scatter((regres.ice.(region_sublist{2})/1e4),regres.NPP.(region_sublist{2}),graphical.SRsize,graphical.shelfcolor,'filled','s')
ylim([0 12.2])
xlim([0 14e4])
xlabel('Annual mean area of open ice-free water on shelf (x10^4 km^2)','Interpreter','tex','FontSize',12)
ylabel('Total Annual NPP on shelf (Tg C)','Interpreter','tex','FontSize',12)
hold on
ax1.XColor=graphical.shelfcolor;
ax1.YColor=graphical.shelfcolor;
    % add regression line
        % REGRESSION LINE VARIABLES
graphical.SRIcemax=max(regres.ice.(region_sublist{2})/1e4);
graphical.SRIcemin=min(regres.ice.(region_sublist{2})/1e4);
graphical.xLSR = graphical.SRIcemin:0.02:graphical.SRIcemax;
graphical.mLSR = table2array(regres.lmsmall.(region_sublist{2}).Coefficients(2,1));
graphical.cLSR = table2array(regres.lmsmall.(region_sublist{2}).Coefficients(1,1));
graphical.yLSR = graphical.mLSR * graphical.xLSR + graphical.cLSR;

plot(graphical.xLSR,graphical.yLSR,'Color',graphical.shelfcolor);
graphical.text1=text(3.6,6, 'R^2=0.83, p<0.001','Color',graphical.shelfcolor,'Interpreter','tex','FontSize',12);

% add 2nd axis with open ocean data
ax2=axes(scatterlm);
scatterOO=scatter(regres.icesmall.(region_sublist{3}),regres.NPP.(region_sublist{3}),graphical.OOsize,graphical.opencolor,'filled')
hold on
%'s',
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color='none';
ax2.XColor=graphical.opencolor;
ax2.YColor=graphical.opencolor;
xlabel('Annual mean area of open ice-free water in open ocean (x10^6 km^2)','Interpreter','tex','FontSize',12)
ylabel('Total Annual NPP in open ocean (Tg C)','Interpreter','tex','FontSize',12)
graphical.OOIcemax=max(regres.icesmall.(region_sublist{3}))
graphical.OOIcemin=min(regres.icesmall.(region_sublist{3}))
graphical.xLOO = graphical.OOIcemin:0.02:(graphical.OOIcemax+0.02);
graphical.mLOO = table2array(regres.lmsmall.(region_sublist{3}).Coefficients(2,1));
graphical.cLOO = table2array(regres.lmsmall.(region_sublist{3}).Coefficients(1,1));
graphical.yLOO = graphical.mLOO * graphical.xLOO + graphical.cLOO;
plot(graphical.xLOO,graphical.yLOO,'Color',graphical.opencolor);
lg2=legend([scatterSR, scatterOO],'Shelf','Open Ocean','Location','northwest')

graphical.text2=text(1.4,80, 'R^2=0.54, p<0.001','Color',graphical.opencolor,'Interpreter','tex','FontSize',12);

% calculate 95% confidence intervals

%% Regression separate
pos1=[0.06 0.1 0.40 0.4];
pos2=[0.55 0.1 0.40 0.4];
pos1=[0.06 0.1 0.40 0.6];
pos2=[0.55 0.1 0.40 0.6];

positionlist=[pos1;pos2];
scatterlm=figure(6);
scatterlm.Position=[150 80 1000 500];
% SR
subplot('Position',positionlist(1,:))
plot(regres.lmsmall.(region_sublist{2}))
xlim([0 14])
hold on
scatterSR=scatter((regres.ice.(region_sublist{2})/1e4),regres.NPP.(region_sublist{2}),graphical.SRsize,graphical.shelfcolor,'filled')
xlabel('Annual ice-free area (x10^4 km^2)','Interpreter','tex','FontSize',12)
ylabel('Annual NPP (Tg C)','Interpreter','tex','FontSize',12)
title('Shelf Region','Interpreter','tex','FontSize',12)
plot(graphical.xLSR,graphical.yLSR,'Color','k');
graphical.text1=text(2,11, 'R^2=0.83, p<0.001','Color','k','Interpreter','tex','FontSize',12);
legend('off')
% OO
subplot('Position',positionlist(2,:))
plot(regres.lmsmall.(region_sublist{3}))
% xlim([0 1.8])
hold on
scatterSR=scatter((regres.ice.(region_sublist{3})/1e6),regres.NPP.(region_sublist{3}),graphical.SRsize,graphical.opencolor,'filled')
xlabel('Annual ice-free area (x10^6 km^2)','Interpreter','tex','FontSize',12)
ylabel('Annual NPP(Tg C)','Interpreter','tex','FontSize',12)
title('Open Ocean','Interpreter','tex','FontSize',12)
plot(graphical.xLOO,graphical.yLOO,'Color','k');
graphical.text2=text(1.1,120, 'R^2=0.54, p<0.001','Color','k','Interpreter','tex','FontSize',12);
legend('off')


%% total annual NPP bar chart
graphical.shelfcolor=[0.8 0.4 0];
graphical.opencolor=[0.4 0.2 0.8];

aix=4;
graphical.x=2003:1:2019
    shelf_open_annual(:,1)=OceanProd.(algorithm{aix}).(region_sublist{2}).NPP_tot_TgC_annual(2:18,2);
    shelf_open_annual(:,2)=OceanProd.(algorithm{aix}).(region_sublist{3}).NPP_tot_TgC_annual(2:18,2);
graphical.y2=shelf_open_annual;

totalNPPplot=figure;
totalNPPplot.Position=[150 80 800 700];
b=bar(graphical.x,graphical.y2,'stacked');
  hold on
b(1).FaceColor=graphical.shelfcolor;
b(2).FaceColor=graphical.opencolor;
  xtickangle(45)
  set(gca, 'XTick', [2003:1:2018])
  ylim([0 150])
 ylabel('Annual NPP (Tg C)')
title('Annual NPP on Shelf and in Open Ocean')
legend('Shelf Region','Open Ocean', 'Position',[0.5 0.8 0.1 0.1])
