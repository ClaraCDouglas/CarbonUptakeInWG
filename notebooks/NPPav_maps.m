load('ProcessedCAFEArrays_May22.mat')
load('WAPSHelfOpenJan22.mat')
load('box_lat_lons.mat', 'andrex_box')
load('SouthernBoundary_Orsi1995.mat')

NPPAnRate_mean=mean(AnnualNPPRate_gperyear(:,:,1:18),3,'omitnan');
[NPPAnRate_max,I]=max(AnnualNPPRate_gperyear(:,:,1:18),[],3,'omitnan');

% to get the average for pixels across the whole time, including when they
% were ice covered
NPPrtemp=AnnualNPPRate_gperyear(:,:,1:18);
NPPrtemp(isnan(NPPrtemp))=0;
NPPAnRate_mean2=mean(NPPrtemp(:,:,:),3);
NPPAnRate_mean2(NPPAnRate_mean2==0)=NaN;


% plot:
figure; tiledlayout(3,1)
nt1=nexttile;
pcolor(lon_wg,lat_wg,NPPAnRate_mean2); shading flat
xlim([-70 35]); ylim([-80 -50])
caxis([0 100])
colorbar
geoshow('landareas.shp','facecolor',[0.8 0.8 0.8])
cmap=cmocean('delta');%,250,'pivot',120);
% cmap(1,:)=1;
colormap(nt1,cmap)
hold on
SR_line=plot(ShelfBox(:,1),ShelfBox(:,2),'color',[.21 .44 .31],'linewi',4);%'#80471C'
O_line=plot(OOBox(:,1),OOBox(:,2),'color',[0.44 0.4 0.91],'linewi',4);%,'LineStyle','--')
WAP_line=plot(WAPBox(:,1),WAPBox(:,2),'color',[0.8 0.4 0],'linewi',4);%,'LineStyle','--')
title('2003-2020 Mean Annual NPP (g m^{-2} a^{-1})')

nt2=nexttile;
pcolor(lon_wg,lat_wg,NPPAnRate_max); shading flat %
xlim([-70 35]); ylim([-80 -50])
caxis([0 100])
colorbar
geoshow('landareas.shp','facecolor',[0.8 0.8 0.8])
cmap=cmocean('delta');
% cmap(1,:)=1;
colormap(nt2,cmap)
hold on
SR_line=plot(ShelfBox(:,1),ShelfBox(:,2),'color',[.21 .44 .31],'linewi',4);%'#80471C'
O_line=plot(OOBox(:,1),OOBox(:,2),'color',[0.44 0.4 0.91],'linewi',4);%,'LineStyle','--')
WAP_line=plot(WAPBox(:,1),WAPBox(:,2),'color',[0.8 0.4 0],'linewi',4);%,'LineStyle','--')
title('Max Annual NPP (g m^{-2} a^{-1}) 2003-2020')

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
title('Year when max annual NPP was seen at each location')


figure;
tiledlayout('flow');
for yix=2003:2020
    nexttile;
    pcolor(lon_wg,lat_wg,AnnualNPPRate_gperyear(:,:,yix-2002)); shading flat %
    xlim([-70 35]); ylim([-80 -50])
    caxis([0 100])
    colorbar
    geoshow('landareas.shp','facecolor',[0.8 0.8 0.8])
    cmap=cmocean('delta');
    colormap(cmap)
    hold on
    AND_line=plot(andrex_box(:,1),andrex_box(:,2),'color','m','linewi',2);%'#80471C'
    title(yix)
end
sgtitle('Area Normalised Annual NPP (gC m^{-2} a^{-1})')

figure; 
pcolor(lon_wg,lat_wg,NPPAnRate_mean2); shading flat
xlim([-70 35]); ylim([-80 -50])
caxis([0 100])
colorbar
geoshow('landareas.shp','facecolor',[0.8 0.8 0.8])
cmap=cmocean('delta');%,250,'pivot',120);
% cmap(1,:)=1;
colormap(cmap)
hold on
SR_line=plot(ShelfBox(:,1),ShelfBox(:,2),'color',[.21 .44 .31],'linewi',4);%'#80471C'
O_line=plot(OOBox(:,1),OOBox(:,2),'color',[0.44 0.4 0.91],'linewi',4);%,'LineStyle','--')
WAP_line=plot(WAPBox(:,1),WAPBox(:,2),'color',[0.8 0.4 0],'linewi',4);%,'LineStyle','--')
title('2003-2020 Mean Annual NPP (g m^{-2} a^{-1})')
set(gca,'FontSize',14)

%% or with mmap
LAT_MIN = -80; LAT_MAX = -50; LONG_MIN = -75; LONG_MAX = 35;
xin=lon_wg; yin=lat_wg; zin=NPPAnRate_mean2; zin2=NPPAnRate_mean;
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

figure;
m_pcolor(xin,yin,zin)
shading interp
set(0,'DefaultAxesColor','none')
caxis([0 100]); 
cmap=cmocean('delta');%,250,'pivot',120);
% cmap(1,:)=1;
colormap(cmap)
hold on
% [C,h]=m_contour(xin,yin,zin,[50:50:100],'w');
% [C2,h2]=m_contour(xin,yin,zin,[200:100:400],'w');
% clabel(C,h,'Color','w','FontSize',12)
% clabel(C2,h2,'Color','w','FontSize',12)
% [C,h]=m_contourf(xin,yin,zin,[0:50:150 200:100:600 1000],'k');
% clabel(C,h,'Color','k','FontSize',13)
hh=colorbar;
ylabel(hh,'g C m^-^2 a^-^1','FontSize',12);
hh.Location='southoutside'; % want to make it a bit shorter still
hh.Position = [0.3 0.07 0.4 0.05];
hh.Position = [0.3 0.1 0.4 0.05];

m_plot(sbdy(:,1),sbdy(:,2),'color',[.8 .8 .8],'linewi',1.25,'LineStyle','-')
SR_line=m_plot(ShelfBox(:,1),ShelfBox(:,2),'color',[.21 .44 .31],'linewi',4);%'#80471C'
O_line=m_plot(OOBox(:,1),OOBox(:,2),'color',[0.44 0.4 0.91],'linewi',4);%,'LineStyle','--')
WAP_line=m_plot(WAPBox(:,1),WAPBox(:,2),'color',[0.8 0.4 0],'linewi',4);%,'LineStyle','--')
m_coast('patch',[0.7 0.7 0.7]);
m_grid('box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);
% caxis([0 500]);
set(gca,'FontSize',12)

get(0,'ScreenSize')
set(gcf,'color','w','position',[500 80 1000 800])

set(gca,'color','none')

title('2003-2020 Mean NPP (g m^{-2} a^{-1})')
% title('Mean Annual NPP (g m^{-2} a^{-1})')
set(get(gca,'title'),'Position',[-0.0169,0.25,-5000000000000000])
