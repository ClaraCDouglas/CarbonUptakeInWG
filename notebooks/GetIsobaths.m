%% forget it for now, just edit normal line for region
% edited in excel to remove individual points whole looking at maps in
% matlab to match up contours.
load('Just2000mContour2.mat')
save('Just2000mContour_cleaned21.txt', 'Just2000mContour2', '-ASCII','-append');
save('ANDREX.txt', 'andrex_box', '-ASCII','-append');

% loaded in SR, OO and 2000m isobath contour as txt files (in data/interim)
open_ocean_ANDbox=OpenOceanBoxcleaned21;
shelf_region_ANDbox=ShelfRegionBoxcleaned21;
Just2000mContour=Just2000mContourcleaned21;
save('openshelfisobath_clean21.mat','open_ocean_ANDbox','shelf_region_ANDbox','Just2000mContour');

%% Creating boundary for shelf region vs open ocean using 2000m contour (50 m either side to get better contour line)
load('box_lat_lons.mat')
LONG_MIN = -180; LONG_MAX = 180; LAT_MIN = -90; LAT_MAX = 90;
[Z,LONG,LAT]=m_tbase([LONG_MIN LONG_MAX LAT_MIN LAT_MAX]);

IN_AND_elev=inpolygon(LONG,LAT,andrex_box(:,1),andrex_box(:,2));
findelevation=find(IN_AND_elev==1);

% findshelf=find(Z>-2000 & Z<0);
% findopen=find(Z<-2000);
findcontour=find(Z==-2000);
findcontour20=find(Z<-1980 & Z>-2020);

% shelfmask=zeros(2160,4321);
% shelfmask(findshelf)=1;
% openmask=zeros(2160,4321);
% openmask(findopen)=1;
contourmask=zeros(2160,4321);
contourmask(findcontour)=1;
contourmask20=zeros(2160,4321);
contourmask20(findcontour20)=1;

    eval(['box_elev=LAT;']);
    eval(['box_elev(:)=0;']);
    eval(['box_elev(findelevation)=1;']);
    
%     eval(['mask_elev_shelf=shelfmask.*(box_elev);']);
%     eval(['mask_elev_open=openmask.*(box_elev);']);
    eval(['mask_elev_contour=contourmask.*(box_elev);']);
    eval(['mask_elev_contour20=contourmask20.*(box_elev);']);

find_elev_position = find(mask_elev_contour==1);
find_elev_position20 = find(mask_elev_contour20==1);

    %LAT2000=LAT(find_elev_position);
    %LONG2000=LONG(find_elev_position);
    %CONTOUR2000=[LONG2000 LAT2000];
    %save('CONTOUR2000.txt', 'CONTOUR2000', '-ASCII','-append');
LAT20=LAT(find_elev_position20);
LONG20=LONG(find_elev_position20);
CONTOUR1980_2020=[LONG20 LAT20];

LAT=LAT(find_elev_position);
LONG=LONG(find_elev_position);
CONTOUR2000=[LONG LAT];

save('CONTOUR1950_2050.txt', 'CONTOUR1950_2050', '-ASCII','-append');


%% Have a looksie

LAT_MIN = -80; LAT_MAX = -50; LONG_MIN = -75; LONG_MAX = 40;

Proj_List = {'Stereographic';'Orthographic';'Azimuthal Equal-area';...
    'Azimuthal Equidistant';'Gnomonic';'Satellite';...
    'Albers Equal-Area Conic';'Lambert Conformal Conic';'Mercator';...
    'Miller Cylindrical';'Equidistant Cylindrical';'Oblique Mercator';...
    'Transverse Mercator';'Sinusoidal';'Gall-Peters';'Hammer-Aitoff';...
    'Mollweide';'Robinson';'UTM';};

proj=11;
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


figure
hold on
v=[-2000,-2000]; % this sets level at -2000m isobath 
[CS,H]=m_tbase('contour',v,'edgecolor','k','linewidth',0.8);
% plot region box and SR/OO line
m_plot(andrex_box(:,1),andrex_box(:,2),'color',[.75 .0 .0],'linewi',1.5,'LineStyle','--')
m_plot(CONTOUR1980_2020(:,1),CONTOUR1980_2020(:,2),'color',[1	0.635 0],'linewi',1.5,'LineStyle','--') % or [0.8 0.4 0.8] or 0.8 0.4 0.2
m_plot(CONTOUR2000(:,1),CONTOUR2000(:,2),'color',[1	0.635 0],'linewi',1.5,'LineStyle','--') % or [0.8 0.4 0.8] or 0.8 0.4 0.2

m_plot(test(:,2),test(:,1),'color',[1	0.635 0],'linewi',1.5,'LineStyle','--') % or [0.8 0.4 0.8] or 0.8 0.4 0.2


m_coast('patch',[0.7 0.7 0.7]);
m_grid('box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);
% caxis([0 500]);
set(gca,'FontSize',12)

get(0,'ScreenSize')


%% test test
test=test';
LAT_MIN = -80; LAT_MAX = -60; LONG_MIN = -75; LONG_MAX = 40;

[Z,LONG,LAT]=m_tbase([LONG_MIN LONG_MAX LAT_MIN LAT_MAX]);
find2000Z=find(Z<=-1910&Z>=-2010);
LONG2000=LONG(find2000Z);
LAT2000=LAT(find2000Z);

figure
hold on
m_plot(andrex_box(:,1),andrex_box(:,2),'color',[.75 .0 .0],'linewi',1.5,'LineStyle','--')
m_plot(LONG2000,LAT2000,'color',[1	0.635 0],'linewi',1.5,'LineStyle','--') % or [0.8 0.4 0.8] or 0.8 0.4 0.2


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LAT_MIN = -80; LAT_MAX = -60; LONG_MIN = -75; LONG_MAX = 40;
[Z,LONG,LAT]=m_tbase([LONG_MIN LONG_MAX LAT_MIN LAT_MAX]);
IN_AND_elev=inpolygon(LONG,LAT,andrex_box(:,1),andrex_box(:,2));
findelevation=find(IN_AND_elev==1);
LONG(IN_AND_elev==0)=NaN;
LAT(IN_AND_elev==0)=NaN;
Z(IN_AND_elev==0)=NaN;

find2000Z=find(Z<=-1990&Z>=-2010);
LONG2000=LONG(find2000Z);
LAT2000=LAT(find2000Z);

figure
hold on
% v=[-2000,-2000]; % this sets level at -2000m isobath 
% [CS,H]=m_tbase('contour',v,'edgecolor','k','linewidth',0.8);
plot(andrex_box(:,1),andrex_box(:,2),'color',[.75 .0 .0],'linewi',1.5,'LineStyle','--')
plot(Just2000mContour2(:,1),Just2000mContour2(:,2),'color',[0.5 .60 .30],'linewi',1.5,'LineStyle','--')
plot(LONG2000,LAT2000,'color',[1	0.635 0],'linewi',1.5,'LineStyle','--') % or [0.8 0.4 0.8] or 0.8 0.4 0.2
geoshow('landareas.shp','facecolor',[0.7 0.7 0.7])
xlim([-75 40])
ylim([-80 -54])
grid on

islands=find(LONG2000>-48 & LAT2000>-63);
LONG2000(islands)=[];
LAT2000(islands)=[];


figure
hold on
v=[-2000,-2000]; % this sets level at -2000m isobath 
[CS,H]=m_tbase('contour',v,'edgecolor','k','linewidth',3);
m_plot(andrex_box(:,1),andrex_box(:,2),'color',[.75 .0 .0],'linewi',1.5,'LineStyle','--')
% m_plot(Just2000mContour2(:,1),Just2000mContour2(:,2),'color',[0.5 .60 .30],'linewi',2,'LineStyle','--')
m_plot(Just2000mContourcleaned21(:,1),Just2000mContourcleaned21(:,2),'color',[0.5 .60 .30],'linewi',1.5,'LineStyle','--')

m_plot(LONG2000,LAT2000,'color',[1	0.635 0],'linewi',1,'LineStyle','-') % or [0.8 0.4 0.8] or 0.8 0.4 0.2
% geoshow('landareas.shp','facecolor',[0.7 0.7 0.7])
m_coast('patch',[0.7 0.7 0.7]);
m_grid('xtick',15,'box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);
% caxis([0 500]);
set(gca,'FontSize',12)

get(0,'ScreenSize')

% %%%%%%%%%%

IN_and=inpolygon(bio_ix_table.longitude,bio_ix_table.latitude,andrex_box(:,1),andrex_box(:,2));
bio_ix_table = bio_ix_table(IN_and,:); % further reduce the profiles to just those within the ANDREX box limits






