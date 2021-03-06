% findweddell = index of locations within WG box in 1080x2160 array
% bin = bin number for each cell corresponding to the list in findweddell

% list values in findweddell where bin==ix --> findweddell(bin==ix)
% make array of zeros
    % array(listedvales)==ix
% test

curious=findweddell(bin==1);
bioregionbin=zeros(1080,2160);
bioregionbin(curious)=1;


%% Bathymetry
bathybins=zeros(1080,2160);

for ix=1:1:7
    curious=findweddell(bin==ix);
    bathybins(curious)=ix;
    clear curious
end
bathybins(find(bathybins==0))=NaN;

map = [0.2 0 .3
    0.2 0.1 0.5
    0.1 0.5 0.8
    0.2 0.7 0.6
    0.8 0.7 0.3
    0.9 1 0
    0 0 0];
colors = distinguishable_colors(21,{'w','k','b',[[0,0.266666666666667,0.105882352941176]]});

B=figure('units','normalized','outerposition',[0 0 1 1])
B(1)=pcolor(lon_m,lat_m,bathybins); shading flat; hold on
colormap(map)
hh=colorbar
caxis([1 7]);
hh.Ticks=[1.45:0.85:7];
hh.TickLabels=[-5000:1000:1000];
hh.Label.String = 'Depth bin (value is upper bin limit)';
hh.Label.Rotation = 90;

xlim([-65,47]);
ylim([-80,-50])
title('Bathymetry bins')


B(2)=plot(coast.long,coast.lat,'k'); hold on

for fix = 1:length(gfixs)
    float_name = float_names{fix};
    if ~isempty(biofloat.(float_name).LONGITUDE)
        B(fix+2) = plot(biofloat.(float_name).LONGITUDE,biofloat.(float_name).LATITUDE,'linewidth',3,'color',colors(fix,:));
        hold on
    end
    
end
legend(B(3:end),float_names(gfixs))
Bax=gca
Bax.FontSize=14



%% NPP rates area-weighted
prodbins=zeros(1080,2160);
for ix=1:length(Yedges)
    curious=findweddell(binY==ix);
    prodbins(curious)=ix;
    clear curious
end

map2 = [0 0 0
    0.5 0.75 1
    0 0.5 0];
prodbins(find(prodbins==0))=NaN;

P=figure('units','normalized','outerposition',[0 0 1 1])
P(1)=pcolor(lon_m,lat_m,prodbins); shading flat; hold on
colormap(map2)
pp=colorbar
caxis([1 3]);
pp.Ticks=[1.45:0.65:3];
pp.TickLabels={'None', '0-200','>200'};
pp.Label.String = 'Area-weighted productivity bins (mg m^-^2 day^-^1)';
pp.Label.Rotation = 90;
xlim([-65,47]);
ylim([-80,-50])

P(2)=plot(coast.long,coast.lat,'k'); hold on

for fix = 1:length(gfixs)
    float_name = float_names{fix};
    if ~isempty(biofloat.(float_name).LONGITUDE)
        P(fix+2) = plot(biofloat.(float_name).LONGITUDE,biofloat.(float_name).LATITUDE,'linewidth',3,'color',colors(fix,:));
        hold on
    end
    
end
legend(P(3:end),float_names(gfixs))
Pax=gca
Pax.FontSize=14
title('Productivity bins')

%%
bioregionbin=zeros(1080,2160);

for ix=1:length(Xedges)
    for nix = 1:length(Yedges)
         if ix==1 & nix==1
             number=1;
         elseif ix==2 & nix==1
             number=1;
         elseif ix==3 & nix==1
             number=1;
         elseif ix==4 & nix==1
             number=1;
         elseif ix==5 & nix==1
             number=1;
         elseif ix==6 & nix==1
             number=1;
         elseif ix==7 & nix==1
             number=1;
         elseif ix==1 & nix==2
             number=2;
         elseif ix==2 & nix==2
             number=3;
         elseif ix==3 & nix==2
             number=4;
         elseif ix==4 & nix==2
             number=5;
         elseif ix==5 & nix==2
             number=6;
         elseif ix==6 & nix==2
             number=7;
         elseif ix==7 & nix==2
             number=14;
         elseif ix==1 && nix==3
             number=8;
         elseif ix==2 && nix==3
             number=9;
         elseif ix==3 && nix==3
             number=10;
         elseif ix==4 && nix==3
             number=11;
         elseif ix==5 && nix==3
             number=12;
         elseif ix==6 && nix==3
             number=13;
         elseif ix==7 && nix==3
             number=14;
         end  
        curious=findweddell(binX == ix & binY == nix);
        bioregionbin(curious)=number;
        clear curious
    end
end
bioregionbin(find(bioregionbin==0))=NaN;

% cmap needs 14 levels
    % black for no PP (#)
    % blue gradient for low PP (#2:7)
    % green gradient for high PP (#8:13)
    % and brown for above sea level (#14)
[cmaplow] = setcmappete('rest',49,'rev');
[cmaphigh] = setcmappete('rest',35,'rev');
cmap_bioregions=zeros(14,3);
cmap_bioregions(8:13,:)=cmaphigh(1:9:50,:)
cmap_bioregions(2:7,:)=cmaplow(1:9:50,:)
cmap_bioregions(14,:)=[0.59 0.29 0]

   
R=figure('units','normalized','outerposition',[0 0 1 1])
R(1)=pcolor(lon_m,lat_m,bioregionbin); shading flat; hold on
colormap(cmap_bioregions)
rr=colorbar
caxis([1 14]);
rr.Ticks=[1.45:0.95:14];
rr.TickLabels={'None', '>5km Low PP','4-5km Low PP','3-4km Low PP','2-3km Low PP','1-2km Low PP','0-1km Low PP',...
    '>5km High PP','4-5km High PP','3-4km High PP','2-3km High PP','1-2km High PP','0-1km High PP','Land'};
rr.Label.String = 'Bioregions';
rr.Label.Rotation = 90;
xlim([-65,47]);
ylim([-80,-50])

P(2)=plot(coast.long,coast.lat,'k'); hold on

for fix = 1:length(gfixs)
    float_name = float_names{fix};
    if ~isempty(biofloat.(float_name).LONGITUDE)
        P(fix+2) = plot(biofloat.(float_name).LONGITUDE,biofloat.(float_name).LATITUDE,'linewidth',3,'color',colors(fix,:));
        hold on
    end
    
end
legend(P(3:end),float_names(gfixs))
Pax=gca
Pax.FontSize=14
title('Bioregion bins and floats')    
    
%% m_map    
LAT_MIN = -80; LAT_MAX = -50; LONG_MIN = -75; LONG_MAX = 40;
xin=lon_m; yin=lat_m; zin=bioregionbin; 
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
% caxis([0 400]); %250 %380 %300
colormap(cmap_bioregions)
hold on
hh=colorbar;
m_coast('patch',[0.7 0.7 0.7]);
m_grid('box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);
caxis([1 14]);
hh.Ticks=[1.45:0.95:14];
hh.TickLabels=[1:1:14];
hh.Label.String = 'Bioregion';
hh.Label.Rotation = 90;
hh.Label.Position = [2 8 0.7];
