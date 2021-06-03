% findweddell = index of locations within WG box in 1080x2160 array
% bin = bin number for each cell corresponding to the list in findweddell

% list values in findweddell where bin==ix --> findweddell(bin==ix)
% make array of zeros
    % array(listedvales)==ix
%
curious=findweddell(bin==1);
curiouser=zeros(1080,2160);
curiouser(curious)=1;


% loop
curiouser=zeros(1080,2160);

for ix=1:1:7
    curious=findweddell(bin==ix);
    curiouser(curious)=ix;
    clear curious
end
curiouser(find(curiouser==0))=NaN;

map = [0.2 0 .3
    0.2 0.1 0.5
    0.1 0.5 0.8
    0.2 0.7 0.6
    0.8 0.7 0.3
    0.9 1 0
    0 0 0];

figure
pcolor(lon_m,lat_m,curiouser) % or pcolor
shading flat
colormap(map)



curiouser=zeros(1080,2160);

for ix=1:length(Yedges)
    curious=findweddell(binY==ix);
    curiouser(curious)=ix;
    clear curious
end

map2 = [0 0 0
    0.5 0.75 1
    0 0.5 0];
curiouser(find(curiouser==0))=NaN;
figure
pcolor(lon_m,lat_m,curiouser) % or pcolor
shading flat
colormap(map2)



curiouser=zeros(1080,2160);

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
        curiouser(curious)=number;
        clear curious
    end
end
curiouser(find(curiouser==0))=NaN;


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

figure
pcolor(lon_m,lat_m,curiouser) % or pcolor
shading flat
colormap(cmap_bioregions)
    
    
    
    
LAT_MIN = -80; LAT_MAX = -50; LONG_MIN = -75; LONG_MAX = 40;
xin=lon_m; yin=lat_m; zin=curiouser; 
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
hh.Label.Rotation = 0;
hh.Label.Position = [0.7 10.65 0.7];
