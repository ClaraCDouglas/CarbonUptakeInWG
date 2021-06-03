%% Load data
% Directory
    % C:\Users\Clara\Documents\Graduate Life\Southampton\Project\USB 0710\Project 0908
% load('latlon_m.mat')
% lat_avCUT=lat_m(721:1080,570:1411);
% lon_avCUT=lon_m(721:1080,570:1411);
% 
% %try to load just the section i want instead of loading whole variable
%     %load('VGPM_data.mat', 'VGPM_npp_all')
% 
% exampleObject = matfile('VGPM_data.mat');
% varlist = who(exampleObject)
% [nrows,ncols,ndim] = size(exampleObject,'VGPM_npp_all');
% VGPM_npp_all_CUT= exampleObject.VGPM_npp_all(721:1080,570:1411,:);
% 
% clear lat_m lon_m ncols ndim nrows varlist
% 
% save AverageVGPM_CUT.mat lat_avCUT lon_avCUT VGPM_npp_all_CUT -v7.3;

load('AverageVGPM_CUT.mat');
%% Average daily NPP


temp_vgpm=VGPM_npp_all_CUT;
% where NPP in ice covered areas = 0
findneg=find(temp_vgpm<0);
temp_vgpm(findneg)=0;
VGPM_av_day=nanmean(temp_vgpm,3);
%quick check
% figure; pcolor(lon_avCUT,lat_avCUT,VGPM_av_day); shading flat;
figure; pcolor(lon_avCUT(121:360,:),lat_avCUT(121:360,:),VGPM_av_day(121:360,:)); shading flat;

%setting places permanently 0 to NaN to represent ice/land
findzero=find(VGPM_av_day==0);
VGPM_av_day_no=VGPM_av_day;
VGPM_av_day_no(findzero)=NaN;
%quick check
figure; pcolor(lon_avCUT(121:360,:),lat_avCUT(121:360,:),VGPM_av_day_no(121:360,:)); shading flat;


% where NPP in ice covered areas = NaN, so not included:
    % Average NPP during ice free conditions
temp_vgpm=VGPM_npp_all_CUT;
findneg=find(temp_vgpm<0);
temp_vgpm(findneg)=NaN;
VGPM_av_day_nan=nanmean(temp_vgpm,3);
figure; pcolor(lon_avCUT(121:360,:),lat_avCUT(121:360,:),VGPM_av_day_nan(121:360,:)); shading flat;

%% load the files needed for next step
load('box_lat_lons.mat', 'andrex_box')
load('Just2000mContour2.mat')
load('BlueColourmaps.mat')
load('SouthernBoundary_Orsi1995.mat')
% Using 'VGPM_av_day_nan' for plots for now --> 'nin'

%% Plot the map
LAT_MIN = -80; LAT_MAX = -50; LONG_MIN = -75; LONG_MAX = 40;

xin=lon_avCUT; yin=lat_avCUT; zin=VGPM_av_day_no; nin=VGPM_av_day_nan; 
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


figure
m_pcolor(xin,yin,nin) %zin %nin
shading interp
set(0,'DefaultAxesColor','none')
caxis([0 350]); %250 %380 %300
colormap(summer)
% [cmap] = setcmappete('rest',482,'norm'); %/456/104
% [cmap] = setcmappete('cpt',185,'norm');
% [cmap] = setcmappete('cpt',182,'norm'); % but can't have the white at the top and representing NaN/0...

    %saving the 182 with grey at the top
    % cmapRAINBOW_blacktogrey=colormap(gca);
    % save cmap_Rainbow.mat cmapRAINBOW_blacktogrey -v7.3;

%     load('cmap_Rainbow.mat')
% colormap(cmapRAINBOW_blacktogrey)
% colormap(flipud(ColourmapC))
% colormap(flipud(ColourmapD))
% colormap(flipud(ColourmapE))
% ColourmapF_flip=flipud(ColourmapF);
% colormap((ColourmapF_flip)) % for nin and medin
% ColourmapF_Blues=colormap(gca);
colormap((ColourmapF_Blues)) % for nin and medin

% ColourmapF2=ColourmapF(5:45,:); 
% ColourmapF2=ColourmapF([2:1:2 5:2:21 23:1:47],:); %6:2:10 11:1:44
% colormap((ColourmapF2)) % for zin

hold on
% [C,h]=m_contour(xin,yin,zin,[50:50:100],'w');
% [C2,h2]=m_contour(xin,yin,zin,[200:100:400],'w');
% clabel(C,h,'Color','w','FontSize',12)
% clabel(C2,h2,'Color','w','FontSize',12)
% [C,h]=m_contourf(xin,yin,zin,[0:50:150 200:100:600 1000],'k');
% clabel(C,h,'Color','k','FontSize',13)

hh=colorbar;
ylabel(hh,'mg C m^-^2 day^-^1','FontSize',12);
 hh.Location='southoutside'; % want to make it a bit shorter still
hh.Position = [0.3 0.07 0.4 0.05];
hh.Position = [0.3 0.1 0.4 0.05];

% plot 2000m isobath
       % To draw the contours at one height (k), 
        % specify levels as a two-element row vector [k k]
v=[-2000,-2000]; % this sets level at -2000m isobath 
m_tbase('contour',v,'edgecolor','k','linewidth',0.8);
% plot SB front
% for nin: --> m_plot(sbdy(:,1),sbdy(:,2),'color',[.7 .7 .7],'linewi',1.25,'LineStyle','-')
m_plot(sbdy(:,1),sbdy(:,2),'color',[.8 .8 .8],'linewi',1.25,'LineStyle','-')
% plot region box and SR/OO line
m_plot(andrex_box(:,1),andrex_box(:,2),'color',[.75 .0 .0],'linewi',1.5,'LineStyle','--')
m_plot(Just2000mContour2(:,1),Just2000mContour2(:,2),'color',[1	0.635 0],'linewi',1.5,'LineStyle','--') % or [0.8 0.4 0.8] or 0.8 0.4 0.2

%  [C,h]=m_contour(xin,yin,zin,[80 80],'b','LineStyle','-','linewi',1.2); % 2.5e9

% hold on
% m_tbase('contour',[-6000:1000:000]);
%      colormap(gca,[a]); 


m_coast('patch',[0.7 0.7 0.7]);
m_grid('box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);
% caxis([0 500]);
set(gca,'FontSize',12)

get(0,'ScreenSize')
set(gcf,'color','w','position',[1000 80 850 620])

set(gca,'color','none')

title('Average Net Primary Production [MODIS VGPM] 06/2002 to 02/2019')
% title('Average Net Primary Production [MODIS VGPM] 20080 zeros')
% title({'Average (Mean) Net Primary Production [MODIS VGPM] 06/2002 to 02/2019';'Average daily rate where NPP is considered 0 when ice covered'})
% title({'Average (Mean) Net Primary Production [MODIS VGPM] 06/2002 to 02/2019';'Average daily rate where NPP is classed as NaN' ; '(i.e. omitted when calculating mean) when ice covered'})
% title({'Average (Median) Net Primary Production [MODIS VGPM] 06/2002 to 02/2019';'Median daily rate where NPP is classed as NaN' ; '(i.e. omitted when calculating median) when ice covered' ; '200mgCm^-^2 day^-^1 white contour'})
% title('') %title({'First line';'Second line'})
set(get(gca,'title'),'Position',[-0.0169,0.25,-5000000000000000])

%% Log NPP

VGPMlog_av_day_no=log(VGPM_av_day_no);

lin=VGPMlog_av_day_no;

figure
m_pcolor(xin,yin,zin)
shading interp
m_coast('patch',[0.7 0.7 0.7]);
m_grid('box','on','tickdir','in','xaxisLocation','top', 'fontsize',12);

% colormap(summer)
% [cmap] = setcmappete('rest',482,'norm'); %/456/104
% [cmap] = setcmappete('cpt',185,'norm');
% [cmap] = setcmappete('cpt',182,'norm'); % but can't have the white at the top and representing NaN/0...

    %saving the 182 with grey at the top
    % cmapRAINBOW_blacktogrey=colormap(gca);
    % save cmap_Rainbow.mat cmapRAINBOW_blacktogrey -v7.3;

%     load('cmap_Rainbow.mat')
colormap(cmapRAINBOW_blacktogrey)
colormap(parula)
set(gca,'colorscale','log')

caxis([0 7]);

[cmap] = setcmappete('rest',270,'norm'); %448 maybe rev
%425 norm okay
%
hold on




%% Average monthly total NPP
 exampleObject = matfile('VGPM_data.mat');
varlist = who(exampleObject)
[nrows,ncols,ndim] = size(exampleObject,'VGPM_npp_tot_gC');
VGPM_npp_tot_gC_CUT= exampleObject.VGPM_npp_tot_gC(721:1080,570:1411,:);
VGPM_npp_tot_gC_NaNs_CUT= exampleObject.VGPM_npp_tot_gC_NaNs(721:1080,570:1411,:);

test=VGPM_npp_tot_gC_CUT(:,:,1);
test2=VGPM_npp_tot_gC_NaNs_CUT(:,:,1);

% again, use VGPM_npp_tot_gC_CUT because there may be places where there is
% 0 one year/month and some the next, so need to include the zeros

VGPM_av_tot=nanmean(VGPM_npp_tot_gC_CUT,3);
%quick check
figure; pcolor(lon_avCUT,lat_avCUT,VGPM_av_tot); shading flat;

% so that is the average total NPP per month
    % bet it looks much better calculated with nans hahaha
