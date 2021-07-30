clearvars
desktop = 1; 
if desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg'))
    cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Matlab Add-ins'));
else
    cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
end

load('vgpm_imported.mat', 'VGPM_npp_tot_gC_all')
load('ProcessedData.mat', 'timedec')

%% just 7 months of full satellite coverage (can alter for all 12 months)
months={'Sept';'Oct';'Nov';'Dec';'Jan';'Feb';'Mar'}
test={'a':'q'}
yearno={'a'}
for yix = 2003%:2019
    findyear=find(timedec>yix-0.3 & timedec<yix+0.25); % removes April-August (partial/no coverage)
    
    map_7mototal_NPP(:,:,yix-2002)=nansum(VGPM_npp_tot_gC_all(:,:,findyear),3);
end
map_7mototal_NPP(map_7mototal_NPP==0)=NaN;
for yix = 2003%:2019
    findyear=find(timedec>yix-0.3 & timedec<yix+0.25); % removes April-August (partial/no coverage)

    yearmaps.(yearno{yix-2002})=(VGPM_npp_tot_gC_all(:,:,findyear));
               
end


load('latlon_m.mat')
load('NPP_cumulativeCMAP.mat')

yearmaps.n1=(yearmaps.a(:,:,1))./map_7mototal_NPP(:,:,1);
figure; pcolor(lon_m,lat_m,yearmaps.n1); shading flat; hold on; colorbar
colormap(cmblueice)
title('September 2003')
xlim([-65,47]);
ylim([-80,-50]);
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')

yearmaps.n2=(yearmaps.a(:,:,1)+yearmaps.a(:,:,2))./map_7mototal_NPP(:,:,1);
figure; pcolor(lon_m,lat_m,yearmaps.n2); shading flat; hold on; colorbar
colormap(cmblueice)
title('Sept+Oct 2003')
xlim([-65,47]);
ylim([-80,-50]);
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')

yearmaps.n3=(yearmaps.a(:,:,1)+yearmaps.a(:,:,2)+yearmaps.a(:,:,3))./map_7mototal_NPP(:,:,1);
figure; pcolor(lon_m,lat_m,yearmaps.n3); shading flat; hold on; colorbar
colormap(cmblueice)
title('Sept+Oct+Nov 2003')
xlim([-65,47]);
ylim([-80,-50]);
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')


yearmaps.n4=(yearmaps.a(:,:,1)+yearmaps.a(:,:,2)+yearmaps.a(:,:,3)+yearmaps.a(:,:,4))./map_7mototal_NPP(:,:,1);
figure; pcolor(lon_m,lat_m,yearmaps.n4); shading flat; hold on; colorbar
colormap(cmblueice)
title('Sept+Oct+Nov+Dec 2003')
xlim([-65,47]);
ylim([-80,-50]);
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')

yearmaps.n5=(yearmaps.a(:,:,1)+yearmaps.a(:,:,2)+yearmaps.a(:,:,3)+yearmaps.a(:,:,4)+yearmaps.a(:,:,5))./map_7mototal_NPP(:,:,1);
figure; pcolor(lon_m,lat_m,yearmaps.n5); shading flat; hold on; colorbar
colormap(cmblueice)
title('Sept+Oct+Nov+Dec 2003 +Jan 2004')
xlim([-65,47]);
ylim([-80,-50]);
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')


yearmaps.n6=(yearmaps.a(:,:,1)+yearmaps.a(:,:,2)+yearmaps.a(:,:,3)+ ...
    yearmaps.a(:,:,4)+yearmaps.a(:,:,5)+yearmaps.a(:,:,6))./map_7mototal_NPP(:,:,1);
figure; pcolor(lon_m,lat_m,yearmaps.n6); shading flat; hold on; colorbar
colormap(cmblueice)
title('Sept+Oct+Nov+Dec 2003 +Jan+Feb 2004')
xlim([-65,47]);
ylim([-80,-50]);
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')

yearmaps.n7=(yearmaps.a(:,:,1)+yearmaps.a(:,:,2)+yearmaps.a(:,:,3)+ ...
    yearmaps.a(:,:,4)+yearmaps.a(:,:,5)+yearmaps.a(:,:,6)+yearmaps.a(:,:,7))./map_7mototal_NPP(:,:,1);
figure; pcolor(lon_m,lat_m,yearmaps.n7); shading flat; hold on; colorbar
% cmap=flipud(hot(40));
colormap(cmblueice)
caxis([0 1])
title('Sept+Oct+Nov+Dec 2003 +Jan+Feb+Mar 2004')
xlim([-65,47]);
ylim([-80,-50]);
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')


