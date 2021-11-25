clearvars
desktop = 0;
if desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg'))
    cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
    addpath(genpath('C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Matlab Add-ins'));
else
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG'))
    cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\m_map'));
    addpath(genpath('C:\Users\ccd1n18\Documents\Projects\SetCMap'));
end

load('vgpm_imported.mat', 'VGPM_npp_tot_gC_all')
load('ProcessedData.mat', 'timedec')
load('latlon_m.mat')

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


%% Working out when 50% of NPP has occurred in each pixel

% months={'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';'Jan';'Feb';'Mar';'Apr';'May';'Jun'}
% test={'a':'q'}
% yearno={'a'}
for yix = 2003%:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    NPP_motot_2003=VGPM_npp_tot_gC_all(:,:,findyear); % cut down to study area
    %   NPP_motot_2003=VGPM_npp_tot_gC_all(850:1050,690:1300,findyear); % cut down to study area
end
% lat_ms=lat_m(850:1050,690:1300);
% lon_ms=lon_m(850:1050,690:1300);
NPP_tot_2003=nansum(NPP_motot_2003,3); % total NPP in each pixel for 2003
NPP_tot_2003(NPP_tot_2003==0)=NaN;

NPP_prop_2003=NPP_motot_2003./NPP_tot_2003; % monthly total NPP as a proportion of annual total

NPP_abscumtot_2003=cumsum(NPP_motot_2003,3); % cumulative total of absolute NPP values
NPP_propcum_2003=cumsum(NPP_prop_2003,3); % cumulative total of proportional values
% ^ (use this one to assess when 50% occurs)


% Want to replicate this on a whole region/world scale, and make it efficient
%     [a,b]=find(NPP_propcum_2003(70,400,:)>=0.5);
%     b_month=min(b);
%     [a2,b2]=find(NPP_propcum_2003(70,400,:)==0);
%     b2_month=max(b2);
%     c_speed=b_month-b2_month;
% for this pixel, it takes 2 months of being ice free for phyto to carry out >= 50% of annual total NPP


% getting month of >=50%
output=zeros(size(NPP_propcum_2003(:,:,:)));
output2=zeros(size(NPP_propcum_2003(:,:,1)));
[r,c,v]=ind2sub(size(NPP_propcum_2003(:,:,1)),(find(NPP_propcum_2003(:,:,:)>=0.5)));

for ix=max(v):-1:min(v)
    selectmonth=ix;
    search=find(v==selectmonth);
    for tix=1:length(search)
        output(r(search(tix)),c(search(tix)),ix)=selectmonth;
        output2(r(search(tix)),c(search(tix)))=selectmonth;
    end
end
% figure; pcolor(output(:,:,7)); shading flat
output2(output2==0)=NaN;
figure; pcolor(output2); shading flat;
caxis([5 9])


[cmap] = setcmappete('rest',480,'reg');
[cmap2] = setcmappete('rest',428,'reg');
[cmap3] = setcmappete('rest',475,'rev');
[cmap4] = setcmappete('rest',4,'rev');
[cmap5] = setcmappete('rest',4,'rev');


figure;
tiledlayout(3,2)
ax1 = nexttile;
pcolor(lon_m,lat_m,output2); shading flat; hold on; pp=colorbar;
xlim([-65,37]);
ylim([-80,-52]);
caxis([5 9])
cmap=cmap(1:round(length(cmap))/5:end,:);
colormap(ax1,cmap)
% colormap(cmap)
pp.XTick=[5:1:9];
pp.Ticks=[5.30:0.85:9];
pp.TickLabels={'Nov','Dec','Jan','Feb','Mar'};
pp.TickLength=0;
title('Month where cumulative total NPP reaches \geq 50% - 2003')
% SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
% OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')

% getting month of ice
outputice=zeros(size(NPP_propcum_2003(:,:,:)));
outputice2=zeros(size(NPP_propcum_2003(:,:,1)));
[ri,ci,vi]=ind2sub(size(NPP_propcum_2003(:,:,1)),(find(NPP_propcum_2003(:,:,:)==0)));

for ix=min(vi):1:max(vi)
    selectmonth=ix;
    search=find(vi==selectmonth);
    for tix=1:length(search)
        outputice(ri(search(tix)),ci(search(tix)),ix)=selectmonth;
        outputice2(ri(search(tix)),ci(search(tix)))=selectmonth;
    end
end
figure; pcolor(outputice(:,:,7)); shading flat
outputice2(outputice2==0)=NaN;
figure; pcolor(outputice2); shading flat;

% figure;
ax2 = nexttile;
pcolor(lon_m,lat_m,outputice2); shading flat; hold on; qq=colorbar;
caxis([1 8])
cmap2=cmap2(1:round(length(cmap2))/8:end,:);
colormap(ax2,cmap2)
qq.XTick=[1:1:8];
qq.Ticks=[1.40:0.88:8];
qq.TickLabels={'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'};
qq.TickLength=0;
xlim([-65,37]);
ylim([-80,-52]);
title('Last month of ice coverage/low light conditions')
% SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
% OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')

% all together - number of months to get to 50%
months2grow=output2-outputice2;
figure; pcolor(months2grow); shading flat;

% figure;
ax3 = nexttile;
pcolor(lon_m,lat_m,months2grow); shading flat; hold on; tt=colorbar;
cmap3=cmap3(1:round(length(cmap3))/7:end,:);
colormap(ax3,cmap3)
colormap(cmap3)
tt.XTick=[1:1:7];
tt.Ticks=[1.30:0.85:7];
tt.TickLabels={1:1:7};
xlim([-65,37]);
ylim([-80,-52]);
title({'Number of months to reach 50% of total annual NPP - 2003';'(i.e. Number of months between ice coverage/low light and 50% of total NPP)'})
% SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
% OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')



% total NPP per year map
load('Colormap_Delta_NPP.mat')
figure
pcolor(lon_m,lat_m,NPP_tot_2003); shading flat; hold on; gg=colorbar;
colormap(ColourmapDelta_NPP)
caxis([0 10^10])
xlim([-65,37]);
ylim([-80,-52]);
title('Total annual NPP (g C) - 2003')
% SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
% OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')




% ice regrowth/low light returns
% use NPP_motot_2003, but only last half to get first time 0 appears in pixel
NPP_motot_2003_nans=NPP_motot_2003;
for ixx= 1:1:12
    temp=NPP_motot_2003(:,:,ixx);
    temp(isnan(NPP_tot_2003))=NaN;
    NPP_motot_2003_nans(:,:,ixx)=temp;
    clearvars temp
end
outputhalf=zeros(size(NPP_motot_2003_nans(:,:,9:12)));
outputhalf2=zeros(size(NPP_motot_2003_nans(:,:,1)));
half_2003=NPP_motot_2003_nans(:,:,9:12); % so 1 = Mar, and so on, to 4 = June
[rh,ch,vh]=ind2sub(size(half_2003(:,:,1)),(find(half_2003(:,:,:)==0)));

% getting month of
for ix=max(vh):-1:min(vh)
    selectmonth=ix
    search=find(vh==selectmonth);
    for tix=1:length(search)
        outputhalf(rh(search(tix)),ch(search(tix)),ix)=selectmonth;
        outputhalf2(rh(search(tix)),ch(search(tix)))=selectmonth;
    end
end
% figure; pcolor(outputhalf(:,:,4)); shading flat
outputhalf2(outputhalf2==0)=NaN;
figure;
pcolor(lon_m,lat_m,outputhalf2); shading flat; hold on;
temp=flipud(winter);
shortwinter=temp(1:round(length(temp))/4:end,:);
colormap(shortwinter)
ww=colorbar
ww.XTick=[1:1:4];
ww.Ticks=[1.4:0.8:4];
ww.TickLabels={'Mar','Apr','May','June'};
ww.TickLength=0;

% calculate growing season

endofGS=outputhalf2+8;
% figure;
% pcolor(lon_m,lat_m,endofGS); shading flat; hold on;

GrowingSeason_2003=endofGS-outputice2;
figure;
pcolor(lon_m,lat_m,GrowingSeason_2003); shading flat; hold on;
temp=parula;
GScmap=temp(1:round(length(temp))/10:end,:);
colormap(GScmap);
ss=colorbar;
ss.XTick=[2:1:11];
ss.Ticks=[2.4:0.9:11];
ss.TickLabels={2:1:11};
ss.TickLength=0;
xlim([-65,37]);
ylim([-80,-52]);
title('Length of growing season - 2003 (months)')


% getting month of 100%
output100=zeros(size(NPP_propcum_2003(:,:,:)));
output100_2=zeros(size(NPP_propcum_2003(:,:,1)));
[r100,c100,v100]=ind2sub(size(NPP_propcum_2003(:,:,1)),(find(NPP_propcum_2003(:,:,:)>=0.99)));

for ix=max(v100):-1:min(v100)
    selectmonth=ix
    search=find(v100==selectmonth);
    for tix=1:length(search)
        output100(r100(search(tix)),c100(search(tix)),ix)=selectmonth;
        output100_2(r100(search(tix)),c100(search(tix)))=selectmonth;
    end
end
figure; pcolor(output100(:,:,7)); shading flat
output100_2(output100_2==0)=NaN;
figure; pcolor(output100_2); shading flat;
caxis([5 9])

figure;
pcolor(lon_m,lat_m,output100_2); shading flat; hold on; hh=colorbar;
caxis([7 11])
cmap=cmap(1:round(length(cmap))/5:end,:);
colormap(cmap)
hh.XTick=[7:1:11];
hh.Ticks=[7.30:0.85:11];
hh.TickLabels={'Jan','Feb','Mar','Apr','May'};
hh.TickLength=0;
xlim([-65,37]);
ylim([-80,-52]);
title('Month where cumulative total NPP reaches \geq 99% - 2003')

%% making mega plot
load('Colormap_Delta_NPP.mat')
load('openshelf_coord.mat')
figure;
tiledlayout(2,3)
% total NPP
ax1 = nexttile;
pcolor(lon_m,lat_m,NPP_tot_2003); shading flat; hold on; gg=colorbar;
colormap(ax1,ColourmapDelta_NPP)
caxis([0 10^10])
xlim([-65,37]);
ylim([-80,-52]);
title('Total annual NPP (g C) - 2003')
SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')
% end of ice/winter season
ax2 = nexttile;
pcolor(lon_m,lat_m,outputice2); shading flat; hold on; qq=colorbar;
caxis([1 8])
cmap2=cmap2(1:round(length(cmap2))/8:end,:);
colormap(ax2,cmap2)
qq.XTick=[1:1:8];
qq.Ticks=[1.40:0.88:8];
qq.TickLabels={'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'};
qq.TickLength=0;
xlim([-65,37]);
ylim([-80,-52]);
title('Last month of ice coverage/low light conditions')
% SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
% OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')
% length of growing season
ax3 = nexttile;
pcolor(lon_m,lat_m,GrowingSeason_2003); shading flat; hold on; ss=colorbar;
xlim([-65,37]);
ylim([-80,-52]);
caxis([2 11])
temp=parula;
GScmap=temp(1:round(length(temp))/10:end,:);
colormap(ax3,GScmap);
ss.XTick=[2:1:11];
ss.Ticks=[2.4:0.9:11];
ss.TickLabels={2:1:11};
ss.TickLength=0;
title('Length of growing season - 2003 (months)')
% SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
% OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')
% number of months to reach 50% of annual total NPP
ax4 = nexttile;
pcolor(lon_m,lat_m,months2grow); shading flat; hold on; tt=colorbar;
xlim([-65,37]);
ylim([-80,-52]);
caxis([1 7])
cmap4=cmap4(1:round(length(cmap4))/7:end,:);
cmap4(7,:)=0.3;
colormap(ax4,cmap4)
tt.XTick=[1:1:7];
tt.Ticks=[1.30:0.85:7];
tt.TickLabels={1:1:7};
tt.TickLength=0;
title({'Number of months to reach 50% of total annual NPP - 2003';'(i.e. Number of months between ice coverage/low light and 50% of total NPP)'})
% SR_line=plot(shelf_region_ANDbox(:,1),shelf_region_ANDbox(:,2),'color',[0.8 0.4 0],'linewi',2);%'#80471C'
% OO_line=plot(open_ocean_ANDbox(:,1),open_ocean_ANDbox(:,2),'color',[0.6 0.2 0.8],'linewi',2);%,'LineStyle','--')

% when 50% is reached
ax5 = nexttile;
pcolor(lon_m,lat_m,output2); shading flat; hold on; pp=colorbar;
xlim([-65,37]);
ylim([-80,-52]);
caxis([5 9])
cmap=cmap(1:round(length(cmap))/5:end,:);
colormap(ax5,cmap)
pp.XTick=[5:1:9];
pp.Ticks=[5.30:0.85:9];
pp.TickLabels={'Nov','Dec','Jan','Feb','Mar'};
pp.TickLength=0;
title('Month where cumulative total NPP reaches \geq 50% - 2003')

%% Ratio between # mos and GS

Ratio_2003=months2grow./GrowingSeason_2003;

% when 50% is reached
ax6 = nexttile;
pcolor(lon_m,lat_m,Ratio_2003); shading flat; hold on; pp=colorbar;
xlim([-65,37]);
ylim([-80,-52]);
caxis([0 0.8])
temp=parula;
Rcmap=temp(1:round(length(temp))/40:end,:);
colormap(ax6,Rcmap);
pp.XTick=[0:0.1:0.8];
title({'Number of months to reach 50% / length of growing season';'2003'})



%% getting month of >=75%
output75=zeros(size(NPP_propcum_2003(:,:,:)));
output75_2=zeros(size(NPP_propcum_2003(:,:,1)));
[r,c,v]=ind2sub(size(NPP_propcum_2003(:,:,1)),(find(NPP_propcum_2003(:,:,:)>=0.75)));

for ix=max(v):-1:min(v)
    selectmonth=ix;
    search=find(v==selectmonth);
    for tix=1:length(search)
        output75(r(search(tix)),c(search(tix)),ix)=selectmonth;
        output75_2(r(search(tix)),c(search(tix)))=selectmonth;
    end
end

months2grow75=output75_2-outputice2;

figure;
pcolor(lon_m,lat_m,months2grow75); shading flat; hold on; tt=colorbar;
xlim([-65,37]);
ylim([-80,-52]);
caxis([1 9])
cmap5=cmap5(1:round(length(cmap5))/9:end,:);
% cmap4(7,:)=0.3;
colormap(cmap5)
tt.XTick=[1:1:7];
tt.Ticks=[1.30:0.85:7];
tt.TickLabels={1:1:7};
tt.TickLength=0;
title({'Number of months to reach 75% of total annual NPP - 2003';'(i.e. Number of months between ice coverage/low light and 75% of total NPP)'})

%% line plot of cumulative NPP
% run box logic (MAKE A FUNCTION)

for yix = 2003%:2019
    findyear=find(timedec>yix-0.5 & timedec<yix+0.5);
    NPP_motot_2003=VGPM_npp_tot_gC_all(:,:,findyear);
end
% lat_m=lat_m(850:1050,690:1300);
% lon_m=lon_m(850:1050,690:1300);
NPP_tot_2003=nansum(NPP_motot_2003,3); % total NPP in each pixel for 2003
NPP_tot_2003(NPP_tot_2003==0)=NaN;

NPP_prop_2003=NPP_motot_2003./NPP_tot_2003; % monthly total NPP as a proportion of annual total

NPP_abscumtot_2003=cumsum(NPP_motot_2003,3); % cumulative total of absolute NPP values
NPP_propcum_2003=cumsum(NPP_prop_2003,3); % cumulative total of proportional values
% ^ (use this one to assess when 50% occurs)
clearvars Results
for rix=1:length(region_sublist)
    for mix=1:1:12
        yemp=NPP_propcum_2003(:,:,mix);
        Results.(region_sublist{rix}).pixelcumsum(:,mix)=yemp(temp.(region_sublist{rix}).box_logic);
    end
    Results.(region_sublist{rix}).pixelcumsum=Results.(region_sublist{rix}).pixelcumsum';

    for yix=1:1:12
        Results.(region_sublist{rix}).avcumsum(yix,1)=mean(Results.(region_sublist{rix}).pixelcumsum(yix,:),'omitnan');
    end
end
figure;
for rix=1:length(region_sublist)
    subplot(1,3,rix)
    plot(Results.(region_sublist{rix}).pixelcumsum);%,'Color',[0.8 0.8 0.8]);
    hold on;
    txt=region_sublist(rix); title(txt);
    plot(Results.(region_sublist{rix}).avcumsum,'k','LineWidth',2)
end


figure; plot(Results.(region_sublist{2}).avcumsum,'Color',[0.8 0.4 0],'LineWidth',2); hold on;
plot(Results.(region_sublist{3}).avcumsum,'Color',[0.6 0.2 0.8],'LineWidth',2)

