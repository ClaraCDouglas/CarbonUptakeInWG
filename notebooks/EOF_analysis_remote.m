%% EOF analysis
% process data to only have within WG box
% NPP_years.anom_annual_day_nan = NPP_years.vgpm_annual_day_nan - NPP.vgpm_annual_av_day_nan;

NPP_years.(algorithm{aix}).annual_day_nan_BOX=NPP_years.(algorithm{aix}).annual_day_nan.*temp.Weddell.box_logic;
NPP_years.(algorithm{aix}).annual_day_nan_BOX(NPP_years.(algorithm{aix}).annual_day_nan_BOX==0)=NaN;

% cut down to box around WG box without setting outside as 0
NPP_years.(algorithm{aix}).annual_day_nan_BOX=NPP_years.vgpm_annual_day_nan(850:1050,690:1300,:);
% to check region:
     lat_m_box=lat_m(850:1050,690:1300);
     lon_m_box=lon_m(850:1050,690:1300);
    % figure;
    % pcolor(lon_m_box,lat_m_box,NPP_years.vgpm_annual_day_nan_BOX(:,:,1)); shading flat
years=2003:1:2019;
years=years';

time=1:length(2003:2019);
time=time';
[ev_index,tda,pev,trends] = calc_pigup_EOF2(NPP_years.(algorithm{aix}).annual_day_nan_BOX,time)

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









