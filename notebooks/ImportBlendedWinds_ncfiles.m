clearvars
cd E:\Data\Wind\Blended\UVData
% Select data folder:
filebase = 'E:\Data\Wind\Blended\UVData';
% filedir = append(filebase,'A');
D0 = dir(fullfile(filebase,'*.nc*'));

% test a file
ls; % copy a file to look at details
nc='uv20130112rt.nc';
A=ncinfo(nc);
A.Variables
A.Variables.Name
ncdisp(nc);
timetest=ncread(nc,'time')
attr_infostart = A.Attributes(11).Value
dt = datetime(1978,1,1)+hours(timetest)
% flag values  = -9999
windtest=ncread(nc,'w');

figure; pcolor(windtest); shading flat;

clearvars dt attr_infostart A timetest nc
%% load in all SIC data
D0 = dir(fullfile(filebase,'*.nc*'));
Winds_time_start=NaN(length(D0),3);
for ix=1:length(D0)
    b=struct2cell(D0);
    list=b{1,ix};
    %datetime
    timetemp=ncread(list,'time');
    Winds_time(ix) = datetime(1978,1,1)+hours(timetemp);
    %time in columns
    fileinfo=ncinfo(list);
    attr_infostart = fileinfo.Attributes(11).Value;
    Winds_time_start(ix,1) = str2double(attr_infostart(1:4));
    Winds_time_start(ix,2) = str2double(attr_infostart(6:7));
    Winds_time_start(ix,3) = str2double(attr_infostart(9:10));
end
Winds_timedec_start=decyear(Winds_time_start);

lat_wind=ncread(list,'lat');
lon_wind=ncread(list,'lon');
lat_north=0; % only importing south of equator
lat_south=-100;
lon_west=-10;
lon_east=370;
J=find(lat_wind>lat_south & lat_wind<lat_north);
K=find(lon_wind>lon_west & lon_wind<lon_east);
lat_windS=lat_wind(J);
lon_windS=lon_wind(K);
% lon_windSw=wrapTo180(lon_windS);
% lon_windy=cat(1,lon_windSw(721:end),lon_windSw(1:720));
clearvars lat_wind lon_wind *north *south *east *west

    Speed_med=NaN(359,1440,length([2003:1:2018]));
    Speed_mean=NaN(359,1440,length([2003:1:2018]));
    Speed_max=NaN(359,1440,length([2003:1:2018]));
    Speed_STD=NaN(359,1440,length([2003:1:2018]));
    U_med=NaN(359,1440,length([2003:1:2018]));
    U_mean=NaN(359,1440,length([2003:1:2018]));
    U_max=NaN(359,1440,length([2003:1:2018]));
    U_STD=NaN(359,1440,length([2003:1:2018]));
    V_med=NaN(359,1440,length([2003:1:2018]));
    V_mean=NaN(359,1440,length([2003:1:2018]));
    V_max=NaN(359,1440,length([2003:1:2018]));
    V_STD=NaN(359,1440,length([2003:1:2018]));
    Dir_med=NaN(359,1440,length([2003:1:2018]));
    Dir_mean=NaN(359,1440,length([2003:1:2018]));
    Dir_STD=NaN(359,1440,length([2003:1:2018]));

for yix = 2003%:2018
disp(yix)
findyear=find(Winds_timedec_start>yix-0.5 & Winds_timedec_start<yix+0.5);
Winds_speed=NaN(359,1440,length(findyear));
Winds_u=NaN(359,1440,length(findyear));
Winds_v=NaN(359,1440,length(findyear));
Winds_dir=NaN(359,1440,length(findyear));
for ix=findyear(1)%:findyear(end)
    minus=findyear(1)-1;
    list=b{1,ix};
    speed=ncread(list,'w');
    speed(speed==-9999)=NaN;
    speed=speed';
    C=speed(J,K);  % extract region
    Winds_speed(:,:,ix-minus)=C;
    
    u=ncread(list,'u');
    u(u==-9999)=NaN;
    u=u';
    Cu=u(J,K);
    Winds_u(:,:,ix-minus)=Cu;
    v=ncread(list,'v');
    v(u==-9999)=NaN;  % fill value is  >250 : check ncdisp(nc)
    v=v';
    Cv=v(J,K);
    Winds_v(:,:,ix-minus)=Cv;
    
    theta_gm=atan(Cu./Cv);
    
    [theta_c2p,rho]=cart2pol(Cu,Cv);
    lunch=(rad2deg(theta_c2p).*(-1))+90;
    findnegl=find(lunch<0);
    lunch(findnegl)=lunch(findnegl)+360;
    Winds_dir(:,:,ix-minus)=lunch;
end
    Speed_med(:,:,yix-2002)=median(Winds_speed,3,'omitnan');
    Speed_mean(:,:,yix-2002)=mean(Winds_speed,3,'omitnan');
    Speed_max(:,:,yix-2002)=max(Winds_speed,[],3,'omitnan');
    Speed_STD(:,:,yix-2002)=std(Winds_speed,0,3,'omitnan');
    U_med(:,:,yix-2002)=median(Winds_u,3,'omitnan');
    U_mean(:,:,yix-2002)=mean(Winds_u,3,'omitnan');
    U_max(:,:,yix-2002)=max(Winds_u,[],3,'omitnan');
    U_STD(:,:,yix-2002)=std(Winds_u,0,3,'omitnan');
    V_med(:,:,yix-2002)=median(Winds_v,3,'omitnan');
    V_mean(:,:,yix-2002)=mean(Winds_v,3,'omitnan');
    V_max(:,:,yix-2002)=max(Winds_v,[],3,'omitnan');
    V_STD(:,:,yix-2002)=std(Winds_v,0,3,'omitnan');
    Dir_med(:,:,yix-2002)=median(Winds_dir,3,'omitnan');
    Dir_mean(:,:,yix-2002)=mean(Winds_dir,3,'omitnan');
    Dir_STD(:,:,yix-2002)=std(Winds_dir,0,3,'omitnan');
end
clearvars C* v u speed minus list lunch r theta* rho find* Winds_speed Winds_u Winds_v Winds_dir
figure; set(0,'DefaultFigureWindowStyle','docked')
pcolor(lon_windS,lat_windS,(lunch)); shading flat; colorbar
colormap(hsv);

C_test=cat(2,C(:,721:end),C(:,1:720));
% clearvars -except '*comp' 'lat_w' 'lon_w' 'speed' 'D0' 'filebase' 'filedir' 'b'
figure; pcolor(lon_windS,lat_windS,(direction(:,:,1))); shading flat;
cmaptemp=parula;
cmaptemp=cat(1,cmocean('thermal'),cmocean('-solar'));
colormap(hsv); colorbar
% get(0,'ScreenSize')
% set(gcf,'color','w','position',[500 80 1000 800])

%
%% 8 day averages
SICv4_8day=[];
time_start_8day_SICv4=[];
year=Winds_time_start(:,1);
yearrange=unique(Winds_time_start(:,1));
for yix=1:length(yearrange)
    yearsel=yearrange(yix);
    yearice=Winds_speed(:,:,year==yearsel);
    yeartime=Winds_time_start(year==yearsel,:);
    for ix = 1:8:361
        if ix<360
        ice_conc_8day_temp=mean(yearice(:,:,(ix:1:ix+7)),3,'omitnan');
        SICv4_8day=cat(3,SICv4_8day, ice_conc_8day_temp);
        else
        ice_conc_8day_temp=mean(yearice(:,:,(ix:1:end)),3,'omitnan');
        SICv4_8day=cat(3,SICv4_8day, ice_conc_8day_temp);
        end
        time_8day_temp=yeartime(ix,:);
        time_start_8day_SICv4=cat(1,time_start_8day_SICv4,time_8day_temp);
        clearvars ice_conc_8day_temp time_8day_temp
    end
end

%% reducing to just the times I have NPP data for before saving
desktop = 0;
laptop=1;
if desktop
    cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\PhD\Projects\carbonuptakeinwg\data\processed' % desktop
    load('cafe_8day_imported_eqWG_withNaNs.mat', 'time_start_all')
elseif laptop
    cd 'C:\Users\ccd1n18\Documents\Projects\CarbonUptakeInWG\data\processed' % laptop
    load('cafe_8day_imported_eqWG_withNaNs.mat', 'time_start_all')
end

SICv4_8day(:,:,time_start_8day_SICv4(:,1)<=time_start_all(1,1)&time_start_8day_SICv4(:,2)<time_start_all(1,2),:)=[];
time_start_8day_SICv4(time_start_8day_SICv4(:,1)<=time_start_all(1,1)&time_start_8day_SICv4(:,2)<time_start_all(1,2),:)=[];
% remove the random 3 8-day slices where there isn't NPP data in 2020
time_start_8day_SICv4(834:836,:)=[];
SICv4_8day(:,:,834:836)=[];

%% final collating of data
cd E:\Data\SeaIceNIMBUS
load('SeaIce_8day_20022020.mat', 'g_area', 'g_lat', 'g_lon')

cd E:\Data\SeaIceV4\sidads.colorado.edu\AllFiles
save('SICv4_8day.mat','g_area','g_lat','g_lon','time_start_8day_SICv4','SICv4_8day');
save('SICv4_daily.mat','g_area','g_lat','g_lon','time_start_SICv4','SICv4_0221', '-v7.3');
