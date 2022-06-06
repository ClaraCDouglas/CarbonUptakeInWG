% load('cafe_8day_imported_eqWG_withNaNs.mat')
load('cafe_8day_imported_May22.mat','cafe_npp_all_8day')
load('cafe_8day_imported_May22.mat', 'lat_wg')
load('cafe_8day_imported_May22.mat', 'lon_wg')
load('cafe_8day_imported_May22.mat', 'time_end_all')
load('cafe_8day_imported_May22.mat', 'time_start_all')
lat=(lat_wg(:,1)); lon=(lon_wg(1,:));
timedec8day=decyear(time_start_all);
timedec8day_end=decyear(time_end_all);
timedec8day_mid=mean([timedec8day,timedec8day_end],2);
datenum_start=datenum(time_start_all);
datenum_end=datenum(time_end_all);
datenum_mid=mean([datenum_start,datenum_end],2);
datetime_mid=datetime(datenum_mid,'ConvertFrom','datenum');

%530 onwards for 2014 onwards
cafe_npp_2014onwards=cafe_npp_all_8day(:,:,530:end);
datetime_mid_2014onwards=datenum_mid(530:end);
time_start_all2014=time_start_all(530:end,:);
figure; pcolor(lon,lat,mean(cafe_npp_2014onwards,3,'omitnan')); shading flat
figure; tiledlayout('flow');
for ix=316:319
    nexttile;
    pcolor(lon,lat,cafe_npp_2014onwards(:,:,ix)); shading flat
    title(ix)
end
% 317 and 318 from 2014 to end data is WRONG
cafe_npp_2014onwards(:,:,317:318)=[];
datetime_mid_2014onwards(317:318)=[];
cafe_npp_2014onwards(cafe_npp_2014onwards<0)=NaN;
figure; pcolor(lon,lat,mean(cafe_npp_2014onwards,3,'omitnan')); shading flat

%make nc file
cd F:\Data\NPP\
% nccreatewrite('sample_file.nc','npp',{'lon','lat','datetime_mid_2014onwards'},cafe_npp_2014onwards);

clearvars -except datetime_mid_2014onwards cafe_npp_2014onwards lat lon

filename='npp_2014onwards_nans.nc';
delete(filename)
lon_att.standard_name='lon';
lon_att.long_name='longitude';
lon_att.units='degrees_east';
lon_att.axis='X';
nccreatewrite(filename,'lon',{'lon'},lon,lon_att)
lat_att.standard_name='lat';
lat_att.long_name='latitude';
lat_att.units='degrees_north';
lat_att.axis='Y';
nccreatewrite(filename,'lat',{'lat'},lat,lat_att)
date_att.long_name='datetime';
date_att.units='datetime_dd-mm-yyyy_hh:mm:ss';
date_att.axis='Z';
nccreatewrite(filename,'date',{'date'},datetime_mid_2014onwards,date_att)
npp_att.long_name='Net_Primary_Production';
npp_att.FillValue='-9999';
npp_att.units='mg/m^2/d^2';
nccreatewrite(filename,'npp',{'lat','lon','date'},cafe_npp_2014onwards(:,:,:),npp_att)

nc='npp_2014onwards_nans.nc'
A=ncinfo(nc)
A.Variables
A.Variables.Name
ncdisp(nc)

%% chl
% load('chl_8day_imported_May22.mat')
load('chl_8day_imported_May22.mat','chl_all_8day')
load('chl_8day_imported_May22.mat', 'lat_wg')
load('chl_8day_imported_May22.mat', 'lon_wg')
load('chl_8day_imported_May22.mat', 'time_end_all')
load('chl_8day_imported_May22.mat', 'time_start_all')

lat=(lat_wg(:,1)); lon=(lon_wg(1,:));
timedec8day=decyear(time_start_all);
timedec8day_end=decyear(time_end_all);
timedec8day_mid=mean([timedec8day,timedec8day_end],2);
datenum_start=datenum(time_start_all);
datenum_end=datenum(time_end_all);
datenum_mid=mean([datenum_start,datenum_end],2);
datetime_mid=datetime(datenum_mid,'ConvertFrom','datenum');

%530 onwards for 2014 onwards
chl_2014onwards=chl_all_8day(:,:,530:end);
datetime_mid_2014onwards=datenum_mid(530:end);
time_start_all2014=time_start_all(530:end,:);
figure; pcolor(lon,lat,mean(chl_2014onwards,3,'omitnan')); shading flat
figure; tiledlayout('flow');
for ix=319:324
    nexttile;
    pcolor(lon,lat,chl_2014onwards(:,:,ix)); shading flat
    title(ix)
end
% 320 and 321 from 2014 to end data is WRONG
chl_2014onwards(:,:,320:321)=[];
datetime_mid_2014onwards(320:321)=[];


%make nc file using ADD-On
% cd ../interim/
cd F:\Data\OceanProductivity\Chl\
% nccreatewrite('sample_file.nc','npp',{'lon','lat','datetime_mid_2014onwards'},cafe_npp_2014onwards);

clearvars -except datetime_mid_2014onwards chl_2014onwards lat lon
figure; pcolor(lon,lat,mean(chl_2014onwards,3,'omitnan')); shading flat

chl_2014onwards(chl_2014onwards<0)=NaN;
filename='chl_2014onwards_nans.nc';
delete(filename)
lon_att.standard_name='lon';
lon_att.long_name='longitude';
lon_att.units='degrees_east';
lon_att.axis='X';
nccreatewrite(filename,'lon',{'lon'},lon,lon_att)
lat_att.standard_name='lat';
lat_att.long_name='latitude';
lat_att.units='degrees_north';
lat_att.axis='Y';
nccreatewrite(filename,'lat',{'lat'},lat,lat_att)
date_att.long_name='datetime';
date_att.units='datetime_dd-mm-yyyy_hh:mm:ss';
date_att.axis='Z';
nccreatewrite(filename,'date',{'date'},datetime_mid_2014onwards,date_att)
chl_att.long_name='Chlorophyll_Concentration';
chl_att.FillValue='-9999';
chl_att.units='mg/m^3';
nccreatewrite(filename,'chl',{'lat','lon','date'},chl_2014onwards(:,:,:),chl_att)

nc='chl_2014onwards_nans.nc'
A=ncinfo(nc)
A.Variables
A.Variables.Name
ncdisp(nc)

%% mld
load('mld_8day_imported_May22.mat','mld_all_8day')
load('mld_8day_imported_May22.mat', 'lat_wg')
load('mld_8day_imported_May22.mat', 'lon_wg')
load('mld_8day_imported_May22.mat', 'time_end_all')
load('mld_8day_imported_May22.mat', 'time_start_all')
lat=(lat_wg(:,1)); lon=(lon_wg(1,:));
timedec8day=decyear(time_start_all);
timedec8day_end=decyear(time_end_all);
timedec8day_mid=mean([timedec8day,timedec8day_end],2);
datenum_start=datenum(time_start_all);
datenum_end=datenum(time_end_all);
datenum_mid=mean([datenum_start,datenum_end],2);
datetime_mid=datetime(datenum_mid,'ConvertFrom','datenum');

%553 onwards for 2014 onwards
mld_2014onwards=mld_all_8day(:,:,553:end);
datetime_mid_2014onwards=datenum_mid(553:end);
time_start_all2014=time_start_all(553:end,:);
figure; pcolor(lon,lat,mean(mld_2014onwards,3,'omitnan')); shading flat
figure; tiledlayout('flow');
for ix=319:324
    nexttile;
    pcolor(lon,lat,mld_2014onwards(:,:,ix)); shading flat
    title(ix)
end

%make nc file using ADD-On
% cd ../interim/
cd F:\Data\OceanProductivity\MLD\
% nccreatewrite('sample_file.nc','npp',{'lon','lat','datetime_mid_2014onwards'},cafe_npp_2014onwards);

clearvars -except datetime_mid_2014onwards mld_2014onwards lat lon
figure; pcolor(lon,lat,mean(mld_2014onwards,3,'omitnan')); shading flat

mld_2014onwards(mld_2014onwards<0)=NaN;
filename='mld_2014onwards_nans.nc';
delete(filename)
lon_att.standard_name='lon';
lon_att.long_name='longitude';
lon_att.units='degrees_east';
lon_att.axis='X';
nccreatewrite(filename,'lon',{'lon'},lon,lon_att)
lat_att.standard_name='lat';
lat_att.long_name='latitude';
lat_att.units='degrees_north';
lat_att.axis='Y';
nccreatewrite(filename,'lat',{'lat'},lat,lat_att)
date_att.long_name='datetime';
date_att.units='datetime_dd-mm-yyyy_hh:mm:ss';
date_att.axis='Z';
nccreatewrite(filename,'date',{'date'},datetime_mid_2014onwards,date_att)
mld_att.standard_name='mld';
mld_att.long_name='HYCOM_MLD';
mld_att.FillValue='-9999';
mld_att.units='m';
nccreatewrite(filename,'mld',{'lat','lon','date'},mld_2014onwards(:,:,:),mld_att)

nc='mld_2014onwards_nans.nc'
A=ncinfo(nc)
A.Variables
A.Variables.Name
ncdisp(nc)
