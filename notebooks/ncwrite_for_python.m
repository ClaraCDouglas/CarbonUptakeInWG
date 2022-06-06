load('cafe_8day_imported_eqWG_withNaNs.mat')
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

%make nc file
cd ../interim/
% nccreatewrite('sample_file.nc','npp',{'lon','lat','datetime_mid_2014onwards'},cafe_npp_2014onwards);

clearvars -except datetime_mid_2014onwards cafe_npp_2014onwards lat lon

cafe_npp_2014onwards(cafe_npp_2014onwards<0)=NaN;
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

nc='npp_2014onwards.nc'
A=ncinfo(nc)
A.Variables
A.Variables.Name
ncdisp(nc)

%% chl
load('chl_8day_imported_May22.mat')
lat=(lat_wg(:,1)); lon=(lon_wg(1,:));
timedec8day=decyear(time_start_all);
timedec8day_end=decyear(time_end_all);
timedec8day_mid=mean([timedec8day,timedec8day_end],2);
datenum_start=datenum(time_start_all);
datenum_end=datenum(time_end_all);
datenum_mid=mean([datenum_start,datenum_end],2);
datetime_mid=datetime(datenum_mid,'ConvertFrom','datenum');

%530 onwards for 2014 onwards
chl_2014onwards=chl_all(:,:,530:end);
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
cd ../interim/
cd E:\Data
% nccreatewrite('sample_file.nc','npp',{'lon','lat','datetime_mid_2014onwards'},cafe_npp_2014onwards);

clearvars -except datetime_mid_2014onwards chl_2014onwards lat lon
% figure; pcolor(lon,lat,mean(chl_2014onwards,3,'omitnan')); shading flat

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

nc='chl_2014onwards.nc'
A=ncinfo(nc)
A.Variables
A.Variables.Name
ncdisp(nc)

nccreate('chl_nctest.nc','chl',...
    'Dimensions',{'lat',length(lat),'lon',length(lon),'date',length(datetime_mid_2014onwards)},...
    'Format','classic')
ncwrite('chl_nctest.nc','chl',chl_2014onwards)

