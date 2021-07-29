%% Load data or calculate from scratch
clearvars
cd E:\Data
% Select data folder:
filebase = 'E:\Data\Chl_a\';
filedir = append(filebase,'A');
D0 = dir(fullfile(filebase,'*.nc*'));
    cd 'E:\Data\Chl_a\';

% test a file
%     ls % copy a file to look at details
%     nc='A20211522021181.L3m_MO_CHL_chlor_a_9km.nc@appkey=886fb1b19604d32206f25f491377e90ae57d87d5'
%     A=ncinfo(nc)
%     A.Variables
%     A.Variables.Name
%     ncdisp(nc)
%     % fill value for chl data is -32767, -999 for lat/lon.
%     % size is 4320x2160
%     testdata=ncread(nc, 'chlor_a'); % example variable name for SSH was Grid_0001
% 
%     testdata(find(testdata)<-900)=NaN;
%     figure; pcolor(testdata); shading flat
%     colorbar; caxis([0 2])
%% For bringing in new data files
%% Importing multiple nc files -- chlorophyll

lat_north=-40;
lat_south=-90;
lon_west=-180;  % note that in this file, longitude is 0 to 360 (not -180 to 180 as in previous files)
lon_east=180;

b=struct2cell(D0);
% chl_all=NaN*ones(2160,4320,length(D0));
chl_south=NaN*ones(600,4320,length(D0));
for ix=1:length(D0)
    list=b{1,ix};
    chl=ncread(list,'chlor_a');
    chl=chl';
    chl(chl<-900)=NaN;  % fill value (missing data) is  <0 : check ncdisp(nc)
    lat_chl=ncread(list,'lat');
    lon_chl=ncread(list,'lon');
    J=find(lat_chl>lat_south & lat_chl<lat_north);
    K=find(lon_chl>lon_west & lon_chl<lon_east);
    chl_south(:,:,ix)=chl(J,K);
%     chl_all(:,:,ix)=chl; % only run on more powerful computer (15.9GB needed)

end
lat_south=lat_chl(J);
lon_south=lon_chl(K);

clearvars -except 'chl*' 'lat_chl' 'lat_south' 'lon_chl' 'lon_south' 'D0' 'filebase' 'filedir' 'b'

%% time start/end data
time_start_all=NaN*ones(length(D0),3);
time_end_all=NaN*ones(length(D0),3);
for tix=1:length(D0)
    %disp(tix)
    list=b{1,tix};
    fileinfo = ncinfo(list);
    attr_infostart = fileinfo.Attributes(11).Value;
    attr_infoend = fileinfo.Attributes(12).Value;
    time_start_all(tix,1) = str2double(attr_infostart(1:4));
    time_start_all(tix,2) = str2double(attr_infostart(6:7));
    time_start_all(tix,3) = str2double(attr_infostart(9:10));
    time_end_all(tix,1) = str2double(attr_infoend(1:4));
    time_end_all(tix,2) = str2double(attr_infoend(6:7));
    time_end_all(tix,3) = str2double(attr_infoend(9:10));
end
clear tix list fileinfo attr_infostart attr_infoend

%% plot check
    figure; pcolor(lon_south,lat_south,chl_south(:,:,1)); shading flat
    colorbar; caxis([0 2])

%% Save variables
save chl_south -regexp ^(?!(chl|chl_all|lat_chl|lon_chl)$). -v7.3
% save chl_all -regexp ^(?!(chl|*_south)$).

%% -- SST
% SST data:
clearvars
filebase = 'E:\Data\SST_MODIS\';
filedir = append(filebase,'A');
D0 = dir(fullfile(filebase,'*.nc*'));
cd 'E:\Data\SST_MODIS\';
    
% test a file
%     ls % copy a file to look at details
%     nc='AQUA_MODIS.20210501_20210531.L3m.MO.SST.sst.9km.nc@appkey=886fb1b19604d32206f25f491377e90ae57d87d5'
%     A=ncinfo(nc)
%     A.Variables
%     A.Variables.Name
%     ncdisp(nc)
%     % fill value for chl data is -32767, -999 for lat/lon.
%     % size is 4320x2160
%     testdata=ncread(nc, 'sst'); % example variable name for SSH was Grid_0001
% 
%     testdata(find(testdata)<-900)=NaN;
%     figure; pcolor(testdata); shading flat
%     colorbar; caxis([-2 45])
%% For bringing in new data files
%% Importing multiple nc files -- chlorophyll

lat_north=-40;
lat_south=-90;
lon_west=-180;  % note that in this file, longitude is 0 to 360 (not -180 to 180 as in previous files)
lon_east=180;

b=struct2cell(D0);
% chl_all=NaN*ones(2160,4320,length(D0));
sst_south=NaN*ones(600,4320,length(D0));
for ix=1:length(D0)
    list=b{1,ix};
    sst=ncread(list,'sst');
    sst=sst';
    sst(sst<-900)=NaN;  % fill value (missing data) is  <0 : check ncdisp(nc)
    lat_sst=ncread(list,'lat');
    lon_sst=ncread(list,'lon');
    J=find(lat_sst>lat_south & lat_sst<lat_north);
    K=find(lon_sst>lon_west & lon_sst<lon_east);
    sst_south(:,:,ix)=sst(J,K);
%     sst_all(:,:,ix)=sst; % only run on more powerful computer (15.9GB needed)

end
lat_south=lat_sst(J);
lon_south=lon_sst(K);

clearvars -except 'sst*' 'lat_sst' 'lat_south' 'lon_sst' 'lon_south' 'D0' 'filebase' 'filedir' 'b'

%% time start/end data
time_start_all=NaN*ones(length(D0),3);
time_end_all=NaN*ones(length(D0),3);
for tix=1:length(D0)
    %disp(tix)
    list=b{1,tix};
    fileinfo = ncinfo(list);
    attr_infostart = fileinfo.Attributes(11).Value;
    attr_infoend = fileinfo.Attributes(12).Value;
    time_start_all(tix,1) = str2double(attr_infostart(1:4));
    time_start_all(tix,2) = str2double(attr_infostart(6:7));
    time_start_all(tix,3) = str2double(attr_infostart(9:10));
    time_end_all(tix,1) = str2double(attr_infoend(1:4));
    time_end_all(tix,2) = str2double(attr_infoend(6:7));
    time_end_all(tix,3) = str2double(attr_infoend(9:10));
end
clear tix list fileinfo attr_infostart attr_infoend

%% plot check
    figure; pcolor(lon_south,lat_south,sst_south(:,:,1)); shading flat
    colorbar; caxis([-2 20])

%% Save variables
save sst_south -regexp ^(?!(sst|sst_all|lat_sst|lon_sst)$). -v7.3
% save sst_all -regexp ^(?!(sst|*_south)$).

