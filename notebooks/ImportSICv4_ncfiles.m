clearvars
cd E:\Data\SeaIceV4\sidads.colorado.edu\AllFiles
% Select data folder:
filebase = 'E:\Data\SeaIceV4\sidads.colorado.edu\AllFiles';
filedir = append(filebase,'A');
D0 = dir(fullfile(filebase,'*.nc*'));

% test a file
ls; % copy a file to look at details
nc='seaice_conc_daily_sh_20211231_f17_v04r00.nc ';
A=ncinfo(nc);
A.Variables;
A.Variables.Name;
ncdisp(nc);
timetest=ncread(nc,'time');
attr_infostart = A.Attributes(39).Value;
dt = datetime(1601,1,1)+days(timetest);
% flag values  = 251,252,253,254,255 (so nan anything over 250)

%% load in all SIC data
b=struct2cell(D0);
% chl_all=NaN*ones(2160,4320,length(D0));
SICv4_0221=NaN(332,316,length(D0));
time_start_SICv4=NaN(length(D0),3);
for ix=1:length(D0)
    %SIC data
    list=b{1,ix};
    SIC=ncread(list,'cdr_seaice_conc');
    SIC(SIC>2.50)=NaN;  % fill value is  >250 : check ncdisp(nc)
    SICv4_0221(:,:,ix)=SIC';
    %datetime
    timetemp=ncread(list,'time');
    datetime_SICv4(ix) = datetime(1601,1,1)+days(timetemp);
    %time in columns
    fileinfo=ncinfo(list);
    attr_infostart = fileinfo.Attributes(39).Value;
    time_start_SICv4(ix,1) = str2double(attr_infostart(1:4));
    time_start_SICv4(ix,2) = str2double(attr_infostart(6:7));
    time_start_SICv4(ix,3) = str2double(attr_infostart(9:10));
end
% clearvars -except '*comp' 'lat_w' 'lon_w' 'speed' 'D0' 'filebase' 'filedir' 'b'
figure; pcolor(SICv4_0221(:,:,1)); shading flat;

%% 8 day averages
SICv4_8day=[];
time_start_8day_SICv4=[];
year=time_start_SICv4(:,1);
yearrange=unique(time_start_SICv4(:,1));
for yix=1:length(yearrange)
    yearsel=yearrange(yix);
    yearice=SICv4_0221(:,:,year==yearsel);
    yeartime=time_start_SICv4(year==yearsel,:);
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
