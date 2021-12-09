%% loading in 2018-2020 data, and adding it to the data from Nov 1978- Dec 2017 from Alex Forran
clearvars
cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\Project\2607 USB\Ice-Alex';
load('nasa_seaice.mat')

cd 'E:\Data\SeaIceNIMBUS\SeaIceMonths';
a=dir('*.bin'); % list the filenames beginning with 'dt'
b=struct2cell(a);

%% monthly data
for n=1:length(a)
    disp(n)
    list=b{1,n};
    filenames=cat(1,filenames,list);
    [temp_ice,header] = read_NASA_SeaIce_1(list);
    ice_conc=cat(3,ice_conc,temp_ice);
end


% datetimeee(1,n)=datetime(b{6,n},'convertfrom','datenum');
for n=1:length(a)
    disp(n)
    list=b{1,n};
    years=list(4:7); years=str2double(years);
    months=list(8:9); months=str2double(months);
    
    year=cat(1,year,years);
    month=cat(1,month,months);
    clearvars years months temp_ice header
end

clearvars list a b n

cd 'E:\Data\SeaIceNIMBUS';
save SeaIce_nasa_1978_2020.mat

%% daily data

% clearvars
cd 'E:\Data\SeaIceNIMBUS';
a=dir('nt_20??????_f1?_v1.1_s.bin'); % list the filenames for daily data only
b=struct2cell(a);

ice_conc=[];
filenames=[];
for n=1:length(a)
    %disp(n)
    list=b{1,n};
    filenames=cat(1,filenames,list);
    [temp_ice,header] = read_NASA_SeaIce_1(list);
    ice_conc=cat(3,ice_conc,temp_ice);
end

% times

% datetimeee(1,n)=datetime(b{6,n},'convertfrom','datenum');
year=[];
month=[];
day=[];
for n=1:length(a)
    disp(n)
    list=b{1,n};
    years=list(4:7); years=str2double(years);
    months=list(8:9); months=str2double(months);
    days=list(10:11); days=str2double(days);
    
    year=cat(1,year,years);
    month=cat(1,month,months);
    day=cat(1,day,days);
    clearvars years months days
end

clearvars list a b n temp_ice header

cd 'E:\Data\SeaIceNIMBUS';
save SeaIce_daily_20022020.mat



%% 8-day averages
% blank final variables to concatonate into
ice_conc_8day=[];
time_start_ice8=[];
time_end_ice8=[];
% list the years covered by dataset
yearrange=unique(year);
% make 3 column time start variable that matches daily seaice data
time_start_allice=[year month day];
for yix=1:length(yearrange)
    yearsel=yearrange(yix)
    yearice=ice_conc(:,:,(find(year==yearsel)));
    timeice=time_start_allice((find(year==yearsel)),:);
    for ix = 1:8:361
        if ix<360
            icesel=yearice(:,:,(ix:1:ix+7));
            timesel=timeice((ix:1:ix+7),:);
            ice_conc_8day_temp=mean(icesel,3);
            ice_conc_8day=cat(3,ice_conc_8day, ice_conc_8day_temp);
            time_start_temp=timesel(1,:);
            time_end_temp=timesel(end,:);
            time_start_ice8=cat(1,time_start_ice8,time_start_temp);
            time_end_ice8=cat(1,time_end_ice8,time_end_temp);
        else
            % ice_conc_8day_temp=mean(yearice(:,:,(ix:1:end)),3);
            % ice_conc_8day=cat(3,ice_conc_8day, ice_conc_8day_temp);
            icesel=yearice(:,:,(ix:1:end));
            timesel=timeice((ix:1:end),:);
            ice_conc_8day_temp=mean(icesel,3);
            ice_conc_8day=cat(3,ice_conc_8day, ice_conc_8day_temp);
            time_start_temp=timesel(1,:);
            time_end_temp=timesel(end,:);
            time_start_ice8=cat(1,time_start_ice8,time_start_temp);
            time_end_ice8=cat(1,time_end_ice8,time_end_temp);
            
        end
        clearvars *temp* *sel
    end
    clearvars yearice timeice yix ix
end

load('SeaIce_nasa_1978_2020.mat', 'g_area')
load('SeaIce_nasa_1978_2020.mat', 'g_lat')
load('SeaIce_nasa_1978_2020.mat', 'g_lon')

% save just 8-day av data (smaller file without the ice_conc variable)
clearvars ice_conc time_start_allice month year day
save SeaIce_8day_20022020.mat -v7.3
