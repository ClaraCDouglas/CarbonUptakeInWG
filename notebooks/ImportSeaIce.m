%% loading in 2018-2020 data, and adding it to the data from Nov 1978- Dec 2017 from Alex Forran
clearvars
cd 'C:\Users\Clara Douglas\OneDrive - University of Southampton\Project\2607 USB\Ice-Alex';
load('nasa_seaice.mat')

cd 'E:\Data\SeaIceNIMBUS\SeaIceMonths';
a=dir('*.bin'); % list the filenames beginning with 'dt' 
b=struct2cell(a);

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
% NEED TO CONNECT TO FOLDER WITH SEA ICE IMPORT FUNCTIONS IN/BRING INTO THIS PROJECT DIRECTORY
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
ice_conc_8day=[];
yearrange=unique(year);
for yix=1:length(yearrange)
    yearsel=yearrange(yix)
    yearice=ice_conc(:,:,find(year==yearsel));
    for ix = 1:8:361
        if ix<360
        ice_conc_8day_temp=mean(yearice(:,:,(ix:1:ix+7)),3);
        ice_conc_8day=cat(3,ice_conc_8day, ice_conc_8day_temp);
        else
        ice_conc_8day_temp=mean(yearice(:,:,(ix:1:end)),3);
        ice_conc_8day=cat(3,ice_conc_8day, ice_conc_8day_temp);
        end
        clearvars ice_conc_8day_temp
    end
end

save SeaIce_daily_8day_20022020.mat -v7.3
