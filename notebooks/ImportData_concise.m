%% Load data or calculate from scratch
clear all
% Select algorithm:
%     algo_choice = 'cbpm'; %not so good
    algo_choice = 'vbpm'; % best?
%     algo_choice = 'eppley'; % lower values than vgpm
%     algo_choice = 'vgpm'; % new

filebase = 'C:\Users\ccd1n18\Documents\Projects\Carbon-Uptake-in-WG_Manuscript\Data\';
filedir = append(filebase,algo_choice);

D0 = dir(fullfile(filedir,'*.hdf'));


% D0 = dir('vgpm*.hdf'); % if D0 above doesn't work, go to folder with data in, and run with correct prefix

% dataload='scratch';
% dataload='saved';

% incmapplots='yes';
% incmapplots='no';

%% For bringing in new data files 
            
%% Importing multiple HDF files
    %% NPP 'algo choice' data
vgpm_npp_all=NaN*ones(1080,2160,length(D0)); 
cd ..\Data\vgpm                    
b=struct2cell(D0);

% update these limits to reflect your region of interest 
    % currently importing whole world
latres=linspace(-90,90,1080);
lonres=linspace(-180,180,2160);
lat_north=100;
lat_south=-100;
lon_west=-200;
lon_east=200;
J=find(latres>lat_south & latres<lat_north);
K=find(lonres>lon_west & lonres<lon_east);
latC=latres(J);
lonC=lonres(K);

% Loop through the multiple files to create one Matlab variable
for iix=1:length(D0)
%disp(n)
    list=b{1,iix};
% scaling equation (for the example npp time series there is no scaling
% applied)
    scale_factor=cell2mat(hdfread(list,'Slope')); 
    offset=cell2mat(hdfread(list,'Intercept'));
    B=hdfread(list,'npp');% import data ** if no scaling equation, and importing whole region, can technically add B straight to "algo"_npp_all
    C=B(J,K);  % extract region
    %C(find(C==fill))=NaN;  %replace missing value % Pete didn't do this
    %with his, so just leave that for now
    C=(C*scale_factor)+offset;
    
    vgpm_npp_all(:,:,iix)=C;
end

fill=cell2mat(hdfread(b{1,1},'Hole Value'))  
clearvars -except 'vgpm_npp_all' 'algo_choice' 'D0' 'filebase' 'filedir' 'b'
% missing data value
%                    vgpmtest=vgpm_npp_all(:,:,1);

        % quick look to check data has imported somewhat (would need to set
        % min color scale to 0, else will show min as hole value of -9999
                    %vgpm_npp_1619=double(vgpm_npp_1619);   
%                     figure;pcolor(lonC,latC,vgpm_npp_all(:,:,6));shading flat 
                    % will be upside down, but will be fine when lat_m and lon_m are made


    %% time start/end data
 time_start_all=NaN*ones(length(D0),6);
 time_end_all=NaN*ones(length(D0),6);
        
for tix=1:length(D0)
%disp(tix)
    list=b{1,tix};
    fileinfo = hdfinfo(list)
    attr_infostart = fileinfo.Attributes(3).Value
    attr_infoend = fileinfo.Attributes(4).Value
    time_start_all(tix,:) = datevec(attr_infostart);
    time_end_all(tix,:) = datevec(attr_infoend);
end

clear tix list fileinfo attr_infostart attr_infoend
    %% lat/lon data
        latrownumber=0:1:1079;
        for aix = 1:1080
            latrow=latrownumber(aix);
            spacing=1/6;
            latdistance=latrow*spacing;
            latpixelnorthedge(aix,1)=90-latdistance;
            latpixelsouthedge(aix,1)=latpixelnorthedge(aix,1)-(1/6);
            latpixelcenter(aix,1)=latpixelnorthedge(aix,1)-(1/12);
        end
        
        lonrownumber=0:1:2160;
        for aix = 1:2160
            lonrow=lonrownumber(aix);
            spacing=1/6;
            londistance=lonrow*spacing;
            lonpixelwestedge(aix,1)=londistance+-180;
            lonpixeleastedge(aix,1)=lonpixelwestedge(aix,1)+(1/6);
            lonpixelcenter(aix,1)=lonpixelwestedge(aix,1)+(1/12);
        end
        
lat_m=NaN*ones(1080,2160);
lon_m=NaN*ones(1080,2160);
for ii = 1:2160
    lat_m(:,ii)=latpixelcenter;
end
for ii = 1:1080
    lon_m(ii,:)=lonpixelcenter;
end
    %% Calculate area of boxes
        earthgeoid = almanac('earth','geoid','km','grs80'); % This is from MATLAB Mapping Toolbox
        disp(['...'])
        area_MODISvgpm_km2=NaN*ones(1080,2160);
        
        for aix = 1:length(latpixelnorthedge)
            lat_n=latpixelnorthedge(aix);
            lat_s=latpixelsouthedge(aix);
            lon_w=lonpixelwestedge(1);
            lon_e=lonpixeleastedge(1);
            area_MODISvgpm_km2(aix,1)=areaint([lat_s,lat_n,lat_n,lat_s,lat_s],[lon_w,lon_w,lon_e,lon_e,lon_w],earthgeoid);
        end
        
        for ii = 2:2160
            area_MODISvgpm_km2(:,ii)=area_MODISvgpm_km2(:,1);
        end
        % convert to m^2
        area_MODISvgpm_m2 = area_MODISvgpm_km2 * 1e6;
        
clear latrownumber lonrownumber aix earthgeoid spacing latpixelcenter latpixelnorthedge latpixelsouthedge lonpixelwestedge ...
lonpixeleastedge latrow lonpixelcenter lonrow lat_n lat_s lon_w lon_e londistance latdistance

    %% Calculate total npp for month in gC
        vgpm_npp_tot_gC_all = vgpm_npp_all;
        vgpm_npp_tot_gC_all(:)=NaN;
        fprintf(['Calculate NPP: (of ',num2str(length(D0)),'):  '])

        for ii = 1:length(D0)
            %fprintf_r(' %i', ii);
            npptp=vgpm_npp_all(:,:,ii);
            findneg=find(npptp<0);
            npptp(findneg)=0;
            vgpm_npp_tot_gC_all(:,:,ii)=(npptp.*area_MODISvgpm_m2.*time_end_all(ii,3))./1000; %Npp (mg C /m2 /day) * area (m2) * number of days / 1000 => gC per pixel in month (/1000 to convert mg to g)
        end

    %% Calculate npp for month in gC - WITH NANS
        vgpm_npp_tot_gC_nans = vgpm_npp_all;
        vgpm_npp_tot_gC_nans(:)=NaN;
        fprintf(['Calculate NPP: (of ',num2str(length(D0)),'):  '])

        for ii = 1:length(D0)
            %fprintf_r(' %i', ii);
            npptp=vgpm_npp_all(:,:,ii);
            findneg=find(npptp<0);
            npptp(findneg)=NaN;
            vgpm_npp_tot_gC_nans(:,:,ii)=(npptp.*area_MODISvgpm_m2.*time_end_all(ii,3))./1000; %Npp (mg C /m2 /day) * area (m2) * number of days / 1000 => gC per pixel in month (/1000 to convert mg to g)
        end
clear ii npptp findneg
%% Save variables
cd ..\..\Workspace_variables
save vgpm_imported -v7.3

        
        