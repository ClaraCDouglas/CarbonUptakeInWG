%% Load data or calculate from scratch
clearvars
desktop = 1;
laptop=0;
% Select algorithm:
%     algo_choice = 'cbpm'; %not so good
%     algo_choice = 'vbpm'; % best?
%     algo_choice = 'eppley'; % lower values than vgpm
    algo_choice = 'cafe'; % new
 
% filebase = 'C:\Users\ccd1n18\Documents\Projects\Carbon-Uptake-in-WG_Manuscript\Data\';
% filedir = append(filebase,algo_choice);
% D0 = dir(fullfile(filedir,'*.hdf'));

% Passport - data holdall
if laptop
    cd 'E:\Data\NPP\NPP_extracted'
    filebase = 'E:\Data\NPP\NPP_extracted';
elseif desktop
    cd 'F:\Data\NPP\NPP_extracted'
    filebase = 'E:\Data\NPP\NPP_extracted';
end
D0 = dir(fullfile(filebase,'*.hdf'));


% D0 = dir('vgpm*.hdf'); % if D0 above doesn't work, go to folder with data in, and run with correct prefix

% dataload='scratch';
% dataload='saved';

% incmapplots='yes';
% incmapplots='no';

%% For bringing in new data files 
            
%% Importing multiple HDF files
    %% NPP 'algo choice' data
% vgpm_npp_all=NaN*ones(1080,2160,length(D0)); 
cafe_npp_all_8day=NaN*ones(1080,1380,length(D0)); %(1080,2160,length(D0)); 

% cd ..\Data\cafe                    
b=struct2cell(D0);
% check attributes/file info
%     temp=D0(1,1);
%     temp=temp.name
%     temp=hdfinfo(temp)
%     temp.SDS.Attributes.Name

% update these limits to reflect your region of interest 
    % whole world dims
latres=linspace(90,-90,2160);
lonres=linspace(-180,180,4320);
    % set limits
lat_north=0; % only importing south of equator
lat_south=-100;
lon_west=-75;
lon_east=40;
J=find(latres>lat_south & latres<lat_north);
K=find(lonres>lon_west & lonres<lon_east);
latC=latres(J);
lonC=lonres(K);

% Loop through the multiple files to create one Matlab variable
for iix=1:length(D0)
    disp(iix)
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

    cafe_npp_all_8day(:,:,iix)=C; %C or B for whole region
end

% figure; pcolor(cafe_npp_all_8day(:,:,863)); shading flat

fill=cell2mat(hdfread(b{1,1},'Hole Value'))  
clearvars -except cafe* algo_choice D0 filebase b fill
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
    fileinfo = hdfinfo(list);
    attr_infostart = fileinfo.Attributes(3).Value; % start time in dd/mm/yyyy hh:mm:ss
    attr_infoend = fileinfo.Attributes(4).Value; % end time in dd/mm/yyyy hh:mm:ss
    time_start_all(tix,:) = datevec(attr_infostart);
    time_end_all(tix,:) = datevec(attr_infoend);
end

clear tix list fileinfo attr_infostart attr_infoend
    %% lat/lon data
        latrownumber=1080:1:2159; %0:1:1079;
        for aix = 1:1080 %1080
            latrow=latrownumber(aix);
            spacing=1/12;
            latdistance=latrow*spacing;
            latpixelnorthedge(aix,1)=90-latdistance;
            latpixelsouthedge(aix,1)=latpixelnorthedge(aix,1)-(1/12);
            latpixelcenter(aix,1)=latpixelnorthedge(aix,1)-(1/24);
        end
        
        lonrownumber=0:1:4320;
        for aix = 1:4320
            lonrow=lonrownumber(aix);
            spacing=1/12;
            londistance=lonrow*spacing;
            lonpixelwestedge(aix,1)=londistance+-180;
            lonpixeleastedge(aix,1)=lonpixelwestedge(aix,1)+(1/12);
            lonpixelcenter(aix,1)=lonpixelwestedge(aix,1)+(1/24);
        end
        
lat_m=NaN*ones(1080,4320);
lon_m=NaN*ones(1080,4320);
for ii = 1:4320
    lat_m(:,ii)=latpixelcenter;
end
for ii = 1:1080
    lon_m(ii,:)=lonpixelcenter;
end

    %% Calculate area of boxes
        earthgeoid = almanac('earth','geoid','km','grs80'); % This is from MATLAB Mapping Toolbox
        disp(['...'])
        area_MODIScafe_km2=NaN*ones(1080,4320);
        
        for aix = 1:length(latpixelnorthedge)
            lat_n=latpixelnorthedge(aix);
            lat_s=latpixelsouthedge(aix);
            lon_w=lonpixelwestedge(1);
            lon_e=lonpixeleastedge(1);
            area_MODIScafe_km2(aix,1)=areaint([lat_s,lat_n,lat_n,lat_s,lat_s],[lon_w,lon_w,lon_e,lon_e,lon_w],earthgeoid);
        end
        
        for ii = 2:4320
            area_MODIScafe_km2(:,ii)=area_MODIScafe_km2(:,1);
        end
        % convert to m^2
        area_MODIScafe_m2 = area_MODIScafe_km2 * 1e6;

% reduce to just wg for computing/matlab size
lon_wg=lon_m(:,1261:2640);
lat_wg=lat_m(:,1261:2640);
area_MODIScafe_km2_wg=area_MODIScafe_km2(:,1261:2640);
area_MODIScafe_m2_wg=area_MODIScafe_m2(:,1261:2640);
clearvars latrownumber lonrownumber aix earthgeoid spacing latpixelcenter latpixelnorthedge latpixelsouthedge lonpixelwestedge ...
lonpixeleastedge latrow lonpixelcenter lonrow lat_n lat_s lon_w lon_e londistance latdistance
clearvars area_MODIScafe_m2 area_MODIScafe_km2 lon_m lat_m

%% sanity check

figure; pcolor(lon_wg,lat_wg,cafe_npp_all_8day(:,:,843)); shading flat

    %% Calculate total npp for 8 day period in gC
        cafe_npp_tot_gC_all_8day = cafe_npp_all_8day;
        cafe_npp_tot_gC_all_8day(:)=NaN;
        fprintf(['Calculate NPP: (of ',num2str(length(D0)),'):  '])

        for ii = 1:length(D0)
            fprintf(' %i', ii);
            npptp=cafe_npp_all_8day(:,:,ii);
            findneg=find(npptp<0);
            npptp(findneg)=0;
%             cafe_npp_tot_gC_all(:,:,ii)=(npptp.*area_MODISvgpm_m2.*(time_end_all(ii,3))./1000; %Npp (mg C /m2 /day) * area (m2) * number of days / 1000 => gC per pixel in month (/1000 to convert mg to g)

            if ~(time_start_all(ii,2)==12)
                cafe_npp_tot_gC_all_8day(:,:,ii)=(npptp.*area_MODISvgpm_m2_wg.*8)./1000; %Npp (mg C /m2 /day) * area (m2) * number of days / 1000 => gC per pixel in month (/1000 to convert mg to g)
            elseif time_start_all(ii,2)==12
                cafe_npp_tot_gC_all_8day(:,:,ii)=(npptp.*area_MODISvgpm_m2_wg.*(time_end_all(ii,3)-time_start_all(ii,3)+1))./1000; %Npp (mg C /m2 /day) * area (m2) * number of days / 1000 => gC per pixel in month (/1000 to convert mg to g)
            else
                error('ERROR in total NPP per 8 day calc')
            end
        end

    %% Calculate npp for 8 day period in gC - WITH NANS
        cafe_npp_tot_gC_nans_8day = cafe_npp_all_8day;
        cafe_npp_tot_gC_nans_8day(:)=NaN;
        fprintf(['Calculate NPP: (of ',num2str(length(D0)),'):  '])

        for ii = 1:length(D0)
            fprintf(' %i', ii);
            npptp=cafe_npp_all_8day(:,:,ii);
            findneg=find(npptp<0);
            npptp(findneg)=NaN;
%             cafe_npp_tot_gC_nans(:,:,ii)=(npptp.*area_MODISvgpm_m2.*time_end_all(ii,3))./1000; %Npp (mg C /m2 /day) * area (m2) * number of days / 1000 => gC per pixel in month (/1000 to convert mg to g)
        
            if ~(time_start_all(ii,2)==12)
                cafe_npp_tot_gC_nans_8day(:,:,ii)=(npptp.*area_MODIScafe_m2_wg.*8)./1000; %Npp (mg C /m2 /day) * area (m2) * number of days / 1000 => gC per pixel in month (/1000 to convert mg to g)
            elseif time_start_all(ii,2)==12
                cafe_npp_tot_gC_nans_8day(:,:,ii)=(npptp.*area_MODIScafe_m2_wg.*(time_end_all(ii,3)-time_start_all(ii,3)+1))./1000; %Npp (mg C /m2 /day) * area (m2) * number of days / 1000 => gC per pixel in month (/1000 to convert mg to g)
            else
                error('ERROR in total NPP per 8 day calc with NaNs')
            end


        end
clear ii npptp findneg

%% sanity check

figure; pcolor(lon_wg,lat_wg,cafe_npp_tot_gC_nans_8day(:,:,843)); shading flat

%% Save variables
% cd ..\..\Workspace_variables
save cafe_8day_imported_eqWG_withNaNs -v7.3

 %723 s to run on laptop       
        