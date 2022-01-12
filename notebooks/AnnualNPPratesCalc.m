leapyear1=time_start_all(:,1);
leapyear2=rem(leapyear1,4)==0;
leapyear3=zeros(length(leapyear1),1);
leapyear3(leapyear2==1)=366;
leapyear3(leapyear3==0)=365;

days_test(:,1)=(timedec8day_end-timedec8day);
days_test(:,2)=days_test(:,1)*365;
days_test(:,3)=days_test(:,1).*leapyear3;

day_chunk=days_test(:,1).*leapyear3;

for aix = 1:length(algorithm)
    for rix = 1:length(region_sublist)
        OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).day_chunk_rates_gm2=(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_av_mgm2_nans.*day_chunk)/1000;
        for yix = 2003:2020
            findyear=find(timedec8day>yix-0.5 & timedec8day<yix+0.5);
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).annual_rate_chunk(yix-2002,1)= ...
                sum(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).day_chunk_rates_gm2(findyear),'omitnan');
            day_chunk_IF=day_chunk(findyear);
            day_chunk_IF(isnan(OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).day_chunk_rates_gm2(findyear)))=NaN;
            OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).days_ice_free(yix-2002,1)=sum(day_chunk_IF,'omitnan')
        end
    end
end




%% attempts 1 and 2 - check over again for validity
    annual_rate=nan(1080,1380,length(setup.yearrange0320));
    for yix = 2003:2020
        findyear=find(timedec8day>yix-0.5 & timedec8day<yix+0.5);
        temp.timestart=time_start_all(findyear,:);
        temp.rate = cafe_npp_all_8day(:,:,findyear);
        temp.rateyr=cafe_npp_all_8day(:,:,findyear);
        temp.rateyr(:)=NaN;
        for ii = 1:length(findyear)
            fprintf(' %i', ii);
            npptp=temp.rate(:,:,ii);
            findneg=find(npptp<0);
            npptp(findneg)=NaN;

            if ~(temp.timestart(ii,2)==12)
                temp.rateyr(:,:,ii)=(npptp.*8); %Npp (mg C /m2 /day) * number of days  => mgC per m2 in 8 days
            elseif temp.timestart(ii,2)==12
                temp.rateyr(:,:,ii)=(npptp.*(temp.timestart(ii,3)-temp.timestart(ii,3)+1)); %Npp (mg C /m2 /day) * number of days  => mgC per m2 in 8 days
            else
                error('ERROR in 8 day NPP rates calc with NaNs')
            end
        end
        temp.rateann=sum(temp.rateyr,3,'omitnan');
        annual_rate(:,:,yix-2002)=temp.rateann;
    end

    rate_8day=cafe_npp_all_8day;
    for ii = 1:length(D0)
    rate_8dayt=cafe_npp_all_8day(:,:,ii);
    findneg=find(rate_8dayt<0);
    rate_8dayt(findneg)=NaN;

        if ~(time_start_all(ii,2)==12)
            rate_8day(:,:,ii)=((rate_8dayt.*8)/1000); %Npp (mg C /m2 /day) * number of days /1000 => gC per m2 in 8 days
        elseif time_start_all(ii,2)==12
            rate_8day(:,:,ii)=((rate_8dayt.*(time_end_all(ii,3)-time_start_all(ii,3)+1))/1000); %Npp (mg C /m2 /day) * number of days /1000 => gC per m2 in 8 days
        else
            error('ERROR in 8 day NPP rates calc with NaNs')
        end
    end
    annual_rate2=nan(1080,1380,length(setup.yearrange0320));
    for yix = 2003:2020
        findyear=find(timedec8day>yix-0.5 & timedec8day<yix+0.5);
        temp.temp=rate_8day(:,:,findyear);
        annual_rate2(:,:,(yix-2002))=sum(temp.temp,3,'omitnan');
    end

    for aix = 1:length(algorithm)
        for rix = 1:length(region_sublist)
            for tix=1:length(setup.yearrange0320)
                % Annual NPP rate (mg m-2 yr-1)
                temp.NPP_annrate=annual_rate(:,:,tix);
                temp.NPP_annrate2=annual_rate2(:,:,tix);
                % temp.NPP_annrate is integrate NPP rate over year 
                    %now summed for pixels within region
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_annrate_gCm2yr1_1(tix,1)=(mean(mean(temp.NPP_annrate(temp.(region_sublist{rix}).box_logic),'omitnan'),'omitnan'))/1000;
                OceanProd_8day.(algorithm{aix}).(region_sublist{rix}).NPP_annrate_gCm2yr1_2(tix,1)=mean(mean(temp.NPP_annrate2(temp.(region_sublist{rix}).box_logic),'omitnan'),'omitnan');
            end
        end
    end


