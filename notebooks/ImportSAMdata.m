% import sam data:
% if in columns for the month, and rows per year:
    % data imported mannually from "newsam.1957.2007.txt" in external
newsam_test=newsam1';
newsam_data=newsam_test(2:end,:);
newsam_reshape=reshape(newsam_data,[],1);

months=1:12;
repeat=length(newsam_reshape)/12;
newsam_reshape(:,2)=repmat(months,1,repeat)

years=newsam(:,1);
years_rep=repelem(years,12);
newsam_reshape(1:length(years_rep),3)=years_rep;
newsam_reshape(769:774,3)=2021;
newsam_reshape(775:end,:)=[];
