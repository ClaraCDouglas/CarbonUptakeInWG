figure;
% tiledlayout(4,2)
tiledlayout(2,2)
for rix = 1:length(region_sublist)
for yix=2003:2021
    Regression.tbl_pixelsYEARsel.(region_sublist{rix})=...
        table(IceFree_pixels_years_COL.(region_sublist{rix})(YEARS_COL.(region_sublist{rix})==yix),...
        AnnualNPPRate_pixels_years_COL.(region_sublist{rix})(YEARS_COL.(region_sublist{rix})==yix),...
        ... %categorical(YEARS_COL.(region_sublist{rix})(YEARS_COL.(region_sublist{rix})==yix)),
        'VariableNames',{'IceFree','AnNPPrate'});
end

rowsrecord=[];
for xix=1:length(AnnualNPPRate_pixels_years.(region_sublist{rix}))
    disp(xix);
   if sum(IceFree_pixels_years.(region_sublist{rix})(xix,:))==0
       rowsrecord=cat(1,rowsrecord,xix);
   end
end

COPYAnnualNPPRate_pixels_years.(region_sublist{rix})=AnnualNPPRate_pixels_years.(region_sublist{rix});
COPYIceFree_pixels_years.(region_sublist{rix})=IceFree_pixels_years.(region_sublist{rix});

COPYAnnualNPPRate_pixels_years.(region_sublist{rix})(rowsrecord,:)=[];
COPYIceFree_pixels_years.(region_sublist{rix})(rowsrecord,:)=[];


copy=COPYAnnualNPPRate_pixels_years.(region_sublist{rix});
[NPP_sorted,indices]=sort(copy,1,'descend');
IFD_sorted=COPYIceFree_pixels_years.(region_sublist{rix})(indices);



nexttile;
pcolor(yearrange0321,NPP_sorted,IFD_sorted);shading flat
colorbar;
xlabel('Year')
ylabel('NPP (g m^-^2 a^-^1)')
title(region_sublist{rix})

% nexttile;
% boxplot(NPP_sorted); title({region_sublist{rix},'NPP'})
% set(gca,'xtick',[1:19],'xticklabel',[2003:1:2021],'FontSize',10)
% nexttile;
% boxplot(IFD_sorted); title({region_sublist{rix},'IFD'})
% set(gca,'xtick',[1:19],'xticklabel',[2003:1:2021],'FontSize',10)
end
sgtitle('Hov plots(?), colormap is IFD')
