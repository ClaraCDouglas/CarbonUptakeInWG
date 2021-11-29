IN_and=inpolygon(g_lon,g_lat,andrex_box(:,1),andrex_box(:,2));
findweddell=find(IN_and==1);

testing=sum(g_area(findweddell));
testice=ice_conc(:,:,1);
testicewed=testice(findweddell);
icefree=find(testicewed<=0.1);
icefreearea=sum(g_area(icefree))
