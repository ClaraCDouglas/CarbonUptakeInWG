timedec_momid_cut=timedec_momid(85:96);
timedec8day_mid_cut=timedec8day_mid(323:368);
motot=OceanProd.cafe.Open.NPP_tot_gC(85:96);
day8tot=OceanProd_8day.cafe.Open.NPP_tot_gC(323:368);

motot_cum=cumsum(motot);
day8tot_cum=cumsum(day8tot);
figure
plot(timedec_momid_cut,motot_cum,'o-')
hold on
plot(timedec8day_mid_cut,day8tot_cum,'o-')
