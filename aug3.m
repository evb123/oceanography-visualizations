bbp7001=combdat.ctd.spikes.ctd_chl;
aa=bbp7001==0;
bbp7001(aa)=nan;
bbp700vc=bbp7001(:);
depthvc=combdat.ctd.zdepth(:);
timevc=combdat.ctd.time(:);

bbp700=combdat.rcf.spikes.rbr_chl;
bb=bbp700==0;
bbp700(bb)=nan;
bbp700vr=bbp700(:);
depthvr=combdat.rcf.zdepth(:);
timevr=combdat.rcf.time(:);

close all
figure
clf
scatter(timevc,depthvc,15,bbp700vc,'filled')
hold on
scatter(timevr,depthvr,15,bbp700vr,'filled')
colormap('cool')
colorbar
datetick('x','mmm/dd', 'keepticks', 'keeplimits')
title('CHL-rbr-ctd Spikes')
xlabel('Time (2017)')
ylabel('Depth (m)')
ylim([-400 -100])
hcb=colorbar;
title(hcb,'chl (ug/l)')