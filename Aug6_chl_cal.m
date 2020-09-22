%% Aug 6 - correcting eco to ctd Aug 2 - ctd correctedregressing...%% July 25, 2018 - this script is a modification of model II regression script sent by Filipa
% the corrected ecopuck chlorophyll is still really off, never approaches
% zero - actually because ctd corrected look a bit high?

%CTD chl correction
cor_chl=combdat.ctd.raw.ctd_chl*1.4905+0.246;

ctd_chl=cor_chl(:);
eco_chl=combdat.ctd.raw.eco_chl(:);

close all
xx=ctd_chl;
y=eco_chl;

%can't have nans
aa=~isnan(xx) & ~isnan(y);

[a,b,sa,sb,r,pv] = rma(xx(aa),y(aa),size(eco_chl,2));
[b,bint,l,ang,~] = rmaregress(xx,y,[2 2]);
P(1)=b(2);
P(2)=b(1);
L(1)=bint(4);
L(2)=bint(3);
U(1)=bint(2);
U(2)=bint(1);
x=[nanmin(xx);nanmax(xx)];
yz=x*L(1)+L(2);
f=fill([x; flipud(x); x(1)], [x*L(1)+L(2); flipud(x*U(1)+U(2)); yz(1)], 'r');
set(f, 'linestyle', 'none', 'facecolor', [1.0000    0.7500    0.7930]);

hold on
t=plot(eco_chl,ctd_chl,'.');
set(t, 'markersize', 7)
xlabel('Ecopuck chl(ug/L)')
ylabel('CTD chl (ug/L)')

hold on;
% Create variable name p1-p3
p=plot(x,x*P(1)+P(2), '--','Color', 'r');
l=plot(x,x*L(1)+L(2), '-','Color', 'r');
u=plot(x,x*U(1)+U(2), '-','Color', 'r');

legend([p, l],{['y=' num2str(P(1)) 'x+' num2str(P(2)) 10 'r^2=' num2str(r) 10 'p=' num2str(pv(2))], [ 10 '95% CI' 10 '.']}, 'location', 'northwest')
title('Model II Regression - ctd chl calibration')

combdat.ctd.raw.ctd_chl_corr=cor_chl;


%% calibrating chla ctd

eco_chl_cal=(0.8805*combdat.ctd.raw.eco_chl)+0.34788;

close all

clf

for i=1:7
    figure(i)
plot(combdat.ctd.raw.eco_chl(:,i),combdat.ctd.zdepth(:,i),'o','markersize',3)
hold on
plot(combdat.ctd.raw.ctd_chl_corr(:,i),combdat.ctd.zdepth(:,i),'o','markersize',3)
hold on
plot(eco_chl_cal(:,i),combdat.ctd.zdepth(:,i))
title('Ecopuck - calibration')
xlabel('Chlorophyll (ug/l)')
ylabel('Depth (m)')
legend({'raw eco chl','raw ctd chl','calibrated eco chl'},'location','south')
end

%% saving ctd ecopuck cal

combdat.ctd.raw.eco_chl_cal=eco_chl_cal;

%% looking at corrected CTD values, do they appear strange? look okay actually

for i=1:7
    figure(i)
plot(combdat.ctd.raw.ctd_chl_corr(:,i),combdat.ctd.zdepth(:,i))
end

%% next, apply correction to RCF eco puck, then calibrate RCF RBR

rcf_eco_chl_cal=(0.8805*combdat.rcf.raw.eco_chl)+0.34788;
close all
for i=10:21
    figure(i)
plot(rcf_eco_chl_cal(:,i),combdat.rcf.zdepth(:,i))
end

aa=rcf_eco_chl_cal<0;
rcf_eco_chl_cal(aa)=nan;

%% linear regression for rcf after calibrating ecopuck rcf

ctd_chl=combdat.rcf.raw.rbr_chl(:);
eco_chl=rcf_eco_chl_cal(:);

close all
xx=ctd_chl;
y=eco_chl;

%can't have nans
aa=~isnan(xx) & ~isnan(y);

[a,b,sa,sb,r,pv] = rma(xx(aa),y(aa),size(eco_chl,2));
[b,bint,l,ang,~] = rmaregress(xx,y,[2 2]);
P(1)=b(2);
P(2)=b(1);
L(1)=bint(4);
L(2)=bint(3);
U(1)=bint(2);
U(2)=bint(1);
x=[nanmin(xx);nanmax(xx)];
yz=x*L(1)+L(2);
f=fill([x; flipud(x); x(1)], [x*L(1)+L(2); flipud(x*U(1)+U(2)); yz(1)], 'r');
set(f, 'linestyle', 'none', 'facecolor', [1.0000    0.7500    0.7930]);

hold on
t=plot(xx(aa),y(aa),'.');
set(t, 'markersize', 7)
ylabel('RCF Ecopuck chl(ug/L)')
xlabel('RBR chl (ug/L)')

hold on;

p=plot(x,x*P(1)+P(2), '--','Color', 'r');
l=plot(x,x*L(1)+L(2), '-','Color', 'r');
u=plot(x,x*U(1)+U(2), '-','Color', 'r');

legend([p, l],{['y=' num2str(P(1)) 'x+' num2str(P(2)) 10 'r^2=' num2str(r) 10 'p=' num2str(pv(2))], [ 10 '95% CI' 10 '.']}, 'location', 'northwest')
title('Model II Regression - RCF CHL calibration')


%% calibrating chla rcf - this looks fantastic

rbr_chl_cal=(combdat.rcf.raw.rbr_chl*0.3327)+0.177;

close all

clf

for i=1
    figure(i)
plot(rcf_eco_chl_cal(:,i),combdat.rcf.zdepth(:,i),'o','markersize',3)
hold on
plot(combdat.rcf.raw.rbr_chl(:,i),combdat.rcf.zdepth(:,i),'o','markersize',3)
hold on
plot(rbr_chl_cal(:,i),combdat.rcf.zdepth(:,i))
title('RBR - calibration')
xlabel('Chlorophyll (ug/l)')
ylabel('Depth (m)')
legend({'rcf eco chl-cal','raw rbr chl','cal rbr chl'},'location','south')
end

%% saving calibrated files

combdat.rcf.raw.eco_chl_cal=rcf_eco_chl_cal;
combdat.rcf.raw.rbr_chl_cal=rbr_chl_cal;


