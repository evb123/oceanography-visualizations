%% Aug - 8 spikes and baseline for bbp700 rcf and ctd

%variable you want to smooth
y=combdat.ctd.raw.bbp700;

%window size
n=15;

z=combdat.ctd.zdepth; %depth

[baseline,spikes] = separate_spikes_median(y,n);

% plotting baseline over ra
close all
for i=1:7 %change for different casts
figure
clf
plot(y(:,i),z(:,i))
hold on
plot(baseline(:,i),z(:,i),'o','markersize',1)
title('ctd bbp700')%change title
xlim([0 0.006])
figure
plot(spikes(:,i),z(:,i))
title('rcf bbp700 spikes')% change title

figure
plot(baseline(:,i),z(:,i))
end

%% Once baseline looks good, find threshold for spikes
%threshold can make minimum to zero and look at average or counts
%can set maximum threshold and convert to nans
%setting minimum threshold to 2x smallest value over 0?
clear aa
spikes_clean=spikes;
aa=spikes_clean<=0;
spikes_clean(aa)=0;
median=nanmedian(spikes_clean(spikes_clean>0));

%min_spike=nanmin(spikes_clean);
%mean_min=nanmean(min_spike);

clear aa
filter=2*median;
aa=spikes_clean<filter; %filter
spikes_clean(aa)=0;

figure
clf
plot(spikes_clean(:,1),z(:,1))

%% saving
% combdat.ctd.spikes.bbp700=spikes_clean;
% combdat.ctd.smooth.bbp700=baseline;

combdat.ctd.spikes.bbp700=spikes_clean;
%combdat.ctd.spikes.bbp700_eq_filter=spikes_clean;
%combdat.ctd.smooth.bbp700=baseline;

%% ploting baselines
close all
for i=1:7
    figure(i)
    plot(combdat.ctd.smooth.bbp700(:,i),combdat.ctd.zdepth(:,i))
end

ctd_5bbp=combdat.ctd.smooth.bbp700(:,5);
kk=nanmin(ctd_5bbp);

bb700_5_int=ctd_5bbp-kk;

figure
plot(bb700_5_int,combdat.ctd.zdepth(:,5))

combdat.ctd.smooth.bbp700(:,5)=bb700_5_int;

%% now plotting RCF baselines to look for negs
%smoothing rcf
%variable you want to smooth
clear y n z
close all
y=combdat.rcf.raw.bbp700;

%window size
n=15;

z=combdat.rcf.zdepth; %depth

[baseline1,spikes1] = separate_spikes_median(y,n);

% plotting baseline over ra
close all
for i=1:21 %change for different casts
figure
clf
plot(y(:,i),z(:,i))
hold on
plot(baseline1(:,i),z(:,i),'o','markersize',1)
title('rcf bbp700')%change title
xlim([0 0.007])
figure
plot(spikes1(:,i),z(:,i))
title('rcf bbp700 spikes')% change title

figure
plot(baseline1(:,i),z(:,i))
end

%% Once baseline looks good, find threshold for spikes
%threshold can make minimum to zero and look at average or counts
%can set maximum threshold and convert to nans
clear aa
spikes_clean1=spikes1;
% aa=spikes_clean1<=0;
% spikes_clean1(aa)=0;
% min_spike1=nanmin(spikes_clean1);
% mean_min1=nanmean(min_spike1);
median1=nanmedian(spikes_clean1(spikes_clean1>0));
clear aa
filter1=2*median1;
aa=spikes_clean1<filter1; %filter
spikes_clean1(aa)=0;

figure
clf
plot(spikes_clean1(:,1),z(:,1))



close all
for i=8
    figure(i)
    plot(combdat.rcf.smooth.bbp700(:,i),combdat.rcf.zdepth(:,i))
end

%% saving rcf
combdat.rcf.spikes.bbp700=spikes_clean1;
%combdat.rcf.smooth.bbp700=baseline;

%% using same threshold

filt_comb=(filter+filter1)/2;

clear aa
spikes_clean2=spikes1;
aa=spikes_clean2<=0;
spikes_clean2(aa)=0;

clear aa

aa=spikes_clean2<filt_comb; %filter
spikes_clean2(aa)=0;

combdat.rcf.spikes.bbp700_eq_filt=spikes_clean2;

clear aa
spikes_clean3=spikes;
aa=spikes_clean3<=0;
spikes_clean3(aa)=0;
clear aa

aa=spikes_clean3<filt_comb; %filter
spikes_clean3(aa)=0;

combdat.ctd.spikes.bbp700_eq_filt=spikes_clean3;

    