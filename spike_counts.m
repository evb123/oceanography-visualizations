%% spikes count analyses 2.0

for i=1:7
    figure
    plot(combdat.ctd.time(:,i),combdat.ctd.zdepth(:,i))
end

%% looking at downcast only
zdepth_dc(5500,7)=nan;
zdepth_uc(5500,7)=nan;

spike_dc(5500,7)=nan;
spike_uc(5500,7)=nan;

spike_eq_dc(5500,7)=nan;
spike_eq_uc(5500,7)=nan;

time_dc(5500,7)=nan;
time_uc(5500,7)=nan;

for i=1:7
zdepth_idx=combdat.ctd.zdepth(:,i);
min_depth=nanmin(zdepth_idx);
k=find(zdepth_idx==min_depth);
downcast=zdepth_idx(1:k);
downcast(isnan(downcast))=[];
zdepth_dc(1:length(downcast),i)=downcast;
upcast=zdepth_idx(k:5500);
upcast(isnan(upcast))=[];
zdepth_uc(1:length(upcast),i)=upcast;
% var to split into downcast/upcast
idx=combdat.ctd.spikes.bbp700(:,i);
idx_dc=idx(1:k);
idx_dc(isnan(idx_dc))=[];
spike_dc(1:length(idx_dc),i)=idx_dc;

idx_uc=idx(k:5500);
idx_uc(isnan(idx_uc))=[];
spike_uc(1:length(idx_uc),i)=idx_uc;

% var to split into downcast/upcast
idx9=combdat.ctd.spikes.bbp700_eq_filt(:,i);
idx_dc=idx9(1:k);
idx_dc(isnan(idx_dc))=[];
spike_eq_dc(1:length(idx_dc),i)=idx_dc;

idx_uc=idx9(k:5500);
idx_uc(isnan(idx_uc))=[];
spike_eq_uc(1:length(idx_uc),i)=idx_uc;
%time
idx4=combdat.ctd.time(:,i);
idx_dc=idx4(1:k);
idx_dc(isnan(idx_dc))=[];
time_dc(1:length(idx_dc),i)=idx_dc;

idx_uc=idx4(k:5500);
idx_uc(isnan(idx_uc))=[];
time_uc(1:length(idx_uc),i)=idx_uc;

clear min_depth zdepth_idx idx downcast upcast idx_dc idx_uc idx9
end

%% check by plotting

aa=zdepth_dc==0;
zdepth_dc(aa)=nan;
spike_dc(aa)=nan;
time_dc(aa)=nan;
aa=zdepth_uc==0;
zdepth_uc(aa)=nan;
spike_uc(aa)=nan;
time_uc(aa)=nan;
%%
for i=1:7
    figure
    plot(spike_eq_dc(:,i),zdepth_dc(:,i))
end

%% saving

spikecount.spike_ctd_dc=spike_dc;
spikecount.zdepth_ctd_dc=zdepth_dc;
spikecount.time_ctd_dc=time_dc;
spikecount.spike_ctd_eq_filt=spike_eq_dc;

%% calculating speed of ctd downcast
%change depth variable & time depending on speed needed
depth=spikecount.zdepth_ctd_dc;
time_idx=spikecount.time_ctd_dc;

%speed calc
distance=nanmin(depth);
time=nanmax(time_idx)-nanmin(time_idx);%[days]
times=time*86400;%convert to seconds to get speed in [m/s]
speed=distance./-times;

mean_speed=nanmean(speed);
std=std(speed);
max_speed=mean_speed+2*std;
min_speed=mean_speed-2*std;

%% saving
spikecount.ctd_speed_dc=speed;


%% counting spikes
%matrix with spike source
% a=spikecount.spike_ctd_dc;
% aa=a==0;
% a(aa)=nan;
% 
% count_column=sum(~isnan(a),1); %1 for column, 2 for row counts
% %count_row=sum(~isnan(a),2); % cant do row because depths different

%% binning count_row
% first, find matrix/vector in depth range, then use nnz to find non-zeros
clear n k i j aa idx2 id3
n=75;%depth bin
k=0:-n:-1050;
count_bin=nan(length(k),1);
count_bin_eq_filt=nan(length(k),1);
for i=1:length(k)-1
    j=i+1;
aa=spikecount.zdepth_ctd_dc<=k(i)&spikecount.zdepth_ctd_dc>k(j);
idx2=spikecount.spike_ctd_dc(aa);
count_bin(i)=nnz(idx2)/nnz(aa);

idx3=spikecount.spike_ctd_eq_filt(aa);
count_bin_eq_filt(i)=nnz(idx3)/nnz(aa);
end

count_depth=k';

%% plotting
close all
figure
scatter(count_bin,count_depth-(n/2))
title('75 m binning, counting ctd spikes')
xlabel('Spike frequency')
ylabel('Depth bin (m)')

figure
scatter(count_bin_eq_filt,count_depth-(n/2))
title('75 m binning, counting ctd spikes')
xlabel('Spike frequency')
ylabel('Depth bin (m)')

%% saving
spikecount.ctd_freq_75=count_bin;
spikecount.ctd_freq_75_eq_filt=count_bin_eq_filt;
spikecount.ctd_depth_75=count_depth;

%% normalizing (w/ 75 m bins)
% 
% count_bin_norm=count_bin;
% count_bin_norm(1:7)=count_bin_norm(1:7)/7;
% count_bin_norm(8:15)=count_bin_norm(8:15)/6;
% 
% %% plot
% 
% figure
% scatter(count_bin_norm,count_depth-(n/2))
% title('75 m binning, normalized ctd spikes')
% xlabel('#spikes/#downcasts')
% ylabel('Depth (m)')
% 
% %% saving
% spikecount.norm.ctd_75=count_bin_norm;
% spikecount.norm.ctd_75d=count_depth;


%% 
ctd_s=count_bin;
ctd_s1=ctd_s;
ctd_s1(1)=nan;

ctd_3=count_bin_eq_filt;
ctd_s2=ctd_3;
ctd_s2(1)=nan;

depth_event=-(count_depth-(n/2));

%% ctd - specific filter

x=1:0.5:1100;
y1=0.1777   *(x/112.5).^(-0.6406  );
y2=0.1777   *(x/112.5).^(-0.5176  );
y3=0.1777   *(x/112.5).^(-0.7637 );
 
figure
scatter(count_bin,count_depth-(n/2),'k','filled')
title('CTD spike frequency spec. filter')
xlabel('Spike frequency')
ylabel('Depth bin (m)')
legend('ctd','b=0.6406','CI 95%','location','southeast')
grid on
box on
hold on
plot(y1,-x,'k')
xlim([0 0.4])
hold on
plot(y2,-x,'k:')
hold on
plot(y3,-x,'k:')
ylim([-1100 0])

%% %% ctd - equal filter to RCF

x=1:0.5:1100;
y1=0.1491    *(x/112.5).^(-0.7829    );
y2=0.1491    *(x/112.5).^(-0.6445  );
y3=0.1491     *(x/112.5).^(-0.9212 );
 
figure
scatter(count_bin_eq_filt,count_depth-(n/2),'k','filled')
title('CTD spike frequency eq. filter')
xlabel('Spike frequency')
ylabel('Depth bin (m)')

grid on
box on
hold on
plot(y1,-x,'k')
xlim([0 0.4])
hold on
plot(y2,-x,'k:')
hold on
plot(y3,-x,'k:')
ylim([-1100 0])
legend('ctd','b=0.78','CI 95%','location','southeast')