%% spikes count analyses 2.0 RCF

%% looking at downcast only
zdepth_dc(18000,21)=nan;
zdepth_uc(18000,21)=nan;

spike_dc(18000,21)=nan;
spike_uc(18000,21)=nan;

spike_eq_dc(18000,21)=nan;
spike_eq_uc(18000,21)=nan;

time_dc(18000,21)=nan;
time_uc(18000,21)=nan;
depth=combdat.rcf.zdepth;%depth source
var=combdat.rcf.spikes.bbp700; %var to split
var2=combdat.rcf.spikes.bbp700_eq_filt;
time=combdat.rcf.time;

for i=1:21
zdepth_idx=depth(:,i);
min_depth=nanmin(zdepth_idx);
k=find(zdepth_idx==min_depth);
downcast=zdepth_idx(1:k);
downcast(isnan(downcast))=[];
zdepth_dc(1:length(downcast),i)=downcast;
upcast=zdepth_idx(k:18000);
upcast(isnan(upcast))=[];
zdepth_uc(1:length(upcast),i)=upcast;
% var to split into downcast/upcast
idx=var(:,i);
idx_dc=idx(1:k);
idx_dc(isnan(idx_dc))=[];
spike_dc(1:length(idx_dc),i)=idx_dc;

idx_uc=idx(k:18000);
idx_uc(isnan(idx_uc))=[];
spike_uc(1:length(idx_uc),i)=idx_uc;

% var2 to split into downcast/upcast
idx9=var2(:,i);
idx_dc=idx9(1:k);
idx_dc(isnan(idx_dc))=[];
spike_eq_dc(1:length(idx_dc),i)=idx_dc;

%time
idx4=time(:,i);
idx_dc=idx4(1:k);
idx_dc(isnan(idx_dc))=[];
time_dc(1:length(idx_dc),i)=idx_dc;

idx_uc=idx4(k:18000);
idx_uc(isnan(idx_uc))=[];
time_uc(1:length(idx_uc),i)=idx_uc;

clear min_depth zdepth_idx idx downcast upcast idx_dc idx_uc
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
close all
for i=1:21
    figure
    plot(spike_eq_dc(:,i),zdepth_dc(:,i))
end

%% saving

spikecount.spike_rcf_dc=spike_dc;
spikecount.zdepth_rcf_dc=zdepth_dc;
spikecount.time_rcf_dc=time_dc;
spikecount.spike_eq_rcf_dc=spike_eq_dc;

%% calculating speed of rcf downcast
%change depth variable & time depending on speed needed
depth=spikecount.zdepth_rcf_dc;
time_idx=spikecount.time_rcf_dc;

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
spikecount.rcf_speed_dc=speed;

%% Delet column five (outside 2 stdev) and 14 (uneven speed)
for i=1:21
    figure
    plot(time_dc(:,i),zdepth_dc(:,i))
end

%% counting spikes
%matrix with spike source
% a=spikecount.spike_rcf_dc;
% aa=a==0;
% a(aa)=nan;
% spikecount.spike_rcf_dc=a;
% count_column=sum(~isnan(a),1); %1 for column, 2 for row counts
%count_row=sum(~isnan(a),2); % cant do row because depths different

%% binning count_row
% % first, find matrix/vector in depth range, then use nnz to find non-zeros
% n=100;%depth bin
% k=0:-n:-1000;
% count_bin=nan(length(k),1);
% for i=1:length(k)-1
%     j=i+1;
% aa=spikecount.zdepth_rcf_dc<=k(i)&spikecount.zdepth_rcf_dc>k(j);
% idx2=spikecount.spike_rcf_dc(aa);
% count_bin(i)=nnz(idx2);
% end
% 
% count_depth=k';

% %% plotting
% close all
% figure
% scatter(count_bin,count_depth-(n/2))
% title('100 m binning, counting rcf spikes')
% xlabel('Number of spikes')
% ylabel('Depth bin (m)')
% 
% %% saving
% spikecount.count_bin_100_rcf=count_bin;
% spikecount.count_depth_100_rcf=count_depth;

%% splitting into three time events
aa=isnan(spikecount.spike_rcf_dc);
spikecount.spike_rcf_dc(aa)=0;
spikecount.spike_eq_rcf_dc(aa)=0;

event1_s=spikecount.spike_rcf_dc(:,1:5);
event1_z=spikecount.zdepth_rcf_dc(:,1:5);
event1_f=spikecount.spike_eq_rcf_dc(:,1:5);

event2_s=spikecount.spike_rcf_dc(:,6:12);
event2_z=spikecount.zdepth_rcf_dc(:,6:12);
event2_f=spikecount.spike_eq_rcf_dc(:,6:12);

event3_s=spikecount.spike_rcf_dc(:,13:19);
event3_z=spikecount.zdepth_rcf_dc(:,13:19);
event3_f=spikecount.spike_eq_rcf_dc(:,13:19);
%% depth binning separate events and getting frequency

clear k n count_bin1 count_bin2 count_bin3
n=75;%depth bin
k=0:-n:-450;
count_bin1=nan(length(k),1);
count_bin2=nan(length(k),1);
count_bin3=nan(length(k),1);
clearvars aa i j idx idx1 idx2 idx3 idx4
for i=1:length(k)-1
    j=i+1;
    clear aa bb cc idx2 idx3 idx4
aa=event1_z<=k(i)&event1_z>k(j);
idx2=event1_s(aa);
count_bin1(i)=nnz(idx2)/nnz(aa);
%event2
bb=event2_z<=k(i)&event2_z>k(j);
idx3=event2_s(bb);
count_bin2(i)=nnz(idx3)/nnz(bb);
%event3
cc=event3_z<=k(i)&event3_z>k(j);
idx4=event3_s(cc);
count_bin3(i)=nnz(idx4)/nnz(cc);
end

count_depth_event=k';

%% plot

figure
scatter(count_bin1,count_depth_event-(n/2),'c','filled')
hold on
scatter(count_bin2,count_depth_event-(n/2),'b','filled')
hold on
scatter(count_bin3,count_depth_event-(n/2),'m','filled')
title('RCF spike frequency')
xlabel('Spike frequency')
ylabel('Depth (m)')
legend('week 1','week 2','week 3','location','southeast')
grid on
box on

%% for equal filter 

clear k n count_bin1 count_bin2 count_bin3
n=75;%depth bin
k=0:-n:-450;
count_bin1f=nan(length(k),1);
count_bin2f=nan(length(k),1);
count_bin3f=nan(length(k),1);
clearvars aa i j idx idx1 idx2 idx3 idx4
for i=1:length(k)-1
    j=i+1;
    clear aa bb cc idx2 idx3 idx4
aa=event1_z<=k(i)&event1_z>k(j);
idx2=event1_f(aa);
count_bin1f(i)=nnz(idx2)/nnz(aa);
%event2
bb=event2_z<=k(i)&event2_z>k(j);
idx3=event2_f(bb);
count_bin2f(i)=nnz(idx3)/nnz(bb);
%event3
cc=event3_z<=k(i)&event3_z>k(j);
idx4=event3_f(cc);
count_bin3f(i)=nnz(idx4)/nnz(cc);
end

count_depth_event=k';

%% plot

figure
scatter(count_bin1f,count_depth_event-(n/2),'c','filled')
hold on
scatter(count_bin2f,count_depth_event-(n/2),'b','filled')
hold on
scatter(count_bin3f,count_depth_event-(n/2),'m','filled')
title('RCF spike frequency')
xlabel('Spike frequency')
ylabel('Depth (m)')
legend('week 1','week 2','week 3','location','southeast')
grid on
box on

%% normalizing spike values/cast
% count_bin1_norm=count_bin1;
% count_bin1_norm(1:3)=count_bin1_norm(1:3)/6;
% count_bin1_norm(4:6)=count_bin1_norm(4:6)/4;
% count_bin1_norm(7:15)=count_bin1_norm(7:15)/2;
% 
% count_bin2_norm=count_bin2;
% count_bin2_norm(1:3)=count_bin2_norm(1:3)/8;
% count_bin2_norm(4:6)=count_bin2_norm(4:6)/4;
% count_bin2_norm(7:15)=nan;
% 
% count_bin3_norm=count_bin3;
% count_bin3_norm(1:3)=count_bin3_norm(1:3)/7;
% count_bin3_norm(4:6)=count_bin3_norm(4:6)/4;
% count_bin3_norm(7:15)=nan;
%% 
% figure
% scatter(count_bin1_norm,count_depth_event-(n/2),'c')
% hold on
% scatter(count_bin2_norm,count_depth_event-(n/2),'b')
% hold on
% scatter(count_bin3_norm,count_depth_event-(n/2),'m')
% title('RCF spike counts by event(75m)- normalized by cast number')
% xlabel('# spikes/# casts')
% ylabel('Depth (m)')
% legend('obs 1','obs 2','obs 3','location','southeast')
% grid on
% box on

%% saving

spikecount.norm.event1=count_bin1;
spikecount.norm.event2=count_bin2;
spikecount.norm.event3=count_bin3;
spikecount.norm.depth_bin=count_depth_event;
spikecount.norm.event1f=count_bin1f;
spikecount.norm.event2f=count_bin2f;
spikecount.norm.event3f=count_bin3f;

%% curve fitting

spike_event1=spikecount.norm.event1;
spike_event1s=spike_event1;
spike_event1s(1)=nan;
depth_event=-(spikecount.norm.depth_bin-(n/2));

spike_event2=spikecount.norm.event2;
spike_event2s=spike_event2;
spike_event2s(1)=nan;

spike_event3=spikecount.norm.event3;
spike_event3s=spike_event3;
spike_event3s(1)=nan;

%% plotting out of mat file specific filt
x=1:0.5:500;
y1=0.1836    *(x/112.5).^(-1.344    );
y2=0.1796     *(x/112.5).^(-0.8668     );
y3=0.1417    *(x/112.5).^(-0.9197  );
figure
scatter(spikecount.norm.event1,spikecount.norm.depth_bin-(n/2),'c','filled')
hold on
scatter(spikecount.norm.event2,spikecount.norm.depth_bin-(n/2),'b','filled')
hold on
scatter(spikecount.norm.event3,spikecount.norm.depth_bin-(n/2),'m','filled')
title('RCF spike frequency spec. filt')
xlabel('Spike frequency')
ylabel('Depth (m)')

grid on
box on
hold on
plot(y1,-x,'c')
xlim([0 0.4])
hold on
plot(y2,-x,'b')
hold on
plot(y3,-x,'m')
% legend('week 1','week 2','week 3','b=1.543','b=1.088','b=1.018','location','southeast')


%%
y1a=0.1836*(x/112.5).^(-1.042);
y1b=0.1836*(x/112.5).^(-1.765);

y2a=0.1796*(x/112.5).^(-0.7322);
y2b=0.1796 *(x/112.5).^(-1.001);

y3a=0.1417 *(x/112.5).^(-0.8387);
y3b=0.1417*(x/112.5).^(- 1.001);

hold on
plot(y1,-x,'c',y1a,-x,'c:',y1b,-x,'c:')
hold on
plot(y2,-x,'b',y2a,-x,'b:',y2b,-x,'b:')
hold on
plot(y3,-x,'m',y3a,-x,'m:',y3b,-x,'m:')
legend('P3A','P3B','P3C','b=1.543','b=1.088','b=1.018','location','southeast')

% title('Spike attenuation RCF')
% xlabel('Spike frequency')
% ylabel('Depth (m)')
% %legend('week 1','week 2','week 3','location','southeast')
% xlim([0 0.2])

%% now doing the same but combining all data

clear k n count_bin1 count_bin2 count_bin3
n=75;%depth bin
k=0:-n:-450;
count_bin_rcf=nan(length(k),1);

clearvars aa i j idx idx1 idx2 idx3 idx4
for i=1:length(k)-1
    j=i+1;
    clear aa bb cc idx2 idx3 idx4
aa=spikecount.zdepth_rcf_dc<=k(i)&spikecount.zdepth_rcf_dc>k(j);
idx2=spikecount.spike_rcf_dc(aa);
count_bin_rcf(i)=nnz(idx2)/nnz(aa);

end

count_depth_rcf=k';

%% plotting

figure
scatter(count_bin_rcf,count_depth_rcf-(n/2),'k','filled')
title('RCF spike frequency')
xlabel('Spike frequency')
ylabel('Depth (m)')
grid on
box on

%% saving
spikecount.rcf.spike_f=count_bin_rcf;
spikecount.rcf.depth_f=count_depth_rcf;

%% spike attenuation

rcf_s=count_bin_rcf;
rcf_s(1)=nan;

%% plotting attenuation

x=1:0.5:1100;
yall=0.1695  *(x/112.5).^(-1.011  );
yall1=0.1695 *(x/112.5).^(-0.8714);
yall2=0.1695 *(x/112.5).^(-1.151);

figure
scatter(spikecount.rcf.spike_f,spikecount.rcf.depth_f-(n/2),'r','filled')
hold on
plot(yall,-x,'r',yall1,-x,'r:',yall2,-x,'r:')
title('Spike frequency - specific filter')
xlabel('Spike frequency')
ylabel('Depth (m)')
xlim([0 0.4])
ylim([-1100 0])
grid on
box on
hold on
scatter(count_bin_eq_filt,count_depth-(n/2),'k','filled')




hold on
plot(y1,-x,'k')
xlim([0 0.4])
hold on
plot(y2,-x,'k:')
hold on
plot(y3,-x,'k:')
ylim([-1100 0])
legend('RCF','b=1.011','CI 95%','CI 95%','CTD','b=0.6406','CI 95%','CI 95%','location','southeast')

%% now repeating for equal filter

spike_event1f=spikecount.norm.event1f;
spike_event1fs=spike_event1f;
spike_event1fs(1)=nan;
depth_event=-(spikecount.norm.depth_bin-(n/2));

spike_event2f=spikecount.norm.event2f;
spike_event2fs=spike_event2f;
spike_event2fs(1)=nan;

spike_event3f=spikecount.norm.event3f;
spike_event3fs=spike_event3f;
spike_event3fs(1)=nan;

%% plotting out of mat file specific filt
x=1:0.5:500;
y1f=0.2125      *(x/112.5).^(-1.267);
y2f=0.2029       *(x/112.5).^(-0.7126 );
y3f=0.1748     *(x/112.5).^(-0.8677 );
figure
scatter(spikecount.norm.event1f,spikecount.norm.depth_bin-(n/2),'c','filled')
hold on
scatter(spikecount.norm.event2f,spikecount.norm.depth_bin-(n/2),'b','filled')
hold on
scatter(spikecount.norm.event3f,spikecount.norm.depth_bin-(n/2),'m','filled')
title('RCF spike frequency equal filter')
xlabel('Spike frequency')
ylabel('Depth (m)')

grid on
box on
hold on
plot(y1f,-x,'c')
xlim([0 0.4])
hold on
plot(y2f,-x,'b')
hold on
plot(y3f,-x,'m')
% legend('week 1','week 2','week 3','b=1.543','b=1.088','b=1.018','location','southeast')


%%
y1af=0.2125*(x/112.5).^(-0.927 );
y1bf=0.2125*(x/112.5).^(-1.608);

y2af=0.2029  *(x/112.5).^(-0.6006);
y2bf=0.2029   *(x/112.5).^(-0.8246);

y3af=0.1748   *(x/112.5).^(-0.7679);
y3bf=0.1748  *(x/112.5).^(- 0.9676);

hold on
plot(y1f,-x,'c',y1af,-x,'c:',y1bf,-x,'c:')
hold on
plot(y2f,-x,'b',y2af,-x,'b:',y2bf,-x,'b:')
hold on
plot(y3f,-x,'m',y3af,-x,'m:',y3bf,-x,'m:')
legend('P3A','P3B','P3C','b=1.267 ','b=0.7126','b=0.8677','location','southeast')

% title('Spike attenuation RCF')
% xlabel('Spike frequency')
% ylabel('Depth (m)')
% %legend('week 1','week 2','week 3','location','southeast')
% xlim([0 0.2])

%% now doing the same but combining all data

clear k n count_bin1 count_bin2 count_bin3
n=75;%depth bin
k=0:-n:-450;
count_bin_rcf_f=nan(length(k),1);

clearvars aa i j idx idx1 idx2 idx3 idx4
for i=1:length(k)-1
    j=i+1;
    clear aa bb cc idx2 idx3 idx4
aa=spikecount.zdepth_rcf_dc<=k(i)&spikecount.zdepth_rcf_dc>k(j);
idx2=spikecount.spike_eq_rcf_dc(aa);
count_bin_rcf_f(i)=nnz(idx2)/nnz(aa);

end

count_depth_rcf=k';

%% plotting

figure
scatter(count_bin_rcf_f,count_depth_rcf-(n/2),'k','filled')
title('RCF spike frequency')
xlabel('Spike frequency')
ylabel('Depth (m)')
grid on
box on

%% saving
spikecount.rcf.spike_eq_filt=count_bin_rcf_f;


%% spike attenuation

rcf_sf=count_bin_rcf_f;
rcf_sf(1)=nan;

%% plotting attenuation

x=1:0.5:1100;
yallf=0.1973    *(x/112.5).^(-0.9026    );
yall1f=0.1973   *(x/112.5).^(-0.7905);
yall2f=0.1973   *(x/112.5).^(-1.015);

figure
scatter(spikecount.rcf.spike_eq_filt,spikecount.rcf.depth_f-(n/2),'r','filled')
hold on
plot(yallf,-x,'r',yall1f,-x,'r:',yall2f,-x,'r:')
title('Spike frequency - equal filter')
xlabel('Spike frequency')
ylabel('Depth (m)')
xlim([0 0.4])
ylim([-1100 0])
grid on
box on
hold on
scatter(count_bin,count_depth-(n/2),'k','filled')




hold on
plot(y1,-x,'k')
xlim([0 0.4])
hold on
plot(y2,-x,'k:')
hold on
plot(y3,-x,'k:')
ylim([-1100 0])
legend('RCF','b=0.9026','CI 95%','CI 95%','CTD','b=0.7829','CI 95%','CI 95%','location','southeast')



