%% Calculating POC flux - RCF specific threshold

%POC flux at depth z = mean spike height (in backscattering units) x 
%(spike frequency (number spikes/number obs)) x
%(bbp-to-carbon-ratio) x sinking rate (m/d)

clear aa
x=1:0.5:2000; %change according to maximum plot depth
spike_height=spikecount.spike_rcf_dc; %source of spike frequencies
aa=spike_height==0;
spike_height(aa)=nan;
spike_height_mean=nanmean(spike_height(:));

z=750; %choose depth for calculating POC flux
POC_bbp700_ratio = 31000; %change depending on ratio
sinking_rate=100; % change depending on POC sinking rate

spike_f= 0.1695 *(z/112.5)^-1.011;
spike_fa=0.1695 *(z/112.5)^-0.8714;%CI
spike_fb=0.1695 *(z/112.5)^-1.151;%CI

flux=spike_height_mean*spike_f*POC_bbp700_ratio*sinking_rate;
fluxa=spike_height_mean*spike_fa*POC_bbp700_ratio*sinking_rate;
fluxb=spike_height_mean*spike_fb*POC_bbp700_ratio*sinking_rate;

%% calculating POC flux - RCF equal threshold
clear aa
spike_height1=spikecount.spike_eq_rcf_dc;
aa=spike_height1==0;
spike_height1(aa)=nan;

spike_height_mean1=nanmean(spike_height1(:));

spike_f1= 0.1973*(z/112.5)^-0.9026;
spike_fa1=0.1973*(z/112.5)^-0.7905; 
spike_fb1=0.1973*(z/112.5)^-1.015; 

flux1=spike_height_mean1*spike_f1*POC_bbp700_ratio*sinking_rate;
fluxa1=spike_height_mean1*spike_fa1*POC_bbp700_ratio*sinking_rate;
fluxb1=spike_height_mean1*spike_fb1*POC_bbp700_ratio*sinking_rate;

%% plotting rcf flux
%RCF equal
x1=spike_height_mean1*(0.1973  *(x/112.5).^-0.9026)*POC_bbp700_ratio*sinking_rate;
xa1=spike_height_mean1*(0.1973  *(x/112.5).^-0.7905)*POC_bbp700_ratio*sinking_rate;
xb1=spike_height_mean1*(0.1973  *(x/112.5).^-1.015)*POC_bbp700_ratio*sinking_rate;

%RCF specific
x2=spike_height_mean*(0.1695  *(x/112.5).^-1.011)*POC_bbp700_ratio*sinking_rate;
xa2=spike_height_mean*(0.1695  *(x/112.5).^-0.8714)*POC_bbp700_ratio*sinking_rate;
xb2=spike_height_mean*(0.1695  *(x/112.5).^-1.151)*POC_bbp700_ratio*sinking_rate;

figure
plot(x1,-x,'k',xa1,-x,'k:',xb1,-x,'k:')
hold on
plot(x2,-x,'g',xa2,-x,'g:',xb2,-x,'g:')
xlim([0 500])
legend('RCF eq. filt','CI 95%','CI 95% ','RCF spec. filt','CI 95%',' CI 95%','location','southeast')
title('Estimated POC flux - RCF eq. and spec. threshold')
xlabel('mgC m^-^2 d^-^1')
ylabel('Depth (m)')

%% for CTD - specific filter
%% 
clear aa
spike_height2=spikecount.spike_ctd_dc;
aa=spike_height2==0;
spike_height2(aa)=nan;

spike_height_mean2=nanmean(spike_height2(:));

%z=175; %choose depth
spike_f2= 0.1777 *(z/112.5)^-0.6406  ;
spike_fa2=0.1777 *(z/112.5)^-0.5176;
spike_fb2=0.1777 *(z/112.5)^-0.7637;
%POC_bbp700_ratio = 31000;
%sinking_rate=100;

flux2=spike_height_mean2*spike_f2*POC_bbp700_ratio*sinking_rate;
fluxa2=spike_height_mean2*spike_fa2*POC_bbp700_ratio*sinking_rate;
fluxb2=spike_height_mean2*spike_fb2*POC_bbp700_ratio*sinking_rate;

%% for plotting flux -

x=1:0.5:2000;

yflux2=spike_height_mean2*(0.1777 *(x/112.5).^-0.6406 )*POC_bbp700_ratio*sinking_rate;
yfluxa2=spike_height_mean2*(0.1777 *(x/112.5).^-0.5176)*POC_bbp700_ratio*sinking_rate;
yfluxb2=spike_height_mean2*(0.1777 *(x/112.5).^-0.7637)*POC_bbp700_ratio*sinking_rate;

plot(yflux2,-x)
xlim([0 500])
%% CTD equal filter
clear aa
spike_height3=spikecount.spike_ctd_eq_filt ;
aa=spike_height3==0;
spike_height3(aa)=nan;

spike_height_mean3=nanmean(spike_height3(:));
% 
% spike_height=spikecount.spike_rcf_dc;
% aa=isnan(spikecount.zdepth_rcf_dc);
% spike_height(aa)=nan;
% bb=~isnan(spike_height);
% num_obs=nnz(bb);
% num_spikes=nnz(spikecount.spike_rcf_dc);
%spike_f=num_spikes/num_obs;
% z=175; %choose depth
spike_f3= 0.1491   *(z/112.5)^-0.7829;
spike_fa3=0.1491   *(z/112.5)^-0.6445;
spike_fb3=0.1491   *(z/112.5)^-0.9212;

flux3=spike_height_mean3*spike_f3*POC_bbp700_ratio*sinking_rate;
fluxa3=spike_height_mean3*spike_fa3*POC_bbp700_ratio*sinking_rate;
fluxb3=spike_height_mean3*spike_fb3*POC_bbp700_ratio*sinking_rate;

%% for plotting flux CTD

x=1:0.5:2000;

yflux3=spike_height_mean3*(0.1491*(x/112.5).^-0.7829 )*POC_bbp700_ratio*sinking_rate;
yfluxa3=spike_height_mean3*(0.1491*(x/112.5).^-0.6445)*POC_bbp700_ratio*sinking_rate;
yfluxb3=spike_height_mean3*(0.1491*(x/112.5).^-0.9212)*POC_bbp700_ratio*sinking_rate;

figure
plot(yflux3,-x,'k',yfluxa3,-x,'k:',yfluxb3,-x,'k:')
hold on
plot(yflux2,-x,'g',yfluxa2,-x,'g:',yfluxb2,-x,'g:')
xlim([0 500])
legend('CTD eq. filt','CI 95%','CI 95% ','CTD spec. filt','CI 95%',' CI 95%','location','southeast')
title('Estimated POC flux - CTD equal and spec. threshold')

%% Now plotting rcf ctd specific

x=1:0.5:2000;

figure
plot(yflux2,-x,'k',yfluxa2,-x,'k:',yfluxb2,-x,'k:')
hold on
plot(x2,-x,'r',xa2,-x,'r:',xb2,-x,'r:')

legend('CTD spec. threshold','CI 95%','CI 95% ','RCF spec. threshold','CI 95%',' CI 95%','location','southeast')
title('Estimated POC flux - platform-specific threshold')
xlim([0 500])
xlabel('mgC m^-^2 d^-^1')
ylabel('Depth (m)')

%% Now plotting rcf ctd equal threshold

x=1:0.5:2000;

figure
% subplot(2,1,1)
% figure
% plot(flux2,-x,'k',fluxa2,-x,'k:',fluxb2,-x,'k:')
% hold on
% plot(x2,-x,'r',xa2,-x,'r:',xb2,-x,'r:')
% 
% legend('CTD spec. threshold','CI 95%','CI 95% ','RCF spec. threshold','CI 95%',' CI 95%','location','southeast')
% title('Estimated POC flux - platform-specific threshold')
% xlim([0 500])
% xlabel('mgC m^-^2 d^-^1')
% ylabel('Depth (m)')
% hold on
% subplot(2,1,2)
plot(yflux3,-x,'k',yfluxa3,-x,'k:',yfluxb3,-x,'k:')
hold on
plot(x1,-x,'r',xa1,-x,'r:',xb1,-x,'r:')
xlim([0 500])
legend('CTD eq. filt','CI 95%','CI 95% ','RCF eq. filt','CI 95%',' CI 95%','location','southeast')
title('Estimated POC flux - equal threshold')
xlabel('mgC m^-^2 d^-^1')
ylabel('Depth (m)')


