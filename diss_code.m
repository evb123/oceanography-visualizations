%% Estimating POC flux from spike frequencies - different platforms and noise thresholds
% August 28, 2018
% Evelyn Byer

%% Calculating POC flux - RCF specific threshold
%POC flux at depth z = mean spike height (in backscattering units) x 
%(spike frequency (number spikes/number obs)) x
%(bbp-to-carbon-ratio) x sinking rate (m/d)

%need spike frequqency matrix 'spikecount.mat'

%outputs=
    %flux= RCF-specific threshold
    %flux1= RCF-equal threshold
    %flux2= CTD-specific threshold
    %flux3= CTD-equal threshold
    %flux#a/b= 95% CI based on b-values

clear aa
x=1:0.5:2000; %change according to maximum plot depth
spike_height=spikecount.spike_rcf_dc; %source of spike frequencies
aa=spike_height==0;
spike_height(aa)=nan;
spike_height_mean=nanmean(spike_height(:));

z=2000; %choose depth for calculating POC flux
POC_bbp700_ratio = 31000; %change depending on ratio
sinking_rate=100; % change depending on POC sinking rate

spike_f= 0.1695 *(z/112.5)^-1.011;
spike_fa=0.1695 *(z/112.5)^-0.8714;%CI
spike_fb=0.1695 *(z/112.5)^-1.151;%CI

flux=spike_height_mean*spike_f*POC_bbp700_ratio*sinking_rate; %POC flux
fluxa=spike_height_mean*spike_fa*POC_bbp700_ratio*sinking_rate; %CI
fluxb=spike_height_mean*spike_fb*POC_bbp700_ratio*sinking_rate; %CI

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

%% for CTD - specific threshold
clear aa
spike_height2=spikecount.spike_ctd_dc;
aa=spike_height2==0;
spike_height2(aa)=nan;

spike_height_mean2=nanmean(spike_height2(:));

spike_f2= 0.1777 *(z/112.5)^-0.6406  ;
spike_fa2=0.1777 *(z/112.5)^-0.5176;
spike_fb2=0.1777 *(z/112.5)^-0.7637;

flux2=spike_height_mean2*spike_f2*POC_bbp700_ratio*sinking_rate;
fluxa2=spike_height_mean2*spike_fa2*POC_bbp700_ratio*sinking_rate;
fluxb2=spike_height_mean2*spike_fb2*POC_bbp700_ratio*sinking_rate;

%% POC flux estimate from CTD - equal threshold
clear aa
spike_height3=spikecount.spike_ctd_eq_filt ;
aa=spike_height3==0;
spike_height3(aa)=nan;
spike_height_mean3=nanmean(spike_height3(:));

spike_f3= 0.1491*(z/112.5)^-0.7829;
spike_fa3=0.1491*(z/112.5)^-0.6445;
spike_fb3=0.1491*(z/112.5)^-0.9212;

flux3=spike_height_mean3*spike_f3*POC_bbp700_ratio*sinking_rate;
fluxa3=spike_height_mean3*spike_fa3*POC_bbp700_ratio*sinking_rate;
fluxb3=spike_height_mean3*spike_fb3*POC_bbp700_ratio*sinking_rate;