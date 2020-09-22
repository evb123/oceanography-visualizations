1000%% July 31st - calculating density in ctd files from smoothed variables

rho_CT_exact = gsw_rho_CT_exact(newctd.smooth.SA,newctd.smooth.CT,newctd.pressure_nan);

%saving

newctd.smooth.rho_exact=rho_CT_exact;

% looking at density profiles

for i=5:7
    figure
    clf
    plot(newctd.smooth.rho_exact(:,i),newctd.zdepth(:,i))
    ylim([-200 0])
end

% something weird happening in upper 10 m in either upcast or downcast,
% rapid decrease in density with depth...consider splitting upcast and
% downcast

%% splitting density upcasts and downcasts

zdepth_dc=nan(5500,7);
zdepth_uc=nan(5500,7);

rho_dc1=nan(5500,7);
rho_uc1=nan(5500,7);

for i=1:7
    
zdepth_idx=newctd.zdepth(:,i);

min_depth=nanmin(zdepth_idx);

k=find(zdepth_idx==min_depth);
downcast=zdepth_idx(1:k);
% downcast=downcast(:);
downcast(isnan(downcast))=[];

zdepth_dc(1:length(downcast),i)=downcast;

upcast=zdepth_idx(k:5500);
% upcast=upcast(:);
upcast(isnan(upcast))=[];

zdepth_uc(1:length(upcast),i)=upcast;

rho_idx=newctd.smooth.rho_exact(:,i);

rho_dc=rho_idx(1:k);
% rho_dc=rho_dc(:);
rho_dc(isnan(rho_dc))=[];

rho_dc1(1:length(rho_dc),i)=rho_dc;

rho_uc=rho_idx(k:5500);
% rho_uc=rho_uc(:);
rho_uc(isnan(rho_uc))=[];

rho_uc1(1:length(rho_uc),i)=rho_uc;

clear min_depth zdepth_idx rho_idx downcast upcast rho_dc rho_uc

end


zdepth_dc(2501:5500,:)=[];
zdepth_uc(2501:5500,:)=[];
rho_uc1(2501:5500,:)=[];
rho_dc1(2501:5500,:)=[];

%% plotting upcasts and downcasts

for i=7
    figure
    clf
    plot(rho_dc1(:,i),zdepth_dc(:,i))
%     ylim([-200 0])
end

%all downcasts have little blip in density in first ~10m

for i=7
    figure
    clf
    plot(rho_uc1(:,i),zdepth_uc(:,i))
%     ylim([-200 0])
end

%upcasts look a lot better

%% splitting N2

zdepth_dc1=nan(2500,7);
zdepth_uc1=nan(2500,7);

n2_dc1=nan(2500,7);
n2_uc1=nan(2500,7);

for i=1:7
    
zdepth_idx=newctd.zdepth(:,i);

min_depth=nanmin(zdepth_idx);

k=find(zdepth_idx==min_depth);
downcast=zdepth_idx(1:k);
% downcast=downcast(:);
% downcast(isnan(downcast))=[];

zdepth_dc1(1:length(downcast),i)=downcast;

upcast=zdepth_idx(k:2500);
% upcast=upcast(:);
%upcast(isnan(upcast))=[];

zdepth_uc1(1:length(upcast),i)=upcast;

rho_idx=newctd.smooth.N2(:,i);

rho_dc=rho_idx(1:k);
% rho_dc=rho_dc(:);
%rho_dc(isnan(rho_dc))=[];

n2_dc1(1:length(rho_dc),i)=rho_dc;

rho_uc=rho_idx(k:2500);
% rho_uc=rho_uc(:);
%rho_uc(isnan(rho_uc))=[];

n2_uc1(1:length(rho_uc),i)=rho_uc;

clear min_depth zdepth_idx rho_idx downcast upcast rho_dc rho_uc

end

%% looking at N2 downcasts and upcasts


% downcast
max_n2_dc=nanmax(n2_dc1);

for i=1:7
    depth=zdepth_dc1(:,i);
    aa=n2_dc1(:,i)==max_n2_dc(i);
    depth_test_dc(1,i)=depth(aa);
end

%upcast
max_n2_uc=nanmax(n2_uc1);

for i=1:7
    depth_idx=zdepth_uc1(:,i);
    aa=n2_uc1(:,i)==max_n2_uc(i);
    depth_n2_uc(1,i)=depth_idx(aa);
end

%% looking at 7

%upcast
for i=7
    figure
    clf
    plot(rho_uc1(:,i),zdepth_uc(:,i))
%     ylim([-200 0])
hline(depth_n2_uc(i))
xlabel('Density kg/m3')
ylabel('Depth (m)')

end

%downcast
for i=1:7
    figure
    clf
    plot(rho_dc1(:,i),zdepth_dc(:,i))
%     ylim([-200 0])
hline(depth_test_dc(i))
xlabel('Density kg/m3')
ylabel('Depth (m)')
end

%% deleting first 10 m of n2

aa=zdepth_dc1>-11;
n2_dc_clean=n2_dc1;
n2_dc_clean(aa)=nan;

aa=zdepth_uc1>-11;
n2_uc_clean=n2_uc1;
n2_uc_clean(aa)=nan;


% downcast
max_n2_dc3=nanmax(n2_dc_clean);

for i=1:7
    depth=zdepth_dc1(:,i);
    aa=n2_dc_clean(:,i)==max_n2_dc3(i);
    depth_clean_dc(1,i)=depth(aa);
end

%upcast
max_n2_uc3=nanmax(n2_uc_clean);

for i=1:7
    depth_idx=zdepth_uc1(:,i);
    aa=n2_uc_clean(:,i)==max_n2_uc3(i);
    depth_clean_uc(1,i)=depth_idx(aa);
end

%% looking at 7

%upcast
for i=1:7
    figure
    clf
    plot(rho_uc1(:,i),zdepth_uc(:,i))
%     ylim([-200 0])
hline(depth_clean_uc(i))
xlabel('Density kg/m3')
ylabel('Depth (m)')

end

%downcast
for i=1:7
    figure
    clf
    plot(rho_dc1(:,i),zdepth_dc(:,i))
%     ylim([-200 0])
hline(depth_clean_dc(i))
xlabel('Density kg/m3')
ylabel('Depth (m)')
end

%% lets look at n2 with sal and temp also

%first lets try recalculating n2 with smooth pressure

[baseline,spikes] = separate_spikes_median(newctd.pressure_nan,11);

for i=1
    figure
plot(newctd.time(:,i),-baseline(:,i))
end

[SA, in_ocean] = gsw_SA_from_SP(newctd.smooth.sal,baseline,newctd.long,newctd.lat);

CT = gsw_CT_from_t(SA,newctd.smooth.temp,baseline);

[N2, p_mid] = gsw_Nsquared(SA,CT,baseline,newctd.lat);

sdat.SA=SA;
sdat.CT=CT;
sdat.N2=N2;
sdat.pressure=baseline;
sdat.N2_depth_old_dc=depth_clean_dc;
sdat.N2_depth_old_uc=depth_clean_uc;

%% now looking at new N2, making first 11 m nan...not sure if all of this is necessary...make sure index
%and matrix match

aa=newctd.zdepth(1:330,:)>-11;
N2(aa)=nan;
aa(331:5500,:)=0;
N2=sdat.N2;
N2(aa)=nan;
N2=sdat.N2;
N2(5500,:)=nan;
N2(aa)=nan;
aa=N2==Inf;
N2(aa)=nan;
aa=N2==-Inf;
N2(aa)=nan;
sdat.N2=N2;

% might as well recalculate zdepth from smoothed pressure...

%% z from smoothed p

zdepth = gsw_z_from_p(sdat.pressure,newctd.lat);

sdat.zdepth=zdepth;

%% oh sheet lets do rho too

rho = gsw_rho(sdat.SA,sdat.CT,sdat.pressure);

sdat.rho=rho;

for i=1:7
    figure
    plot(rho(:,i),zdepth(:,i))
end

%looks good
%% Scattering N2

zdepth_N2=sdat.zdepth;
aa=isnan(sdat.N2);
zdepth_N2(aa)=nan;
time_N2=newctd.time;
time_N2(aa)=nan;

N2v=sdat.N2(:);
zdepthv=zdepth_N2(:);
timev=time_N2(:);

figure
scatter(timev,zdepthv,15,N2v,'filled')
colormap('cool')
colorbar
datetick('x','mmm/dd', 'keepticks', 'keeplimits')
title('N2 CTD')
xlabel('Time (2017)')
ylabel('Depth (m)')
ylim([-1000 0])
hcb=colorbar;
title(hcb,'N2')

% for graphing one at a time

for i=1:7
    figure
    clf
    scatter(newctd.time(:,i),sdat.zdepth(:,i),20,sdat.N2(:,i),'filled')
    colormap('cool')
colorbar
datetick('x','mmm/dd', 'keepticks')
title('N2 CTD')
xlabel('Time (2017)')
ylabel('Depth (m)')
hcb=colorbar;
title(hcb,'N2')
end

%% now refinding max n2 and depth of max n2...should i smooth n2 first?

zdepth_dc_sdat=nan(2500,7);
zdepth_uc_sdat=nan(2500,7);

n2_dc_sdat=nan(2500,7);
n2_uc_sdat=nan(2500,7);

pmid_dc_sdat=nan(2500,7);
pmid_uc_sdat=nan(2500,7);

rho_dc_sdat=nan(2500,7);
rho_uc_sdat=nan(2500,7);

sal_dc_sdat=nan(2500,7);
sal_uc_sdat=nan(2500,7);

temp_dc_sdat=nan(2500,7);
temp_uc_sdat=nan(2500,7);

for i=1:7

    %change where index comes from
zdepth_idx=sdat.zdepth(:,i);

min_depth=nanmin(zdepth_idx);

k=find(zdepth_idx==min_depth);
downcast=zdepth_idx(1:k);


zdepth_dc_sdat(1:length(downcast),i)=downcast;

upcast=zdepth_idx(k:2500);


zdepth_uc_sdat(1:length(upcast),i)=upcast;

%change where index comes from
idx=sdat.N2(:,i);
idx_dc=idx(1:k);
n2_dc_sdat(1:length(idx_dc),i)=idx_dc;
idx_uc=idx(k:2500);
n2_uc_sdat(1:length(idx_uc),i)=idx_uc;

% %change where index comes from
% idx2=p_mid(:,i);
% idx_dc2=idx2(1:k);
% pmid_dc_sdat(1:length(idx_dc2),i)=idx_dc2;
% idx_uc2=idx2(k:2500);
% pmid_uc_sdat(1:length(idx_uc2),i)=idx_uc2;

%change where index comes from
idx3=sdat.rho(:,i);
idx_dc3=idx3(1:k);
rho_dc_sdat(1:length(idx_dc3),i)=idx_dc3;
idx_uc3=idx3(k:2500);
rho_uc_sdat(1:length(idx_uc3),i)=idx_uc3;

idx4=sdat.SA(:,i);
idx_dc4=idx4(1:k);
sal_dc_sdat(1:length(idx_dc4),i)=idx_dc4;
idx_uc4=idx4(k:2500);
sal_uc_sdat(1:length(idx_uc4),i)=idx_uc4;

%change where index comes from
idx3=sdat.rho(:,i);
idx_dc3=idx3(1:k);
rho_dc_sdat(1:length(idx_dc3),i)=idx_dc3;
idx_uc3=idx3(k:2500);
rho_uc_sdat(1:length(idx_uc3),i)=idx_uc3;

idx5=sdat.CT(:,i);
idx_dc5=idx5(1:k);
temp_dc_sdat(1:length(idx_dc5),i)=idx_dc5;
idx_uc5=idx5(k:2500);
temp_uc_sdat(1:length(idx_uc5),i)=idx_uc5;
end

% sdat.downcast.N2=n2_dc_sdat;
% sdat.upcast.N2=n2_uc_sdat;
% sdat.downcast.zdepth=zdepth_dc_sdat;
% sdat.upcast.zdepth=zdepth_uc_sdat;
% sdat.downcast.pmid=pmid_dc_sdat;
% sdat.upcast.pmid=pmid_uc_sdat;
% sdat.downcast.rho=rho_dc_sdat;
% sdat.upcast.rho=rho_uc_sdat;
sdat.downcast.SA=sal_dc_sdat;
sdat.upcast.SA=sal_uc_sdat;
sdat.downcast.CT=temp_dc_sdat;
sdat.upcast.CT=temp_uc_sdat;


%% looking at N2 downcasts and upcasts


% downcast
max_n2_dc=nanmax(n2_dc_sdat);

for i=1:7
    depth=sdat.downcast.pmid(:,i);
    aa=n2_dc_sdat(:,i)==max_n2_dc(i);
    pmid_dc(1,i)=depth(aa);
end

%upcast
max_n2_uc=nanmax(n2_uc_sdat);

for i=1:7
    depth_idx=sdat.upcast.pmid(:,i);
    aa=n2_uc_sdat(:,i)==max_n2_uc(i);
    pmid_uc(1,i)=depth_idx(aa);
end

sdat.downcast.pmid_max=pmid_dc;
sdat.upcast.pmid_max=pmid_uc;

%% convert pmid to z

zdepth_pmid_dc = gsw_z_from_p(sdat.downcast.pmid_max,newctd.lat(1,:));

zdepth_pmid_uc = gsw_z_from_p(sdat.upcast.pmid_max,newctd.lat(1,:));

sdat.downcast.n2_max_z=zdepth_pmid_dc;
sdat.upcast.n2_max_z=zdepth_pmid_uc;

% now lets graph to look at rho etc. found upcast and downcast of rho
% above...shit, should do for sal and temp......done

%upcast
for i=1:7
    figure
    clf
    plot(sdat.upcast.rho(:,i),sdat.upcast.zdepth(:,i))
%     ylim([-200 0])
hline(sdat.upcast.n2_max_z(i))
xlabel('Density kg/m3')
ylabel('Depth (m)')

end
%% downcasts look better
%downcast
close all
for i=1:7
    figure(i)
    title('Max N2 depth downcasts')
    subplot(1,3,1);
    
    plot(sdat.downcast.rho(:,i),sdat.downcast.zdepth(:,i))
hline(sdat.downcast.n2_max_z(i))
xlabel('Density kg/m3')
ylabel('Depth (m)')
xlim([1026.5 1030.5])

subplot(1,3,2);

 plot(sdat.downcast.CT(:,i),sdat.downcast.zdepth(:,i))
hline(sdat.downcast.n2_max_z(i))
xlabel('Conservative Temperature (C)')
ylabel('Depth (m)')

subplot(1,3,3);

 plot(sdat.downcast.SA(:,i),sdat.downcast.zdepth(:,i))
hline(sdat.downcast.n2_max_z(i))
xlabel('Absolute Salinity (g/kg)')
ylabel('Depth (m)')
xlim([33.75 35])

end


