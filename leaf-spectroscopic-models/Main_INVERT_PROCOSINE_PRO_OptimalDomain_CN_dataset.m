% prospect: invert model to estimate leaf parameters
clear;clc;

T = readtable('..\dataset.csv');

Ref_meas = table2array(T(2:height(T),51:2151)); % get the leaf reflectance from 400 to 2500nm
Ref_meas=Ref_meas.';
length=size(Ref_meas,2); % number of samples

df_mod = csvread('..\prior.csv',0,0);
N_mod = df_mod(:,1);
Bspec_mod= df_mod(:,10);
%---------- get the measured traits
traits = table2array(T(2:height(T),2152:2156));
Chl_mes = str2double(traits(:,1));
Car_mes = str2double(traits(:,2));
LMA_mes = str2double(traits(:,3));
EWT_mes = str2double(traits(:,4));
N_mes = str2double(traits(:,5));
%---------- Set the range of leaf parameters
%-----------[N, Cab, Car, Anth, Cbrown, Cw, Cm, Prot, CBC]
P0=[1.5 40 10 0.1 0 0.01 0 0.001 0.009 0.2];
LB=[0.5 40 10 0.1 0 0.001 0 0.0001 0.0009 -0.2]; 
UB=[3.5 40 10 0.1 0 0.03 0 0.003 0.027 0.6]; 
%---------- Optimal domains for Protein: 2100-2139nm and 2160-2179nm (Féret et al. 2021) 
wl_range = [1700:1739,1760:1779];
Ref_meas_sub=Ref_meas(wl_range,:);

for i=1:length
P0(1)=N_mod(i); % use leaf structure as prior information
LB(1)=N_mod(i);
UB(1)=N_mod(i);
P0(10)=Bspec_mod(i);
LB(10)=Bspec_mod(i);
UB(10)=Bspec_mod(i);
[sol(i,:),fval(i),exitflag(i),output(i)]=fminsearchbnd('chi2P5B_PROCOSINE',P0,LB,UB,[],Ref_meas_sub(:,i),wl_range);
end

%---------- get the estimated traits
N_est = sol(:,8)*1000/4.43;
Pro_est = sol(:,8)*1000;
%---------- Optimal domains for CBC: (Féret et al. 2021) 
wl_range = [1480:1499,1560:1579,1760:1799,2040:2059,2120:2139,2160:2239,2260:2279,2340:2359,2380:2399];
wl_range = wl_range - 400;
Ref_meas_sub=Ref_meas(wl_range,:);

for i=1:length
P0(1)=N_mod(i); % use leaf structure as prior information
LB(1)=N_mod(i);
UB(1)=N_mod(i);
P0(10)=Bspec_mod(i);
LB(10)=Bspec_mod(i);
UB(10)=Bspec_mod(i);
[sol(i,:),fval(i),exitflag(i),output(i)]=fminsearchbnd('chi2P5B_PROCOSINE',P0,LB,UB,[],Ref_meas_sub(:,i),wl_range);
end

CBC_est = sol(:,9)*1000;
LMA_est = Pro_est + CBC_est;
