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
traits = table2array(T(2:height(T),2152:2153));
LMA_mes = str2double(traits(:,1));
EWT_mes = str2double(traits(:,2));

%---------- (1) Optimal domains for EWT and LMA: 1700-2400nm (Feret  et al. 2019)
wl_range = 1300:2000;
Ref_meas_sub=Ref_meas(wl_range,:);

%---------- Set the range of leaf parameters
%-----------[N, Cab, Car, Anth, Cbrown, Cw, Cm, Prot, CBC, Bspec]
P0=[1.5 40 10 0.1 0 0.01 0.01 0 0 0.2];
LB=[0.5 40 10 0.1 0 0.001 0.001 0 0 -0.2]; 
UB=[3.5 40 10 0.1 0 0.06 0.03 0 0 0.6]; 

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
EWT_est = sol(:,6)*1000;
LMA_est = sol(:,7)*1000;
