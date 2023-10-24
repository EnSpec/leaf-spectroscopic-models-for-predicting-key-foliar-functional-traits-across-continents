% prospect: invert model to estimate leaf parameters
clear;clc;

T = readtable('..\dataset.csv');

Ref_meas = table2array(T(2:height(T),51:2151)); % get the leaf reflectance from 400 to 2500nm
Ref_meas=Ref_meas.';
length=size(Ref_meas,2); % number of samples
%---------- get the measured traits
traits = table2array(T(2:height(T),2152:2156));
Chl_mes = str2double(traits(:,1));
Car_mes = str2double(traits(:,2));
LMA_mes = str2double(traits(:,3));
EWT_mes = str2double(traits(:,4));
N_mes = str2double(traits(:,5));

%----------  Set the range of leaf parameters
%-----------[N, Cab, Car, Anth, Cbrown, Cw, Cm, Prot, CBC, Bspec]
P0=[1.5 40 10 0.1 0 0.01 0 0.001 0.009 0.2];
LB=[0.5 0.5 0.5 0.0 0 0.001 0 0.0001 0.0009 -0.2]; 
UB=[3.5 120 25 20 0 0.06 0 0.003 0.027 0.6]; 

%---------- Set the spectral domain used in model inversion
%---------- Full domain: 400-2500nm
wl_range = 1:2101;
Ref_meas=Ref_meas(wl_range,:);

for i=1:length
[sol(i,:),fval(i),exitflag(i),output(i)]=fminsearchbnd('chi2P5B_PROCOSINE',P0,LB,UB,[],Ref_meas(:,i),wl_range);
end

%---------- get the estimated traits
Chl_est = sol(:,2);
Car_est = sol(:,3);
EWT_est = sol(:,6)*1000;
LMA_est = sol(:,8)*1000 + sol(:,9)*1000;
N_est = sol(:,8)*1000/4.43;
