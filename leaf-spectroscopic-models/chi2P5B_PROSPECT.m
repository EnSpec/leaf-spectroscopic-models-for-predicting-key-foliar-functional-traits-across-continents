% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2P5B_PROSPECT(x,rmes,wl_range)

N=x(1);
Cab=x(2);
Car=x(3);
Ant=x(4);
Brown=x(5);
Cw=x(6);
Cm=x(7);
Prot=x(8);
CBC=x(9);


RT=prospect_PRO_v3(N,Cab,Car,Ant,Brown,Cw,Cm,Prot,CBC);

chi2=sqrt(sum((RT(wl_range,2)-rmes).^2));
