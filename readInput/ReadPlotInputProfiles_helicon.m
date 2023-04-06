close all;
clear all;
clc;

% Physical constants:
% =========================================================================
e_c = 1.6020e-19;
ep_0= 8.85E-12;
k_B = 1.3806e-23;
m_p = 1.6726e-27;

m_e = 9.1094e-31;

m_i = 25*m_e;
mu0 = 4*pi*1e-7;
c   = 299792458;
E_0 = m_p*c^2;

n_0=3.2E26;
% T_e=1.27E5;
v_i=0.1*c;
% 
 T_e=v_i^2*m_i/k_B

wpe=sqrt(n_0*e_c^2/(ep_0*m_e))
di=c/wpe

l_d=sqrt((ep_0*k_B*T_e)/(n_0*e_c^2))


%%
x=ncread('profilesHelicon.nc','x');
z=ncread('profilesHelicon.nc','z');

% x=ncread('profilesHelicon.nc','x');
% z=ncread('profilesHelicon.nc','z');

ne=ncread('profilesHelicon.nc','ne');
figure; imagesc(z,x,ne);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Density')
colorbar;



te=ncread('profilesHelicon.nc','te');
figure; imagesc(z,x,te);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Te')
colorbar;


br=ncread('profilesHelicon.nc','br');
figure;imagesc(z,x,br);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input br')

bz=ncread('profilesHelicon.nc','bz');
figure;imagesc(z,x,bz);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Bz')



