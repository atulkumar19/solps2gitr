close all;
clear all;
clc;

file = 'profilesProtoMPEX_base.nc';
% file1='/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/protoMPEX/parametricScan/no_diffision_flag2/densityScan/te8to8ne1e18to1e19/input/profilesProtoMPEX.nc';
x0 = ncread(file,'x');
z0 = ncread(file,'z');
ni0 = ncread(file,'ni');
ne0 = ncread(file,'ne');
ti0 = ncread(file,'ti');
te0 = ncread(file,'te');
vr0 = ncread(file,'vr');
vt0 = ncread(file,'vt');
vz0 = ncread(file,'vz');
br0 = ncread(file,'br');
bt0 = ncread(file,'bt');
bz0 = ncread(file,'bz');
% vz1= ncread(file1,'vz');

% Physical constants:
% =========================================================================
e_c = 1.6020e-19;
k_B = 1.3806e-23;
m_p = 1.6726e-27;
m_e = 9.1094e-31;
mu0 = 4*pi*1e-7;
c   = 299792458;
E_0 = m_p*c^2;


% Calculate Psi from Bz
% ----------------------
dx=x0(2)-x0(1);

for ii=1:numel(x0)
    for jj= 1:numel(z0)
        psi(ii,jj)  = 2*pi.*sum((bz0 (1:ii,jj).*x0(1:ii))).*dx;
    end
end

% Normalize Psi with the value of Psi at helicon location
% -------------------------------------------------------
[m,n]=size(psi);
% x=1:m;y=1:n;
% psi_new=(interp2(psi,0.0575/(max(x0)-min(x0))*m,1.745/(max(z0)-min(z0))*n));
psi_new=(interp2(psi,1.745/(max(z0)-min(z0))*n,0.0625/(max(x0)-min(x0))*m))


% Nornalized Psi
% --------------
psi=psi./psi_new;

figure;imagesc(x0(1:25:end-1),z0(1:25:end-1),psi(1:25:end-1,1:25:end-1)')

%% Construct a temperature profile using Psi_N for ProtoMPEX
% ==========================================================
te_max=8;
te_min=8;
te= (psi<1).*(te_max-te_min).*(1-psi).^1.75 + te_min;
figure;pcolor(x0(1:25:end-1),z0(1:25:end-1),te(1:25:end-1,1:25:end-1)');
figure;imagesc(x0(1:25:end-1),z0(1:25:end-1),ne0(1:25:end-1,1:25:end-1)');

%% Construct density profile using Psi_N for ProtoMPEX
% ==========================================================
ne_max=1.0E18;
ne_min=1.0E20;
ne= (psi<1).*(ne_max-ne_min).*(1-psi).^1.75 + ne_min;
figure;imagesc(x0(1:25:end-1),z0(1:25:end-1),ne(1:25:end-1,1:25:end-1)');
figure;imagesc(x0(1:25:end-1),z0(1:25:end-1),ne0(1:25:end-1,1:25:end-1)');



%% Construct a Mach profile
% =======================================

MachNum=vz0./(sqrt((8/(2*pi))*k_B.*11604.*te0./(2*m_p)));
figure; plot(z0,MachNum(1,:))



%% Update the netcdf file for GITR
% ================================

%%Modify variable here
% -------------------
ti0 = te;
te0 = te;
% vz0 = MachNum.*(sqrt((8/(2*pi))*k_B.*11604.*te./(2*m_p)));
% ne0 = ne;
% vz0=vz1;
figure;pcolor(z0(1:25:end-1),x0(1:25:end-1),vz0(1:25:end-1,1:25:end-1));

nR = length(x0);
nZ = length(z0);
ncid = netcdf.create('profilesProtoMPEX.nc','NC_WRITE');

dimR = netcdf.defDim(ncid,'nX',nR);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);
gridZnc = netcdf.defVar(ncid,'z','float',dimZ);
Ne2Dnc = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Ni2Dnc = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
Te2Dnc = netcdf.defVar(ncid,'te','float',[dimR dimZ]);
Ti2Dnc = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);
vrnc = netcdf.defVar(ncid,'vr','float',[dimR dimZ]);
vtnc = netcdf.defVar(ncid,'vt','float',[dimR dimZ]);
vznc = netcdf.defVar(ncid,'vz','float',[dimR dimZ]);
brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,x0);
netcdf.putVar(ncid,gridZnc,z0);
netcdf.putVar(ncid,Ne2Dnc,ne0);
netcdf.putVar(ncid,Ni2Dnc,ni0);
netcdf.putVar(ncid,Te2Dnc,te0);
netcdf.putVar(ncid,Ti2Dnc,ti0);

netcdf.putVar(ncid,vrnc,vr0);
netcdf.putVar(ncid,vtnc,vt0);
netcdf.putVar(ncid,vznc,vz0);

netcdf.putVar(ncid,brnc,br0);
netcdf.putVar(ncid,btnc,bt0);
netcdf.putVar(ncid,bznc,bz0);

netcdf.close(ncid);
% ti1 = ncread('profilesHelicon_new.nc','ti');