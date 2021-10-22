%% The output file is then used to compute the EOFs.
clear all;
close all;
clc

%% --> Paths and Grid;
fpath  = '/home/hajsong/PDAF/tutorial_global_oce_latlon_v1.16_omi_hj/ini_ens/';
inpath = '/home/hajsong/PDAF/tutorial_global_oce_latlon_v1.16_omi_hj/run/';
%% opath  = '/net/wave3/ganeshgopal/WAVE5_DATA/GOM_DART/EOF_GOM/NEWEOF/';
opath  = '/home/hajsong/PDAF/tutorial_global_oce_latlon_v1.16_omi_hj/ini_ens/';

%eval(['load ' fpath 'GRID.mat']);
hfacc = rdmds([inpath, 'hFacC']);
hfacs = rdmds([inpath, 'hFacS']);
hfacw = rdmds([inpath, 'hFacW']);
eval(['load ' fpath 'FMT.mat']);

mskc = hfacc * 1;
mskc(mskc<1) = 0;
msks = hfacs * 1;
msks(msks<1) = 0;
mskw = hfacw * 1;
mskw(mskw<1) = 0;
%% --> Read fields and Store in 1D vector.from 5 year model tun from 2004-2008

tf = 120;
dt = 30;
frame = [dt*(1:tf)];
nt = length(frame);

for i = 1:nt
f = frame(i);

%% Salt
s  = rdmds([inpath 'Stave'],f);
iz = find(mskc==1); 
x  = s(iz);
n  = length(x);
disp(['dimS = ' num2str(length(iz))]);

%% Temp
t  = rdmds([inpath 'Ttave'],f);
iz = find(mskc==1);
x(n+1:n+length(iz)) = t(iz);
n  = length(x);
disp(['dimT = ' num2str(length(iz))]);

%% Vel-U
u  = rdmds([inpath 'uVeltave'],f);
iz = find(mskw==1);
x(n+1:n+length(iz)) = u(iz);
n  = length(x);
disp(['dimU = ' num2str(length(iz))]);

%% Vel-V
v  = rdmds([inpath 'vVeltave'],f);
iz = find(msks==1);
x(n+1:n+length(iz)) = v(iz);
n  = length(x);
disp(['dimV = ' num2str(length(iz))]);

%% SSH
h  = rdmds([inpath 'ETAtave'],f);
iz = find(mskc(:,:,1)==1);
x(n+1:n+length(iz)) = h(iz);
disp(['dimH = ' num2str(length(iz))]);

%% Write output
wrslice([opath 'eof_vecs.bin'],x,i,fmt,Ieee);
disp(['dimension = ' num2str(length(x))]);
end  

display(['Number of vectors = ' num2str(nt)])
