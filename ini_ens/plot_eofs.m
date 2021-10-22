%%  ****************************************************************
%%  *                                                              *
%%  *             Plot first few eofs for each variable            *
%%  *                                                              *
%%  ****************************************************************


clear all;
close all;
clc

%% Set fields dimensions.

dimstate = 3661614;
dims = 917225; 
dimt = 917225; 
dimu = 898205; 
dimv = 895215; 
dimh = 33744;
     
if dimstate ~= dims+dimt+dimu+dimv+dimh
  disp('DIMENSIONS DO NOT MATCH');
end


%% Set number of EOFs.

neofs = 10;         


%% Set input and output paths.
fpath  = '/net/wave3/ganeshgopal/GOM/FORCINGS_2004_2008/';
inpath = '/net/wave3/ganeshgopal/GOM/Tmp_NCEP_QSCAT_2004_2008/';
dpath  = '/net/wave3/ganeshgopal/WAVE5_DATA/GOM_DART/EOF_GOM/NEWEOF/';
opath = '/net/wave3/ganeshgopal/WAVE5_DATA/GOM_DART/EOF_GOM/NEWEOF/Figs/';

eval(['load ' fpath 'GRID.mat']);
eval(['load ' fpath 'FMT.mat']);



%% Create masks for each variable.
%% Read fields.
tf = 257;
dt = 672;
frame = [dt*(1:tf)];
nt = length(frame);

for i = 1
f = frame(i);
masks = rdmds([inpath 'S'],f);
maskt = rdmds([inpath 'T'],f);
masku = rdmds([inpath 'U'],f);
maskv = rdmds([inpath 'V'],f);
maskh = rdmds([inpath 'Eta'],f);
%% Find land points.
is = find(masks~=0);
it = find(maskt~=0);
iu = find(masku~=0);
iv = find(maskv~=0);
ih = find(maskh~=0);
end

%% Set extension for output files names.
for i = 1:neofs
  ind = length(num2str(i));
  tistr(i,:) = '0000000000';
  tistr(i,10-ind+1:10) = num2str(i);
end


%% Convert EOFs into 3D and 2D seperate variables.
for i = 1:neofs
  %% Read EOFs
  eof = rdslice([dpath 'eof_basis.bin'],[dimstate 1],i,'real*4');
  %% Save EOFs for each variable
  eofs = eof(1:dims);
  eoft = eof(dims+1:dims+dimt);
  eofu = eof(dims+dimt+1:dims+dimt+dimu);
  eofv = eof(dims+dimt+dimu+1:dims+dimt+dimu+dimv);
  eofh = eof(dims+dimt+dimu+dimv+1:dims+dimt+dimu+dimv+dimh);
  %% Initialize 3D and 2D arrays to store the EOFs for each variable.
  eofs3d = zeros(nxc,nyc,nzc);
  eoft3d = zeros(nxc,nyc,nzc);
  eofu3d = zeros(nxc,nyc,nzc);
  eofv3d = zeros(nxc,nyc,nzc);
  eofh2d = zeros(nxc,nyc);
  %% Make 3D and 2D fields.
  eofs3d(is) = eofs;
  eoft3d(it) = eoft;
  eofu3d(iu) = eofu;
  eofv3d(iv) = eofv;
  eofh2d(ih) = eofh;
  %5 Plot.
  eofs3d(eofs3d==0) = nan;
  eoft3d(eoft3d==0) = nan;
  eofu3d(eofu3d==0) = nan;
  eofv3d(eofv3d==0) = nan;
  eofh2d(eofh2d==0) = nan;
  figure
  set(gca,'fontsize',12)
  subplot(3,2,1);
  imagesc(xc,yc,eofs3d(:,:,1)');
  axis xy; colorbar;
  title(['Salinity: EOF-' num2str(i)]);
  subplot(3,2,2);
  imagesc(xc,yc,eoft3d(:,:,1)');
  axis xy; colorbar;
  title(['Temperature: EOF-' num2str(i)]);
  subplot(3,2,3);
  imagesc(xc,yc,eofu3d(:,:,1)');
  axis xy; colorbar;
  title(['Zonal Vel: EOF-' num2str(i)]);
  subplot(3,2,4);
  imagesc(xc,yc,eofv3d(:,:,1)');
  axis xy; colorbar;
  title(['Meridional Vel: EOF-' num2str(i)]);
  subplot(3,2,5);
  imagesc(xc,yc,eofh2d');
  axis xy; colorbar;
  title(['Sea Surface Height: EOF-' num2str(i)]);
  eval(['print -dpsc  ' opath 'plot_gom_eof_' num2str(i)]);
end


