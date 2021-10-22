%%  ****************************************************************
%%  *                                                              *
%%  *   Convert Initial Ensemble into DART format (3D variables)   *
%%  *                                                              *
%%  ****************************************************************


clear all;

%% Set fields dimensions.

dimstate = 107376;
dims = 27405; 
dimt = 27405; 
dimu = 25459; 
dimv = 24792; 
dimh = 2315;
     
if dimstate ~= dims+dimt+dimu+dimv+dimh
  disp('DIMENSIONS DO NOT MATCH');
end


%% Set number of EOFs.

neofs = 10;         


%% Set input and output paths.
fpath  = '/home/hajsong/PDAF/tutorial_global_oce_latlon_v1.16_omi_hj/ini_ens/';
inpath = '/home/hajsong/PDAF/tutorial_global_oce_latlon_v1.16_omi_hj/run/';
dpath  = fpath;
opath  = '/home/hajsong/Research/DART_LLC90/NEWEOF/';

eval(['load ' fpath 'FMT.mat']);
hfacc = rdmds([inpath, 'hFacC']);
hfacs = rdmds([inpath, 'hFacS']);
hfacw = rdmds([inpath, 'hFacW']);
xc = rdmds([inpath, 'XC'])';
yc = rdmds([inpath, 'YC'])';
[nxc, nyc, nzc] = size(hfacc);
 
%% Create masks for each variable.
%% Read fields.

mskc = hfacc * 1;
mskc(mskc<1) = 0;
msks = hfacs * 1;
msks(msks<1) = 0;
mskw = hfacw * 1;
mskw(mskw<1) = 0;

is = find(mskc==1);
it = find(mskc==1);
iu = find(mskw==1);
iv = find(msks==1);
ih = find(mskc(:,:,1)==1);

%% Set extension for output files names.
for i = 1:neofs
  ind = length(num2str(i));
  tistr(i,:) = '0000000000';
  tistr(i,10-ind+1:10) = num2str(i);
end


%% Convert EOFs into 3D and 2D seperate variables.
for i = 1:neofs
  disp(i);
  %% Read EOFs
  eof = rdslice([dpath 'ini_ens.bin'],[dimstate 1],i,'real*4');
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
  %% Save variables.
  wrslice([dpath 'S.'   tistr(i,:) '.data'],eofs3d,1,'real*4');
  wrslice([dpath 'T.'   tistr(i,:) '.data'],eoft3d,1,'real*4');
  wrslice([dpath 'U.'   tistr(i,:) '.data'],eofu3d,1,'real*4');
  wrslice([dpath 'V.'   tistr(i,:) '.data'],eofv3d,1,'real*4');
  wrslice([dpath 'Eta.' tistr(i,:) '.data'],eofh2d,1,'real*4');
  %% Plot.
  if i < 10
  eofs3d(eofs3d==0) = nan;
  eoft3d(eoft3d==0) = nan;
  eofu3d(eofu3d==0) = nan;
  eofv3d(eofv3d==0) = nan;
  eofh2d(eofh2d==0) = nan;
  figure
  set(gca,'fontsize',12)
  subplot(3,2,1);
  pcolor(xc,yc,eofs3d(:,:,1)');
  shading flat; colorbar;
  title(['Salinity: EOF-' num2str(i)]);
  subplot(3,2,2);
  pcolor(xc,yc,eoft3d(:,:,1)');
  shading flat; colorbar;
  title(['Temperature: EOF-' num2str(i)]);
  subplot(3,2,3);
  pcolor(xc,yc,eofu3d(:,:,1)');
  shading flat; colorbar;
  title(['Zonal Vel: EOF-' num2str(i)]);
  subplot(3,2,4);
  pcolor(xc,yc,eofv3d(:,:,1)');
  shading flat; colorbar;
  title(['Meridional Vel: EOF-' num2str(i)]);
  subplot(3,2,5);
  pcolor(xc,yc,eofh2d');
  shading flat; colorbar;
  title(['Sea Surface Height: EOF-' num2str(i)]);
%  eval(['print -dpsc  ' opath 'plot_gom_ens_' num2str(i)]);
  end
end


