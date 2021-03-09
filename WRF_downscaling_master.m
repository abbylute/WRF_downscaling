%% WRF Downscaling Master Script

% parameters
outTR = 4; % output temporal resolution (hrs)
outSRf = 210; % fine output spatial resolution (m)
outSRc = 1050; % coarse output spatial resolution (m)
chunksize = outSRf*400; % chunk size without buffer (m). should be a multiple of outSR.
buffer = 16000; % size of buffer around each chunk (m)
    % buffer should allow > than (window-1)/2 WRF grid cells to be
    % available on all sides for lapse rate calculations
window = 7; % number of WRF grid cells to use for lapse rate downscaling
era = 'CTRL'; % time period
elev_dif_thres = 50; % maximum elevation difference across outSRf resolution cells within each outSRc cell (m)
solar_dif_thres = 10; % as above, but for solar radiation from matlab (%)


% directories
%mdir = '/home/abby/'; % master directory on thunder
mdir = '/mnt/ceph/alute/'; % master directory on tesla
wrfhdir = [mdir,'DATA/WRF/WUS_',num2str(outTR),'hr/']; % hourly WRF data
wrfmdir = [mdir,'DATA/WRF/monthly/']; % monthly WRF data
prismppt = [mdir,'DATA/PRISM/Monthly_ppt_4km_Oct2000_Sep2013.mat'];
outdir = [mdir,'DATA/WRF/downscaled/WUS/']; % output directory for downscaled data
inDEM = [wrfhdir,'lon_lat_hgt_trimmed.mat']; % filename for 4km WRF DEM
outDEMf = [mdir,'DATA/Mapping/WUS_NED_',num2str(outSRf),'m.mat']; % filename for fine resolution output DEM
outDEMftif = [mdir,'DATA/Mapping/WUS_NED_',num2str(outSRf),'m.tif'];
outDEMc = [mdir,'DATA/Mapping/WUS_NED_',num2str(outSRc),'m.mat']; % filename for coarse resolution output DEM
outDEMctif = [mdir,'DATA/Mapping/WUS_NED_',num2str(outSRc),'m.tif'];
us_latlon = [mdir,'DATA/Mapping/US_latlon.mat']; % location of US latlon
addpath([mdir,'DATA/WRF/downscaled/Code/'])
%pathtoR = '/opt/modules/devel/R/3.6.0/lib64/R/bin/Rscript'; % location of R program on thunder. in R: file.path(R.home("bin"), "R")
pathtoR = '/opt/modules/devel/R/3.5.1/lib64/R/bin/Rscript'; % location on tesla
solartcRscript = [mdir,'DATA/WRF/downscaled/Code/get_solar_terrain_corrections.R']; % location of R solar terrain correction script
solarparamdir = [outdir,'solar_param_files/']; % file to store parameters temporarily for solar downscaling
reggmttz = -7; % GMT timezone of the full domain

%% Define spatial chunks

% at finer resolution:
%define_spatial_chunks(outDEMf, outSRf, chunksize, buffer, outdir, us_latlon)

% at medium resolution:
%define_spatial_chunks(outDEMc, outSRc, chunksize, buffer, outdir, us_latlon)


%% For each spatial chunk,

finechunks = matfile([outdir,'chunks/chunk_coordinates_',num2str(outSRf),'m.mat']);
finechunks = finechunks.chunk_coords;
in_us = finechunks.in_us;
nchunk = size(finechunks.st_col,2);
clear finechunks
ch_to_run = find(in_us>0);

%chdone = table2array(readtable('/mnt/ceph/alute/DATA/WRF/downscaled/WUS/CTRL/ACSWDNB/chunks_done.txt'));
%f = ismember(1:nchunk,chdone);
%ch_to_run = find(in_us>0 & ~f');
%nchunk = size(ch_to_run,1);
%in_us = in_us(~f);


parfor chu = 1:sum(in_us)
    ch = ch_to_run(chu);
  %ch=42; %ch=291;%NCASC; %ch=62;%Olympics %ch=45;%GNP
    
        % define where to model at what resolution
        [outlonf, outlatf, outlonc, outlatc] = pick_modeling_locations(ch, outDEMf, outDEMc, outSRf, outSRc, elev_dif_thres, solar_dif_thres, outdir);


        % downscale outTR hourly datasets dependent on spatial scale

        % only downscale at fine res
        tic
        if (isempty(outlonc) && ~isempty(outlonf)) 
            downscale_WRF_lapse_rates(ch, outSRf, inDEM, outDEMf, outDEMftif, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era, 'fine', solarparamdir, pathtoR, solartcRscript, reggmttz)

        % only downscale at coarse res
        elseif (~isempty(outlonc) && isempty(outlonf)) 
            downscale_WRF_lapse_rates(ch, outSRc, inDEM, outDEMc, outDEMctif, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era, 'coarse', solarparamdir, pathtoR, solartcRscript, reggmttz)    

        % downscale at both res    
        elseif (~isempty(outlonc) && ~isempty(outlonf)) 
            downscale_WRF_lapse_rates(ch, outSRf, inDEM, outDEMf, outDEMftif, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era, 'fine', solarparamdir, pathtoR, solartcRscript, reggmttz)
            downscale_WRF_lapse_rates(ch, outSRc, inDEM, outDEMc, outDEMctif, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era, 'coarse', solarparamdir, pathtoR, solartcRscript, reggmttz)

        else
            warning(['chunk ',num2str(ch),' did not have any points to downscale'])
        end
        toc
end








