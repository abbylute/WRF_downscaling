%% WRF Downscaling Master Script

% parameters
outTR = 4; % output temporal resolution (hrs)
outSR = 210; % output spatial resolution (m)
chunksize = outSR*500; % chunk size without buffer (m). should be a multiple of outSR.
buffer = 16000; % size of buffer around each chunk (m)
    % buffer should allow > than (window-1)/2 WRF grid cells to be
    % available on all sides for lapse rate calculations
window = 7; % number of WRF grid cells to use for lapse rate downscaling

% directories
mdir = '/Volumes/WDPassport/'; % master directory
wrfhdir = ; % hourly WRF data
wrfmdir = [mdir,'DATA/WRF/Monthly_means/']; % monthly WRF data
outdir = [mdir,'DATA/WRF/Downscaled/']; % output directory for downscaled data
inDEM = ; % filename for 4km WRF DEM
outDEM = [mdir,'DATA/DEM/NED/',num2str(outSR),'m/WUS_NED_',num2str(outSR),'m.mat']; % filename for fine resolution DEM
addpath([outdir,'Code/'])


%% Define spatial chunks

define_spatial_chunks(outDEM, outSR, chunksize, buffer, outdir)


%% For each spatial chunk,

load([outdir,'chunk_coordinates.mat']);
nchunk = size(chunk_coords.st_col,2);

for ch = 1:nchunk
% - calculate monthly lapse rates for LW, PSFC, Q2, WIND, PPT

% - calculate outTR hourly lapse rates for temperature

% - downscale outTR hourly datasets for all variables except solar, save

    downscale_WRF_lapse_rates(ch, outDEM, outSR, outTR, window, outdir, wrfhdir, wrfmdir)

    
% - in R, run solar terrain correction script, save

% - downscale solar using terrain corrections, save

% - bias correct downscaled ppt using PRISM. (OR DO THIS BEFORE DOWNSCALING?)

end

% OUTPUT
% outTR x outSR datasets for each chunk for each variable








